#pragma once

/**
@file job_context.h
@brief Internal per-job state shared by the engine: RNG, logger, limits and I/O policy.

A JobContext is owned by the caller's stack frame (run_simulation_full() or
run_simulation_memory()). ScopedJob publishes a pointer to it in a thread_local
slot for the duration of the job and restores the previous pointer on exit,
including exceptional exit, so concurrent jobs on different threads and nested
jobs on one thread never share RNG, noise, logger, counters or limits.

This header is internal plumbing; library users only need library.h.
**/

#include <cstdint>
#include <limits>
#include <ostream>
#include <random>
#include <vector>

#include "KarstNSim/errors.h"

namespace KarstNSim {
	namespace detail {

		class JobContext {
		public:
			static constexpr std::uint64_t kUnlimited = std::numeric_limits<std::uint64_t>::max();
			//! Number of work units between two cancellation polls.
			static constexpr std::uint64_t kPollInterval = 4096;

			std::mt19937 rng; //!< The job RNG. Reseeded at the start of every iteration.
			std::vector<std::uint32_t> seed; //!< Seed of the current iteration.
			std::ostream* out = nullptr; //!< Informational log; null means quiet.
			std::ostream* err = nullptr; //!< Diagnostic log; null means quiet.
			bool allow_filesystem = true; //!< False for the in-memory API: any file read or write by the engine throws std::logic_error.

			std::uint64_t work_limit = kUnlimited;
			std::uint64_t max_points = kUnlimited;
			std::uint64_t max_edges = kUnlimited;
			bool (*is_cancelled)(void*) = nullptr;
			void* user = nullptr;

			/*!
			\brief Adds deterministic work units, enforces work_limit and polls cancellation periodically.
			*/
			void tick(std::uint64_t units = 1) {
				work_done_ = units > kUnlimited - work_done_ ? kUnlimited : work_done_ + units;
				if (work_done_ > work_limit) {
					fail_work_limit();
				}
				if (work_done_ >= next_poll_) {
					poll();
				}
			}

			/*! \brief Polls cancellation immediately. */
			void checkpoint();

			/*! \brief Throws PointLimitExceeded if count > max_points. */
			void admit_points(std::uint64_t count, const char* stage) const;

			/*! \brief Throws EdgeLimitExceeded if slots > max_edges. */
			void admit_edges(std::uint64_t slots, const char* stage) const;

			std::uint64_t work_done() const noexcept { return work_done_; }

		private:
			void poll();
			[[noreturn]] void fail_work_limit() const;

			std::uint64_t work_done_ = 0;
			std::uint64_t next_poll_ = kPollInterval;
		};

		/*! \brief Innermost job active on this thread, or null. */
		JobContext* current_job() noexcept;

		/*! \brief Innermost job active on this thread; throws std::logic_error outside a job. */
		JobContext& require_job(const char* caller);

		/*!
		\brief Installs a stack-owned JobContext as the current job of this thread.
		The previous job pointer is restored by the destructor.
		*/
		class ScopedJob {
		public:
			explicit ScopedJob(JobContext& context) noexcept;
			~ScopedJob();
			ScopedJob(const ScopedJob&) = delete;
			ScopedJob& operator=(const ScopedJob&) = delete;
		private:
			JobContext* previous_;
		};

		/*! \brief Informational log of the current job; a per-thread discarding stream when quiet or outside a job. */
		std::ostream& log_out();

		/*! \brief Diagnostic log of the current job; a per-thread discarding stream when quiet or outside a job. */
		std::ostream& log_err();

		/*! \brief JobContext::tick() on the current job; no-op outside a job. */
		inline void tick(std::uint64_t units = 1) {
			if (JobContext* job = current_job()) {
				job->tick(units);
			}
		}

		/*! \brief JobContext::checkpoint() on the current job; no-op outside a job. */
		void checkpoint();

		/*! \brief JobContext::admit_points() on the current job; no-op outside a job. */
		void admit_points(std::uint64_t count, const char* stage);

		/*! \brief JobContext::admit_edges() on the current job; no-op outside a job. */
		void admit_edges(std::uint64_t slots, const char* stage);

		/*! \brief True outside a job and for jobs that allow filesystem access. */
		bool filesystem_enabled() noexcept;

		/*! \brief Throws std::logic_error when the current job forbids filesystem access. */
		void require_filesystem(const char* what);
	}
}
