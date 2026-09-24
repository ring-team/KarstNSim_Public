#include "KarstNSim/job_context.h"

#include <stdexcept>
#include <string>

namespace KarstNSim {

	const char* to_string(ErrorCode code) noexcept {
		switch (code) {
		case ErrorCode::InvalidInput: return "invalid_input";
		case ErrorCode::Cancelled: return "cancelled";
		case ErrorCode::WorkLimitExceeded: return "work_limit_exceeded";
		case ErrorCode::PointLimitExceeded: return "point_limit_exceeded";
		case ErrorCode::EdgeLimitExceeded: return "edge_limit_exceeded";
		case ErrorCode::NoRoute: return "no_route";
		case ErrorCode::SimulationFailed: return "simulation_failed";
		}
		return "unknown";
	}

	namespace detail {

		namespace {
			thread_local JobContext* current_job_ = nullptr;

			std::ostream& discard_stream() {
				// One per thread: formatting state and error bits are never shared across jobs.
				thread_local std::ostream stream(nullptr);
				return stream;
			}
		}

		void JobContext::checkpoint() {
			next_poll_ = work_done_;
			poll();
		}

		void JobContext::poll() {
			next_poll_ = work_done_ > kUnlimited - kPollInterval ? kUnlimited : work_done_ + kPollInterval;
			if (is_cancelled != nullptr && is_cancelled(user)) {
				throw CancelledError("KarstNSim job cancelled after " +
					std::to_string(work_done_) + " work units.");
			}
		}

		void JobContext::fail_work_limit() const {
			throw LimitExceededError(ErrorCode::WorkLimitExceeded,
				"KarstNSim work limit exceeded: " + std::to_string(work_done_) +
				" work units > work_limit " + std::to_string(work_limit) + ".",
				work_limit, work_done_);
		}

		void JobContext::admit_points(std::uint64_t count, const char* stage) const {
			if (count > max_points) {
				throw LimitExceededError(ErrorCode::PointLimitExceeded,
					std::string("KarstNSim point limit exceeded during ") + stage + ": " +
					std::to_string(count) + " points > max_points " + std::to_string(max_points) + ".",
					max_points, count);
			}
		}

		void JobContext::admit_edges(std::uint64_t slots, const char* stage) const {
			if (slots > max_edges) {
				throw LimitExceededError(ErrorCode::EdgeLimitExceeded,
					std::string("KarstNSim edge limit exceeded during ") + stage + ": " +
					std::to_string(slots) + " directed edge slots > max_edges " +
					std::to_string(max_edges) + ".",
					max_edges, slots);
			}
		}

		JobContext* current_job() noexcept {
			return current_job_;
		}

		JobContext& require_job(const char* caller) {
			if (current_job_ == nullptr) {
				throw std::logic_error(std::string(caller) +
					" requires an active KarstNSim job (use run_simulation_full or run_simulation_memory).");
			}
			return *current_job_;
		}

		ScopedJob::ScopedJob(JobContext& context) noexcept
			: previous_(current_job_) {
			current_job_ = &context;
		}

		ScopedJob::~ScopedJob() {
			current_job_ = previous_;
		}

		std::ostream& log_out() {
			JobContext* job = current_job_;
			return job != nullptr && job->out != nullptr ? *job->out : discard_stream();
		}

		std::ostream& log_err() {
			JobContext* job = current_job_;
			return job != nullptr && job->err != nullptr ? *job->err : discard_stream();
		}

		void checkpoint() {
			if (JobContext* job = current_job_) {
				job->checkpoint();
			}
		}

		void admit_points(std::uint64_t count, const char* stage) {
			if (JobContext* job = current_job_) {
				job->admit_points(count, stage);
			}
		}

		void admit_edges(std::uint64_t slots, const char* stage) {
			if (JobContext* job = current_job_) {
				job->admit_edges(slots, stage);
			}
		}

		bool filesystem_enabled() noexcept {
			return current_job_ == nullptr || current_job_->allow_filesystem;
		}

		void require_filesystem(const char* what) {
			if (!filesystem_enabled()) {
				throw std::logic_error(std::string("KarstNSim in-memory job attempted filesystem access: ") + what);
			}
		}
	}
}
