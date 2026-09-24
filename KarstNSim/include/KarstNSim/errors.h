#pragma once

/**
@file errors.h
@brief Typed errors reported by the in-memory KarstNSim library API.

Every error thrown by KarstNSim::run_simulation_memory() derives from
KarstNSim::Error (itself a std::runtime_error) and carries an ErrorCode.
std::bad_alloc is not translated and propagates unchanged.
**/

#include <cstdint>
#include <stdexcept>
#include <string>

namespace KarstNSim {

	/*!
	\brief Stable numeric codes for library failures.
	*/
	enum class ErrorCode : int {
		InvalidInput = 1, //!< Parameters, geometry, connectivity or options are invalid, or a filesystem-only feature was requested.
		Cancelled = 2, //!< RunOptions::is_cancelled returned true at a checkpoint.
		WorkLimitExceeded = 3, //!< The deterministic work counter exceeded RunOptions::work_limit.
		PointLimitExceeded = 4, //!< The sampling cloud (or a surface point set) would exceed RunOptions::max_points.
		EdgeLimitExceeded = 5, //!< The cost graph would allocate more directed edge slots than RunOptions::max_edges.
		NoRoute = 6, //!< No inlet-to-outlet path could be computed for an iteration.
		SimulationFailed = 7 //!< Any other engine failure after the inputs were accepted.
	};

	/*!
	\brief Returns a short identifier for an error code, for example "invalid_input".
	*/
	const char* to_string(ErrorCode code) noexcept;

	/*!
	\brief Base class of all library errors.
	*/
	class Error : public std::runtime_error {
	public:
		Error(ErrorCode code, const std::string& message)
			: std::runtime_error(message), code_(code) {
		}
		ErrorCode code() const noexcept { return code_; }
	private:
		ErrorCode code_;
	};

	/*! \brief ErrorCode::InvalidInput. */
	class InvalidInputError : public Error {
	public:
		explicit InvalidInputError(const std::string& message)
			: Error(ErrorCode::InvalidInput, message) {
		}
	};

	/*! \brief ErrorCode::Cancelled. */
	class CancelledError : public Error {
	public:
		explicit CancelledError(const std::string& message)
			: Error(ErrorCode::Cancelled, message) {
		}
	};

	/*!
	\brief ErrorCode::WorkLimitExceeded, PointLimitExceeded or EdgeLimitExceeded.
	\details limit() is the configured limit and requested() the count that would have exceeded it.
	*/
	class LimitExceededError : public Error {
	public:
		LimitExceededError(ErrorCode code, const std::string& message,
			std::uint64_t limit, std::uint64_t requested)
			: Error(code, message), limit_(limit), requested_(requested) {
		}
		std::uint64_t limit() const noexcept { return limit_; }
		std::uint64_t requested() const noexcept { return requested_; }
	private:
		std::uint64_t limit_;
		std::uint64_t requested_;
	};

	/*!
	\brief ErrorCode::NoRoute. iteration() is the zero-based iteration without any route.
	*/
	class NoRouteError : public Error {
	public:
		NoRouteError(const std::string& message, int iteration)
			: Error(ErrorCode::NoRoute, message), iteration_(iteration) {
		}
		int iteration() const noexcept { return iteration_; }
	private:
		int iteration_;
	};
}
