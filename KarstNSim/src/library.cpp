#include "KarstNSim/library.h"

#include <ios>
#include <new>
#include <string>
#include <utility>

namespace KarstNSim {

	namespace {

		// Restores the caller's stream formatting, which engine messages modify.
		class LogFormatGuard {
		public:
			explicit LogFormatGuard(std::ostream* stream)
				: stream_(stream) {
				if (stream_ != nullptr) {
					flags_ = stream_->flags();
					precision_ = stream_->precision();
					fill_ = stream_->fill();
				}
			}
			~LogFormatGuard() {
				if (stream_ != nullptr) {
					stream_->flags(flags_);
					stream_->precision(precision_);
					stream_->fill(fill_);
				}
			}
			LogFormatGuard(const LogFormatGuard&) = delete;
			LogFormatGuard& operator=(const LogFormatGuard&) = delete;
		private:
			std::ostream* stream_;
			std::ios_base::fmtflags flags_{};
			std::streamsize precision_ = 0;
			char fill_ = ' ';
		};

		void reject_filesystem_flag(bool enabled, const char* flag, const char* files) {
			if (enabled) {
				throw InvalidInputError(std::string("[library] '") + flag + "' writes " + files +
					" and is not supported by run_simulation_memory; set it to false.");
			}
		}

		void validate_memory_request(const ParamsSource& parameters) {
			reject_filesystem_flag(parameters.create_vset_sampling, "create_vset_sampling",
				"<name>_pts.txt and <name>_s_pts.txt");
			reject_filesystem_flag(parameters.create_nghb_graph, "create_nghb_graph",
				"<name>_nghb_graph.txt");
			reject_filesystem_flag(parameters.create_nghb_graph_property, "create_nghb_graph_property",
				"<name>_nghb_graph.txt properties");
			reject_filesystem_flag(parameters.create_grid, "create_grid", "<name>_box.txt");

			if (!parameters.sections_simulation_only && parameters.use_user_connectivity_matrix &&
				parameters.connectivity_matrix.empty()) {
				throw InvalidInputError(
					"[library] 'use_user_connectivity_matrix' is true but ParamsSource::connectivity_matrix "
					"is empty. run_simulation_memory never reads connectivity_matrix.txt: supply the matrix "
					"in memory, or set use_user_connectivity_matrix to false for automatic connectivity.");
			}

			validate_simulation_parameters(parameters);
		}

		// Runs `step`, translating non-library exceptions: invalid input before the engine
		// starts, SimulationFailed afterwards. Library errors and std::bad_alloc pass through.
		template <typename Step>
		auto translate_errors(ErrorCode code, Step&& step) -> decltype(step()) {
			try {
				return step();
			}
			catch (const Error&) {
				throw;
			}
			catch (const std::bad_alloc&) {
				throw;
			}
			catch (const std::exception& error) {
				if (code == ErrorCode::InvalidInput) {
					throw InvalidInputError(error.what());
				}
				throw Error(code, error.what());
			}
		}
	}

	std::vector<KarstNetworkResult> run_simulation_memory(ParamsSource parameters, const RunOptions& options) {
		detail::JobContext job;
		struct WorkReport {
			const detail::JobContext& job;
			std::uint64_t* output;
			~WorkReport() { if (output) *output = job.work_done(); }
		} work_report{job, options.work_used};
		job.out = options.log;
		job.err = options.log;
		job.allow_filesystem = false;
		job.work_limit = options.work_limit;
		job.max_points = options.max_points;
		job.max_edges = options.max_edges;
		job.is_cancelled = options.is_cancelled;
		job.user = options.user;

		LogFormatGuard log_format(options.log);
		detail::ScopedJob job_scope(job);

		translate_errors(ErrorCode::InvalidInput, [&] { validate_memory_request(parameters); });
		job.checkpoint();

		std::vector<KarstNetworkResult> results;
		// Grow only as iterations finish. A very large requested iteration count
		// must not reserve unbounded memory before work/cancellation checks apply.

		for (int i = 0; i < parameters.number_of_iterations; i++) {
			const unsigned int used_seed = detail::iteration_seed(parameters, i);

			detail::log_out() << "========== SIMULATION " << i
				<< " STARTED WITH SEED " << used_seed
				<< " ==========\n\n";

			initializeRng(std::vector<std::uint32_t>{ used_seed });

			GeologicalParameters params;
			std::vector<KeyPoint> keypts;
			const std::string sim_name_iter = detail::iteration_name(parameters, i);

			KarsticNetwork karst(sim_name_iter, parameters.domain, params, keypts, parameters.surf_wat_table);
			translate_errors(ErrorCode::InvalidInput, [&] {
				detail::configure_network(karst, parameters, sim_name_iter, /*in_memory=*/true);
			});
			job.checkpoint();

			std::optional<KarstNetworkResult> result = translate_errors(ErrorCode::SimulationFailed, [&] {
				return detail::run_network(karst, parameters);
			});
			if (!result.has_value()) {
				throw NoRouteError("[library] No path found between inlets and outlets for iteration " +
					std::to_string(i) + " (seed " + std::to_string(used_seed) + ").", i);
			}
			results.push_back(std::move(*result));
			job.checkpoint();
		}

		return results;
	}
}
