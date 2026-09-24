/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.

***************************************************************/

#pragma once

/*!
\file packed_adjacency.h
\brief Compact fixed-shape directed adjacency storage for the cost graph.

Targets are stored in one contiguous int32 array (row-major, `cols()` slots per row) and
edge costs in one contiguous float array (edge-major, `channels()` floats per edge).
A slot whose target is -1 is unused. Views and references returned by this class are
lightweight (no allocation) and remain valid only until the next reset() or destruction.
*/

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace KarstNSim {

	class PackedAdjacency
	{
	public:
		static constexpr std::int32_t no_target = -1; //!< Target value of an unused slot.

		/*!
		\brief Mutable view over the cost channels of one edge.
		*/
		class WeightView
		{
		public:
			WeightView() = default;
			WeightView(float* data, std::size_t size) : data_(size ? data : nullptr), size_(size) {}

			float* begin() const { return data_; }
			float* end() const { return size_ ? data_ + size_ : data_; }
			float* data() const { return data_; }
			std::size_t size() const { return size_; }
			bool empty() const { return size_ == 0; }
			float& operator[](std::size_t i) const { assert(i < size_); return data_[i]; }
			float& at(std::size_t i) const {
				if (i >= size_) throw std::out_of_range("PackedAdjacency::WeightView::at: channel out of range");
				return data_[i];
			}

		private:
			float* data_ = nullptr;
			std::size_t size_ = 0;
		};

		/*!
		\brief Read-only view over the cost channels of one edge.
		*/
		class ConstWeightView
		{
		public:
			ConstWeightView() = default;
			ConstWeightView(const float* data, std::size_t size) : data_(size ? data : nullptr), size_(size) {}
			ConstWeightView(WeightView v) : data_(v.data()), size_(v.size()) {}

			const float* begin() const { return data_; }
			const float* end() const { return size_ ? data_ + size_ : data_; }
			const float* data() const { return data_; }
			std::size_t size() const { return size_; }
			bool empty() const { return size_ == 0; }
			const float& operator[](std::size_t i) const { assert(i < size_); return data_[i]; }
			const float& at(std::size_t i) const {
				if (i >= size_) throw std::out_of_range("PackedAdjacency::ConstWeightView::at: channel out of range");
				return data_[i];
			}

		private:
			const float* data_ = nullptr;
			std::size_t size_ = 0;
		};

		//! Mutable reference to one slot.
		struct EdgeRef
		{
			std::int32_t& target;
			WeightView weight;
		};

		//! Read-only reference to one slot.
		struct ConstEdgeRef
		{
			const std::int32_t& target;
			ConstWeightView weight;
		};

		//! Mutable proxy over the slots of one row.
		class RowRef
		{
		public:
			RowRef(PackedAdjacency& owner, std::size_t row) : owner_(&owner), row_(row) {}
			std::size_t size() const { return owner_->cols(); }
			EdgeRef operator[](std::size_t slot) const { return (*owner_)(row_, slot); }

		private:
			PackedAdjacency* owner_;
			std::size_t row_;
		};

		//! Read-only proxy over the slots of one row.
		class ConstRowRef
		{
		public:
			ConstRowRef(const PackedAdjacency& owner, std::size_t row) : owner_(&owner), row_(row) {}
			std::size_t size() const { return owner_->cols(); }
			ConstEdgeRef operator[](std::size_t slot) const { return (*owner_)(row_, slot); }

		private:
			const PackedAdjacency* owner_;
			std::size_t row_;
		};

		PackedAdjacency() = default;
		PackedAdjacency(const PackedAdjacency&) = default;
		PackedAdjacency& operator=(const PackedAdjacency&) = default;
		PackedAdjacency(PackedAdjacency&& other) noexcept { swap(other); }
		PackedAdjacency& operator=(PackedAdjacency&& other) noexcept {
			if (this != &other) {
				PackedAdjacency tmp(std::move(other));
				swap(tmp);
			}
			return *this;
		}

		void swap(PackedAdjacency& other) noexcept {
			std::swap(rows_, other.rows_);
			std::swap(cols_, other.cols_);
			std::swap(channels_, other.channels_);
			targets_.swap(other.targets_);
			weights_.swap(other.weights_);
		}

		/*!
		\brief Replaces the whole graph with `rows` x `columns` unused slots of `channels` zero costs.
		Throws std::length_error if the shape cannot be represented. On failure the graph is unchanged.
		Rows must fit in int32 because targets are row indices.
		*/
		void reset(std::size_t rows, std::size_t columns, std::size_t channels) {
			const std::size_t max_size = std::numeric_limits<std::size_t>::max();
			if (rows > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
				throw std::length_error("PackedAdjacency::reset: row count exceeds int32 target range");
			if (columns != 0 && rows > max_size / columns)
				throw std::length_error("PackedAdjacency::reset: rows * columns overflows");
			const std::size_t edges = rows * columns;
			if (channels != 0 && edges > max_size / channels)
				throw std::length_error("PackedAdjacency::reset: edges * channels overflows");
			const std::size_t weight_count = edges * channels;
			if (edges > std::vector<std::int32_t>().max_size() || weight_count > std::vector<float>().max_size())
				throw std::length_error("PackedAdjacency::reset: graph too large");

			std::vector<std::int32_t> targets(edges, no_target);
			std::vector<float> weights(weight_count, 0.0f);
			targets_.swap(targets);
			weights_.swap(weights);
			rows_ = rows;
			cols_ = columns;
			channels_ = channels;
		}

		std::size_t size() const { return rows_; }         //!< Number of rows (source nodes).
		std::size_t cols() const { return cols_; }         //!< Number of slots per row.
		std::size_t channels() const { return channels_; } //!< Number of cost channels per edge.

		//! Bytes actually reserved by the target and cost arrays.
		std::size_t storage_bytes() const {
			return targets_.capacity() * sizeof(std::int32_t) + weights_.capacity() * sizeof(float);
		}

		//! Unchecked slot access (assertions only).
		EdgeRef operator()(std::size_t row, std::size_t slot) {
			assert(row < rows_ && slot < cols_);
			const std::size_t e = row * cols_ + slot;
			return EdgeRef{ targets_[e], WeightView(weight_ptr(e), channels_) };
		}
		ConstEdgeRef operator()(std::size_t row, std::size_t slot) const {
			assert(row < rows_ && slot < cols_);
			const std::size_t e = row * cols_ + slot;
			return ConstEdgeRef{ targets_[e], ConstWeightView(weight_ptr(e), channels_) };
		}

		RowRef operator[](std::size_t row) { assert(row < rows_); return RowRef(*this, row); }
		ConstRowRef operator[](std::size_t row) const { assert(row < rows_); return ConstRowRef(*this, row); }
		RowRef row(std::size_t r) { return (*this)[r]; }
		ConstRowRef row(std::size_t r) const { return (*this)[r]; }

		//! Bounds-checked slot access. Throws std::out_of_range.
		EdgeRef at(std::size_t row, std::size_t slot) { check_slot(row, slot, "at"); return (*this)(row, slot); }
		ConstEdgeRef at(std::size_t row, std::size_t slot) const { check_slot(row, slot, "at"); return (*this)(row, slot); }

		/*!
		\brief Writes the target and all cost channels of one slot.
		`target` must be -1 (unused) or a valid row index; `weights.size()` must equal channels().
		Throws std::out_of_range / std::invalid_argument; nothing is written on failure.
		*/
		void set(std::size_t row, std::size_t slot, int target, const std::vector<float>& weights) {
			check_slot(row, slot, "set");
			check_channels(weights, "set");
			if (target < no_target || (target >= 0 && static_cast<std::size_t>(target) >= rows_))
				throw std::out_of_range("PackedAdjacency::set: target out of range");
			const std::size_t e = row * cols_ + slot;
			targets_[e] = static_cast<std::int32_t>(target);
			copy_weights(e, weights);
		}

		//! Writes all cost channels of one slot, keeping its target. Same validation as set().
		void set_weights(std::size_t row, std::size_t slot, const std::vector<float>& weights) {
			check_slot(row, slot, "set_weights");
			check_channels(weights, "set_weights");
			copy_weights(row * cols_ + slot, weights);
		}

	private:
		float* weight_ptr(std::size_t edge) { return channels_ ? weights_.data() + edge * channels_ : nullptr; }
		const float* weight_ptr(std::size_t edge) const { return channels_ ? weights_.data() + edge * channels_ : nullptr; }

		void copy_weights(std::size_t edge, const std::vector<float>& weights) {
			float* dst = weight_ptr(edge);
			for (std::size_t c = 0; c < channels_; ++c)
				dst[c] = weights[c];
		}

		void check_slot(std::size_t row, std::size_t slot, const char* where) const {
			if (row >= rows_ || slot >= cols_)
				throw std::out_of_range(std::string("PackedAdjacency::") + where + ": slot out of range");
		}

		void check_channels(const std::vector<float>& weights, const char* where) const {
			if (weights.size() != channels_)
				throw std::invalid_argument(std::string("PackedAdjacency::") + where + ": channel count mismatch");
		}

		std::size_t rows_ = 0;
		std::size_t cols_ = 0;
		std::size_t channels_ = 0;
		std::vector<std::int32_t> targets_;
		std::vector<float> weights_;
	};

} // namespace KarstNSim
