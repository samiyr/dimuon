#ifndef DM_SIDIS_DISTRIBUTION_MATRIX_H
#define DM_SIDIS_DISTRIBUTION_MATRIX_H

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <ranges>
#include <span>
#include <type_traits>

#include "dimuon/common/constants.hpp"
#include "dimuon/common/flavor.hpp"

#include "dimuon/interpolation/linear_interpolation_basis.hpp"

#include "dimuon/pdf/is_pdf_interface.hpp"

#include "dimuon/utility/shims.hpp"
#include "dimuon/utility/sparse_vector.hpp"

namespace dm {

namespace sidis {
	/// The values of a parton distribution, a PDF or FF whose values are numbers (see `is_scalar_pdf_interface`), at \p column_count momentum fractions:
	/// a matrix with a row for each flavor and a column for each momentum fraction. The values of one flavor are its row (see `rows`).
	template <std::size_t column_count>
	class parton_distribution_matrix {
	public:
		/// A row of the matrix: the values of a flavor at each momentum fraction
		using row_t = std::array<double, column_count>;
		/// Type alias for the values of all flavors, indexed by `flavor::to_index` (see the all-flavor `xfxQ2` of the PDF interfaces)
		using column_t = std::array<double, constants::total_flavors>;

		/// Constructs the matrix whose columns are \p column_values, the values of all flavors at each momentum fraction
		explicit parton_distribution_matrix(const std::same_as<column_t> auto &...column_values) noexcept
			requires (sizeof...(column_values) == column_count)
			: values{column_values...} {}

		/// Returns the values of \p flavor, its row, as a matrix with a row for each component of the vector of values, of which there is one:
		/// a sparse vector of the rows with the single row 0, or no rows if the values are zero (see `single_element_sparse_vector`)
		auto rows(const flavor::type flavor) const noexcept {
			row_t row;

			for (const auto &[column_idx, column] : std::views::enumerate(values)) {
				row[column_idx] = column[flavor::to_index(flavor)];
			}

			return single_element_sparse_vector(0, row);
		}

	private:
		/// The values of all flavors in each column
		const std::array<column_t, column_count> values;
	};

	/// The values of the interpolation kernels, or of other functions in the place of a PDF or FF whose values are vectors (see `is_vector_pdf_interface`),
	/// at \p column_count momentum fractions: a matrix with a row for each kernel and a column for each momentum fraction, of which only the rows of the non-zero kernels
	/// are kept. The kernels are the same for all flavors, so where the matrix of a parton distribution has a row for each flavor (see `parton_distribution_matrix`),
	/// this one has a row for each kernel, and the values of every flavor are the whole matrix (see `rows`).
	template <std::size_t column_count, std::size_t max_rows = linear_interpolation_basis::maximum_kernels_per_point * column_count>
	class kernel_distribution_matrix {
	public:
		/// A row of the matrix: a kernel and its values at each momentum fraction
		struct row_t {
			kernel_index index;
			std::array<double, column_count> values;
		};
		using column_t = vector_kernel_interface::values;

		/// Constructs the matrix whose columns are \p column_values, sparse vectors over the kernels (see `is_vector_pdf_interface`):
		/// a row for each kernel that is non-zero in any column
		explicit kernel_distribution_matrix(const std::same_as<column_t> auto &...column_values)
			requires (sizeof...(column_values) == column_count)
			: stored_rows{make_rows(column_values...)} {}

		/// Returns the values of every flavor: the whole matrix, as a range of its stored rows
		std::span<const row_t> rows(const flavor::type) const noexcept { return std::span{stored_rows}; }

	private:
		/// The rows of the non-zero kernels, in no particular order. At most two kernels are non-zero in each column (see `vector_kernel_interface::values`),
		/// so a matrix has at most `2 * column_count` rows (more throw `std::bad_alloc`).
		using storage_t = shims::inplace_vector<row_t, max_rows>;

		/// The rows of the non-zero kernels
		const storage_t stored_rows;

		/// Returns the rows of the matrix whose columns are \p column_values: a row for each kernel that is non-zero in any column
		static storage_t make_rows(const std::same_as<column_t> auto &...column_values) {
			storage_t rows_;
			
			// Each column contains the values of the non-zero kernels at some momentum fraction. The kernel indices depend on the momentum fraction
			// and are not ordered. Because of these facts, this loop has to search for each kernel index to assign a value.
			for (const auto &[column_idx, column] : std::views::enumerate(std::array{column_values...})) {
				for (const auto &[kernel_idx, value] : column) {

					const auto location = std::ranges::find(rows_, kernel_idx, &row_t::index);

					if (location != rows_.end()) {
						location->values[column_idx] = value;
					} else {
						rows_.emplace_back(kernel_idx, decltype(row_t::values){}).values[column_idx] = value;
					}
				}
			}

			return rows_;
		}
	};

	namespace distribution_matrix {
		/// The values of a PDF or FF of type \p distribution_t at \p columns momentum fractions, as a matrix: a row for each flavor for a parton distribution,
		/// whose values are numbers (`parton_distribution_matrix`), and a row for each kernel for the interpolation kernels, whose values are vectors
		/// that are the same for all flavors (`kernel_distribution_matrix`)
		template <is_pdf_interface distribution_t, std::size_t columns>
		using type = std::conditional_t<is_vector_pdf_interface<distribution_t>, kernel_distribution_matrix<columns>, parton_distribution_matrix<columns>>;

		/// @brief Returns the values of a PDF or FF at \p momentum_fraction for all flavors, or zero if \p momentum_fraction is above one
		/// (i.e. where PDFs and FFs vanish). This gives the lower bounds of the convolution integrals, e.g. `z / zeta <= 1` for `zeta >= z`.
		/// @param[in] distribution PDF or FF
		/// @param[in] momentum_fraction Momentum fraction
		/// @param[in] scale Momentum scale (GeV^2)
		auto compute_values(const is_pdf_interface auto &distribution, const double momentum_fraction, const double scale) {
			using values_t = decltype(distribution.xfxQ2(momentum_fraction, scale));
			return momentum_fraction <= 1.0 ? distribution.xfxQ2(momentum_fraction, scale) : values_t{};
		}

	}

}

}

#endif
