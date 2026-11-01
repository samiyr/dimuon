#ifndef DM_PDF_IS_PDF_INTERFACE_H
#define DM_PDF_IS_PDF_INTERFACE_H

#include <array>
#include <concepts>
#include <string>
#include <cstddef>
#include <utility>

#include "dimuon/common/constants.hpp"
#include "dimuon/common/flavor.hpp"

#include "dimuon/pdf/interface/lha_interface.hpp"
#include "dimuon/pdf/interface/functional_form_interface.hpp"
#include "dimuon/pdf/interface/vector_kernel_interface.hpp"

#include "dimuon/utility/collections.hpp"

namespace dm {

/// A concept of a PDF (or FF) whose values are numbers, such as a PDF set (`lha_interface`)
template <typename pdf_t>
concept is_scalar_pdf_interface = requires(pdf_t pdf, const flavor::type flavor, const double x, const double Q2) {
	{ pdf.xfxQ2(flavor, x, Q2) } -> std::convertible_to<double>;
	{ pdf.xfxQ2(x, Q2) } -> std::convertible_to<std::array<double, constants::total_flavors>>;
	{ pdf.alpha_s(Q2) } -> std::convertible_to<double>;
	{ pdf.Q2_min() } -> std::convertible_to<double>;
	{ pdf.Q2_max() } -> std::convertible_to<double>;
	{ pdf.quark_mass(flavor) } -> std::convertible_to<double>;
	{ pdf.get_set_name() } -> std::convertible_to<std::string>;
	{ pdf.get_set_member() } -> std::convertible_to<std::size_t>;
	{ pdf.dimension() } -> std::convertible_to<std::size_t>;
};

/// A concept of a PDF (or FF) whose values are vectors
template <typename pdf_t>
concept is_vector_pdf_interface = requires(pdf_t pdf, const flavor::type flavor, const double x, const double Q2) {
	{ pdf.xfxQ2(flavor, x, Q2) } -> collections::range_convertible_to<std::pair<std::size_t, double>>;
	{ pdf.xfxQ2(x, Q2) } -> collections::range_convertible_to<std::pair<std::size_t, double>>;
	{ pdf.alpha_s(Q2) } -> std::convertible_to<double>;
	{ pdf.dimension() } -> std::convertible_to<std::size_t>;
};

/// A concept of a PDF (or FF) whose values are either numbers (`is_scalar_pdf_interface`) or vectors (`is_vector_pdf_interface`)
template <typename pdf_t>
concept is_pdf_interface = is_scalar_pdf_interface<pdf_t> || is_vector_pdf_interface<pdf_t>;

static_assert(is_scalar_pdf_interface<lha_interface<>>, "Type `lha_interface` doesn't satisfy concept `is_scalar_pdf_interface`.");
static_assert(
	is_scalar_pdf_interface<functional_form_interface<decltype([](const flavor::type, const double, const double) { return 0.0; })>>,
	"Type `functional_form_interface` doesn't satisfy concept `is_scalar_pdf_interface`."
);
static_assert(is_vector_pdf_interface<vector_kernel_interface>, "Type `vector_kernel_interface` doesn't satisfy concept `is_vector_pdf_interface`.");

}

#endif