#ifndef DM_PDF_INTERFACE_FUNCTIONAL_FORM_INTERFACE_H
#define DM_PDF_INTERFACE_FUNCTIONAL_FORM_INTERFACE_H

#include <array>
#include <limits>
#include <string>
#include <type_traits>

#include "dimuon/common/constants.hpp"
#include "dimuon/common/flavor.hpp"

namespace dm {

template <typename function_t>
class functional_form_interface {
public:
	std::string get_set_name() const { return name; }
	std::size_t get_set_member() const noexcept { return 0; }

	constexpr explicit functional_form_interface(function_t function_, std::string name_ = "functional_form") noexcept(std::is_nothrow_move_constructible_v<function_t>)
	: function{std::move(function_)}, name{std::move(name_)} { }

	constexpr double xfxQ2(const flavor::type flavor, const double x, const double Q2) const noexcept(noexcept(function(flavor, x, Q2))) {
		return function(flavor, x, Q2);
	}
	/// Computes the values of all flavors at once, indexed by `flavor::to_index`
	std::array<double, constants::total_flavors> xfxQ2(const double x, const double Q2) const noexcept(noexcept(function(flavor::g, x, Q2))) {
		std::array<double, constants::total_flavors> values;
		for (const flavor::type flavor : flavor::all_flavors) {
			values[flavor::to_index(flavor)] = function(flavor, x, Q2);
		}
		return values;
	}
	constexpr double alpha_s(const double) const noexcept {
		return 1.0;
	}
	constexpr double quark_mass(const flavor::type) const noexcept {
		return 0.0;
	}
	constexpr double Q2_min() const noexcept {
		return 0.0;
	}
	constexpr double Q2_max() const noexcept {
		return std::numeric_limits<double>::max();
	}
	
	std::size_t dimension() const noexcept { return 1; }

private:
	function_t function;
	std::string name;
};

}

#endif