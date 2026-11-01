#ifndef DM_COMMON_PERTURBATIVE_QUANTITY_H
#define DM_COMMON_PERTURBATIVE_QUANTITY_H

#include "dimuon/common/perturbative_order.hpp"
#include "dimuon/utility/collections.hpp"
#include "dimuon/utility/math.hpp"
#include <array>
#include <cstddef>
#include <format>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace dm {

/// A quantity with LO, NLO and NNLO contributions, with support for basic arithmetic operations and printing
class perturbative_quantity {
public:
	perturbative_quantity() : perturbative_quantity{0.0} {}

	explicit perturbative_quantity(
		double alpha_s_,
		double lo_correction = 0.0,
		double nlo_correction = 0.0,
		double nnlo_correction = 0.0
	) : alpha_s{alpha_s_}, corrections{lo_correction, nlo_correction, nnlo_correction} {}

	explicit perturbative_quantity(
		const perturbative_order pto,
		double value
	) : perturbative_quantity{} {
		raw_correction(pto) = value;
	}

	/// @brief Calculates the perturbative coefficient based on the given perturbative order
	/// @param[in] pto Perturbative order
	/// @return `(alpha_s / (2 pi))^n` where `n` is 0 for LO, 1 for NLO, and 2 for NNLO
	/// @note Passing a value for \p pto that has not been defined in the enum \p perturbative_order results in undefined behaviour
	double pto_coefficient(const perturbative_order pto) const noexcept {
		switch(pto) {
			case perturbative_order::lo: return 1.0;
			case perturbative_order::nlo: return alpha_s / (2.0 * std::numbers::pi);
			case perturbative_order::nnlo: return math::power<2>(alpha_s / (2.0 * std::numbers::pi));
		}
		std::unreachable();
	}

	/// The correction of the order \p pto, without multiplication by the prefactor given by `pto_coefficient`
	double correction(const perturbative_order pto) const noexcept {
		return corrections[static_cast<std::size_t>(std::to_underlying(pto))];
	}

	/// @brief Gives access to the underlying data for the given perturbative order
	/// @param[in] pto The perturbative order
	/// @return Reference to the stored value, which can be read or written into
	double &raw_correction(const perturbative_order pto) noexcept {
		return corrections[static_cast<std::size_t>(std::to_underlying(pto))];
	}

	/// The leading-order value
	double lo_correction() const noexcept { return corrections[0]; }
	/// The next-to-leading-order correction, without multiplication by the prefactor (~ `alpha_s`) given by `pto_coefficient`
	double nlo_correction() const noexcept { return corrections[1]; }
	/// The next-to-next-to-leading-order correction, without multiplication by the prefactor (~ `alpha_s^2`) given by `pto_coefficient`
	double nnlo_correction() const noexcept { return corrections[2]; }

	/// The next-to-leading-order correction, with multiplication by the prefactor (~ `alpha_s`) given by `pto_coefficient`
	double nlo_correction_with_alpha_s() const noexcept { return pto_coefficient(perturbative_order::nlo) * nlo_correction(); }
	/// The next-to-next-to-leading-order correction, with multiplication by the prefactor (~ `alpha_s^2`) given by `pto_coefficient`
	double nnlo_correction_with_alpha_s() const noexcept { return pto_coefficient(perturbative_order::nnlo) * nnlo_correction(); }

	/// The leading-order value
	double lo() const noexcept { return lo_correction(); }
	/// The next-to-leading-order value (i.e. LO + NLO correction)
	double nlo() const noexcept { return lo_correction() + nlo_correction_with_alpha_s(); }
	/// The next-to-next-to-leading-order value (i.e. LO + NLO correction + NNLO correction)
	double nnlo() const noexcept { return lo_correction() + nlo_correction_with_alpha_s() + nnlo_correction_with_alpha_s(); }

private:
	/// The alpha_s value
	double alpha_s;
	/// The corrections of each order, indexed by the order (see `correction`).
	/// The arithmetic operators act on them element by element (see the operators of `std::array` in `collections`).
	std::array<double, 3> corrections;

	/// @throws std::runtime_error if the `alpha_s` value of \p rhs differs from this one, in which case the two cannot be combined
	void check_same_alpha_s(const perturbative_quantity &rhs) const {
		if (alpha_s != rhs.alpha_s) {
			throw std::runtime_error(std::format("Cannot combine perturbative objects with non-equal alpha_s values {} and {}", alpha_s, rhs.alpha_s));
		}
	}

public:
	perturbative_quantity &operator+=(const perturbative_quantity &rhs) {
		check_same_alpha_s(rhs);
		using collections::operator+=;
		corrections += rhs.corrections;
		return *this;
	}
	perturbative_quantity &operator-=(const perturbative_quantity &rhs) {
		check_same_alpha_s(rhs);
		using collections::operator-=;
		corrections -= rhs.corrections;
		return *this;
	}
	perturbative_quantity &operator*=(const perturbative_quantity &rhs) {
		check_same_alpha_s(rhs);
		using collections::operator*=;
		corrections *= rhs.corrections;
		return *this;
	}
	perturbative_quantity &operator/=(const perturbative_quantity &rhs) {
		check_same_alpha_s(rhs);
		using collections::operator/=;
		corrections /= rhs.corrections;
		return *this;
	}

	perturbative_quantity &operator+=(const double rhs) noexcept {
		using collections::operator+=;
		corrections += rhs;
		return *this;
	}
	perturbative_quantity &operator-=(const double rhs) noexcept {
		using collections::operator-=;
		corrections -= rhs;
		return *this;
	}
	perturbative_quantity &operator*=(const double rhs) noexcept {
		using collections::operator*=;
		corrections *= rhs;
		return *this;
	}
	perturbative_quantity &operator/=(const double rhs) noexcept {
		using collections::operator/=;
		corrections /= rhs;
		return *this;
	}

	friend perturbative_quantity operator+(perturbative_quantity lhs, const perturbative_quantity &rhs) { lhs += rhs; return lhs; }
	friend perturbative_quantity operator-(perturbative_quantity lhs, const perturbative_quantity &rhs) { lhs -= rhs; return lhs; }
	friend perturbative_quantity operator*(perturbative_quantity lhs, const perturbative_quantity &rhs) { lhs *= rhs; return lhs; }
	friend perturbative_quantity operator/(perturbative_quantity lhs, const perturbative_quantity &rhs) { lhs /= rhs; return lhs; }

	friend perturbative_quantity operator*(perturbative_quantity lhs, const double rhs) noexcept { lhs *= rhs; return lhs; }
	friend perturbative_quantity operator*(const double lhs, perturbative_quantity rhs) noexcept { rhs *= lhs; return rhs; }

	friend std::ostream& operator<<(std::ostream &os, const perturbative_quantity &r) {
		return os << r.lo() << ", " << r.nlo() << ", " << r.nnlo();
	}
};

}

#endif