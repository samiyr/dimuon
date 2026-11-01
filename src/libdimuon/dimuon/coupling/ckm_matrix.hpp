#ifndef DM_COUPLING_CKM_MATRIX_H
#define DM_COUPLING_CKM_MATRIX_H

#include <array>
#include <cstddef>
#include <cmath>
#include <utility>

#include "dimuon/common/flavor.hpp"
#include "dimuon/common/constants.hpp"

namespace dm {

/// Convenient access to the squared CKM matrix elements
struct ckm_matrix {
	/// Constructs the CKM matrix from the squared elements, which are to be given in the order
	/// V_ud, V_us, V_ub, V_cd, V_cs, V_cb, V_td, V_ts, V_tb
	constexpr explicit ckm_matrix(const std::array<double, 9> &ckm) noexcept : ckm_full_matrix{
		{
			/* d: 				(d) 					u 						(s) 					c 						(b) 					t 				*/
								0,					ckm[0],						0,					ckm[3],						0,					ckm[6]	,
			/* u: 				d 						(u)						s 						(c)						b 						(t) 			*/
							ckm[0],						0,					ckm[1],						0,					ckm[2],						0	,
			/* s:				(d) 					u 						(s) 					c 						(b) 					t 				*/
								0,					ckm[1],						0,					ckm[4],						0,					ckm[7]	,
			/* c: 				d 						(u)						s 						(c)						b 						(t) 			*/
							ckm[3],						0,					ckm[4],						0,					ckm[5],						0	,
			/* b: 				(d) 					u 						(s) 					c 						(b) 					t 				*/
								0,					ckm[2],						0,					ckm[5],						0,					ckm[8]	,
			/* t: 				d 						(u)						s 						(c)						b 						(t) 			*/
							ckm[6],						0,					ckm[7],						0,					ckm[8],						0		
		}
	} {}

	/// Default CKM matrix, with the values taken from the PDG global fit
	constexpr static ckm_matrix global_fit() noexcept {
		using namespace constants;
		return ckm_matrix{ {
			V_ud_sq, V_us_sq, V_ub_sq, 
			V_cd_sq, V_cs_sq, V_cb_sq, 
			V_td_sq, V_ts_sq, V_tb_sq
		} };
	}

	/// Returns the squared CKM matrix element corresponding to the two flavors. The order of the flavors doesn't matter, and they can also be antiflavors
	/// (i.e. with a negative sign). If an invalid flavor combination is given, such as dd or ds, returns zero. If flavors out of the range [-6, -1] + [1, 6]
	/// are provided, the result is undefined behaviour.
	double operator()(const flavor::type flavor1, const flavor::type flavor2) const noexcept {
		const std::size_t f1 = static_cast<std::size_t>(std::abs(std::to_underlying(flavor1)));
		const std::size_t f2 = static_cast<std::size_t>(std::abs(std::to_underlying(flavor2)));
		const std::size_t row = f1 - 1;
		const std::size_t column = f2 - 1;
		constexpr std::size_t column_count = 6; 
		return ckm_full_matrix[row * column_count + column];
	}

private:
	std::array<double, 36> ckm_full_matrix;
};

}

#endif