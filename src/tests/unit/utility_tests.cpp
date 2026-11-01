#include <array>
#include <cmath>
#include <filesystem>
#include <numbers>
#include <ranges>
#include <set>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <dimuon/common/constants.hpp>
#include <dimuon/sidis/vectorization.hpp>
#include <dimuon/utility/io.hpp>
#include <dimuon/utility/math.hpp>
#include <dimuon/utility/sparse_vector.hpp>
#include <dimuon/utility/tensor.hpp>
#include <dimuon/utility/vectorization.hpp>

#include "support/assertions.hpp"

using namespace dm;

namespace {
	using std::numbers::pi, std::numbers::ln2;

	std::array<double, 8> make_dense(const sparse_vector<double> &v) {
		std::array<double, 8> values{};
		for (const auto &[index, value] : v) { values[index] = value; }
		return values;
	}

	sparse_vector<double> make_sparse(const std::vector<std::pair<std::size_t, double>> &elements) {
		sparse_vector<double> v;
		for (const auto &[index, value] : elements) { v[index] = value; }
		return v;
	}

	namespace sv = sidis::vectorization;
	const std::array pdf_values{sv::distribution_value::pdf, sv::distribution_value::pdf_hat};
	const std::array ff_values{sv::distribution_value::ff, sv::distribution_value::ff_hat, sv::distribution_value::ff_hat_p};
}

TEST(math, power) {
	const double x = 1.37;
	const auto check = [x]<int exponent>() { EXPECT_REL_NEAR(math::power<exponent>(x), std::pow(x, exponent), 4e-15) << "exponent " << exponent; };
	[&]<int... n>(std::integer_sequence<int, n...>) { (check.template operator()<n - 6>(), ...); }(std::make_integer_sequence<int, 17>{});
}

TEST(math, polylog_values) {
	const double zeta3 = constants::zeta3;

	EXPECT_EQ(math::polylogT<2>(0.0), 0.0);
	EXPECT_REL_NEAR(math::polylogT<2>(-1.0), -math::power<2>(pi) / 12.0, 2e-15);
	EXPECT_REL_NEAR(math::polylogT<2>(0.5), math::power<2>(pi) / 12.0 - math::power<2>(ln2) / 2.0, 2e-15);

	EXPECT_EQ(math::polylogT<3>(0.0), 0.0);
	EXPECT_REL_NEAR(math::polylogT<3>(-1.0), -3.0 * zeta3 / 4.0, 2e-15);
	EXPECT_REL_NEAR(math::polylogT<3>(0.5), (-2.0 * math::power<2>(pi) * ln2 + 4.0 * math::power<3>(ln2) + 21.0 * zeta3) / 24.0, 2e-15);
}

TEST(math, polylog_reflection) {
	// Li2(x) + Li2(1 - x) = pi^2 / 6 - log(x) log(1 - x)
	for (const double x : {0.1, 0.3, 0.75}) {
		EXPECT_REL_NEAR(math::polylogT<2>(x) + math::polylogT<2>(1.0 - x), math::power<2>(pi) / 6.0 - std::log(x) * std::log(1.0 - x), 2e-15) << "x = " << x;
	}
}

TEST(math, polylog_truncation) {
	// polylogT vanishes from one onwards
	EXPECT_EQ(math::polylogT<2>(1.0), 0.0);
	EXPECT_EQ(math::polylogT<3>(2.0), 0.0);
}

TEST(math, space_generation) {
	// linear_space and log_space take a step (here one that divides the range, so that max is included) and chebyshev_space a number of points
	const double min = 1.0;
	const double max = 1000.0;
	for (const std::vector<double> &space : {math::linear_space(min, max, 1.0), math::log_space(min, max, 1.0), math::chebyshev_space(min, max, 10)}) {
		EXPECT_FALSE(space.empty());
		EXPECT_EQ(space.front(), min);
		EXPECT_EQ(space.back(), max);
		EXPECT_TRUE(std::ranges::is_sorted(space));
	}
}

TEST(sparse_vector, arithmetic) {
	const double a0 = 2.0;
	const double a3 = 1.5;
	const double a7 = -1.0;

	const double b3 = 0.5;
	const double b5 = 4.0;

	const sparse_vector<double> a = make_sparse({{7, a7}, {0, a0}, {3, a3}});
	const sparse_vector<double> b = make_sparse({{3, b3}, {5, b5}});
	
	EXPECT_EQ(make_dense(a + b), (std::array{a0, 0.0, 0.0, a3 + b3, 0.0, b5, 0.0, a7}));
	EXPECT_EQ(make_dense(a - b), (std::array{a0, 0.0, 0.0, a3 - b3, 0.0, -b5, 0.0, a7}));
	EXPECT_EQ(make_dense(a * b), (std::array{0.0, 0.0, 0.0, a3 * b3, 0.0, 0.0, 0.0, 0.0}));
	EXPECT_EQ(make_dense(2.5 * a), (std::array{2.5 * a0, 0.0, 0.0, 2.5 * a3, 0.0, 0.0, 0.0, 2.5 * a7}));
	EXPECT_EQ(make_dense(a * 2.5), (std::array{2.5 * a0, 0.0, 0.0, 2.5 * a3, 0.0, 0.0, 0.0, 2.5 * a7}));
}

TEST(sparse_vector, ordering) {
	const std::size_t first = 0;
	const std::size_t second = 3;
	const std::size_t third = 7;

	const sparse_vector<double> v = make_sparse({{third, -1.0}, {first, 2.0}, {second, 3.0}});
	
	EXPECT_EQ(v | std::views::keys | std::ranges::to<std::vector>(), (std::vector{first, second, third}));
}

TEST(sparse_vector, single_element) {
	const std::size_t index = 4;
	const double value = 2.5;
	const auto element = single_element_sparse_vector(index, value);

	ASSERT_EQ(std::ranges::distance(element), 1);
	EXPECT_EQ(*std::ranges::begin(element), (std::pair{index, value}));

	EXPECT_TRUE(std::ranges::empty(single_element_sparse_vector(index, 0.0)));
	EXPECT_TRUE(std::ranges::empty(single_element_sparse_vector(1, std::array<double, 2>{})));
}

TEST(tensor, extents) {
	const std::size_t n0 = 2;
	const std::size_t n1 = 3;
	const std::size_t n2 = 4;

	const tensor<double, 3> t{tensor<double, 3>::extents_t{n0, n1, n2}, 0.0};

	EXPECT_EQ(t.extent(0), n0);
	EXPECT_EQ(t.extent(1), n1);
	EXPECT_EQ(t.extent(2), n2);
}

TEST(tensor, copy) {
	const double original_value = 1.0;
	const double copy_value = 5.0;

	const tensor<double, 2> original{tensor<double, 2>::extents_t{2, 2}, original_value};
	tensor<double, 2> copy = original;
	copy[1, 0] = copy_value;

	EXPECT_EQ((original[1, 0]), original_value);
	EXPECT_EQ((original.view()[1, 0]), original_value);
	EXPECT_EQ((copy[1, 0]), copy_value);
	EXPECT_EQ((copy.view()[1, 0]), copy_value);
}

TEST(vectorization, masking_vectors) {
	// The product of the masking vectors of a PDF value and an FF value is one in a single lane, a different lane for each product
	std::set<std::size_t> lanes;
	for (const sv::distribution_value a : pdf_values) {
		for (const sv::distribution_value b : ff_values) {
			const sv::distribution_product product = sv::make_masking_vector(a) * sv::make_masking_vector(b);
			std::vector<std::size_t> ones;
			for (std::size_t lane = 0; lane < sv::lane_count; ++lane) {
				EXPECT_TRUE(product[lane] == 0.0 || product[lane] == 1.0);
				if (product[lane] == 1.0) { ones.push_back(lane); }
			}
			ASSERT_EQ(ones.size(), 1);
			lanes.insert(ones.front());
		}
	}
	EXPECT_EQ(lanes.size(), sv::lane_count);
}

namespace {
	/// Creates temporary files with `io::create_temporary_file` and removes them in TearDown, which runs also when the test fails
	class io_temporary_file : public ::testing::Test {
	protected:
		std::filesystem::path create_temporary_file(const std::string &extension) {
			paths.push_back(io::create_temporary_file(extension));
			return paths.back();
		}

		void TearDown() override {
			for (const std::filesystem::path &path : paths) {
				std::error_code error;
				std::filesystem::remove(path, error);
			}
		}

	private:
		std::vector<std::filesystem::path> paths;
	};
}

TEST_F(io_temporary_file, temp_file_creation) {
	const std::string extension = ".csv";
	const std::filesystem::path first = create_temporary_file(extension);
	const std::filesystem::path second = create_temporary_file(extension);

	EXPECT_NE(first, second);
	
	for (const std::filesystem::path &path : {first, second}) {
		EXPECT_TRUE(std::filesystem::is_regular_file(path));
		EXPECT_EQ(std::filesystem::file_size(path), 0);
		EXPECT_EQ(path.extension(), extension);
		EXPECT_TRUE(std::filesystem::equivalent(path.parent_path(), std::filesystem::temp_directory_path()));
	}
}
