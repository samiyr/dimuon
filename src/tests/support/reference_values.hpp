#ifndef DM_TESTS_SUPPORT_REFERENCE_VALUES_H
#define DM_TESTS_SUPPORT_REFERENCE_VALUES_H

#include <filesystem>
#include <format>
#include <fstream>
#include <istream>
#include <iterator>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace dm::testing {
	struct reference_point {
		std::vector<double> inputs;
		double output;
	};

	inline std::istream &operator>>(std::istream &stream, reference_point &point) {
		std::string line;
		if (!std::getline(stream, line)) { return stream; }
		std::istringstream values(line);
		point.inputs.assign(std::istream_iterator<double>{values}, {});
		point.output = point.inputs.back();
		point.inputs.pop_back();
		return stream;
	}

	inline std::filesystem::path data_directory() { return DIMUON_TEST_DATA_DIR; }

	inline std::optional<std::vector<reference_point>> read_reference_points(const std::filesystem::path &file) {
		std::ifstream stream(data_directory() / file);
		if (!stream) { return std::nullopt; }
		while (stream.peek() == '#') { stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); }
		return std::vector<reference_point>(std::istream_iterator<reference_point>{stream}, {});
	}

	inline std::string missing_reference_message(const std::filesystem::path &file) {
		return std::format("No reference values: '{}' does not exist (see tests/data/README.md)", (data_directory() / file).string());
	}
}

#endif
