#ifndef DM_ARGPARSER_ARG_PARSER_H
#define DM_ARGPARSER_ARG_PARSER_H

#include <vector>
#include <string>
#include <optional>
#include <algorithm>
#include <ranges>
#include <format>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include <dimuon/common/flavor.hpp>

struct arg_parser {
	arg_parser(int argc, char **argv) : arguments(argv + 1, argv + argc) {}

	std::optional<std::string> get_argument(const std::string name) const {
		read_flags.insert("--" + name);

		auto iterator = std::find(arguments.cbegin(), arguments.cend(), "--" + name);
		if (iterator == arguments.cend()) { 
			return std::nullopt; 
		} else {
			iterator++;
			if (iterator == arguments.cend() || iterator->starts_with("--")) {
				throw std::runtime_error(std::format("Argument --{} has no value", name));
			}
			return *iterator;
		}
	}

	std::string get_argument(const std::string name, const std::string default_value) const {
		return get_argument(name).value_or(default_value);
	}

	std::string get_required_argument(const std::string name) const {
		const std::optional<std::string> argument = get_argument(name);
		if (!argument) {
			throw std::runtime_error(std::format("Missing required argument --{}", name));
		}
		return *argument;
	}

	bool get_flag(const std::string name) const {
		read_flags.insert("--" + name);

		auto iterator = std::find(arguments.cbegin(), arguments.cend(), "--" + name);
		return iterator != arguments.cend();
	}

	std::optional<int> get_int_argument(const std::string name) const {
		return get_argument(name).transform([](const std::string &s) { return std::stoi(s); });
	}

	int get_int_argument(const std::string name, const int default_value) const {
		return get_int_argument(name).value_or(default_value);
	}
	
	std::optional<unsigned int> get_unsigned_int_argument(const std::string name) const {
		return get_argument(name).transform([](const std::string &s) { return static_cast<unsigned int>(std::stoul(s)); });
	}

	unsigned int get_unsigned_int_argument(const std::string name, const unsigned int default_value) const {
		return get_unsigned_int_argument(name).value_or(default_value);
	}

	long long int get_long_long_int_argument(const std::string name, const long long int default_value) const {
		return get_argument(name).transform([](const std::string &s) { return std::stoll(s); }).value_or(default_value);
	}

	std::optional<double> get_double_argument(const std::string name) const {
		return get_argument(name).transform([](const std::string &s) { return std::stod(s); });
	}

	double get_double_argument(const std::string name, const double default_value) const {
		return get_double_argument(name).value_or(default_value);
	}

	double get_required_double_argument(const std::string name) const {
		return std::stod(get_required_argument(name));
	}

	template <typename value_t>
	value_t get_choice_argument(const std::string name, const std::unordered_map<std::string, value_t> &choices, const value_t default_value) const {
		return get_argument(name).transform([&](const std::string &key) { return choice(name, key, choices); }).value_or(default_value);
	}

	template <typename value_t>
	value_t get_required_choice_argument(const std::string name, const std::unordered_map<std::string, value_t> &choices) const {
		return choice(name, get_required_argument(name), choices);
	}

	std::vector<std::string> get_string_list(const std::string name, const std::vector<std::string> default_value) const {
		return get_argument(name).transform(split_list).value_or(default_value);
	}

	std::vector<int> get_int_list(const std::string name, const std::vector<int> default_value) const {
		const std::vector<std::string> string_default = default_value | std::views::transform([](const int i) { return std::to_string(i); }) | std::ranges::to<std::vector>();
		const std::vector<std::string> string_list = get_string_list(name, string_default);

		return string_list | std::views::transform([](const std::string &s) { return std::stoi(s); }) | std::ranges::to<std::vector>();
	}

	std::vector<double> get_required_double_list(const std::string name) const {
		return split_list(get_required_argument(name)) | std::views::transform([](const std::string &s) { return std::stod(s); }) | std::ranges::to<std::vector>();
	}

	dm::flavor::flavor_vector get_flavor_list(const std::string name, const dm::flavor::flavor_vector &default_value) const {
		const std::vector<int> converted_default_value = default_value | std::views::transform(std::to_underlying<dm::flavor::type>) | std::ranges::to<std::vector>();
		const std::vector<int> flavor_int_list = get_int_list(name, converted_default_value);

		return flavor_int_list | std::views::transform([](const int f) { return static_cast<dm::flavor::type>(f); }) | std::ranges::to<std::vector>();
	}

	void check_arguments() const {
		const std::vector<std::string> flags = arguments
			| std::views::filter([](const std::string &argument) { return argument.starts_with("--"); })
			| std::ranges::to<std::vector>();

		for (const std::string &flag : flags) {
			if (std::ranges::count(flags, flag) > 1) {
				throw std::runtime_error(std::format("Argument {} is given more than once", flag));
			}
			if (!read_flags.contains(flag)) {
				throw std::runtime_error(std::format("Argument {} is unknown, or not used in this run", flag));
			}
		}
	}

	private:
	std::vector<std::string> arguments;
	mutable std::unordered_set<std::string> read_flags;

	static std::vector<std::string> split_list(const std::string &input) {
		std::vector<std::string> result;
		std::string current_string;
		for (const char character : input) {
			if (character == ',') {
				result.push_back(current_string);
				current_string.clear();
			} else {
				current_string.push_back(character);
			}
		}
		result.push_back(current_string);
		return result;
	}

	template <typename value_t>
	static value_t choice(const std::string &name, const std::string &key, const std::unordered_map<std::string, value_t> &choices) {
		const auto iterator = choices.find(key);
		if (iterator == choices.end()) {
			throw std::runtime_error(std::format("Invalid value '{}' of argument --{}", key, name));
		}
		return iterator->second;
	}
};

#endif