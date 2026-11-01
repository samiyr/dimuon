#ifndef DM_UTILITY_IO_H
#define DM_UTILITY_IO_H

#include <format>
#include <filesystem>
#include <concepts>
#include <ostream>
#include <streambuf>
#include <string>
#include <cerrno>
#include <stdlib.h>
#include <system_error>
#include <unistd.h>

namespace dm {

namespace io {
	/// @brief Create the directory structure necessary for \p path if it doesn't exist already
	/// @return Success or failure as a bool
	inline bool create_directory_tree(const std::filesystem::path &path) {
		const auto directory = path.parent_path();
		return std::filesystem::create_directories(directory);
	}

	/// Creates a new, empty file with a random name and the extension \p extension (e.g. ".csv") in the temporary directory
	/// ($TMPDIR, or /tmp if it is not set), and returns its path.
	/// The file is created atomically with `mkstemps`, so that concurrent processes never get the same file.
	inline std::filesystem::path create_temporary_file(const std::string &extension) {
		std::string path = (std::filesystem::temp_directory_path() / ("XXXXXX" + extension)).string();
		const int descriptor = mkstemps(path.data(), static_cast<int>(extension.size()));
		if (descriptor == -1) {
			throw std::filesystem::filesystem_error("Could not create a temporary file", path, std::error_code(errno, std::generic_category()));
		}
		close(descriptor);
		return path;
	}

	/// Creates a string from the integer \p i and pads it with zeroes to width \p width
	inline std::string leading_zeroes(const std::integral auto i, const unsigned int width) {
		return std::format("{:0{}}", i, width);
	}

	/// No-op output stream (does nothing with the data it receives)
	struct null_stream : public std::ostream {
		struct null_buffer : public std::streambuf {
			int overflow( int c ) noexcept { return c; }
		} buffer;
		null_stream() : std::ostream{&buffer} {}
	};

}

}

#endif