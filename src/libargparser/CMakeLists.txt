add_library(argparser INTERFACE)
target_compile_features(argparser INTERFACE cxx_std_23)

target_sources(argparser INTERFACE FILE_SET HEADERS)
add_subdirectory(argparser)

target_link_libraries(argparser INTERFACE dimuon::dimuon)
