include(CMakeFindDependencyMacro)

# LHAPDF and Cuba have no CMake packages of their own. Their find modules are installed next to this file.
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
find_dependency(GSL)
find_dependency(LHAPDF)
find_dependency(Cuba)
list(POP_FRONT CMAKE_MODULE_PATH)

include(${CMAKE_CURRENT_LIST_DIR}/dimuonTargets.cmake)
