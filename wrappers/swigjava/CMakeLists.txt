cmake_minimum_required (VERSION 2.8)
set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/Modules)
include(${CMAKE_CURRENT_SOURCE_DIR}/Modules/CMakeUtilityFunctions.cmake)

project (adiosjava)

add_subdirectory(src)
add_subdirectory(test)

include(CTest)
