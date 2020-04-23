cmake_minimum_required(VERSION 3.17)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_BUILD_TYPE Debug)
set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")

project(Polymers)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src)