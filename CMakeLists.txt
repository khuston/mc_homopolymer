cmake_minimum_required(VERSION 3.10)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_BUILD_TYPE Debug)

project(mc_chain)

add_executable(mc_chain mc_chain.cpp loggers.cpp)

target_link_libraries(mc_chain gsl gslcblas)