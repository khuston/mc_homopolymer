cmake_minimum_required(VERSION 3.17)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_BUILD_TYPE Debug)

project(Polymers)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/python)

#find_file (MCChainSourcePath mc_chain.cpp)
#find_file (LoggersSourcePath loggers.cpp)

#add_executable(mc_chain ${MCChainSourcePath} ${LoggersSourcePath})
