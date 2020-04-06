find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

include_directories(${Python3_INCLUDE_DIRS})
include_directories(${Polymers_SOURCE_DIR}/include/mc_homopolymer)

link_directories(${Polymers_SOURCE_DIR}/src)

add_library(PyPolymers python_bindings.cpp)

target_link_libraries(PyPolymers gsl gslcblas)