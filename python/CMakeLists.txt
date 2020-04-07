find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

include_directories(${Python3_INCLUDE_DIRS})
include_directories(${Polymers_SOURCE_DIR}/include/mc_homopolymer)

#link_directories(${Polymers_SOURCE_DIR}/src)

add_subdirectory(pybind11)
pybind11_add_module(PyPolymers python_bindings.cpp)

target_link_libraries(PyPolymers PRIVATE gsl gslcblas)
target_link_libraries(PyPolymers PRIVATE Polymers)