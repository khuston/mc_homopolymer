include_directories(${Polymers_SOURCE_DIR}/include/mc_homopolymer)

add_library(Polymers mc_chain.cpp loggers.cpp)

target_link_libraries(Polymers gsl gslcblas)