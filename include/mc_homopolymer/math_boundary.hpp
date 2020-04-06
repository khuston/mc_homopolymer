#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

class Matrix {
    Matrix() = delete;
    Matrix(const& Matrix other);
    Matrix(const& rows, const& columns);
    
}