#include <mc_chain.hpp>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <loggers.hpp>
#if __cplusplus > 201703L
// C++20 code
#include <format>
#endif

namespace Polymers
{
MonteCarloChain::MonteCarloChain(const int& count, const double& epsilon, const double& sigma)
: logger(Loggers::LoggerFactory::get_null_logger())
{
    // eps = 0.1;
    // s = 1.0;
    V_prefactor = 2.0 * M_PI * epsilon;
    const1 = sqrt(2.0) * pow(sigma, 3.0) / 3.0;
    const2 = sigma * 0.61 / sqrt(2.0);
    crankshaft_log2max = (int)(std::log((double)(count - 2)) / std::log(2.0));

    positions = gsl_matrix_alloc(count, dimensions);
    v = gsl_vector_alloc(dimensions);
    B = gsl_matrix_alloc(count, dimensions);
    a = gsl_vector_alloc(dimensions);
    Q = gsl_matrix_alloc(dimensions, dimensions);
    axis = gsl_vector_alloc(dimensions);
    y = gsl_vector_alloc(count);
    bl = gsl_vector_alloc(count - 1);

    const gsl_rng_type *T;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    init_matrix(positions);

    int frame_number = 0;
}

MonteCarloChain::~MonteCarloChain()
{
    gsl_matrix_free(positions);
    gsl_vector_free(v);
    gsl_rng_free(r);
    gsl_matrix_free(B);
    gsl_vector_free(y);
    gsl_vector_free(bl);
}

void MonteCarloChain::set_logger(const Loggers::Logger& logger)
{
    this->logger = logger;
}

void MonteCarloChain::log(const Loggers::LogLevel& message_level, const std::string& message)
{
    this->logger.log(message_level, message);
}


void MonteCarloChain::print_matrix(gsl_matrix *m)
{
    int i, j;
    for (i = 0; i < m->size1; i++)
        for (j = 0; j < m->size2; j++)
            printf("m(%d,%d) = %g\n", i, j, gsl_matrix_get(m, i, j));
}

void MonteCarloChain::print_vector(gsl_vector *v)
{
    int i;
    for (i = 0; i < v->size; i++)
        printf("v(%d) = %g\n", i, gsl_vector_get(v, i));
}

void MonteCarloChain::set_random_vector(gsl_vector *u, double scale)
{
    int j;
    for (j = 0; j < dimensions; j++)
        gsl_vector_set(u, j, gsl_ran_gaussian(r, 1));
    gsl_vector_scale(u, scale / gsl_blas_dnrm2(u));
}

void MonteCarloChain::set_random_unit_vector(gsl_vector *u)
{
    set_random_vector(u, 1.0);
}

void MonteCarloChain::bead_diffusion(gsl_matrix *m, int i)
{
    set_random_vector(v, 0.1);
    gsl_vector_view row = gsl_matrix_row(m, i);
    gsl_vector_add(&row.vector, v);
}

int MonteCarloChain::end_rotate(gsl_matrix *m, int rotation_type)
{
    int i;
    int direction = gsl_rng_uniform_int(r, 2);
    // Rotation type 1 produces a random orientation
    if (rotation_type == 1)
    {
        set_random_unit_vector(v);

        if (direction == 0)
            i = 0;
        else if (direction == 1)
            i = m->size1 - 1;
        else
            return -2; // Unexpected value

        gsl_vector_view row = gsl_matrix_row(m, i);
        gsl_vector_add(&row.vector, v);

        return i;
    }
    return -1; // Unrecognized rotation type
}

void MonteCarloChain::slither(gsl_matrix *m)
{
    int i, j, k, N, offset;
    int direction = gsl_rng_uniform_int(r, 2);
    N = m->size1;
    if (direction == 0)
    {
        k = 0;
        offset = -1;
        /* Shift all bead positions up 1, effectively destroying the (n-1)th bead.
         * A new end bead will be created at the 0th position. */
        for (i = N - 1; i > 0; i--)
            for (j = 0; j < 3; j++)
                gsl_matrix_set(m, i, j, gsl_matrix_get(m, i + offset, j));
    }
    else if (direction == 1)
    {
        k = N - 1;
        offset = 1;
        /* Shift all bead positions down 1, effectively destroying the 0th bead.
         * A new end bead will be created at the (n-1)th position. */
        for (i = 0; i < N - 1; i++)
            for (j = 0; j < 3; j++)
                gsl_matrix_set(m, i, j, gsl_matrix_get(m, i + offset, j));
    }

    /* Set new end bead to same position as adjacent bead.
     * Then displace it by a unit vector. */
    for (j = 0; j < 3; j++)
        gsl_matrix_set(m, k, j, gsl_matrix_get(m, k - offset, j));

    set_random_unit_vector(v);
    gsl_vector_view row = gsl_matrix_row(m, k);
    gsl_vector_add(&row.vector, v);
}

void MonteCarloChain::rotation_from_axis_angle(gsl_matrix *m, gsl_vector *u, double theta)
{
    /* vector u must be normalized. There is no check here. */
    double cost = cos(theta);
    double sint = sin(theta);
    double ux = gsl_vector_get(u, 0);
    double uy = gsl_vector_get(u, 1);
    double uz = gsl_vector_get(u, 2);
    gsl_matrix_set(m, 0, 0, cost + ux * ux * (1.0 - cost));
    gsl_matrix_set(m, 0, 1, ux * uy * (1.0 - cost) - uz * sint);
    gsl_matrix_set(m, 0, 2, ux * uz * (1.0 - cost) + uy * sint);

    gsl_matrix_set(m, 1, 0, uy * ux * (1.0 - cost) + uz * sint);
    gsl_matrix_set(m, 1, 1, cost + uy * uy * (1.0 - cost));
    gsl_matrix_set(m, 1, 2, uy * uz * (1.0 - cost) - ux * sint);

    gsl_matrix_set(m, 2, 0, uz * ux * (1.0 - cost) - uy * sint);
    gsl_matrix_set(m, 2, 1, uz * uy * (1.0 - cost) - ux * sint);
    gsl_matrix_set(m, 2, 2, cost + uz * uz * (1.0 - cost));
}

void MonteCarloChain::rotate_around_points(gsl_matrix *m, gsl_vector *a0, gsl_vector *ar, double theta)
{
    int i;
    gsl_vector_view row;
    gsl_vector_memcpy(axis, ar);
    gsl_vector_sub(axis, a0);
    gsl_vector_scale(axis, 1.0 / gsl_blas_dnrm2(axis));
    rotation_from_axis_angle(Q, axis, theta);

    gsl_matrix_view B_sub = gsl_matrix_submatrix(B, 0, 0, m->size1, 3);
    gsl_matrix_memcpy(&B_sub.matrix, m);
    for (i = 0; i < m->size1; i++)
        row = gsl_matrix_row(&B_sub.matrix, i);
    gsl_vector_sub(&row.vector, a0);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &B_sub.matrix, Q, 0.0, m);
    for (i = 0; i < m->size1; i++)
        row = gsl_matrix_row(m, i);
    gsl_vector_add(&row.vector, a0);
}

int MonteCarloChain::crankshaft(gsl_matrix *m)
{
    int N = m->size1;
    int n = pow(2.0, (double)(gsl_rng_uniform_int(r, crankshaft_log2max) + 1.0));
    int k = gsl_rng_uniform_int(r, N - n);
    double theta = 2 * M_PI * gsl_rng_uniform(r);

    gsl_vector_view a0;
    gsl_vector_view ar;

    if (k == 0)
        a0 = gsl_matrix_row(m, k + n + 1);
    else
        a0 = gsl_matrix_row(m, k - 1);

    if (k + n == N)
        ar = gsl_matrix_row(m, k - 2);
    else
        ar = gsl_matrix_row(m, k + n);

    gsl_matrix_view m_sub = gsl_matrix_submatrix(m, k, 0, n, 3);
    rotate_around_points(&m_sub.matrix, &a0.vector, &ar.vector, theta);

    return k;
}

void MonteCarloChain::bond_lengths(gsl_vector *b, gsl_matrix *m)
{
    int i, j;

    for (i = 0; i < b->size; i++)
    {
        for (j = 0; j < 3; j++)
        {
            gsl_matrix_set(B, i, j, gsl_matrix_get(m, i + 1, j) - gsl_matrix_get(m, i, j));
        }
        gsl_vector_view row = gsl_matrix_row(B, i);
        gsl_vector_set(b, i, gsl_blas_dnrm2(&row.vector));
    }
}

void MonteCarloChain::pivot_z(gsl_matrix *m)
{
    int i, j;
    int N = m->size1;
    int k = gsl_rng_uniform_int(r, N - 2);
    int direction = gsl_rng_uniform_int(r, 2);
    const int num_symm_ops = 2;
    int symmetry_operation = gsl_rng_uniform_int(r, num_symm_ops);
    double theta;
    int offset = 0;
    gsl_vector_view a1;
    gsl_matrix_view m_sub;

    k += 1;

    if (symmetry_operation == 0) // reflection
    {
    }
    else if (symmetry_operation == 1) // rotation
    {
        if (direction == 0)
        {
            offset = -1;
            m_sub = gsl_matrix_submatrix(m, k, 0, N - k, 3);
        }
        else if (direction == 1)
        {
            offset = 1;
            m_sub = gsl_matrix_submatrix(m, 0, 0, k, 3);
        }

        theta = 2.0 * M_PI * gsl_rng_uniform(r);
        a1 = gsl_matrix_row(m, k);
        gsl_vector_memcpy(v, &a1.vector);
        gsl_vector_set(v, 2, gsl_vector_get(v, 2) + 1.0);
        rotate_around_points(&m_sub.matrix, &a1.vector, v, theta);
    }
}

int MonteCarloChain::pivot(gsl_matrix *m)
{
    int i, j;
    int N = m->size1;
    int k = gsl_rng_uniform_int(r, N - 2);
    int direction = gsl_rng_uniform_int(r, 2);
    const int num_symm_ops = 2;
    int symmetry_operation = gsl_rng_uniform_int(r, num_symm_ops);
    double amap2;
    double theta;
    int offset = 0;
    gsl_vector_view row;
    gsl_vector_view a1;
    gsl_vector_view a2;
    gsl_matrix_view m_sub;

    k += 1;

    if (symmetry_operation == 0) // reflection
    {
        // B' = -B + 2B[0] + 2(a - a')/|a - a'| (x) (B - a) . (a - a')/|a - a'|
        gsl_matrix_memcpy(B, m);

        // a = B[k]
        row = gsl_matrix_row(m, k);
        gsl_vector_memcpy(a, &row.vector);

        // B -> (B - a)
        if (direction == 0)
        {
            for (i = k; i < N; i++)
                for (j = 0; j < 3; j++)
                    gsl_matrix_set(B, i, j, gsl_matrix_get(B, i, j) - gsl_vector_get(a, j));

            offset = -1;
        }
        else if (direction == 1)
        {
            for (i = 0; i < k; i++)
                for (j = 0; j < 3; j++)
                    gsl_matrix_set(B, i, j, gsl_matrix_get(B, i, j) - gsl_vector_get(a, j));

            offset = 1;
        }
        // a -> (a - a')  | a' = B[k-1]
        row = gsl_matrix_row(m, k + offset);
        gsl_vector_sub(a, &row.vector);

        // |a - a'|^2
        gsl_blas_ddot(a, a, &amap2);

        gsl_blas_dgemv(CblasNoTrans, 1.0 / amap2, B, a, 0.0, y);

        if (direction == 0)
        {
            for (i = N - 1; i >= k; i--)
                for (j = 0; j < 3; j++)
                    gsl_matrix_set(m, i, j,
                                   -gsl_matrix_get(m, i, j) + 2 * gsl_matrix_get(m, k, j) +
                                       2 * gsl_vector_get(a, j) * gsl_vector_get(y, i));
        }
        else if (direction == 1)
        {
            for (i = 0; i < k; i++)
                for (j = 0; j < 3; j++)
                    gsl_matrix_set(m, i, j,
                                   -gsl_matrix_get(m, i, j) + 2 * gsl_matrix_get(m, k, j) +
                                       2 * gsl_vector_get(a, j) * gsl_vector_get(y, i));
        }
    }
    else if (symmetry_operation == 1) // rotation
    {
        if (direction == 0)
        {
            offset = -1;
            m_sub = gsl_matrix_submatrix(m, k, 0, N - k, 3);
        }
        else if (direction == 1)
        {
            offset = 1;
            m_sub = gsl_matrix_submatrix(m, 0, 0, k, 3);
        }
        // a -> (a - a')  | a' = B[k-1]
        a1 = gsl_matrix_row(m, k);
        a2 = gsl_matrix_row(m, k + offset);
        theta = 2.0 * M_PI * gsl_rng_uniform(r);
        rotate_around_points(&m_sub.matrix, &a1.vector, &a2.vector, theta);
    }
    else
        return -1;

    return 0;
}

void MonteCarloChain::init_matrix(gsl_matrix *m)
{
    int i;

    for (i = 0; i < m->size1; i++)
    {
        gsl_matrix_set(m, i, 0, 0.0);
        gsl_matrix_set(m, i, 1, 1.0 * sin(i / 2.0));
        gsl_matrix_set(m, i, 2, 1.0 + 1.0 * i);
    }
}

int MonteCarloChain::any_adsorbed(gsl_matrix *m)
{
    int i;
    double x;
    int N = m->size1;
    gsl_vector_view column = gsl_matrix_column(m, 2);
    for (i = 0; i < N; i++)
    {
        x = gsl_vector_get(&column.vector, i);
        if (x <= 2.0)
            return 1;
    }
    return 0;
}

double MonteCarloChain::evaluate_V_total(gsl_matrix *m)
{
    int i, N;
    double r2_bond;
    double x;
    // Calculate bond lengths
    // Calculate bond energy
    // Calculate bead-wall interaction
    // Calculate
    N = m->size1;
    gsl_vector_view bl_sub = gsl_vector_subvector(bl, 0, N - 1);
    bond_lengths(&bl_sub.vector, m);
    gsl_vector_add_constant(&bl_sub.vector, -1.0);
    gsl_blas_ddot(&bl_sub.vector, &bl_sub.vector, &r2_bond);

    gsl_vector_view column = gsl_matrix_column(m, 2);
    double V_wall = 0;
    for (i = 0; i < N; i++)
    {
        x = gsl_vector_get(&column.vector, i);
        if (x <= 0)
            V_wall += 1.0 / 0.0;
        else
            V_wall += V_prefactor * (twofifths * pow(s / x, 10) - pow(s / x, 4) - const1 / pow(x + const2, 3));
    }

    double V = V_wall + 100.0 * r2_bond;

    return V;
}

void MonteCarloChain::fix_bead0(gsl_matrix *m)
{
    gsl_matrix_set(m, 0, 0, 0.0);
    gsl_matrix_set(m, 0, 1, 0.0);
    gsl_matrix_set(m, 0, 2, 0.0);
}

void MonteCarloChain::run(const int& nsteps, const int& write_stride, const bool& tethered)
{
    int i, k;
    gsl_matrix *m_new;
    gsl_matrix_view m_partial;
    gsl_matrix_view m_partial_new;
    double V_new;
    const int number_of_move_types = 5;
    int move_type;
    int N = positions->size1;

    m_new = gsl_matrix_alloc(N, positions->size2);

    for (i = 0; i < nsteps; i++)
    {
        gsl_matrix_memcpy(m_new, positions);

        if (i % (N + 1) == 0)
        {
            move_type = gsl_rng_uniform_int(r, number_of_move_types);
            // printf ("debug: move_type=%d\n", move_type);
            if (move_type == 0)
            {
                pivot(m_new);
                if (tethered)
                    fix_bead0(m_new);
                V = evaluate_V_total(positions);
                V_new = evaluate_V_total(m_new);
            }
            if (move_type == 1)
            {
                pivot_z(m_new);
                if (tethered)
                    fix_bead0(m_new);
                V = evaluate_V_total(positions);
                V_new = evaluate_V_total(m_new);
            }
            else if (move_type == 2)
            {
                crankshaft(m_new);
                if (tethered)
                    fix_bead0(m_new);

                V = evaluate_V_total(positions);
                V_new = evaluate_V_total(m_new);
            }
            else if (move_type == 3)
            {
                slither(m_new);
                if (tethered)
                    fix_bead0(m_new);
                V = evaluate_V_total(positions);
                V_new = evaluate_V_total(m_new);
            }
            else if (move_type == 4)
            {
                k = end_rotate(m_new, 1);
                if (tethered)
                    fix_bead0(m_new);
                if (k == 0)
                {
                    m_partial = gsl_matrix_submatrix(positions, k, 0, 2, 3);
                    m_partial_new = gsl_matrix_submatrix(m_new, k, 0, 2, 3);
                }
                else
                {
                    m_partial = gsl_matrix_submatrix(positions, k - 1, 0, 2, 3);
                    m_partial_new = gsl_matrix_submatrix(m_new, k - 1, 0, 2, 3);
                }
                V = evaluate_V_total(&m_partial.matrix);
                V_new = evaluate_V_total(&m_partial_new.matrix);
            }
        }
        else
        {
            k = (i % (N + 1)) - 1;
            bead_diffusion(m_new, k);
            if (tethered)
                fix_bead0(m_new);
            if (k == 0)
            {
                m_partial = gsl_matrix_submatrix(positions, k, 0, 2, 3);
                m_partial_new = gsl_matrix_submatrix(m_new, k, 0, 2, 3);
            }
            else if (k == N - 1)
            {
                m_partial = gsl_matrix_submatrix(positions, k - 1, 0, 2, 3);
                m_partial_new = gsl_matrix_submatrix(m_new, k - 1, 0, 2, 3);
            }
            else
            {
                m_partial = gsl_matrix_submatrix(positions, k - 1, 0, 3, 3);
                m_partial_new = gsl_matrix_submatrix(m_new, k - 1, 0, 3, 3);
            }
            V = evaluate_V_total(&m_partial.matrix);
            V_new = evaluate_V_total(&m_partial_new.matrix);
        }
        if (any_adsorbed(m_new) & (gsl_rng_uniform(r) < pow(M_E, V - V_new)))
        {
            gsl_matrix_memcpy(positions, m_new);
        }
    }
    gsl_matrix_free(m_new);
}
} // namespace MCChain

int main(int argc, char **argv)
{
    int opt;
    int N;
    double epsilon;
    double sigma;

    N = 5;
    epsilon = 1;
    sigma = 1;

    Polymers::MonteCarloChain chain = Polymers::MonteCarloChain(N, epsilon, sigma);

    std::unique_ptr<Loggers::Logger> logger = Loggers::LoggerFactory::make_stdout_logger(Loggers::LogLevel::Debug);

    chain.set_logger(*logger);
    chain.run(10000 * N, 1000 * N, false);

    return 0;
}