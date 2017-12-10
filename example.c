#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static gsl_rng * r;

static gsl_vector * v; // working vector (3)
static gsl_matrix * B; // working matrix (N, 3)
static gsl_matrix * Q; // working matrix (3, 3)
static gsl_vector * axis; // working vector (3)
static gsl_vector * a; // working vector (3)
static gsl_vector * y; // working vector (N)
static gsl_vector * bl; // working matrix (N-1)
static double V; // working double


void print_matrix (gsl_matrix * m)
{
  int i, j;
  for (i = 0; i < m->size1; i++)
    for (j = 0; j < m->size2; j++)
      printf ("m(%d,%d) = %g\n", i, j,
              gsl_matrix_get (m, i, j));
}


void print_vector (gsl_vector * v)
{
  int i;
  for (i = 0; i < v->size; i++)
    printf ("v(%d) = %g\n", i,
            gsl_vector_get (v, i));
}


void set_random_vector (gsl_vector * u, double scale)
{
    int j;
    for (j = 0; j < 3; j++)
        gsl_vector_set(u, j, gsl_ran_gaussian(r, 1));
    gsl_vector_scale(u, scale/gsl_blas_dnrm2(u));
}


void set_random_unit_vector (gsl_vector * u) {
    set_random_vector(u, 1.0);
}


void bead_diffusion (gsl_matrix * m, int i)
{
    set_random_vector(v, 0.1);
    gsl_vector_view row = gsl_matrix_row(m, i);
    gsl_vector_add(&row.vector, v);
}


int end_rotate (gsl_matrix * m, int rotation_type)
{
    int i;
    int direction = gsl_rng_uniform_int(r, 2);
    // Rotation type 1 produces a random rotation
    if (rotation_type == 1)
    {
        set_random_unit_vector(v);

        if (direction == 0)
            i = 0;
        else if (direction == 1)
            i = m->size1-1;
        else return -2; // Unexpected value

        gsl_vector_view row = gsl_matrix_row(m, i);
        gsl_vector_add(&row.vector, v);

      return 0;
    }
  return -1; // Unrecognized rotation type
}


void slither (gsl_matrix * m)
{
    int i, j, N;
    int direction = gsl_rng_uniform_int(r, 2);
    N = m->size1;
    if (direction == 0)
    {
        /* Shift all bead positions up 1, effectively destroying the (n-1)th bead.
         * A new end bead will be created at the 0th position. */
        for (i = N-1; i > 0; i--)
            for (j = 0; j < 3; j++)
                gsl_matrix_set(m, i, j, gsl_matrix_get(m, i-1, j));

        /* Set new end bead to same position as adjacent bead.
         * Then displace it by a unit vector. */
        for (j = 0; j < 3; j++)
            gsl_matrix_set(m, 0, j, gsl_matrix_get(m, 1, j));

        set_random_unit_vector(v);
        gsl_vector_view row = gsl_matrix_row(m, 0);
        gsl_vector_add(&row.vector, v);
    }
    else if (direction == 1)
    {
        /* Shift all bead positions down 1, effectively destroying the 0th bead.
         * A new end bead will be created at the (n-1)th position. */
        for (i = 0; i < N-1; i++)
            for (j = 0; j < 3; j++)
                gsl_matrix_set(m, i, j, gsl_matrix_get(m, i+1, j));

        /* Set new end bead to same position as adjacent bead.
         * Then displace it by a unit vector. */
        for (j = 0; j < 3; j++)
            gsl_matrix_set(m, N-1, j, gsl_matrix_get(m, N-2, j));

        set_random_unit_vector(v);
        gsl_vector_view row = gsl_matrix_row(m, N-1);
        gsl_vector_add(&row.vector, v);

    }

}


void rotation_from_axis_angle (gsl_matrix * m, gsl_vector * u, double theta)
{
    /* vector u must be normalized. There is no check here. */
    double cost = cos (theta);
    double sint = sin (theta);
    double ux = gsl_vector_get(u, 0);
    double uy = gsl_vector_get(u, 1);
    double uz = gsl_vector_get(u, 2);
    gsl_matrix_set (m, 0, 0, cost + ux*ux*(1.0 - cost));
    gsl_matrix_set (m, 0, 1, ux*uy*(1.0 - cost) - uz*sint);
    gsl_matrix_set (m, 0, 2, ux*uz*(1.0 - cost) + uy*sint);

    gsl_matrix_set (m, 1, 0, uy*ux*(1.0 - cost) + uz*sint);
    gsl_matrix_set (m, 1, 1, cost + uy*uy*(1.0 - cost));
    gsl_matrix_set (m, 1, 2, uy*uz*(1.0 - cost) - ux*sint);

    gsl_matrix_set (m, 2, 0, uz*ux*(1.0 - cost) - uy*sint);
    gsl_matrix_set (m, 2, 1, uz*uy*(1.0 - cost) - ux*sint);
    gsl_matrix_set (m, 2, 2, cost + uz*uz*(1.0 - cost));
}


void rotate_around_points (gsl_matrix * m, gsl_vector * a0, gsl_vector * ar, double theta)
{
    int i;
    gsl_vector_view row;
    gsl_vector_memcpy(axis, ar);
    gsl_vector_sub(axis, a0);
    gsl_vector_scale(axis, 1.0/gsl_blas_dnrm2(axis));
    rotation_from_axis_angle (Q, axis, theta);
        
    gsl_matrix_view B_sub = gsl_matrix_submatrix(B, 0, 0, m->size1, 3);
    gsl_matrix_memcpy(&B_sub.matrix, m);
    for (i = 0; i < m->size1; i++)
        row = gsl_matrix_row(&B_sub.matrix, i);
        gsl_vector_sub(&row.vector, a0);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &B_sub.matrix, Q, 0.0, &B_sub.matrix);
    for (i = 0; i < m->size1; i++)
        row = gsl_matrix_row(&B_sub.matrix, i);
        gsl_vector_add(&row.vector, a0);
    gsl_matrix_swap(m, &B_sub.matrix);
}


void crankshaft (gsl_matrix * m)
{
    int N = m->size1;
    int k = gsl_rng_uniform_int(r, N-1);
    double theta = 2*M_PI*gsl_rng_uniform(r);
    
    gsl_vector_view a0;
    gsl_vector_view ar;

    if (k == 0)
        a0 = gsl_matrix_row(m, k+2);
    else
        a0 = gsl_matrix_row(m, k-1);

    if (k == 0)
        ar = gsl_matrix_row(m, k+3);
    else if (k == N-2)
        ar = gsl_matrix_row(m, k-2);
    else
        ar = gsl_matrix_row(m, k+2);
    
    gsl_matrix_view m_sub = gsl_matrix_submatrix(m, k, 0, 2, 3);
    rotate_around_points(&m_sub.matrix, &a0.vector, &ar.vector, theta);

    //gsl_vector_memcpy(axis, &ar.vector);
    //gsl_vector_sub(axis, &a0.vector);
    //gsl_vector_scale(axis, 1.0/gsl_blas_dnrm2(axis));
    //rotation_from_axis_angle (Q, axis, theta);
    //    
    //gsl_vector_view a1 = gsl_matrix_row(m, k);
    //gsl_vector_memcpy(v, &a1.vector);
    //gsl_vector_sub(v, &a0.vector);
    //gsl_blas_dgemv(CblasNoTrans, 1.0, Q, v, 0.0, v);
    //gsl_vector_add(v, &a0.vector);
    //gsl_vector_swap(&a1.vector, v);

    //gsl_vector_view a2 = gsl_matrix_row(m, k+1);
    //gsl_vector_memcpy(v, &a2.vector);
    //gsl_vector_sub(v, &a0.vector);
    //gsl_blas_dgemv(CblasNoTrans, 1.0, Q, v, 0.0, v);
    //gsl_vector_add(v, &a0.vector);
    //gsl_vector_swap(&a2.vector, v);
}


void bond_lengths (gsl_matrix * m, gsl_vector * b)
{
    int i, j;

    for (i = 1; i < b->size; i++)
    {
        for (j = 0; j < 3; j++)
            gsl_matrix_set(B, i, j, gsl_matrix_get(m, i, j) - gsl_matrix_get(m, i-1, j));
        gsl_vector_view row = gsl_matrix_row (B, i);
        gsl_vector_set(b, i, gsl_blas_dnrm2(&row.vector));
    }
}


int pivot (gsl_matrix * m)
{
    int i, j;
    int N = m->size1;
    int k = gsl_rng_uniform_int(r, N-2);
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
    //if (direction == 1)
    //    k = m->size1 - k;

    //printf("symmetry_op == %d\n", symmetry_operation);
    //printf("direction == %d\n", direction);
    //printf("k == %d\n", k);

    if (symmetry_operation == 0) // reflection
    {
        //a = chosen one
        //ap = behind chosen one
        //b = anything to be rotated

        // B' = -B + 2B[0] + 2(a - a')/|a - a'| (x) (B - a) . (a - a')/|a - a'|
        
        gsl_matrix_memcpy (B, m);

        // a = B[k]
        row = gsl_matrix_row(m, k);
        gsl_vector_memcpy (a, &row.vector);

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
        row = gsl_matrix_row(m, k+offset);
        gsl_vector_sub(a, &row.vector);
        
        // |a - a'|^2
        gsl_blas_ddot(a, a, &amap2);

        gsl_blas_dgemv(CblasNoTrans, 1.0/amap2, B, a, 0.0, y);

        if (direction == 0)
        {
            for (i = N-1; i >= k; i--)
                for (j = 0; j < 3; j++)
                    gsl_matrix_set (m, i, j, -gsl_matrix_get(m, i, j) + 2*gsl_matrix_get(m, k, j) + 2*gsl_vector_get(a, j)*gsl_vector_get(y, i));
        }
        else if (direction == 1)
        {
            for (i = 0; i < k; i++)
                for (j = 0; j < 3; j++)
                    gsl_matrix_set (m, i, j, -gsl_matrix_get(m, i, j) + 2*gsl_matrix_get(m, k, j) + 2*gsl_vector_get(a, j)*gsl_vector_get(y, i));
        }
    }
    else if (symmetry_operation == 1) // rotation
    {
        if (direction == 0)
        {
            offset = -1;
            m_sub = gsl_matrix_submatrix(m, k, 0, N-k, 3);
        }
        else if (direction == 1)
        {
            offset = 1;
            m_sub = gsl_matrix_submatrix(m, 0, 0, k, 3);
        }
        // a -> (a - a')  | a' = B[k-1]
        a1 = gsl_matrix_row(m, k);
        a2 = gsl_matrix_row(m, k+offset);
        theta = 2.0*M_PI*gsl_rng_uniform(r);
        rotate_around_points(&m_sub.matrix, &a1.vector, &a2.vector, theta);
    }
    else return -1;

  return 0;
}


void init_matrix (gsl_matrix * m)
{
  int i;

  for (i = 0; i < m->size1; i++)
  {
    gsl_matrix_set (m, i, 0, 1.0*i);
    gsl_matrix_set (m, i, 1, 1.0*sin(i/2.0));
    gsl_matrix_set (m, i, 2, 0.0);
  }
}


double evaluate_V (gsl_matrix * m)
{
    // Calculate bond lengths
    // Calculate bond energy
    // Calculate bead-wall interaction
    // Calculate 
    return 0.0;
}


void init_pdb (FILE * f)
{
    fprintf(f, "AUTHOR    GENERATED BY Kyle Huston\n");
}


void mc_run (gsl_matrix * m, int nsteps, int write_stride)
{
    int i;
    gsl_matrix * m_new;
    double V_new;
    const int number_of_move_types = 4;
    int move_type;
    int N = m->size1;
    FILE * f = fopen("output.pdb", "w");

    init_pdb(f);

    m_new = gsl_matrix_alloc (N, m->size2);

    V = evaluate_V (m);

    for (i = 0; i < nsteps; i++)
    {
        move_type = gsl_rng_uniform_int(r, number_of_move_types);

        gsl_matrix_memcpy(m_new, m);
        
        if (i % (N+1) == 0)
        {
            if (move_type == 0)
                pivot (m_new);
            else if (move_type == 1)
                crankshaft (m_new);
            else if (move_type == 2)
                slither (m_new);
            else if (move_type == 3)
                end_rotate (m_new, 1);
        }
        else
        {
            bead_diffusion (m_new, (i % (N+1))-1);
        }

        V_new = evaluate_V(m_new);

        if (gsl_rng_uniform(r) < pow(M_E, V - V_new))
        {
            gsl_matrix_memcpy(m, m_new);
            V = V_new;
        }

        if (i % write_stride == 0)
        {
            write_frame(m, "frame_%5d.pdb");
        }
        
    }

    fclose(f);
    gsl_matrix_free (m_new);
}



int main (void)
{
  gsl_matrix * m = gsl_matrix_alloc (4, 3);

  v = gsl_vector_alloc(3); // allocate global working vector
  B = gsl_matrix_alloc (m->size1, 3); // allocate global working matrix
  a = gsl_vector_alloc (3); // allocate global working vector
  Q = gsl_matrix_alloc (3, 3); // allocate global working matrix
  axis = gsl_vector_alloc (3); // allocate global working vector
  y = gsl_vector_alloc (m->size1); // allocate global working vector
  bl = gsl_vector_alloc (m->size1-1); // allocate global working vector
  const gsl_rng_type * T;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  init_matrix(m);
  mc_run(m, 1000000, 1000);

  gsl_matrix_free (m);
  gsl_vector_free(v);
  gsl_rng_free (r);
  gsl_matrix_free (B);
  gsl_vector_free(y);
  gsl_vector_free(bl);

  return 0;
}
