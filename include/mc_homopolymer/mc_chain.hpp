#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <loggers.hpp>
#include <memory>
#include <vector>

#ifndef MCCHAIN_H
#define MCCHAIN_H

namespace Polymers
{
class MonteCarloChain
{
  public:  
    MonteCarloChain() = delete;
    MonteCarloChain(const int& count, const double& epsilon, const double& sigma);
    ~MonteCarloChain();    
    void set_logger(const Loggers::Logger& logger); // todo: Separate logging into decorator.
    void log(const Loggers::LogLevel& message_level, const std::string& message);
    void run(const int& nsteps, const int& write_stride, const bool& tethered);

  private:
    static constexpr double twofifths = 2.0 / 5.0;
    static const int dimensions = 3;

    gsl_rng *r;
    gsl_vector *v;    // working vector (3)
    gsl_matrix *B;    // working matrix (N, 3)
    gsl_matrix *Q;    // working matrix (3, 3)
    gsl_vector *axis; // working vector (3)
    gsl_vector *a;    // working vector (3)
    gsl_vector *y;    // working vector (N)
    gsl_vector *bl;   // working matrix (N-1)
    gsl_matrix *positions;
    double V; // working double
    double s;   // polymer-wall interaction length
    double eps; // polymer-wall interaction energy
    double V_prefactor;
    double const1;
    double const2;
    int crankshaft_log2max;

    Loggers::Logger& logger;

    void print_matrix(gsl_matrix *m);
    void print_vector(gsl_vector *v);
    void set_random_vector(gsl_vector *u, double scale);
    void set_random_unit_vector(gsl_vector *u);
    void bead_diffusion(gsl_matrix *m, int i);
    int end_rotate(gsl_matrix *m, int rotation_type);
    void slither(gsl_matrix *m);
    void rotation_from_axis_angle(gsl_matrix *m, gsl_vector *u, double theta);
    void rotate_around_points(gsl_matrix *m, gsl_vector *a0, gsl_vector *ar, double theta);
    int crankshaft(gsl_matrix *m);
    void bond_lengths(gsl_vector *b, gsl_matrix *m);
    void pivot_z(gsl_matrix *m);
    int pivot(gsl_matrix *m);
    void init_matrix(gsl_matrix *m);
    int any_adsorbed(gsl_matrix *m);
    double evaluate_V_total(gsl_matrix *m);
    void fix_bead0(gsl_matrix *m);
};
} // namespace MCChain

#endif