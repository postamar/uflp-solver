// (C) 2009-2011, by Paul-Virak Khuong, under the MIT License in ../LICENSE.txt.

#ifndef __BUNDLE_SOLVER_H__
#define __BUNDLE_SOLVER_H__

#include <stdint.h>

typedef enum {
        kEINorm = 0,
        kEIAbort,
        kEILoopNow,
        kEIContAnyway
} bundle_step_status_t;

typedef enum {
        kOK = 0,
        kUnbounded,
        kUnfeasible,
        kAbort,
        kMaxIter,
        kError
} bundle_status_t;

struct bundle_callbacks
{
        void (*SetGi)(struct bundle_callbacks *,
                      double * subgradient, unsigned name);
        double (*Fi)(struct bundle_callbacks *,
                     const double * Lam_val, const unsigned * Lam_idx,
                     unsigned Lam_count);
        int (*GetGi)(struct bundle_callbacks *,
                     double * Eps, const unsigned ** Indices);
        bundle_step_status_t (*EveryIteration)(struct bundle_callbacks *);
};

struct bundle_sparse_vector
{
        const double   * values;
        const unsigned * indices;
        unsigned count;
};

struct bundle_solver;

struct bundle_config {
        size_t iteration_count;
        size_t inactive_delay;
        size_t bundle_size;
        double t_enlarge_factor;
        double t_reduce_factor;
        double serious_step_factor;
        double medium_step_factor;
        double null_step_factor;

        double t_optimal;
        double eps_lin;

        double t_max;
        double t_min;
        double t_init;

        double t_strategy_1;
        double t_strategy_2;
        double t_strategy_eps;

        size_t pricing_warmup;
        size_t pricing_period;
        size_t pricing_max_age;
};

void 
bundle_config_init(struct bundle_config * OUT_config);
void 
bundle_config_parse(struct bundle_config * config, const char * config_string);

typedef struct bundle_solver bundle_solver_t;

bundle_solver_t * 
bundle_solver_create(const char * config_string, unsigned nvar);
bundle_solver_t * 
bundle_solver_create_from_config(const struct bundle_config * config,
                                 unsigned nvar);
void
bundle_solver_destroy(bundle_solver_t *);

void 
bundle_set_lambda(bundle_solver_t *, const double * lambda);
void
bundle_set_uc_range(bundle_solver_t *, unsigned start, unsigned end);
void
bundle_set_uc_indices(bundle_solver_t *, const unsigned * indices);

void
bundle_get_solution(bundle_solver_t *, struct bundle_sparse_vector * OUT_vector);
void
bundle_write_solution(bundle_solver_t *, double * OUT_lambda);
double 
bundle_get_Fi(bundle_solver_t *);
double 
bundle_get_best_Fi(bundle_solver_t *);

int 
bundle_current_solution_is_last(bundle_solver_t *);
bundle_status_t
bundle_get_status(bundle_solver_t *);

bundle_status_t
bundle_solve_with_callbacks(bundle_solver_t *, struct bundle_callbacks * callbacks);

void 
bundle_start_solving(bundle_solver_t *);
void 
bundle_stop_solving(bundle_solver_t *);
int 
bundle_solve_done_p(bundle_solver_t *);

void
bundle_get_current_lambda(bundle_solver_t *, struct bundle_sparse_vector * OUT_vector);
double * 
bundle_get_gi_dest(bundle_solver_t *, unsigned * name);
int 
bundle_register_values(bundle_solver_t *, double value,
                       const double * subgradient,
                       const unsigned * indices,
                       double eps);

#endif
