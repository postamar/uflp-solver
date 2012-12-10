// (C) 2009-2011, by Paul-Virak Khuong, under the MIT License in ../LICENSE.txt.

#ifndef __BUNDLE_SOLVER_HPP__
#define __BUNDLE_SOLVER_HPP__
#include <iostream>

extern "C" {
#include "bundle_solver.h"
};

namespace bundle {
#if 0
}
#endif

typedef bundle_status_t status_t;
typedef bundle_step_status_t step_status_t;

struct callbacks : bundle_callbacks
{
        callbacks ()
        {
                SetGi          = 0;
                Fi             = 0;
                GetGi          = 0;
                EveryIteration = 0;
        };

        callbacks (const bundle_callbacks &tmp)
                : bundle_callbacks(tmp)
        {};
};

struct sparse_vector : bundle_sparse_vector
{
        sparse_vector ()
        {
                values  = 0;
                indices = 0;
                count   = 0;
        };

        sparse_vector (const bundle_sparse_vector &tmp)
                : bundle_sparse_vector(tmp)
        {};
};

struct solver;

struct config : bundle_config {
        config ();
        config (const bundle_config &tmp);
        explicit config (std::istream &str);

#  define DEF(name, type)                               \
        config with_##name (type value) const;

        DEF(iteration_count, size_t)

        DEF(inactive_delay, size_t)
        DEF(bundle_size, size_t)

        DEF(t_enlarge_factor, double)
        DEF(t_reduce_factor, double)

        DEF(serious_step_factor, double)
        DEF(medium_step_factor, double)
        DEF(null_step_factor, double)

        DEF(t_optimal, double)
        DEF(eps_lin, double)

        DEF(t_max, double)
        DEF(t_min, double)
        DEF(t_init, double)

        DEF(t_strategy_1, double)
        DEF(t_strategy_2, double)
        DEF(t_strategy_eps, double)

        DEF(pricing_warmup, size_t)
        DEF(pricing_period, size_t)
        DEF(pricing_max_age, size_t)
#undef DEF

        void write_config(std::ostream &str) const;
        solver make_solver (size_t nv) const;
};

struct solver
{
        solver(std::istream &iStrm, const unsigned NV);
        ~solver();

        // no copy con
        solver(const solver&);

        void set_lambda(const double * lambda);
        void set_uc(unsigned start, unsigned end);
        void set_uc(const unsigned * indices);

        sparse_vector get_solution();
        void get_solution(double * lambda);
        double get_Fi(void);
        double get_best_Fi(void);

        bool current_solution_is_last(void);
        status_t get_bundle_status(void);

        status_t solve_with_callbacks(callbacks *);

        void start_solving();
        void stop_solving();
        bool solve_done_p();
        sparse_vector get_current_lambda();
        double * get_gi_dest(unsigned * name);
        bool register_values(double value,
                             const double * subgradient,
                             const unsigned * indices,
                             double eps);

        struct solver_impl * const impl;
};
#if 0
{
#endif
};
#endif
