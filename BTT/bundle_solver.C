// (C) 2009-2011, by Paul-Virak Khuong, under the MIT License in ../LICENSE.txt.

#include "bundle_solver.hpp"
// #include <pthread.h>
#include <cfloat>
#include <climits>
#include <cstdio>
#include <sstream>
#include <string>
#include "OPTvect.h"
#include "QPBundle.h"

namespace bundle {
#if 0
}
#endif

using namespace OPTtypes_di_unipi_it;


struct solver_impl : public Bundle_di_unipi_it::QPBundle
{
public:
        solver_impl(std::istream &iStrm, const unsigned NV);
        ~solver_impl();

        // no copy con
        solver_impl(const solver_impl&);

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
protected:
        void SetGi(double * SubG, const unsigned Name);
        double Fi(const double * Lam1, const unsigned * LamB, unsigned LamBd);
        int GetGi(double &Eps1, const unsigned * &SGBse);
        EIStatus EveryIteration(void);
#if DO_AGGR
        void AggrPrimalSol(const double * Mlt, const unsigned * NmSt,
                           const unsigned Dm, const unsigned NwNm);
#endif
private:
        //  W: transition by worker
        //  M:     "      "  main thread
        //
        //                          W
        //                   .----------------.
        //                  /                  '.
        //           W     v     M              |
        // starting --> waiting --> got_values -'
        //                |             |
        //                |M            |W
        //                v      W      v   M (on join)
        //              killed  --->  done ----> dead
        //
        enum loop_state_t {
                dead = 0, // can only be true if the worker thread is dead
                starting,
                waiting_for_values,
                got_values,
                done,
                killed
        };

        bool wait_for_state (loop_state_t waiting_state);
        bool set_state (loop_state_t new_state);
        static void * worker_routine(void * self);
        void worker_method();

        callbacks * callbacks_ptr;
/*
        pthread_t worker;

        // message "queue"
        pthread_mutex_t mutex;
        pthread_cond_t condvar;
*/
        // message
        double * volatile SubG;
        volatile unsigned Name;
        volatile double Fi_val;
        volatile sparse_vector Lambda;
        volatile sparse_vector Gi;
        volatile double GiEps;

        // queue state
        volatile loop_state_t loop_state;
};

solver_impl::solver_impl(istream &iStrm , cIndex NV)
        : QPBundle(&iStrm , NV),
          callbacks_ptr(0),
          SubG(0),
          Name(-1U),
          Fi_val(-HpINF),
          GiEps(nan("")),
          loop_state(dead)
{
#define check_enum(ENUM)                                                \
        assert((long)::ENUM == (long)Bundle_di_unipi_it::Bundle::ENUM)

        check_enum(kOK);
        check_enum(kUnbounded);
        check_enum(kUnbounded);
        check_enum(kAbort);
        check_enum(kMaxIter);
        check_enum(kError);

        check_enum(kEINorm);
        check_enum(kEIAbort);
        check_enum(kEILoopNow);
        check_enum(kEIContAnyway);
#undef check_enum
/*
        pthread_mutex_init(&mutex, NULL);
        pthread_cond_init(&condvar, NULL);
*/
        assert(InINF == -1U);
        assert(HpINF == DBL_MAX);
}

solver_impl::~solver_impl()
{
	/*
        pthread_mutex_destroy(&mutex);
        pthread_cond_destroy(&condvar);
		*/
}

void solver_impl::set_lambda(const double * lambda)
{ SetLambda(lambda); }

void solver_impl::set_uc(unsigned start, unsigned end)
{ SetUC(start, end, 0); }
void solver_impl::set_uc(const unsigned * indices)
{ SetUC(0, 0, indices); }

sparse_vector solver_impl::get_solution()
{
        sparse_vector ret;
        ret.values = ReadSol(ret.indices, ret.count);
        return ret;
}
void solver_impl::get_solution (double * lambda) { ReadSol(lambda); }
double solver_impl::get_Fi () { return ReadFiVal(); }
double solver_impl::get_best_Fi () { return ReadBestFiVal(); }

bool solver_impl::current_solution_is_last() { return CurrentIsLast(); }
status_t solver_impl::get_bundle_status ()
{ return (status_t)GetBStatus(); }

status_t solver_impl::solve_with_callbacks (callbacks * callbacks_ptr)
{
		RemoveItems();
        this->callbacks_ptr = callbacks_ptr;
        status_t ret = (status_t)Solve();
        this->callbacks_ptr = 0;
        return ret;
}
bool solver_impl::wait_for_state(loop_state_t waiting)
{
        assert(0 == callbacks_ptr);
        assert(dead != loop_state);
 
        bool ret = (loop_state != waiting);
 /*      
        if (ret
            && (loop_state != killed)
            && (loop_state != done)) {
                pthread_mutex_lock(&mutex);
                while ((loop_state != waiting)
                       && (loop_state != killed)
                       && (loop_state != done))
                        pthread_cond_wait(&condvar, &mutex);
                ret = loop_state != waiting;
                pthread_mutex_unlock(&mutex);
        }
*/
        return ret;
}

bool solver_impl::set_state(loop_state_t new_state)
{
        assert(0 == callbacks_ptr);
        assert(dead != loop_state);

        bool ret = 1;
/*
        pthread_mutex_lock(&mutex);
        if ((loop_state != done) && (loop_state != killed)) {
                loop_state = new_state;
                ret = 0;
                pthread_cond_signal(&condvar);
        }
        pthread_mutex_unlock(&mutex);
		*/
        return ret;
}

void solver_impl::start_solving()
{
	/*
        assert(0 == callbacks_ptr);
        assert(dead == loop_state);
        loop_state = starting;
        pthread_create(&worker, NULL, worker_routine, (void*)this);
        wait_for_state(waiting_for_values);
		*/
}

void solver_impl::stop_solving()
{
	/*
        assert(0 == callbacks_ptr);
        if (dead == loop_state) return;

        set_state(killed);

        pthread_join(worker, 0);
        loop_state = dead;
		*/
}

bool solver_impl::solve_done_p() 
{
	/*
        if (got_values == loop_state)
                wait_for_state(waiting_for_values);

        if (!((dead == loop_state)
              || (done == loop_state)
              || (killed == loop_state)))
                return 0;

        if (loop_state != dead) {
                pthread_join(worker, 0);
                loop_state = dead;
        }
*/
        return 1;
}

sparse_vector solver_impl::get_current_lambda ()
{
        assert(0 == callbacks_ptr);

        if (wait_for_state(waiting_for_values)) {
                Lambda.values = NULL;
                Lambda.indices = NULL;
                Lambda.count = 0;
        }

        sparse_vector ret;
        ret.values = Lambda.values;
        ret.indices = Lambda.indices;
        ret.count = Lambda.count;
        return ret;
}

double * solver_impl::get_gi_dest(unsigned * name)
{
        assert(0 == callbacks_ptr);

        if (wait_for_state(waiting_for_values)) {
                if (name) *name = -1U;
                return 0;
        }

        if (name) *name = Name;
        return SubG;
}

bool solver_impl::register_values(double value,
                                  const double * subgradient,
                                  const unsigned * indices,
                                  double eps)
{
        assert(0 == callbacks_ptr);
        if (wait_for_state(waiting_for_values))
                return 1;

        Fi_val = value;
        Gi.values = subgradient;
        Gi.indices = indices;
        GiEps = eps;

        return (set_state(got_values)
                || wait_for_state(waiting_for_values));
}

void * solver_impl::worker_routine (void * self)
{
        ((solver_impl*)self)->worker_method();
        return 0;
}

void solver_impl::worker_method(void)
{
        assert(starting == loop_state);
        Solve();
        set_state(done);
        SubG = NULL;
        Name = -1U;
        Fi_val = nan("");
        Lambda.values  = Gi.values  = 0;
        Lambda.indices = Gi.indices = 0;
        Lambda.count   = Gi.count   = 0;
        GiEps = nan("");
}

void solver_impl::SetGi(SgRow SubG, cIndex Name)
{
        if (callbacks_ptr)
                return callbacks_ptr->SetGi(callbacks_ptr, SubG, Name);

        this->SubG = SubG;
        this->Name = Name;
}

HpNum solver_impl::Fi(cLMRow Lam1, cIndex_Set LamB, Index LamBd)
{ 
        if (callbacks_ptr)
                return callbacks_ptr->Fi(callbacks_ptr, Lam1, LamB, LamBd); 

        Lambda.values  = Lam1;
        Lambda.indices = LamB;
        Lambda.count   = LamBd;

        if (set_state(waiting_for_values)
            || wait_for_state(got_values))
                return -HpINF;

        return Fi_val;
}

SIndex solver_impl::GetGi(HpNum &Eps1, cIndex_Set &SGBse)
{
        if (callbacks_ptr)
                return callbacks_ptr->GetGi(callbacks_ptr, &Eps1, &SGBse); 

        assert(got_values == loop_state);

        Eps1 = GiEps;
        SGBse = Gi.indices;
        if (Gi.values != SubG) {
                if (Gi.indices)
                        for (size_t i = 0; Gi.indices[i] < -1U; i++)
                                SubG[i] = Gi.values[i];
                else
                        for (size_t i = 0; i < NumVar; i++)
                                SubG[i] = Gi.values[i];
        }

#if SPRBL_FI != 0
# error "SPRBL_FI != 0 not supported."
#endif

        return 0;
}

solver_impl::EIStatus solver_impl::EveryIteration(void)
{
        if (callbacks_ptr) {
                if (callbacks_ptr->EveryIteration)
                        return (EIStatus)
                                callbacks_ptr->EveryIteration(callbacks_ptr);

                return kEINorm;
        }

        if (killed == loop_state)
                return kEIAbort;

        assert((starting == loop_state) || (got_values == loop_state));
        return kEINorm;
}

#if( DO_AGGR )
void 
solver_impl::AggrPrimalSol(cHpRow Mlt, cIndex_Set NmSt,
                           cIndex Dm, cIndex NwNm)
{
        (void)Mlt;
        (void)NmSt;
        (void)Dm;
        (void)NwNm;
}
#endif

// implementation hiding stubs
#define DEF(ret, name, args, vars)               \
        ret solver::name args { return impl->name vars; }

DEF(void, set_lambda, (const double * lambda), (lambda))
DEF(void, set_uc, (unsigned start, unsigned end), (start, end))
DEF(void, set_uc, (const unsigned * indices), (indices))

DEF(sparse_vector, get_solution, (void), ())
DEF(void, get_solution, (double * lambda), (lambda))
DEF(double, get_Fi, (void), ())
DEF(double, get_best_Fi, (void), ())

DEF(bool, current_solution_is_last, (void), ())
DEF(status_t, get_bundle_status, (void), ())

DEF(status_t, solve_with_callbacks,
    (callbacks * callbacks_ptr), 
    (callbacks_ptr))

DEF(void, start_solving, (void), ())
DEF(void, stop_solving, (void), ())
DEF(bool, solve_done_p, (void), ())
DEF(sparse_vector, get_current_lambda, (void), ())
DEF(double *, get_gi_dest, (unsigned * name), (name))
DEF(bool, register_values,
    (double value, 
     const double * subgradient,
     const unsigned * indices,
     double eps),
    (value, subgradient, indices, eps))
#undef DEF

solver::solver (std::istream &iStrm, const unsigned NV)
      : impl(new solver_impl(iStrm, NV))
{}

solver::~solver()
{ delete impl; }

config::config ()
{
        iteration_count = 100;

        inactive_delay = 10;
        bundle_size = 10;

        t_enlarge_factor = 10;
        t_reduce_factor = 0.1;

        serious_step_factor = 0.1;
        medium_step_factor = 0.3;
        null_step_factor = 3.0;

        t_optimal = 100.0;
        eps_lin = 1e-6;

        t_max = 1e6;
        t_min = 1e-6;
        t_init = 1.0;

        t_strategy_1 = 0.0;
        t_strategy_2 = 0.0;
        t_strategy_eps = 0.0;

        pricing_warmup = 30;
        pricing_period = 10;
        pricing_max_age = 5;
}

config::config (const bundle_config &tmp)
        : bundle_config(tmp)
{}

config::config (std::istream &str)
{
        str >> iteration_count

            >> inactive_delay
            >> bundle_size

            >> t_enlarge_factor
            >> t_reduce_factor
                
            >> serious_step_factor
            >> medium_step_factor
            >> null_step_factor
                
            >> t_optimal
            >> eps_lin

            >> t_max
            >> t_min
            >> t_init
                
            >> t_strategy_1
            >> t_strategy_2
            >> t_strategy_eps

            >> pricing_warmup
            >> pricing_period
            >> pricing_max_age;
}

#define DEF(name, type)                                               \
        config config::with_##name                                      \
        (type value) const                                              \
        {                                                               \
                config ret(*this);                                      \
                ret.name = value;                                       \
                return ret;                                             \
        }

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

void config::write_config(std::ostream &str) const
{
        str << iteration_count << std::endl

            << inactive_delay  << std::endl
            << bundle_size     << std::endl
                
            << t_enlarge_factor << std::endl
            << t_reduce_factor  << std::endl

            << serious_step_factor << std::endl
            << medium_step_factor  << std::endl
            << null_step_factor    << std::endl

            << t_optimal           << std::endl
            << eps_lin             << std::endl

            << t_max               << std::endl
            << t_min               << std::endl
            << t_init              << std::endl

            << t_strategy_1        << std::endl
            << t_strategy_2        << std::endl
            << t_strategy_eps      << std::endl

            << pricing_warmup      << std::endl
            << pricing_period      << std::endl
            << pricing_max_age     << std::endl;
}

solver config::make_solver (size_t nv) const
{
        std::stringstream str;
        write_config(str);
        
        return solver(str, nv);
}

#if 0
{
#endif
};

extern "C" {

// C wrappers
void bundle_config_init (struct bundle_config *OUT_config)
{
        *OUT_config = bundle::config();
}

void
bundle_config_parse (struct bundle_config * OUT_config, const char * string)
{
        std::stringstream str;
        str << string;
        *OUT_config = bundle::config(str);
}

bundle_solver_t *
bundle_solver_create (const char * config_string, unsigned nvar)
{
        std::stringstream stream;
        stream << config_string;
        return (bundle_solver_t*)new bundle::solver(stream, nvar);
}

bundle_solver_t *
bundle_solver_create_from_config(const bundle_config * config_struct,
                                 unsigned nvar)
{
        std::stringstream stream;
        ((bundle::config*)config_struct)->write_config(stream);
        return (bundle_solver_t*)new bundle::solver(stream, nvar);
}

void bundle_solver_destroy(bundle_solver_t * bundle)
{
        delete (bundle::solver*)bundle;
}

void bundle_set_lambda (bundle_solver_t * bundle, const double * lambda)
{
        ((bundle::solver*)bundle)->set_lambda(lambda);
}

void bundle_set_uc_range (bundle_solver_t * bundle, 
                          unsigned start, unsigned end)
{
        ((bundle::solver*)bundle)->set_uc(start, end);
}

void bundle_set_uc_indices(bundle_solver_t * bundle, const unsigned * indices)
{
        ((bundle::solver*)bundle)->set_uc(indices);
}

void
bundle_get_solution(bundle_solver_t * bundle, 
                    struct bundle_sparse_vector *OUT_vector)
{
        *OUT_vector = ((bundle::solver*)bundle)->get_solution();
}

void 
bundle_write_solution(bundle_solver_t * bundle, double * OUT_lambda)
{
        ((bundle::solver*)bundle)->get_solution(OUT_lambda);
}

double 
bundle_get_Fi(bundle_solver_t * bundle)
{
        return ((bundle::solver*)bundle)->get_Fi();
}

double 
bundle_get_best_Fi(bundle_solver_t * bundle)
{
        return ((bundle::solver*)bundle)->get_best_Fi();
}

int 
bundle_current_solution_is_last(bundle_solver_t * bundle)
{
        return ((bundle::solver*)bundle)->current_solution_is_last();
}

bundle_status_t
bundle_get_status (bundle_solver_t * bundle)
{
        return ((bundle::solver*)bundle)->get_bundle_status();
}

bundle_status_t
bundle_solve_with_callbacks (bundle_solver_t * bundle,
                             struct bundle_callbacks * callbacks_ptr)
{
        return ((bundle::solver*)bundle)->
                solve_with_callbacks((bundle::callbacks*)callbacks_ptr);
}

void bundle_start_solving(bundle_solver_t * bundle)
{
        ((bundle::solver*)bundle)->start_solving();
}

void bundle_stop_solving(bundle_solver_t * bundle)
{
        ((bundle::solver*)bundle)->stop_solving();
}

int bundle_solve_done_p(bundle_solver_t * bundle)
{
        return ((bundle::solver*)bundle)->solve_done_p();
}

void 
bundle_get_current_lambda(bundle_solver_t * bundle,
                          struct bundle_sparse_vector * OUT_vector)
{
        *OUT_vector = ((bundle::solver*)bundle)->get_current_lambda();
}

double * bundle_get_gi_dest(bundle_solver_t * bundle, unsigned * name)
{
        return ((bundle::solver*)bundle)->get_gi_dest(name);
}

int bundle_register_values(bundle_solver_t * bundle,
                           double value,
                           const double * subgradient,
                           const unsigned * indices,
                           double eps)
{
        return ((bundle::solver*)bundle)
                ->register_values(value, subgradient, indices, eps);
}

}
