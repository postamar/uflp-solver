// (C) 2011, Marius Posta (postamar@iro.umontreal.ca)
// Check LICENSE.txt for the legal blah-blah.

#ifndef INSTANCE_H
#define INSTANCE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "RngStream.h"


double ddot (int n, const double *x, const double *y); // BLAS ddot_ wrapper

static const double bound_eps = 1e-3; // numerical error tolerance

typedef struct
{
	int i;
	double c;
} ic_pair_t;

int ic_pair_cmp (const void *a, const void *b);
void quickselect (RngStream rng, ic_pair_t *queue, int n, int k); // selection procedure


enum instance_format_t {inst_fmt_orlib, inst_fmt_cap}; // used to parse instance data
enum node_selection_t {search_bf, search_df}; // search strategy parameter values

typedef struct { // used to specify default search parameters
	char name[256];
	double primalbias;
	char strategy[256];
	int cachesize;
	int guiding;
	int bundlesize;
	int rootnodeiter;
	int subnodeiter;
	double time;
	double upperbound;
	char implbranch[255];
	char log[256]; 
	char output[256];
} search_param_t; 


typedef struct { // used to store search parameters, and instance data (size and costs)
	int subg_cache_size; // TODO: as it is, this parameter cannot be changed on the fly. fix this.

	int request_period;
	int bundle_max_size;
	int n_root_node_max_iter;
	int n_node_max_iter_phase1;
	int n_node_max_iter_phase2;
	double max_utime;
	double ub;
	double primal_bias;
	int fix;
	int log;
	int printsol;
	enum node_selection_t strategy;
	char name[64];

	enum instance_format_t fmt; // used to read the instance costs
	int n_words_in_bit_vector; // n-bit-vector size (in words)

	int n; // number of locations
	int m; // number of clients
	double *f, **c; // opening and service costs
	int **inc; // ordinals for service costs in increasing order, by client

} instance_t;


void instance_alloc (instance_t *inst, int n, int m);
instance_t *instance_copy (instance_t *inst); // returns deep copy of *inst
void instance_sync (instance_t *master, instance_t *slave); // copies costs from *master to *slave
void instance_destroy (instance_t *inst);
int instance_read_header (instance_t *inst, int first_read); // reads instance header from stdin
int instance_read (instance_t *inst); // reads instance costs from stdin
double instance_reset (instance_t *inst); // recomputes optimum upper bound and client orderings by cost
void instance_set_parameter (instance_t *inst, const char name[256], const char val[256]); // sets a search parameter 
int instance_parse_args (instance_t *inst, int argc, char **argv, const search_param_t defarg); // sets the search parameters as ordered in the command line (or else using the default values stored in defarg)




// ** Shared data **
// this data structure implements message passing in our method
//
enum search_state_t {state_waiting, state_running, state_solved, state_timeout, state_finish};

typedef struct
{
	instance_t *inst;
	int n_solved;

	double exec_time;
	double max_utime;

	int *guiding_x;
	int *improving_partial_x;
	int *best_x;
	double best_z;
	enum search_state_t search_state;

	double global_lb;
	int n_node_eval;
	int n_lag_eval;
	int n_moves;
} shared_t;

shared_t *shared_create (instance_t *inst);
void shared_destroy (shared_t *shared);
void shared_reset (shared_t *shared, instance_t *inst, double ub);
void shared_output_results (shared_t *shared);




// Bit vector related definitions
typedef unsigned long long word_t;

inline int get_bit (word_t *w, int i)
{
	int wi = (int) i / (8 * sizeof(word_t));
	int oi = i % (8 * sizeof(word_t));
	return (w[wi] & (((word_t) 1) << oi)) ? 1 : 0;
}


inline void set_bit (word_t *w, int i, int v)
{
	int wi = (int) i / (8 * sizeof(word_t));
	int oi = i % (8 * sizeof(word_t));
	if (v)
		w[wi] |= ((word_t) 1) << oi;
	else
		w[wi] &= ~(((word_t) 1) << oi);
}


#endif

