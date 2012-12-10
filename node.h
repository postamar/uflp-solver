// (C) 2011, Marius Posta (postamar@iro.umontreal.ca)
// Check LICENSE.txt for the legal blah-blah.

#ifndef NODE_H
#define NODE_H

// Include the correct BTT wrapper header,
#include "instance.h"
#ifdef __cplusplus
// compiling in c++ (e.g. g++)
#include "BTT//bundle_solver.hpp"
#else
// compiling in c99 (e.g. clang)
#include "BTT//bundle_solver.h"
#endif


// Adjustable array data structure
typedef struct
{
	void **e;
	int n;
	int alloc_size;
	int alloc_size_increment;
} adj_array_t;

adj_array_t *adj_array_create (int alloc_size_increment);
void adj_array_destroy (adj_array_t *adj_array);
void adj_array_push (adj_array_t *adj_array, void *e);
void *adj_array_pop (adj_array_t *adj_array);

// Branch-and-bound search tree node data structure
typedef struct node_t node_t;
struct node_t
{
	node_t *left, *right; // subnodes
	int min_open, max_open; // minimum and maximum number of open locations
	word_t *fix_mask; // n-bit-vector of fixed locations, mask
	word_t *fix_val; // n-bit-vector of fixed locations, value
	double *mul; // best multipliers
	double lb, ub; // lower and upper bounds
	int branch_i; // branching location
	int id; // node number (low is old, high is recent) 
};

node_t *node_alloc (instance_t *inst);
void node_set (instance_t *inst, node_t *parent, node_t *node); // in-place copy
void node_destroy (node_t *node); 
void node_create_children (instance_t *inst, node_t *parent); // creates parent->left and parent->right

// Heap functions
void heap_sift_down (adj_array_t *h_node, int root, enum node_selection_t strategy);
void heap_sift_up (adj_array_t *h_node, int current, enum node_selection_t strategy);


// Node evaluation data structures
typedef struct node_eval_t node_eval_t;

// Subgradient cache element data structure
typedef struct
{
	long min_open, max_open, n_x1;
	unsigned long x_hash, fixed_x_hash, subg_hash;
	word_t *x;
	word_t *fixed_x;
	double *subg;
	long double eps0;
} hyperplane_t;

// Argument to BTT callback functions
typedef struct
{
	struct bundle_callbacks callbacks;
	node_t *node;
	node_eval_t *node_eval;
} bundle_callback_data_t;


struct node_eval_t // stores everything necessary for node evaluation
{
	instance_t *inst;

	// statistical counters
	int node_counter;
	int n_node_eval;
	int n_lag_eval;

	// bundle-related stuff
	char *bundle_cfg;
	bundle_callback_data_t cbdata;
	bundle_solver_t *bundle;
	double *bundle_mul;
	double *bundle_mul_alloc;
	double *bundle_subg;
	int bundle_max_gi;
	int bundle_n_gi;
	int bundle_n_fi;

	// lag. relax. evaluation related stuff
	double z;
	double *rc;
	int *x;
	int *int_subg;
	int *x1;
	int n_x1;
	int n_fixed_0, n_fixed_1;
	ic_pair_t *queue;
	double *min_c;

	// impl. br. related stuff
	int *partition;
	int *I1;
	int *Is;
	int *I0;
	double z_b0, z_b1;

	// subgradient cache and bloom filter related stuff
	hyperplane_t current_hyperplane;
	hyperplane_t *cache;
	ic_pair_t *selection;
	int n_selected;
	int fresh_idx;
	unsigned short *bloom;
	int bloom_size;
	unsigned long *hash_salt;
	int n_hash_salt;

	// misc.
	RngStream rng;
};


node_eval_t *node_eval_create (instance_t *inst);
void node_eval_reset (node_eval_t *node_eval, double approx_z); 
void node_eval_destroy (node_eval_t *node_eval);

// Lagragian relaxation evaluation
int node_eval_lagrangian_fixed (node_eval_t *node_eval, node_t *node, double *mul);
double node_eval_lagrangian (node_eval_t *node_eval, node_t *node, double *mul);

// Bundle search & BTT callbacks
int node_eval_improve_mul (node_eval_t *node_eval, node_t *node, int max_subg, double ub);
double bundle_eval (struct bundle_callbacks *callbacks_ptr, const double *Lam_val, const unsigned *Lam_idx, unsigned Lam_count);
int bundle_get_subg (struct bundle_callbacks *callbacks_ptr, double *Eps, const unsigned **Indices);
void bundle_set_subg (struct bundle_callbacks *callbacks_ptr, double *subgradient, unsigned );

// Subgradient cache & Bloom filter functions
void select_hyperplanes_in_cache (node_eval_t *node_eval, node_t *node);
unsigned long get_hash_x (word_t *x, int n_words);
unsigned long get_hash_subg (int *subg, int m);
int check_bloom_filter (node_eval_t *node_eval);
void add_bloom_filter (node_eval_t *node_eval, short val);

// SSE intrinsics
double clamped_sum (double * x, double * y, size_t count);
void initialize_count (int * z, size_t count, int value);
void clamped_count (int * z, double * x, double * y, size_t count);

#endif
