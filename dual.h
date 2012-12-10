// (C) 2011, Marius Posta (postamar@iro.umontreal.ca)
// Check LICENSE.txt for the legal blah-blah.

#ifndef DUAL_H
#define DUAL_H

#include "instance.h"
#include "node.h"


// States of execution for the dual worker thread
enum dual_state_t {dual_state_undef, dual_state_wait, dual_state_run, dual_state_solved};
enum dual_detailed_state_t {dual_state_asleep, dual_state_prep_phase1, dual_state_run_phase1, dual_state_prep_phase2, dual_state_run_phase2};



// Data structure for the dual worker thread
typedef struct {
	instance_t *inst;

	// state variables
	enum dual_state_t state;
	enum dual_detailed_state_t detailed_state;
	int phase; // set to 1 during preprocessing, to 2 during resolution

	// shared data
	int *guiding_x;
	int *improving_partial_x;
	int *best_x;
	double best_z;

	// lower bounds
	double lb; // global
	double *lb_k; // for each k, for which sum_i x_i == k

	// data structure for evaluating bounds on nodes
	node_eval_t *node_eval;

	// open nodes are stored in a min-heap
	adj_array_t *open_node_heap;
	int n_opened_nodes;
	int improving_partial_counter;

	// bogus node used to store global min_open and max_open
	// and also for its data structures for computing the improv. part. sol.
	node_t *global;
} dual_t;


dual_t *dual_create (instance_t *inst);
void dual_destroy (dual_t *dual);

void dual_reset (dual_t *dual); // reinitializes everything depending on instance costs (but not size)
void dual_new_upper_bound (dual_t *dual); // updates the open node queue
void dual_run (dual_t *dual); // performs one iteration of the dual process, managing state

void dual_prep_phase1 (dual_t *dual); // initializes preprocessing
void dual_prep_phase2 (dual_t *dual); // initializes resolution
void dual_iteration (dual_t *dual); // evaluates a node of the b-and-b tree

node_t *dual_node_tighten_bounds (dual_t *dual, node_t *node); // performs implicit branching
void dual_update_improving_partial (dual_t *dual); // updates improv. part. sol.

double dual_ascent (instance_t *inst, double *mul); // applies dual ascent heuristic

node_t *dual_node_fix_open (dual_t *dual, node_t *node, int k); // fixes location k to open in node
node_t *dual_node_fix_closed (dual_t *dual, node_t *node, int k); // fixes location k to close in node
node_t *dual_node_increase_min_open (dual_t *dual, node_t *node, int n); // increases min_open by n in node
node_t *dual_node_decrease_max_open (dual_t *dual, node_t *node, int n); // decreases max_open by n in node

int dual_node_goldengorin (dual_t *dual, node_t *node); // Goldengorin et al. preprocessing rules

#endif

