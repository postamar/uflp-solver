// (C) 2011, Marius Posta (postamar@iro.umontreal.ca)
// Check LICENSE.txt for the legal blah-blah.


#ifndef PRIMAL_H
#define PRIMAL_H

#include "instance.h"


// Data structure for storing a solution moving in the 1OPT neighborhood
typedef struct
{
	double total_cost; 
	double ub;

	int *x;
	double *flip_gain; // array mapping each location to total cost differences following 1OPT move

	// for each client, all open locations are stored in a heap, sorted by service cost
	int heap_size; // number of open locations
	int **heap; // locations heap
	int **heap_inv; // for each client, array mapping each location to index in heap, or -1 if not in heap (i.e. location is not open)

	int *besti1; // array mapping each client to cheapest location (-1 if all locations are closed)
	int *besti2; // array mapping each client to second cheapest location (-1 if only one location is open)
} solution_t;

solution_t *solution_alloc (instance_t *inst);
void solution_reset (instance_t *inst, solution_t *s, int *x); // replaces current solution with 'x'
void solution_flip (instance_t *inst, solution_t *s, int w); // performs 1OPT move on location w
void check_solution (instance_t *inst, solution_t *s); // used for debugging



// States of execution for the primal process
enum primal_state_t {primal_state_undef, primal_state_wait, primal_state_run};

// Data structure for the primal worker process
typedef struct {
	instance_t *inst;
	enum primal_state_t state;

	// shared data
	int *guiding_x;
	int *improving_partial_x;
	int *best_x;
	double best_z;

	// current solution
	solution_t *sol;

	// tabu states
	double *tabu;
	double tabu_length;

	// statistics
	int n_moves;
	int n_moves_at_last_improvement;

	// a strong random number generator used for move selection
	RngStream rng;
} primal_t;


primal_t *primal_create (instance_t *inst);
void primal_destroy (primal_t *primal);
void primal_reset (primal_t *primal); // prepares resolution of new instance with different costs (but same size)
void primal_new_upper_bound (primal_t *primal); // resets tabu search state with new best solution
void primal_run (primal_t *primal); // performs one iteration of primal process



#endif

