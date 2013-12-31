// (C) 2011, Marius Posta (postamar@iro.umontreal.ca)
// Check LICENSE.txt for the legal blah-blah.



#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "primal.h"
#include "dual.h"


void solve_instance (primal_t *primal, dual_t *dual, shared_t *shared);


int main (int argc, char **argv)
{
	instance_t *inst = (instance_t*) calloc(1, sizeof(instance_t));
	shared_t *shared;
	primal_t *primal;
	dual_t *dual;
	search_param_t default_param = {"noname", 1.0, "bfs", 127, 250, 40, 2000, 50, INFINITY, INFINITY, "on", "off", "short"};

	// Parse command line arguments
	if (instance_parse_args(inst, argc, argv, default_param)) {
		free(inst);
		return 1;
	}

	// Read parameter settings and instance header
	if (!instance_read_header(inst, 1)) {
		free(inst);
		return 1;
	}

	// Create data structures 
	instance_alloc(inst, inst->n, inst->m);
	shared = shared_create(instance_copy(inst));
	primal = primal_create(inst);
	dual = dual_create(inst);
	
	// Main loop
	do {
		// Read instance data (i.e. opening costs and service costs)
		if (!instance_read(inst)) {
			fprintf(stderr, "ERROR: failed to read instance data.\n");
			primal_destroy(primal);
			dual_destroy(dual);
			shared_destroy(shared);
			instance_destroy(inst);
			return 1;
		}
		shared_reset(shared, inst, instance_reset(inst));
		
		// Solve instance
		solve_instance(primal, dual, shared);

		// Output computation results
		shared_output_results(shared);

		// Read parameter settings and instance header 
	} while (instance_read_header(inst, 0));

	// Clean up
	shared->search_state = state_finish;  
	primal_destroy(primal);
	dual_destroy(dual);
	shared_destroy(shared);
	instance_destroy(inst);
	return 0;
}



// Performs one iteration of the primal process
void primal_tick (primal_t *primal, shared_t *shared)
{
	int new_upper_bound; // flag
	const size_t vsize = shared->inst->n * sizeof(int);

	if (shared->search_state != state_running) {
		primal->state = primal_state_wait;
		return;
	}

	assert(primal->state != primal_state_wait);

	//% Algorithm 1 step 4.
	// Update shared computation data (message passing)
	new_upper_bound = 0;
	shared->n_moves = primal->n_moves;
	memcpy(primal->guiding_x, shared->guiding_x, vsize);
	memcpy(primal->improving_partial_x, shared->improving_partial_x, vsize);
	if (!(primal->best_z == INFINITY && shared->best_z == INFINITY)) {
		if (primal->best_z <= shared->best_z - bound_eps) {
			new_upper_bound = 1;
			shared->best_z = primal->best_z;
			memcpy(shared->best_x, primal->best_x, vsize);
		} else if (shared->best_z < primal->best_z - bound_eps) {
			new_upper_bound = 1;
			primal->best_z = shared->best_z;
			memcpy(primal->best_x, shared->best_x, vsize);
		}
	}

	//% Algorithm 1 step 4 (cont'd) followed by steps 1 to 3.
	if (new_upper_bound) 
		primal_new_upper_bound(primal);
	primal_run(primal);
}

// Performs one iteration of the dual process
void dual_tick (dual_t *dual, shared_t *shared)
{
	int run, new_upper_bound; // flags
	const size_t vsize = shared->inst->n * sizeof(int);

	if (shared->search_state != state_running) {
		dual->state = dual_state_wait;
		return;
	}

	assert(dual->state != dual_state_wait);
	run = new_upper_bound = 0;

	// Update shared computation data (message passing)
	shared->global_lb = dual->lb;
	shared->n_node_eval = dual->node_eval->n_node_eval;
	shared->n_lag_eval = dual->node_eval->n_lag_eval;
	memcpy(shared->guiding_x, dual->guiding_x, vsize);
	memcpy(shared->improving_partial_x, dual->improving_partial_x, vsize);
	if (!(dual->best_z == INFINITY && shared->best_z == INFINITY)) {
		if (dual->best_z <= shared->best_z - bound_eps) {
			new_upper_bound = 1;
			shared->best_z = dual->best_z;
			memcpy(shared->best_x, dual->best_x, vsize);
		} else if (shared->best_z <= dual->best_z - bound_eps) {
			new_upper_bound = 1;
			dual->best_z = shared->best_z;
			memcpy(dual->best_x, shared->best_x, vsize);
		}
	}

	// Check for termination
	if (dual->state == dual_state_solved) 
		shared->search_state = state_solved;
	else 
		run = 1;

	// Do the work
	if (new_upper_bound)
		dual_new_upper_bound(dual);

	if (run) 
		dual_run(dual);
}


// Returns 1 if time limit is reached, 2 if instance is solved, 0 otherwise.
inline int check_done (const struct rusage *begin, shared_t *shared)
{
	struct rusage current_usage;

	getrusage(RUSAGE_SELF, &current_usage);
	shared->exec_time = (current_usage.ru_utime.tv_sec - begin->ru_utime.tv_sec);
	shared->exec_time += (current_usage.ru_utime.tv_usec - begin->ru_utime.tv_usec) * 1e-6;

	if (shared->exec_time >= shared->inst->max_utime) {
		shared->search_state = state_timeout;
		return 1;
	} else if (shared->search_state == state_solved) {
		++shared->n_solved;
		return 2;
	}
	return 0;
}


// Solves the current instance
void solve_instance (primal_t *primal, dual_t *dual, shared_t *shared)
{
	struct rusage init_usage_solve;

	// process scheduling variables
	double balance_factor = 2.0, balance_factor_limit = 1000.0;
	int n_consec_primal = 0, n_consec_dual = 0, n_consec_primal_max = 10, n_consec_dual_max = 100;

	if (shared->search_state == state_finish) {
		primal->state = primal_state_undef;
		dual->state = dual_state_undef;
		return;
	}

	// Copy instance data to each process
	instance_sync(shared->inst, primal->inst); 
	instance_sync(shared->inst, dual->inst); 
	getrusage(RUSAGE_SELF, &init_usage_solve);

	// Initialize the search
	primal_reset(primal);
	dual_reset(dual);

	// Perform primal and dual iterations according to scheduling policy
	while (!check_done(&init_usage_solve, shared)) {
		if (n_consec_primal < n_consec_primal_max 
				&& shared->inst->primal_bias > 0.0
				&& (n_consec_dual >= n_consec_dual_max
					|| (double) shared->n_moves * balance_factor < (double) shared->n_lag_eval)) {
			primal_tick(primal, shared);
			++n_consec_primal;
			n_consec_dual = 0;
		} else {
			dual_tick(dual, shared);
			n_consec_primal = 0;
			++n_consec_dual;
		}
		if (balance_factor > balance_factor_limit || shared->inst->primal_bias == 0.0) 
			balance_factor = balance_factor_limit;
		else 
			balance_factor += 1.0 / (balance_factor * shared->inst->primal_bias);
	}
}

