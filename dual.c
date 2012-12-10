// (C) 2011, Marius Posta (postamar@iro.umontreal.ca)
// Check LICENSE.txt for the legal blah-blah.

#include "dual.h"


// performs an iteration of the dual process
void dual_run (dual_t *dual)
{
	switch(dual->detailed_state) {
		case dual_state_prep_phase1:
			dual_prep_phase1(dual); 
			break;
		case dual_state_prep_phase2:
			dual_prep_phase2(dual); 
			break;
		case dual_state_run_phase1:
		case dual_state_run_phase2:
			dual_iteration(dual);
			break;
		default:
			assert(0);
	}
}


void dual_prep_phase1 (dual_t *dual)
{
	instance_t *inst = dual->inst;
	int k;
	node_t *root;

	assert(dual->open_node_heap->n == 0);
	assert(dual->n_opened_nodes == 0);

	root = node_alloc(inst);
	node_set(dual->inst, dual->global, root);
	adj_array_push(dual->open_node_heap, root);	//% Algorithm 2 step 0

	for (k = 0; k <= dual->inst->n; k++)
		dual->lb_k[k] = INFINITY;

	dual->detailed_state = dual_state_run_phase1;
	dual_new_upper_bound(dual);

	if (inst->log) 
		fprintf(stderr, "%8s\t[ %5s , %5s ]\t%15s\t%15s\t%15s\t%12s\t%12s\n", "#open", "nmin", "nmax", "node lb", "global lb", "global ub", "#node eval", "#lag. eval");
}


void dual_prep_phase2 (dual_t *dual)
{
	instance_t *inst = dual->inst;
	node_t *root;

	dual->phase = 2;
	assert(dual->open_node_heap->n == 0);

	root = node_alloc(inst);
	node_set(dual->inst, dual->global, root); 
	adj_array_push(dual->open_node_heap, root); //% Algorithm 2 step 0

	dual->detailed_state = dual_state_run_phase2;
	dual_new_upper_bound(dual);

	if (inst->log) 
		fprintf(stderr, "%8s\t[ %5s , %5s ]\t%15s\t%15s\t%15s\t%12s\t%12s\n", "#open", "nmin", "nmax", "node lb", "global lb", "global ub", "#node eval", "#lag. eval");
}


// Evaluates and separates a node of the branch and bound search tree
void dual_iteration (dual_t *dual)
{
	instance_t *inst = dual->inst;
	node_eval_t *node_eval = dual->node_eval;
	node_t *node;
	int max_subg = (dual->phase == 1) ? inst->n_node_max_iter_phase1 : inst->n_node_max_iter_phase2;

	// ** Node removal and selection **
	//% Algorithm 2 step 4(d)
	if (dual->open_node_heap->n == 0) {		// no more open nodes in the queue?
		if (dual->phase == 1) { 
			dual->detailed_state = dual_state_prep_phase2; // preprocessing is over
		} else {
			dual->lb = dual->best_z;
			dual->detailed_state = dual_state_asleep;
			dual->state = dual_state_solved; // search is over
		}
		return; // in either case, nothing left to do right now
	}
	//% Algorithm 2 step 4(b)
	dual_update_improving_partial(dual); 

	//% Algorithm 2 step 1
	node = (node_t*) dual->open_node_heap->e[0]; // select the best node, i.e. the node at the top of the heap
	if (inst->strategy == search_bf) 
		dual->lb = node->lb; // if the node search strategy is BFS, then node lb is global

	// Remove the node from the top of the heap
	if (dual->open_node_heap->n == 1) {
		dual->open_node_heap->e[0] = NULL;
		dual->open_node_heap->n = 0;
	} else {
		dual->open_node_heap->e[0] = (node_t*) adj_array_pop(dual->open_node_heap);
		heap_sift_down(dual->open_node_heap, 0, inst->strategy);
	}

	if (dual->lb > dual->best_z - bound_eps) {
		while (dual->open_node_heap->n > 0)
			node_destroy((node_t*) adj_array_pop(dual->open_node_heap));
		if (inst->log)
			fprintf(stderr, "%8i\t[ %-5i , %5i ]\t      cutoff\t%12.2f\t%12.2f\t%12i\t%12i\n", dual->open_node_heap->n, node->min_open, node->max_open, dual->lb, dual->best_z, node_eval->n_node_eval, node_eval->n_lag_eval);
		node_destroy(node);
		return;
	}

	// ** Node evaluation **

	// Goldengorin et al. preprocessing - disabled 
/*	if (dual_node_goldengorin(dual, node)) {
		if (node->ub < dual->best_z - bound_eps) {
			dual->best_z = node->ub;
			memcpy(dual->best_x, dual->node_eval->x, inst->n * sizeof(int));
		}
	}
	else */ 

	// Perform bundle search to improve the node multipliers
	//% Algorithm 2 step 2
	if (node_eval_improve_mul(dual->node_eval, node, max_subg, dual->best_z)) {
		// A better primal solution was found!
		//% Algorithm 4 step (a)i. - step (a)ii. is implicit
		dual->best_z = node->ub;
		memcpy(dual->best_x, dual->node_eval->x, inst->n * sizeof(int));
	} 

	//% Algorithm 2 step 3(a)
	if (node->lb > dual->best_z - bound_eps) { // fathom the node
		if (inst->log)
			fprintf(stderr, "%8i\t[ %-5i , %5i ]\t      cutoff\t%12.2f\t%12.2f\t%12i\t%12i\n", dual->open_node_heap->n, node->min_open, node->max_open, dual->lb, dual->best_z, node_eval->n_node_eval, node_eval->n_lag_eval);
		node_destroy(node);
		return;
	}

	//% Algorithm 2 step 3(b)

	// Apply implicit branching
	node = dual_node_tighten_bounds(dual, node); 
	memcpy(dual->guiding_x, node_eval->x, inst->n * sizeof(int));
	if (node == NULL) { // implicit branching fathoms the node
		// note: not sure that this ever happens in the current implementation
		if (inst->log)
			fprintf(stderr, "%8i\t[       ,       ]\t      cutoff\t%12.2f\t%12.2f\t%12i\t%12i\n", dual->open_node_heap->n, dual->lb, dual->best_z, node_eval->n_node_eval, node_eval->n_lag_eval);
		return;
	}

	if (inst->log) 
		fprintf(stderr, "%8i\t[ %-5i , %5i ]\t%12.2f\t%12.2f\t%12.2f\t%12i\t%12i\n", dual->open_node_heap->n, node->min_open, node->max_open, node->lb, dual->lb, dual->best_z, node_eval->n_node_eval, node_eval->n_lag_eval);

	//% Algorithm 3 - case \underline{n} = \bar{n}
	if (dual->phase == 1 && node->min_open == node->max_open) { 
		// this node is terminal in the preprocessing tree
		dual->lb_k[node->min_open] = node->lb;
		node_destroy(node);
		return;
	}
	//% Algorithm 4 - case p \in \{0,1\}^n
	if (dual->phase == 2 && node->branch_i == -1) { 
		// all locations are fixed, cannot branch further 
		node_destroy(node);
		return;
	}


	// ** Node separation ** 
	
	node_create_children(inst, node);
	dual->n_opened_nodes += 2;
	node->right->id = ++node_eval->node_counter;
	node->left->id = ++node_eval->node_counter;

	if (dual->phase == 1) {
		//% Algorithm 3 - case otherwise
		node->left->max_open = lround(floor((double) (node->min_open + node->max_open) / 2.0));
		node->right->min_open = node->left->max_open + 1;
	} else {
		//% Algorithm 4 - case otherwise
		set_bit(node->left->fix_mask, node->branch_i, 1);
		set_bit(node->left->fix_val, node->branch_i, 0);
		assert(node->left->max_open == node->max_open);
		if (node->max_open == inst->n - node_eval->n_fixed_0)
			node->left = dual_node_decrease_max_open(dual, node->left, 1);

		set_bit(node->right->fix_mask, node->branch_i, 1);
		set_bit(node->right->fix_val, node->branch_i, 1);
		assert(node->right->max_open == node->max_open);
		if (node->min_open == node_eval->n_fixed_1)
			node->right = dual_node_increase_min_open(dual, node->right, 1);
	}

	// Add each subnode to queue (or fathom it)
	//% Algorithm 2 step 3(c)
	if (node->right != NULL) {
		node->right->lb = node_eval->z_b1;
		if (node->right->lb <= dual->best_z - bound_eps) {
			adj_array_push(dual->open_node_heap, node->right);
			heap_sift_up(dual->open_node_heap, dual->open_node_heap->n - 1, inst->strategy);
		} else {
			node_destroy(node->right);
			node->right = NULL;
		}
	}
	if (node->left != NULL) {
		node->left->lb = node_eval->z_b0;
		if (node->left->lb <= dual->best_z - bound_eps) {
			adj_array_push(dual->open_node_heap, node->left);
			heap_sift_up(dual->open_node_heap, dual->open_node_heap->n - 1, inst->strategy);
		} else {
			node_destroy(node->left);
			node->left = NULL;
		}
	}

	// The current node has been evaluated
	node_destroy(node);
}


// Updates the search state after the upper bound has been improved
//% Algorithm 2 step 4(a)
void dual_new_upper_bound (dual_t *dual)
{
	instance_t *inst = dual->inst;
	node_t *node;
	int diff, n = dual->open_node_heap->n;
	int k;

	// Identify min_open and max_open using bounds from preprocessing
	if (dual->phase == 2) {
		for (; dual->global->min_open < dual->inst->n; dual->global->min_open++)
			if (dual->lb_k[dual->global->min_open] <= dual->best_z - bound_eps)
				break;

		for (; dual->global->max_open > 1; dual->global->max_open--)
			if (dual->lb_k[dual->global->max_open] <= dual->best_z - bound_eps)
				break;
	}


	// Reset the global lower bound, we shall recompute it
	dual->lb = INFINITY;

	// Update the open node queue
	for (k = 0; k < n; k++) {
		node = (node_t*) dual->open_node_heap->e[k];
		assert(node != NULL);
		if (node->lb <= dual->best_z - bound_eps) { 
			// The node remains open
			if (dual->lb > node->lb)
				dual->lb = node->lb; // update the global lower bound
			if (node->ub > dual->best_z) 
				node->ub = dual->best_z; // update the local upper bound (this flags the node as being up to date)
			diff = dual->global->min_open - node->min_open;
			if (dual->phase == 2 && diff > 0) // update min_open (ImBranchP)
				node = dual_node_increase_min_open(dual, node, diff);
			diff = node->max_open - dual->global->max_open;
			if (dual->phase == 2 && diff > 0) // update max_open (ImBranchP)
				node = dual_node_decrease_max_open(dual, node, diff);
		} else {
			// The node is fathomed
			node_destroy((node_t*) dual->open_node_heap->e[k]);
			dual->open_node_heap->e[k] = NULL;
		}
	}

	// Heapify remaining open nodes
	for (dual->open_node_heap->n = k = 0; k < n; k++) {
		node = (node_t*) dual->open_node_heap->e[k];
		dual->open_node_heap->e[k] = NULL;
		if (node != NULL) {
			assert(node->lb < dual->best_z - bound_eps);
			dual->open_node_heap->e[dual->open_node_heap->n++] = node;
		}
	}
	for (k = lround(floor((dual->open_node_heap->n - 2) / 2)); k >= 0; --k) 
		heap_sift_down(dual->open_node_heap, k, inst->strategy);


	// Update improving partial solution
	dual->improving_partial_counter = inst->request_period;
	dual_update_improving_partial(dual);
}



void dual_update_improving_partial (dual_t *dual) 
{
	instance_t *inst = dual->inst;
	node_t *node;
	int wi, i, k;
	word_t fix0, fix1, fix_mask, fix_val;

	if (dual->phase == 1) 
		return; // do not perform during preprocessing
	if (dual->improving_partial_counter++ < inst->request_period)
		return; // do not perform too often
	dual->improving_partial_counter = 0;

	// Bit vector magic
	for (wi = 0; wi < inst->n_words_in_bit_vector; wi++) { // iterate over each word
		fix0 = fix1 = ~0; // initialized to 11...1
		for (k = 0; k < dual->open_node_heap->n; k++) { // iterate over each open node
			node = (node_t*) dual->open_node_heap->e[k];
			assert(node != NULL);
			fix1 &= node->fix_mask[wi] & node->fix_val[wi]; // AND all locations fixed to 1
			fix0 &= node->fix_mask[wi] & (~node->fix_val[wi]); // AND all locations fixed to 0
		}
		// Right now, a bit in fix1 is set to 1 if the corresponding 
		// location is forced open in all open nodes of the queue,
		// and likewise for fix0 and locations forced closed.
		// We now generate the improving partial solution mask and value words.
		fix_mask = fix0 ^ fix1; // XOR (OR would be OK too) to generate the mask
		fix_val = fix1 & fix_mask; // zero all bits not covered by the mask
		dual->global->fix_mask[wi] |= fix_mask; // '=' instead of '|=' would be OK too
		dual->global->fix_val[wi] |= fix_val;
	}

	// Set and check the improving partial solution.
	// Once a component is fixed it should never change
	for (i = 0; i < inst->n; i++)
		if (dual->improving_partial_x[i] == -1) {
			if (get_bit(dual->global->fix_mask, i)) 
				dual->improving_partial_x[i] = get_bit(dual->global->fix_val, i);
		} else {
			assert(get_bit(dual->global->fix_mask, i));
			assert(get_bit(dual->global->fix_val, i) == dual->improving_partial_x[i]);
		}
}


// Performs Erlenkotter's dual ascent heuristic on the 'global' lagrangian dual
//% Algorithm 2 step 0
double dual_ascent (instance_t *inst, double *mul)
{
	double *slack = (double*) calloc(inst->n, sizeof(double));
	int *k = (int*) calloc(inst->m, sizeof(int));
	int i, j, flag = 1;
	double delta, z = 0.0;
	
	// initialize dual ascent
	for (j = 0; j < inst->m; j++) {
		if (mul[j] < inst->c[inst->inc[j][0]][j])
			mul[j] = inst->c[inst->inc[j][0]][j];
		z += mul[j];
		while (k[j] < inst->n && mul[j] >= inst->c[inst->inc[j][k[j]]][j])
			++k[j];
	}
	for (i = 0; i < inst->n; i++) 
		slack[i] = inst->f[i];

	// dual ascent
	while (flag) {
		flag = 0;
		for (j = 0; j < inst->m; j++) {
			for (i = 0, delta = INFINITY; i < inst->n; i++)
				if (mul[j] >= inst->c[i][j] && slack[i] < delta)
					delta = slack[i];
			if (k[j] < inst->n && delta > inst->c[inst->inc[j][k[j]]][j] - mul[j]) {
				delta = inst->c[inst->inc[j][k[j]++]][j] - mul[j];
				flag = 1;
			}
			for (i = 0; i < inst->n; i++) 
				if (mul[j] >= inst->c[i][j])
					slack[i] -= delta;
			mul[j] += delta;
			z += delta;
		}
	}

	free(slack);
	free(k);

	return z;
}


// Allocates memory for the data structures of the dual process
dual_t *dual_create (instance_t *inst)
{
	int i;
	dual_t *dual = (dual_t*) malloc(sizeof(dual_t));

	dual->inst = instance_copy(inst);
	dual->state = dual_state_undef;
	dual->detailed_state = dual_state_asleep;

	dual->guiding_x = (int*) calloc(inst->n, sizeof(int));
	dual->improving_partial_x = (int*) calloc(inst->n, sizeof(int));
	dual->best_x = (int*) calloc(inst->n, sizeof(int));
	dual->best_z = INFINITY;

	for (i = 0; i < inst->n; i++) {
		dual->best_x[i] = dual->guiding_x[i] = 0;
		dual->improving_partial_x[i] = -1;
	}

	dual->lb_k = (double*) calloc(inst->n + 1, sizeof(double));
	dual->node_eval = node_eval_create(inst);
	dual->open_node_heap = adj_array_create(1024);
	dual->global = node_alloc(dual->inst);

	return dual;
}


// Frees memory allocated to the data structures of the dual process
void dual_destroy (dual_t *dual)
{
	int i;
	for (i = 0; i < dual->open_node_heap->n; i++) {
		assert(dual->open_node_heap->e[i] != NULL);
		node_destroy((node_t*) dual->open_node_heap->e[i]);
	}
	for (; i < dual->open_node_heap->alloc_size; i++) 
		assert(dual->open_node_heap->e[i] == NULL);

	adj_array_destroy(dual->open_node_heap);

	node_destroy(dual->global);

	node_eval_destroy(dual->node_eval);

	free(dual->lb_k);

	free(dual->best_x);
	free(dual->improving_partial_x);
	free(dual->guiding_x);
	instance_destroy(dual->inst);

	free(dual);
}


// Resets the dual process search state prior to new optimization,
// assuming that instance size remains the same
void dual_reset (dual_t *dual)
{
	int i, k, no_ascent_flag;

	dual->state = dual_state_run;
	dual->detailed_state = dual_state_asleep;
	dual->phase = 1;
	dual->lb = 0.0;
	dual->best_z = INFINITY;

	for (k = 0; k <= dual->inst->n; k++)
		dual->lb_k[k] = 0.0;

	for (i = 0; i < dual->inst->n; i++) {
		dual->best_x[i] = dual->guiding_x[i] = 0;
		dual->improving_partial_x[i] = -1;
	}

	dual->global->lb = 0.0;
	dual->global->ub = INFINITY;
	dual->global->min_open = 1;
	dual->global->max_open = dual->inst->n;
	for (k = 0; k < dual->inst->n_words_in_bit_vector; k++)
		dual->global->fix_mask[k] = dual->global->fix_val[k] = (word_t) 0;

	for (no_ascent_flag = i = 0; i < dual->inst->n && !no_ascent_flag; i++)
		no_ascent_flag = (dual->global->mul[i] != 0.0);

	if (no_ascent_flag) 
		dual->global->lb = dual->lb = 0.0;
	else 
		dual->global->lb = dual->lb = dual_ascent(dual->inst, dual->global->mul);

	node_eval_reset(dual->node_eval, dual->lb);

	while (dual->open_node_heap->n > 0)
		adj_array_pop(dual->open_node_heap);

	dual->n_opened_nodes = 0;
	dual->improving_partial_counter = 0;

	if (dual->inst->fix)
		dual->detailed_state = dual_state_prep_phase1;
	else
		dual->detailed_state = dual_state_prep_phase2;
}



