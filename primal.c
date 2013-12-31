// (C) 2011, Marius Posta (postamar@iro.umontreal.ca)
// Check LICENSE.txt for the legal blah-blah.


#include "primal.h"


// Resets tabu search to new best solution
void primal_new_upper_bound (primal_t *primal)
{
	int i;

	solution_reset(primal->inst, primal->sol, primal->best_x);

	for (i = 0; i < primal->inst->n; i++)
		if (primal->tabu[i] != INFINITY)
			primal->tabu[i] = 0.0;
}


// Resets primal process prior to solving new instance,
// with different costs but same size
//% Algorithm 1 step 0
void primal_reset (primal_t *primal)
{
	instance_t *inst = primal->inst;
	solution_t *s = primal->sol;
	int i;

	primal->tabu_length = 10.0; // initial tabu tenure
	primal->n_moves_at_last_improvement = 0; // reset counter
	primal->n_moves = 0; // reset counter
	for (i = 0; i < inst->n; i++) {
		primal->guiding_x[i] = 1; // bogus guiding solution
		primal->improving_partial_x[i] = -1; // reset improv. part. sol
		primal->tabu[i] = 0.0; // nothing is tabu
	}
	
	solution_reset(inst, s, primal->best_x); // reset with old best-known solution (chances are, it's not too bad)
	primal->best_z = s->total_cost; // reset upper bound

	primal->state = primal_state_run; // set state as ready to run
}

// Updates tabu tenure after 1OPT move on location w
void tabu_tenure (primal_t *primal, int w, int improving_move)
{
	primal->tabu[w] = primal->tabu_length + ((double) primal->n_moves);

	if (improving_move) {
		if (primal->tabu_length > 2.0)
			--primal->tabu_length;
	} else if (primal->tabu_length < 10.0)
		++primal->tabu_length;

	++primal->n_moves;
}


// Performs one iteration of the primal process
void primal_run (primal_t *primal)
{
	instance_t *inst = primal->inst;
	solution_t *s = primal->sol;
	int i, w, n, aspiration, n_best, n_free, n_diver;
	double best_gain;

	check_solution(inst, s);

	//% Algorithm 1 Step 4 (cont'd)
	// Perform moves to conform with the improving partial solution
	for (i = 0; i < inst->n; i++) 
		if (primal->improving_partial_x[i] != -1 && primal->tabu[i] != INFINITY) {
			primal->tabu[i] = INFINITY;
			if (s->x[i] != primal->improving_partial_x[i]) {
				solution_flip(inst, s, i);
				check_solution(inst, s);
				++primal->n_moves;
				if (s->total_cost < primal->best_z) {
					primal->n_moves_at_last_improvement = primal->n_moves;
					primal->best_z = s->total_cost;
					memcpy(primal->best_x, s->x, inst->n * sizeof(int));
				}
			}
		}
	// Perturb tabu state using the guiding solution
	if (primal->n_moves_at_last_improvement + inst->request_period < primal->n_moves) {
		primal->n_moves_at_last_improvement = primal->n_moves;
		for (i = 0; i < inst->n; i++) 
			if (primal->tabu[i] != INFINITY && s->x[i] == primal->guiding_x[i])
				primal->tabu[i] = ((double) primal->n_moves) + primal->tabu_length; 
	}

	// Analyze search state
	//% Algorithm 1 step 1
	for (aspiration = 0, best_gain = -INFINITY, n_best = n_free = n_diver = 0, i = 0; i < inst->n; i++) 
		if (primal->tabu[i] == INFINITY) {
			assert(s->x[i] == primal->improving_partial_x[i]);
		} else {
			++n_free; // count the number of unfixed locations
			if ((primal->n_moves == 0 || (s->total_cost - s->flip_gain[i] <= primal->best_z - bound_eps)) && s->flip_gain[i] != INFINITY) 
				aspiration = 1; // aspiration criterion is satisfied
			if ((double) primal->n_moves > primal->tabu[i] || aspiration) {
				++n_diver; // count the number of non-tabu locations
				if (s->flip_gain[i] > best_gain + bound_eps && best_gain != INFINITY)
					n_best = 1, best_gain = s->flip_gain[i]; // store the best location to perform a move on 
				else if (fabs(s->flip_gain[i] - best_gain) < 2.0 * bound_eps || (s->flip_gain[i] == INFINITY && best_gain == INFINITY))
					++n_best; // count the number of best locations
			}
		}

	// Select the 1OPT move to perform
	w = -1; 
	if (best_gain <= bound_eps) { // if there exist no improving moves
		if (n_free == 0) { // if there exist no moves
			//% Algorithm 1 step 1(d)
			return;
		}
		if (n_diver == 0) { // if there exist no non-tabu moves
			// select one location at random
			//% Algorithm 1 step 1(c)
			n = RngStream_RandInt(primal->rng, 0, n_free - 1);
			for (i = 0; i < inst->n && w == -1; i++)
				if (primal->tabu[i] != INFINITY && n-- <= 0) 
					w = i;
		} else { // there exist non-tabu moves
			// select one location at random
			//% Algorithm 1 step 1(b)
			n = RngStream_RandInt(primal->rng, 0, n_diver - 1);
			for (i = 0; i < inst->n && w == -1; i++)
				if (((double) primal->n_moves > primal->tabu[i] || aspiration) && n-- <= 0) 
					w = i;
		}
	} else { // there exists at least one improving move
		// select one of the best at random
		//% Algorithm 1 step 1(a)
		assert(n_best > 0);
		n = RngStream_RandInt(primal->rng, 0, n_best - 1);
		for (w = -1, i = 0; i < inst->n && w == -1; i++) 
			if (((double) primal->n_moves > primal->tabu[i] || aspiration) 
					&& (fabs(s->flip_gain[i] - best_gain) < 2.0 * bound_eps || (s->flip_gain[i] == INFINITY && best_gain == INFINITY))
					&& n-- <= 0) 
				w = i;
	}

	// check that in all cases we have found something
	assert(w >= 0); 
	assert(w < inst->n);
	assert(primal->tabu[w] != INFINITY);

	// Update search state
	//% Algorithm 1 step 2
	tabu_tenure(primal, w, (s->flip_gain[w] > bound_eps)); 
	//% Algorithm 1 step 3(a,b,c)
	solution_flip(inst, s, w);
	check_solution(inst, s);
	//% Algorithm 1 step 3(d)
	if (s->total_cost < primal->best_z) {
		// Update best known solution
		primal->n_moves_at_last_improvement = primal->n_moves;
		primal->best_z = s->total_cost;
		memcpy(primal->best_x, s->x, inst->n * sizeof(int));
	}
}


primal_t *primal_create (instance_t *inst)
{
	primal_t *primal = (primal_t*) malloc(sizeof(primal_t));

	primal->inst = instance_copy(inst);
	primal->state = primal_state_undef;

	primal->guiding_x = (int*) calloc(inst->n, sizeof(int));
	primal->improving_partial_x = (int*) calloc(inst->n, sizeof(int));
	primal->best_x = (int*) calloc(inst->n, sizeof(int));
	primal->best_z = INFINITY;

	primal->sol = solution_alloc(inst);
	primal->tabu = (double*) calloc(inst->n, sizeof(double));
	primal->n_moves = 0;
	primal->n_moves_at_last_improvement = 0;

	primal->rng = RngStream_CreateStream("");

	return primal;
}


void primal_destroy (primal_t *primal)
{
	instance_destroy(primal->inst);
	RngStream_DeleteStream(primal->rng);
	free(primal->best_x);
	free(primal->improving_partial_x);
	free(primal->guiding_x);
	free(primal);
}




// ** Location heap operations **


void sift_down(double **c, int j, int *heap, int *heap_inv, int i, int size)
{
	int s, t;

	while (i * 2 + 1 < size) {
		s = i * 2 + 1;
		
		if (s+1 < size && c[heap[s]][j] > c[heap[s+1]][j])
			++s;
		if (c[heap[i]][j] <= c[heap[s]][j])
			return;

		heap_inv[heap[s]] = i;
		heap_inv[heap[i]] = s;
		t = heap[i];
		heap[i] = heap[s];
		heap[s] = t;
		i = s;
	}
}

void sift_up(double **c, int j, int *heap, int *heap_inv, int i)
{
	int p, t;
	while (i > 0) {
		p = (i - 1) / 2;
		if (c[heap[p]][j] <= c[heap[i]][j])
			return;
		heap_inv[heap[p]] = i;
		heap_inv[heap[i]] = p;
		t = heap[p];
		heap[p] = heap[i];
		heap[i] = t;
		i = p;
	}
}



solution_t *solution_alloc (instance_t *inst)
{
	solution_t *s = (solution_t*) malloc(sizeof(solution_t));
	int j;

	s->x = (int*) calloc(inst->n, sizeof(int));
	s->flip_gain = (double*) calloc(inst->n, sizeof(double));
	s->heap = (int**) calloc(inst->m, sizeof(int*));
	s->heap_inv = (int**) calloc(inst->m, sizeof(int*));
	s->besti1 = (int*) calloc(inst->m, sizeof(int));
	s->besti2 = (int*) calloc(inst->m, sizeof(int));

	for (j = 0; j < inst->m; j++) {
		// allocate j-th heap
		s->heap[j] = (int*) calloc(inst->n, sizeof(int));
		s->heap_inv[j] = (int*) calloc(inst->n, sizeof(int));
	}

	return s;
}


// Performs 1OPT move on location w in O(m log n) time
void solution_flip(instance_t *inst, solution_t *s, int w)
{
	int i, j, w2, i2;

	// Update costs
	s->total_cost -= s->flip_gain[w]; 
	s->flip_gain[w] *= -1.0;

	// Deal with special case where the solution was 'all locations closed'
	if (s->heap_size == 0) { 
		s->x[w] = 1;
		for (j = 0; j < inst->m; j++) {
			// the cheapest location in each heap is the only open location
			s->heap[j][0] = w;
			s->heap_inv[j][w] = 0;
			s->besti1[j] = w;
		}
		for (i = 0; i < inst->n; i++) 
			if (i != w) {
				// recompute cost differences for all remaining closed locations
				s->flip_gain[i] = -inst->f[i];
				for (j = 0; j < inst->m; j++) {
					if (inst->c[i][j] < inst->c[w][j])
						s->flip_gain[i] += inst->c[w][j] - inst->c[i][j];
				}
			}
		s->heap_size = 1;
		return;
	}

	// Deal with the special case where we close the only remaining open solution
	if (s->heap_size == 1 && s->x[w]) {
		s->x[w] = 0;
		for (j = 0; j < inst->m; j++) 
			s->besti1[j] = s->heap_inv[j][w] = s->heap[j][0] = -1; //remove w from the heap
		for (i = 0; i < inst->n; i++) 
			if (i != w) {
				// compute fictional cost differences
				s->flip_gain[i] = inst->ub - inst->f[i];
				for (j = 0; j < inst->m; j++) 
					s->flip_gain[i] -= inst->c[i][j];
			}
		s->heap_size = 0;
		return;
	}

	// Deal with the general cases
	if (s->x[w]) { // case where we close location w
		s->x[w] = 0;
		--s->heap_size;
		for (j = 0; j < inst->m; j++) {
			// carefully remove it from each heap
			i = s->heap_inv[j][w];
			assert(i != -1);
			assert(s->heap[j][i] == w);
			w2 = s->heap[j][s->heap_size];
			assert(s->heap_inv[j][w2] == s->heap_size);
			s->heap_inv[j][w] = -1;
			if (w != w2) {
				s->heap_inv[j][w2] = i;
				s->heap[j][i] = w2;
				sift_up(inst->c, j, s->heap[j], s->heap_inv[j], i);
				sift_down(inst->c, j, s->heap[j], s->heap_inv[j], s->heap_inv[j][w2], s->heap_size);
			}
		}
	} else { // case where we open location w
		s->x[w] = 1;
		for (j = 0; j < inst->m; j++) {
			// add it to each heap
			s->heap[j][s->heap_size] = w;
			s->heap_inv[j][w] = s->heap_size;
			sift_up(inst->c, j, s->heap[j], s->heap_inv[j], s->heap_size);
		}
		++s->heap_size;
	}

	// Update move costs of all closed locations
	for (j = 0; j < inst->m; j++) {
		if (inst->c[s->heap[j][0]][j] <= inst->c[s->besti1[j]][j]) { 
			// the client j can use a cheaper location than previously
			for (i = 0; i < inst->n; i++) 
				if (i != w && !s->x[i]) {
					if (inst->c[i][j] < inst->c[s->heap[j][0]][j]) 
						s->flip_gain[i] -= inst->c[s->besti1[j]][j] - inst->c[s->heap[j][0]][j];
					else if (inst->c[i][j] < inst->c[s->besti1[j]][j]) 
						s->flip_gain[i] -= inst->c[s->besti1[j]][j] - inst->c[i][j];
				}
		} else {
			for (i = 0; i < inst->n; i++) 
				if (i != w && !s->x[i]) {
					if (inst->c[i][j] < inst->c[s->besti1[j]][j]) 
						s->flip_gain[i] += inst->c[s->heap[j][0]][j] - inst->c[s->besti1[j]][j] ;
					else if (inst->c[i][j] < inst->c[s->heap[j][0]][j]) 
						s->flip_gain[i] += inst->c[s->heap[j][0]][j] - inst->c[i][j];
				}
		}
	}

	if (s->heap_size == 1) { // Special case where one location is open
		// This special case only occurs when we close one of two open locations
		// Recompute all move costs and two cheapest locations from scratch
		i = s->heap[0][0];
		assert(s->x[i]);
		s->flip_gain[i] = -inst->ub + inst->f[i];
		for (j = 0; j < inst->m; j++) {
			s->flip_gain[i] += inst->c[i][j];
			s->besti2[j] = -1;
			s->besti1[j] = s->heap[j][0];
		}
	} else if (s->heap_size == 2 && s->x[w]) { // Special case where two locations are open
		// This special case only occurs when only one location was open an we open another 
		// Recompute all move costs and two cheapest locations from scratch
		i = (s->heap[0][0] == w) ? s->heap[0][1] : s->heap[0][0];
		assert(s->x[i]);
		assert(i != w);
		s->flip_gain[i] = inst->f[i];
		for (j = 0; j < inst->m; j++) {
			assert(s->besti2[j] == -1);
			if (s->heap[j][0] == i)
				s->flip_gain[i] += inst->c[i][j] - inst->c[s->heap[j][1]][j];
			s->besti2[j] = s->heap[j][1];
			s->besti1[j] = s->heap[j][0];
		}
	} else {
		// Update all move costs and the two cheapest locations 
		for (j = 0; j < inst->m; j++) {
			i2 = (s->heap_size == 2 || inst->c[s->heap[j][1]][j] <= inst->c[s->heap[j][2]][j]) ? 1 : 2;
			if (s->heap[j][0] != s->besti1[j] || inst->c[s->heap[j][i2]][j] != inst->c[s->besti2[j]][j]) {
				i = s->besti1[j];
				if (i != w && s->x[i]) 
					s->flip_gain[i] += inst->c[s->besti2[j]][j] - inst->c[i][j];
				i = s->heap[j][0];
				if (i != w && s->x[i])
					s->flip_gain[i] -= inst->c[s->heap[j][i2]][j] - inst->c[i][j];
			}
			s->besti2[j] = s->heap[j][i2];
			s->besti1[j] = s->heap[j][0];
		}
	}
}



void check_solution(instance_t *inst, solution_t *s) 
{
	int i, j;
	double g, t;

	return; // Delete this line if assert checks are needed

	// check heaps
	assert(s->heap_size >= 0);
	assert(s->heap_size <= inst->n);
	for (j = 0; j < inst->m; j++) {
		for (i = 0; i < inst->n; i++) 
			if (s->heap_inv[j][i] == -1) {
				assert(!s->x[i]);
			} else {
				assert(s->heap_inv[j][i] >= 0);
				assert(s->heap_inv[j][i] < s->heap_size);
				assert(s->heap[j][s->heap_inv[j][i]] == i);
			}
		for (i = 0; i < s->heap_size; i++) {
			assert(s->x[s->heap[j][i]]);
			if (i * 2 + 1 < s->heap_size)
				assert(inst->c[s->heap[j][i]][j] <= inst->c[s->heap[j][i*2+1]][j]);
			if (i * 2 + 2 < s->heap_size)
				assert(inst->c[s->heap[j][i]][j] <= inst->c[s->heap[j][i*2+2]][j]);
		}

		assert(s->besti1[j] == -1 || s->besti1[j] == s->heap[j][0]);
		if (s->heap_size == 0) {
			assert(s->besti1[j] == s->besti2[j] && s->besti1[j] == -1);
		} else if (s->heap_size == 1) {
			assert(s->besti2[j] == -1);
		} else if (s->heap_size == 2) {
			assert(s->besti2[j] == s->heap[j][1]);
		} else {
			assert(s->besti2[j] == s->heap[j][1] || s->besti2[j] == s->heap[j][2]);
			if (s->besti2[j] == s->heap[j][2]) 
				assert(inst->c[s->heap[j][2]][j] <= inst->c[s->heap[j][1]][j]);
			else
				assert(inst->c[s->heap[j][1]][j] <= inst->c[s->heap[j][2]][j]);
		}
	}

	// check flip costs
	if (s->heap_size == 0) {
		for (i = 0; i < inst->n; i++) {
			assert(!s->x[i]);
			for (g = inst->ub - inst->f[i], j = 0; j < inst->m; j++) 
				g -= inst->c[i][j];
			assert(fabs(g - s->flip_gain[i]) < bound_eps);
		}
	} else {
		for (i = 0; i < inst->n; i++) 
			if (s->x[i] && s->heap_size == 1) {
				for (g = -inst->ub + inst->f[i], j = 0; j < inst->m; j++) 
					g += inst->c[i][j];
				assert(fabs(g - s->flip_gain[i]) < bound_eps);
			} else if (s->x[i]) {
				assert(s->heap_size >= 2);
				for (g = inst->f[i], j = 0; j < inst->m; j++) 
					if (s->besti1[j] == i)
						g -= inst->c[s->besti2[j]][j] - inst->c[i][j];
				assert(fabs(g - s->flip_gain[i]) < bound_eps);
			} else {
				assert(s->x[i] == 0 && s->heap_size >= 1);
				for (g = -inst->f[i], j = 0; j < inst->m; j++) 
					if (inst->c[s->besti1[j]][j] > inst->c[i][j])
						g += inst->c[s->besti1[j]][j] - inst->c[i][j];
				assert(fabs(g - s->flip_gain[i]) < bound_eps);
			}
	}

	// check solution cost
	if (s->heap_size == 0) {
		assert(fabs(s->total_cost - inst->ub) < bound_eps);
	} else {
		t = 0.0;
		for (i = 0; i < inst->n; i++)
			if (s->x[i])
				t += inst->f[i];
		for (j = 0; j < inst->m; j++) {
			for (g = 1e50, i = 0; i < inst->n; i++)
				if (s->x[i] && inst->c[i][j] < g)
					g = inst->c[i][j];
			t += g;
		}
		assert(fabs(t - s->total_cost) < bound_eps);
	}
}


// Rebuild the solution_t data structure using 'x' 
void solution_reset (instance_t *inst, solution_t *s, int *x)
{
	int i, j, k, l;

	// Empty all heaps 
	for (j = 0; j < inst->m; j++) {
		s->besti1[j] = s->besti2[j] = -1;
		for (i = 0; i < inst->n; i++) 
			s->heap[j][i] = s->heap_inv[j][i] = -1;
	}

	// Compute total opening costs and heap size
	s->heap_size = 0;
	s->total_cost = 0.0;
	if (NULL != x) {
		for (s->heap_size = i = 0; i < inst->n; i++) {
			s->x[i] = x[i];
			if (x[i] == 1) {
				++s->heap_size;
				s->total_cost += inst->f[i];
			}
		}
	}

	// Populate heaps, compute move costs and total solution cost
	if (NULL == x || s->heap_size == 0) { // Special case where all locations are closed
		s->total_cost = inst->ub;
		for (i = 0; i < inst->n; i++) {
			s->x[i] = 0;
			s->flip_gain[i] = s->total_cost - inst->f[i];
			for (j = 0; j < inst->m; j++) 
				s->flip_gain[i] -= inst->c[i][j];
		}
	} else { // General case
		// Populate heaps by simply adding all open locations in increasing service costs,
		// as this verifies the heap property.
		for (j = 0; j < inst->m; j++) {
			for (i = k = l = 0; k < inst->n && l < s->heap_size; k++) {
				i = inst->inc[j][k];
				if (x[i] == 1) {
					s->heap[j][l] = i;
					s->heap_inv[j][i] = l;
					if (l == 0) {
						s->besti1[j] = i;
						s->total_cost += inst->c[i][j];
					} else if (l == 1)
						s->besti2[j] = i;
					++l;
				}
			}
			assert(l == s->heap_size);
		}

		// Compute move costs
		for (i = 0; i < inst->n; i++) {
			if (s->x[i] && s->heap_size == 1) {
				for (s->flip_gain[i] = - inst->ub + inst->f[i], j = 0; j < inst->m; j++) 
					s->flip_gain[i] += inst->c[i][j];
			} else if (s->x[i]) {
				assert(s->heap_size >= 2);
				for (s->flip_gain[i] = inst->f[i], j = 0; j < inst->m; j++) 
					if (s->besti1[j] == i)
						s->flip_gain[i] -= inst->c[s->besti2[j]][j] - inst->c[i][j];
			} else {
				assert(s->x[i] == 0 && s->heap_size >= 1);
				for (s->flip_gain[i] = -inst->f[i], j = 0; j < inst->m; j++) 
					if (inst->c[s->besti1[j]][j] > inst->c[i][j])
						s->flip_gain[i] += inst->c[s->besti1[j]][j] - inst->c[i][j];
			}
		}
	}

	check_solution(inst, s);
}





