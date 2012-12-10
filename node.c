// (C) 2011, Marius Posta (postamar@iro.umontreal.ca)
// Check LICENSE.txt for the legal blah-blah.

#include "node.h"


// BTT callback to evaluate the lagrangian dual at the given trial point
double bundle_eval (struct bundle_callbacks *callbacks_ptr, const double *Lam_val, const unsigned *Lam_idx, unsigned Lam_count)
{
	bundle_callback_data_t *callback_data = (bundle_callback_data_t*) callbacks_ptr;
	node_t *node = callback_data->node;
	node_eval_t *node_eval = callback_data->node_eval;
	instance_t *inst = node_eval->inst;
	unsigned k;

	// get multipliers and store them in dense vector format
	if (Lam_idx == NULL) {
		assert(Lam_count == (unsigned) inst->m);
		node_eval->bundle_mul = (double*) Lam_val;
	} else {
		node_eval->bundle_mul = node_eval->bundle_mul_alloc;
		memset(node_eval->bundle_mul, 0, inst->m * sizeof(double));
		for (k = 0; k < (unsigned) Lam_count; k++)
			node_eval->bundle_mul[Lam_idx[k]] = Lam_val[k];
	}

	// evaluate lagrangian relaxation
	++node_eval->bundle_n_fi;
	node_eval_lagrangian(node_eval, node, node_eval->bundle_mul);

	// store best multipliers
	if (node_eval->z > node->lb) {
		memcpy(node->mul, node_eval->bundle_mul, inst->m * sizeof(double));
		node->lb = node_eval->z;
	}

	// cutoff tests
	if (node->lb >= node->ub - bound_eps) 
		return -DBL_MAX;
	if (node_eval->bundle_n_gi > node_eval->bundle_max_gi)
		return -DBL_MAX;

	return -node_eval->z; // return opposite because BTT minimizes
}


// BTT callback to generate one or several subgradients at the current trial point
int bundle_get_subg (struct bundle_callbacks *callbacks_ptr, double *Eps, const unsigned **Indices)
{
	bundle_callback_data_t *callback_data = (bundle_callback_data_t*) callbacks_ptr;
	node_eval_t *node_eval = callback_data->node_eval;
	instance_t *inst = node_eval->inst;
	node_t *node = callback_data->node;
	hyperplane_t htemp, *h;
	ic_pair_t temp;
	int i, j, k;
	int queue_free_begin = node->min_open - node_eval->n_fixed_1;
	int queue_free_end = node->max_open - node_eval->n_fixed_1;

	++node_eval->bundle_n_gi; // increment bundle update counter

	if (node_eval->bundle_n_fi == 1 && inst->subg_cache_size) {
		// This is the first trial point and we are using a cache.
		if (node_eval->bundle_n_gi == 1) { 
			// Begin by retrieving all other subgradients from the cache.
			select_hyperplanes_in_cache(node_eval, node); 
		} else {
			// Return one of the many subgradients retrieved from the cache.
			assert(node_eval->n_selected > 0);
			h = &node_eval->cache[node_eval->selection[--node_eval->n_selected].i];
			*Eps = node_eval->selection[node_eval->n_selected].c + bound_eps;
			memcpy(node_eval->bundle_subg, h->subg, inst->m * sizeof(double));
			return -node_eval->n_selected; // 0 means no more, < 0 means call this callback again
		}
	} else 
		node_eval->n_selected = 0;

	*Eps = 0.0;
	*Indices = NULL;

	// Compute the subgradient at the current multipliers
	h = &node_eval->current_hyperplane;
	initialize_count(node_eval->int_subg, inst->m, 1); 
	for (k = 0; k < node_eval->n_x1; k++) 
		clamped_count(node_eval->int_subg, inst->c[node_eval->x1[k]], node_eval->bundle_mul, inst->m);
	for (j = 0; j < inst->m; j++) 
		node_eval->bundle_subg[j] = h->subg[j] = -node_eval->int_subg[j];

	// Store it in the cache, if we maintain one.
	if (inst->subg_cache_size) {
		h->eps0 = node_eval->z + ddot(inst->m, h->subg, node_eval->bundle_mul);
		h->subg_hash = get_hash_subg(node_eval->int_subg, inst->m);

		// First, compute the 'naive' key,
		h->min_open = node->min_open;
		h->max_open = node->max_open;
		h->n_x1 = node_eval->n_x1;
		for (k = 0; k < inst->n_words_in_bit_vector; k++) 
			h->x[k] = 0;
		for (k = 0; k < node_eval->n_x1; k++) 
			set_bit(h->x, node_eval->x1[k], 1);
		h->x_hash = get_hash_x(h->x, inst->n_words_in_bit_vector);
		memcpy(h->fixed_x, node->fix_mask, inst->n_words_in_bit_vector * sizeof(word_t));
		// then improve it by relaxing the fixed locations
		for (i = 0; i < inst->n; i++) 
			if (get_bit(node->fix_mask, i)) {
				temp.c = node_eval->rc[i];
				temp.i = i;
				if ((get_bit(node->fix_val, i) 
							&& -1 == ic_pair_cmp(&temp, &node_eval->queue[queue_free_begin]))
						|| (!get_bit(node->fix_val, i) 
							&& 1 == ic_pair_cmp(&temp, &node_eval->queue[queue_free_end]))) 
					set_bit(h->fixed_x, i, 0);
			}
		h->fixed_x_hash = get_hash_x(h->fixed_x, inst->n_words_in_bit_vector);

		// Check Bloom filter 
		if (0 == check_bloom_filter(node_eval)) {
			// Add to Bloom filter
			add_bloom_filter(node_eval, 1);
			if (node_eval->fresh_idx == inst->subg_cache_size - 1)
				node_eval->fresh_idx = 0;
			// Add to cache
			htemp = node_eval->cache[node_eval->fresh_idx];
			node_eval->cache[node_eval->fresh_idx++] = node_eval->current_hyperplane;
			node_eval->current_hyperplane = htemp;
			// Update Bloom filter if element was removed from cache
			if (node_eval->current_hyperplane.eps0 != INFINITY)
				add_bloom_filter(node_eval, -1);
		}
	}

	// 0 means no more subgradients available, < 0 means call this callback again
	return -node_eval->n_selected; 
}


// BTT callback which provides the address to which we should write the current subgradient to.
void bundle_set_subg (struct bundle_callbacks *callbacks_ptr, double *subgradient, unsigned u)
{
	bundle_callback_data_t *callback_data = (bundle_callback_data_t*) callbacks_ptr;

	(void) u;
	callback_data->node_eval->bundle_subg = subgradient;
}


// Retrieves a selection of valid subgradients from the cache
void select_hyperplanes_in_cache (node_eval_t *node_eval, node_t *node) 
{
	instance_t *inst = node_eval->inst;
	hyperplane_t *h;
	int i, j, k;
	double eps = 0.0;

	node_eval->n_selected = 0;

	// Visit the whole cache
	for (k = 0, i = node_eval->fresh_idx - 1; i != node_eval->fresh_idx; --i, k++) {
		assert(k < inst->subg_cache_size);
		if (i == -1) {
			i = inst->subg_cache_size - 1;
			if (i == node_eval->fresh_idx) 
				break;
			else
				continue;
		}

		h = &node_eval->cache[i];
		if (h->eps0 == INFINITY) 
			break;

		// Perform tests on the key
		if (h->min_open <= node->min_open 
				&& node->min_open <= h->n_x1 
				&& h->n_x1 <= node->max_open
				&& node->max_open <= h->max_open) {
			for (j = 0; j < inst->n_words_in_bit_vector; j++) {
				if (((h->x[j] ^ node->fix_val[j]) & node->fix_mask[j])
						|| (node->fix_mask[j] | h->fixed_x[j]) != node->fix_mask[j])
					break;
			}
			if (j == inst->n_words_in_bit_vector) {
				// Tests succeeded
				eps = h->eps0 - node_eval->z - ddot(inst->m, h->subg, node_eval->bundle_mul);
				assert(eps > -bound_eps);
				assert(node_eval->n_selected >= 0);
				assert(node_eval->n_selected < inst->subg_cache_size);

				node_eval->selection[node_eval->n_selected].i = i;
				node_eval->selection[node_eval->n_selected++].c = eps;
			}
		}
	}

	// Select only the best elements if we have too many
	if (node_eval->n_selected >= inst->bundle_max_size) {
		quickselect(node_eval->rng, node_eval->selection, node_eval->n_selected, inst->bundle_max_size - 1);
		node_eval->n_selected = inst->bundle_max_size - 1;
	}
}


// Applies the bundle search to optimize the multipliers for the node,
// returns 1 if a new improving primal solution was found, 0 otherwise.
int node_eval_improve_mul (node_eval_t *node_eval, node_t *node, int max_subg, double ub)
{
	instance_t *inst = node_eval->inst;
	double z = 0.0;
	int i, j, k;

	++node_eval->n_node_eval;

	// bundle-search for the best multipliers
	node->lb = 0.0;
	node->ub = ub;
	node_eval->n_selected = node_eval->bundle_n_fi = node_eval->bundle_n_gi = 0;
	node_eval->bundle_max_gi = max_subg;
	bundle_set_lambda(node_eval->bundle, node->mul);	

	node_eval->cbdata.node = node;
	bundle_solve_with_callbacks(node_eval->bundle, (struct bundle_callbacks*) &node_eval->cbdata);

	// compute primal cost
	node_eval_lagrangian(node_eval, node, node->mul);
//	assert(fabs(node_eval->z - node->lb) < bound_eps);

	for (j = 0; j < inst->m; j++)
		node_eval->min_c[j] = INFINITY;

	for (k = 0; k < node_eval->n_x1; k++) {
		i = node_eval->x1[k];
		z += inst->f[i];
		for (j = 0; j < inst->m; j++) 
			if (inst->c[i][j] < node_eval->min_c[j]) 
				node_eval->min_c[j] = inst->c[i][j]; 
	}

	for (j = 0; j < inst->m; j++)
		z += node_eval->min_c[j];

	if (z < node->ub) {
		node->ub = z;
		return 1;
	}

	return 0;
}


// Computes reduced costs and identifies unfixed locations in lagragian relaxation
//% Subsection 3.1
int node_eval_lagrangian_fixed (node_eval_t *node_eval, node_t *node, double *mul)
{
	instance_t *inst = node_eval->inst;
	int i, j, nfree, wi, oi;
	word_t w;

	++node_eval->n_lag_eval;

	node_eval->z = 0.0;
	node_eval->n_fixed_0 = node_eval->n_fixed_1 = 0;

	for (j = 0; j < inst->m; j++) 
		node_eval->z += mul[j];

	for (nfree = i = wi = oi = 0; i < inst->n; i++, oi++) {
		node_eval->rc[i] = inst->f[i] + clamped_sum(inst->c[i], mul, inst->m);

		if (oi == 8 * sizeof(word_t)) 
			oi = 0, ++wi;
		w = ((word_t) 1) << oi;
		if (node->fix_mask[wi] & w) {
			if ((node->fix_val[wi] & node->fix_mask[wi]) & w) {
				node_eval->x[i] = 1;
				node_eval->x1[node_eval->n_fixed_1++] = i;
				node_eval->z += node_eval->rc[i];
			} else {
				node_eval->x[i] = 0;
				++node_eval->n_fixed_0;
			}
		} else {
			node_eval->x[i] = -1;
			node_eval->queue[nfree].i = i;
			node_eval->queue[nfree++].c = node_eval->rc[i];
		}
	}
	
	node_eval->n_x1 = node_eval->n_fixed_1;
	assert(nfree + node_eval->n_fixed_0 + node_eval->n_fixed_1 == inst->n);

	if ((node_eval->n_fixed_1 > node->max_open) 
			|| (node->min_open + node_eval->n_fixed_0 > inst->n)) {
		node_eval->z = INFINITY;
		return 0;
	}
	return 1;
}


// Solves the lagrangian relaxation
//% Subsection 3.1
double node_eval_lagrangian (node_eval_t *node_eval, node_t *node, double *mul)
{
	instance_t *inst = node_eval->inst;
	int ni1, ni0, nfree;
	int i, k;

	if (!node_eval_lagrangian_fixed(node_eval, node, mul))
		return node_eval->z; // INFINITY

	nfree = inst->n - (node_eval->n_fixed_1 + node_eval->n_fixed_0); // number of unfixed locations

	// Identify I1
	ni1 = node->min_open - node_eval->n_fixed_1;
	if (ni1 > 0) {
		quickselect(node_eval->rng, node_eval->queue, nfree, ni1);
		for (k = 0; k < ni1; k++) {
			i = node_eval->queue[k].i;
			node_eval->x[i] = 1;
			node_eval->x1[node_eval->n_x1++] = i;
			node_eval->z += node_eval->rc[i];
		}
	} else 
		ni1 = 0;

	// Identify I0
	ni0 = inst->n - node->max_open - node_eval->n_fixed_0;
	if (ni0 > 0) {
		quickselect(node_eval->rng, &node_eval->queue[ni1], nfree - ni1, nfree - ni1 - ni0);
		for (k = nfree - 1; k >= nfree - ni0; --k) {
			i = node_eval->queue[k].i;
			node_eval->x[i] = 0;
		}
	} else
		ni0 = 0;

	// Find optimal states of locations in I*
	for (k = ni1; k < nfree - ni0; k++) {
		i = node_eval->queue[k].i;
		node_eval->x[i] = 0;
		if (node_eval->rc[i] < 0.0) {
			node_eval->x[i] = 1;
			node_eval->x1[node_eval->n_x1++] = i;
			node_eval->z += node_eval->rc[i];
		}
	}

	return node_eval->z;
}



