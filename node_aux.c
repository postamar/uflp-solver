// (C) 2011, Marius Posta (postamar@iro.umontreal.ca)
// Check LICENSE.txt for the legal blah-blah.


#include "node.h"


// ** Adjustable array ** 

adj_array_t *adj_array_create (int alloc_size_increment)
{
	int i;
	adj_array_t *adj_array = (adj_array_t*) malloc(sizeof(adj_array_t));

	adj_array->e = (void**) calloc(alloc_size_increment, sizeof(void*));
	adj_array->alloc_size = adj_array->alloc_size_increment = alloc_size_increment;
	adj_array->n = 0;
	for (i = 0; i < alloc_size_increment; i++)
		adj_array->e[i] = NULL;

	return adj_array;
}


void adj_array_destroy (adj_array_t *adj_array)
{
	free(adj_array->e);
	free(adj_array);
}


// Stack functions

void adj_array_push (adj_array_t *adj_array, void *e)
{
	int i;

	if (adj_array->n == adj_array->alloc_size) {
		adj_array->alloc_size += adj_array->alloc_size_increment;
		adj_array->e = (void**) realloc(adj_array->e, adj_array->alloc_size * sizeof(void*));
		for (i = adj_array->n; i < adj_array->alloc_size; i++)
			adj_array->e[i] = NULL;
	}
	for (i = adj_array->n; i < adj_array->alloc_size; i++)
		assert(adj_array->e[i] == NULL);


	assert(adj_array->e[adj_array->n] == NULL);
	adj_array->e[adj_array->n++] = e;
}


void *adj_array_pop (adj_array_t *adj_array)
{
	void *e = NULL;

	assert(adj_array->n > 0);
	
	e = adj_array->e[--adj_array->n];
	adj_array->e[adj_array->n] = NULL;

	return e;
}


// Heap functions

void heap_sift_down (adj_array_t *h_node, int root, enum node_selection_t strategy) 
{
	int child, swap, n = h_node->n;
	node_t *temp, **e = (node_t**) h_node->e;

	while (root * 2 + 1 < n) {
		assert(e[root] != NULL);
		child = root * 2 + 1;
		swap = root;
		switch (strategy) {
			case search_bf:
			if (e[child] != NULL && e[root]->lb > e[child]->lb)
				swap = child;
			if (child + 1 < n && e[child+1] != NULL && e[swap]->lb > e[child+1]->lb)
				swap = child + 1;
			break;
			case search_df:
			if (e[child] != NULL && e[root]->id < e[child]->id)
				swap = child;
			if (child + 1 < n && e[child+1] != NULL && e[swap]->id < e[child+1]->id)
				swap = child + 1;
			break;
			default:
			assert(0);
		}
		if (swap == root)
			return;

		temp = e[root];
		e[root] = e[swap];
		e[swap] = temp;
		root = swap;
	}
}


void heap_sift_up (adj_array_t *h_node, int current, enum node_selection_t strategy) 
{
	int parent;
	node_t *temp, **e = (node_t**) h_node->e;

	while (current > 0) {
		assert(e[current] != NULL);
		parent = lround(floor((current - 1) / 2));
		assert(e[parent] != NULL);

		switch (strategy) {
			case search_bf:
			if (e[parent]->lb <= e[current]->lb)
				return;
			break;
			case search_df:
			if (e[parent]->id > e[current]->id)
				return;
			break;
			default:
			assert(0);
		}

		temp = e[current];
		e[current] = e[parent];
		e[parent] = temp;
		current = parent;
	}
}



// ** Node definitions **

node_t *node_alloc (instance_t *inst)
{
	node_t *node = (node_t*) malloc(sizeof(node_t));

	node->left = node->right = NULL;
	node->fix_mask = (word_t*) calloc(inst->n_words_in_bit_vector, sizeof(word_t));
	node->fix_val = (word_t*) calloc(inst->n_words_in_bit_vector, sizeof(word_t));
	node->mul = (double*) calloc(inst->m, sizeof(double));
	node->lb = 0.0;
	node->ub = INFINITY;
	node->min_open = 1;
	node->max_open = inst->n;
	node->branch_i = -1;
	node->id = 0;

	return node;
}


void node_set (instance_t *inst, node_t *parent, node_t *node)
{
	node->branch_i = -1;

	node->left = node->right = NULL;
	node->min_open = parent->min_open;
	node->max_open = parent->max_open;
	node->lb = parent->lb;
	node->ub = parent->ub;
	memcpy(node->fix_mask, parent->fix_mask, inst->n_words_in_bit_vector * sizeof(word_t));
	memcpy(node->fix_val, parent->fix_val, inst->n_words_in_bit_vector * sizeof(word_t));
	memcpy(node->mul, parent->mul, inst->m * sizeof(double));
}

void node_create_children (instance_t *inst, node_t *parent)
{
	parent->left = node_alloc(inst);
	parent->right = node_alloc(inst);
	node_set(inst, parent, parent->left);
	node_set(inst, parent, parent->right);
}


void node_destroy (node_t *node)
{
	if (node->mul != NULL)
		free(node->mul);
	free(node->fix_val);
	free(node->fix_mask);
	free(node);
}






// Node evaluation data structures

node_eval_t *node_eval_create (instance_t *inst)
{
	int i;
	node_eval_t *node_eval = (node_eval_t*) malloc(sizeof(node_eval_t));

	node_eval->n_node_eval = node_eval->n_lag_eval = 0;
	node_eval->inst = inst;
	node_eval->rng = RngStream_CreateStream("");

	node_eval->bundle_cfg = (char*) calloc(2048, sizeof(char));
	node_eval->bundle = NULL;

	node_eval->x = (int*) calloc(inst->n, sizeof(int));
	node_eval->x1 = (int*) calloc(inst->n, sizeof(int));
	node_eval->rc = (double*) calloc(inst->n, sizeof(double));
	node_eval->int_subg = (int*) calloc(inst->m, sizeof(int));
	node_eval->min_c = (double*) calloc(inst->m, sizeof(double));
	node_eval->bundle_mul_alloc = (double*) calloc(inst->m, sizeof(double));
	node_eval->queue = (ic_pair_t*) calloc(inst->n, sizeof(ic_pair_t));

	node_eval->cbdata.callbacks.SetGi = bundle_set_subg;
	node_eval->cbdata.callbacks.Fi = bundle_eval;
	node_eval->cbdata.callbacks.GetGi = bundle_get_subg;
	node_eval->cbdata.callbacks.EveryIteration = NULL;
	node_eval->cbdata.node = NULL;
	node_eval->cbdata.node_eval = node_eval;

	node_eval->cache = (hyperplane_t*) calloc(inst->subg_cache_size, sizeof(hyperplane_t));
	for (i = 0; i < inst->subg_cache_size; i++) {
		node_eval->cache[i].subg = (double*) calloc(inst->m, sizeof(double));
		node_eval->cache[i].x = (word_t*) calloc(inst->n_words_in_bit_vector, sizeof(word_t));
		node_eval->cache[i].fixed_x = (word_t*) calloc(inst->n_words_in_bit_vector, sizeof(word_t));
	}
	node_eval->current_hyperplane.subg = (double*) calloc(inst->m, sizeof(double));
	node_eval->current_hyperplane.x = (word_t*) calloc(inst->n_words_in_bit_vector, sizeof(word_t));
	node_eval->current_hyperplane.fixed_x = (word_t*) calloc(inst->n_words_in_bit_vector, sizeof(word_t));
	node_eval->selection = (ic_pair_t*) calloc(inst->subg_cache_size, sizeof(ic_pair_t));

	node_eval->bloom_size = 635291;
	node_eval->n_hash_salt = 7;
	node_eval->hash_salt = (unsigned long*) calloc(node_eval->n_hash_salt, sizeof(unsigned long));
	for (i = 0; i < node_eval->n_hash_salt; i++)
		node_eval->hash_salt[i] = RngStream_RandInt(node_eval->rng, 0, node_eval->bloom_size - 1) - i;
	node_eval->bloom = (unsigned short*) calloc(node_eval->bloom_size, sizeof(unsigned short));


	node_eval->partition = (int*) calloc(inst->n, sizeof(int));
	node_eval->I1 = (int*) calloc(inst->n, sizeof(int));
	node_eval->Is = (int*) calloc(inst->n, sizeof(int));
	node_eval->I0 = (int*) calloc(inst->n, sizeof(int));

	return node_eval;
}


void node_eval_destroy (node_eval_t *node_eval)
{
	int i;

	free(node_eval->partition);
	free(node_eval->I1);
	free(node_eval->Is);
	free(node_eval->I0);

	free(node_eval->bloom);
	free(node_eval->hash_salt);

	free(node_eval->selection);
	free(node_eval->current_hyperplane.fixed_x);
	free(node_eval->current_hyperplane.x);
	free(node_eval->current_hyperplane.subg);
	for (i = 0; i < node_eval->inst->subg_cache_size; i++) {
		free(node_eval->cache[i].fixed_x);
		free(node_eval->cache[i].x);
		free(node_eval->cache[i].subg);
	}
	free(node_eval->cache);

	free(node_eval->queue);
	free(node_eval->bundle_mul_alloc);
	free(node_eval->min_c);
	free(node_eval->rc);
	free(node_eval->int_subg);
	free(node_eval->x1);
	free(node_eval->x);

	bundle_solver_destroy(node_eval->bundle);
	free(node_eval->bundle_cfg);
	RngStream_DeleteStream(node_eval->rng);
	free(node_eval);
}


// Returns a config string to initialize a BTT bundle object with
char* generate_config_string  (char *config_string, int bundle_max_size, int bundle_max_iter, double approx_z_opt)
{
    int inactive_delay = bundle_max_size;
    int bundle_size = bundle_max_size;

    float t_enlarge_factor = 10;
    float t_reduce_factor = 0.1;
    float serious_step_factor = 0.2;
    float medium_step_factor = 0.3;
    float null_step_factor = 2.9;

    float t_optimal = 50.0;
    float eps_lin = 0.02 / approx_z_opt;

    float t_max = 1e5;
    float t_min = 1e-5;
    float t_init = 0.1;

    float t_strategy_1 = 0.0;
    float t_strategy_2 = 0.0;
    float t_strategy_eps = 0.0;

    int pricing_warmup = 0;
    int pricing_period = 0;
    int pricing_max_age = 0;

	sprintf(config_string, "%i\n%i\n%i\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%i\n%i\n%i\n",
        bundle_max_iter,
        inactive_delay,
        bundle_size,
        t_enlarge_factor,
        t_reduce_factor,
        serious_step_factor,
        medium_step_factor,
        null_step_factor,
        t_optimal,
        eps_lin,
        t_max,
        t_min,
        t_init,
        t_strategy_1,
        t_strategy_2,
        t_strategy_eps,
        pricing_warmup,
        pricing_period,
        pricing_max_age);

	return config_string;
}


void node_eval_reset (node_eval_t *node_eval, double approx_z)
{
	instance_t *inst = node_eval->inst;
	int i;

	// Reset bundle object
	generate_config_string(node_eval->bundle_cfg, inst->bundle_max_size, inst->n_root_node_max_iter, approx_z);
	if (node_eval->bundle != NULL)
		bundle_solver_destroy(node_eval->bundle);
	node_eval->bundle = bundle_solver_create(node_eval->bundle_cfg, inst->m);
	bundle_set_uc_range(node_eval->bundle, 0, inst->m); // tell BTT that multipliers are not >= 0

	// Reset subgradient cache
	node_eval->node_counter = 0;
	node_eval->fresh_idx = 0;
	for (i = 0; i < inst->subg_cache_size; i++) 
		node_eval->cache[i].eps0 = INFINITY;

	// Reset Bloom filter
	for (i = 0; i < node_eval->bloom_size; i++)
		node_eval->bloom[i] = 0;
}



// ** Bloom filter functions ** 


unsigned long hash (unsigned long *key, int n_key, unsigned long hash0)
{
	int i;

	for (i = 0; i < n_key; i++)
        hash0 = key[i] + (hash0 << 6) + (hash0 << 16) - hash0;

	return hash0;
}

unsigned long get_hash_x (word_t *x, int n_words)
{
	return hash((unsigned long *) x, n_words * sizeof(word_t) / sizeof(unsigned long), 0);
}

unsigned long get_hash_subg (int *subg, int m)
{
	int off = m % (sizeof(unsigned long) / sizeof(int));
	unsigned long hash0;
	memcpy(&hash0, subg, off * sizeof(int));
	return hash((unsigned long*) &subg[off], (int) (m * sizeof(int) / sizeof(unsigned long)), hash0); 
}

unsigned long hyperplane_hash (hyperplane_t h, unsigned long hash0) 
{
	unsigned long hh = hash0;

	hh = hash((unsigned long*) &h.min_open, 1, hh);
	hh = hash((unsigned long*) &h.max_open, 1, hh);
	hh = hash((unsigned long*) &h.n_x1, 1, hh);
	hh = hash(&h.x_hash, 1, hh);
	hh = hash(&h.fixed_x_hash, 1, hh);
	hh = hash(&h.subg_hash, 1, hh);

	return hh;
}

// Returns 0 if current subgradient is certainly not already in Bloom filter.
int check_bloom_filter (node_eval_t *node_eval)
{
	unsigned long hh;
	int k;

	for (k = 0; k < node_eval->n_hash_salt; k++) {
		hh = hyperplane_hash(node_eval->current_hyperplane, node_eval->hash_salt[k]);
		if (0 == node_eval->bloom[hh % node_eval->bloom_size])
			return 0;
	}
	return 1;
}

// Adds (val == +1) or removes (val == -1) current subgradient to/from Bloom filter
void add_bloom_filter (node_eval_t *node_eval, short val)
{
	unsigned long hh;
	int k;

	for (k = 0; k < node_eval->n_hash_salt; k++) {
		hh = hyperplane_hash(node_eval->current_hyperplane, node_eval->hash_salt[k]);
		node_eval->bloom[hh % node_eval->bloom_size] += val;
	}
}



// SSE intrinsics for fast floating-point arithmetic
//
// (C) Paul-Virak Khuong, 2011 (pvk@pvk.ca)
// Thanks Paul!
#include <xmmintrin.h>

typedef double v2df __attribute__ ((vector_size (16)));
typedef float  v4sf __attribute__ ((vector_size (16)));
typedef size_t v2ui __attribute__ ((vector_size (16)));
typedef int    v4si __attribute__ ((vector_size (16)));

// sum min{x_i - y_i, 0}
double clamped_sum (double * x, double * y, size_t count)
{
        v2df * xx = (v2df *)x,
              * yy = (v2df *)y;

        size_t vector_count = count/2;

        const v2df zero = {0.0, 0.0};
        v2df acc = zero;
        for (size_t i = 0; i < vector_count; i++)
                acc += __builtin_ia32_minpd((xx[i] - yy[i]), zero);
        
        union { double scalar[2]; v2df vector; } final;
		final.vector = acc;

        double partial = final.scalar[0] + final.scalar[1];
        if (count & 1) {
                size_t last = count-1;
                double diff = x[last] - y[last];
                return partial + ((diff < 0) ? diff : 0);
        } else 
			return partial;
}

void initialize_count (int * z, size_t count, int value)
{
        v4si * zz = (v4si *)z;
        v4si init = {value, value, value, value};
        size_t vector_count = count/4;
        for (size_t i = 0; i < vector_count; i++)
                zz[i] = init;
 
        size_t last = count & ~3UL;
        switch (count & 3) {
        case 3:
                z[last+2] = value;
        case 2:
                z[last+1] = value;
        case 1:
                z[last] = value;
        case 0:;
        }
}

// z_i += x_i < y_i
void clamped_count (int * z, double * x, double * y, size_t count)
{
        v2df * xx = (v2df *)x,
              * yy = (v2df *)y;
        v4si * zz = (v4si *)z;
        size_t vector_count = count/4;

        for (size_t i = 0; i < vector_count; i++) {
                v2df d1 = (v2df)_mm_cmplt_pd(xx[2*i], yy[2*i]);
                v2df d2 = (v2df)_mm_cmplt_pd(xx[2*i+1], yy[2*i+1]);
                // [d1_0 d1_2 d2_0, d2_2]
                v4si d = (v4si)_mm_shuffle_ps((v4sf)d1, (v4sf)d2, 0|2<<2|0<<4|2<<6);
                                                     
                zz[i] += d;
        }
        
        size_t last = count & ~3UL;
        switch (count & 3) {
        case 3:
                z[last+2] -= (x[last+2] < y[last+2]) ? 1 : 0;
        case 2:
                z[last+1] -= (x[last+1] < y[last+1]) ? 1 : 0;
        case 1:
                z[last]   -= (x[last] < y[last]) ? 1 : 0;
        case 0:;
        }
}
