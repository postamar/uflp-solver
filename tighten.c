// (C) 2011, Marius Posta (postamar@iro.umontreal.ca)
// Check LICENSE.txt for the legal blah-blah.

#include "dual.h"


// ** Node separation functions used in implicit branching ** 

node_t *dual_node_increase_min_open (dual_t *dual, node_t *node, int n)
{
	(void) dual;
	node->min_open = node->min_open + n;

	return node;
}

node_t *dual_node_decrease_max_open (dual_t *dual, node_t *node, int n)
{
	(void) dual;
	node->max_open = node->max_open - n;

	return node;
}


node_t *dual_node_fix_open (dual_t *dual, node_t *node, int k)
{
	set_bit(node->fix_mask, k, 1);
	set_bit(node->fix_val, k, 1);
	++dual->node_eval->n_fixed_1;

	if (dual->node_eval->n_fixed_1 > node->min_open) {
		assert(dual->node_eval->n_fixed_1 == node->min_open + 1);
		return dual_node_increase_min_open(dual, node, 1);
	}
	return node;
}

node_t *dual_node_fix_closed (dual_t *dual, node_t *node, int k)
{
	set_bit(node->fix_mask, k, 1);
	set_bit(node->fix_val, k, 0);
	++dual->node_eval->n_fixed_0;

	if (dual->node_eval->n_fixed_0 > dual->inst->n - node->max_open) {
		assert(dual->node_eval->n_fixed_0 == dual->inst->n - (node->max_open + 1));
		return dual_node_decrease_max_open(dual, node, 1);
	}
	return node;
}


// ** Reoptimization functions **


// Returns the value of the lagrangian dual if location 'i' in I* was closed
double dual_bound_fix0Is (node_eval_t *node_eval, int i, int min0)
	// min0 is the cheapest location in I0, in terms of reduced costs
{
	double *rc = node_eval->rc;
	double z = node_eval->z - ((rc[i] < 0) ? rc[i] : 0);

	if (min0 != -1 && rc[min0] < 0)
		z += rc[min0];
	return z;

}

// Returns the value of the lagrangian dual if location 'i' in I* was opened
double dual_bound_fix1Is (node_eval_t *node_eval, int i, int max1)
	// max1 is the dearest location in I1, in terms of reduced costs
{
	double *rc = node_eval->rc;
	double z = node_eval->z + ((rc[i] > 0) ? rc[i] : 0);

	if (max1 != -1 && rc[max1] > 0)
		z -= rc[max1];
	return z;
}

// Returns the value of the lagrangian dual if location 'i' in I1 was closed 
double dual_bound_fix0I1 (node_eval_t *node_eval, int i, int min0, int mins)
	// mins is the cheapest location in I*, in terms of reduced costs
{
	double *rc = node_eval->rc;
	double z = node_eval->z - rc[i];

	if (min0 == -1 && mins == -1)
		return INFINITY;

	if (mins == -1) {
		z += rc[min0];
	} else {
		z += (rc[mins] > 0) ? rc[mins] : 0;
		if (min0 != -1) 
			z += (rc[min0] < 0) ? rc[min0] : 0;
	}
	return z;
}

// Returns the value of the lagrangian dual if location 'i' in I0 was opened 
double dual_bound_fix1I0 (node_eval_t *node_eval, int i, int max1, int maxs)
	// maxs is the dearest location in I*, in terms of reduced costs
{
	double *rc = node_eval->rc;
	double z = node_eval->z + rc[i];

	if (max1 == -1 && maxs == -1)
		return INFINITY;

	if (maxs == -1) {
		z -= rc[max1];
	} else {
		z -= (rc[maxs] < 0) ? rc[maxs] : 0;
		if (max1 != -1) 
			z -= (rc[max1] > 0) ? rc[max1] : 0;
	}
	return z;
}

// Implementation of implicit branching
// It additionally selects the branching location for [case otherwise] in Algorithm 4
//% Algorithms 3 and 4 - Subsection 3.3  
node_t *dual_node_tighten_bounds (dual_t *dual, node_t *node)
{
	node_eval_t *node_eval = dual->node_eval;
	instance_t *inst = node_eval->inst;
	ic_pair_t *queue = node_eval->queue;
	int *lrp = node_eval->partition; // array mapping each location to the part it is in (-1 for F0 or F1, 0 for I0, 1 for I1, 2 for I*)
	int *I1 = node_eval->I1;
	int *Is = node_eval->Is;
	int *I0 = node_eval->I0;
	double *rc = node_eval->rc; // reduced costs
	int nI0, nI1, nIs, nu; // respective sizes of I0, I1, I* and their union
	int i, k;
	double ub = dual->best_z - bound_eps;
	int min0, max1, mins, maxs;
	int min0k, max1k, minsk, maxsk; // indices of the above in I0, I1, I*, I*, respectively
	double z, zmaxmax, zminmin, z0, z1, rci;

	// Initialization
	node->branch_i = -1;
	assert(node->min_open <= node->max_open);
	assert(node->min_open >= node_eval->n_fixed_1);
	assert(node->max_open <= inst->n - node_eval->n_fixed_0);

	// Compute part sizes
	nu = inst->n - (node_eval->n_fixed_0 + node_eval->n_fixed_1);
	nI1 = node->min_open - node_eval->n_fixed_1;
	nI0 = (inst->n - node->max_open) - node_eval->n_fixed_0;
	nIs = nu - (nI1 + nI0);

	// Generate partition using node_eval->queue and associated data
	for (i = 0; i < inst->n; i++)
		lrp[i] = -1;
	min0k = max1k = minsk = maxsk = -1;
	min0 = max1 = mins = maxs = -1;
	for (k = 0; k < nI1; k++) {
		// the first nI1 elements in the queue form I1
		I1[k] = i = queue[k].i;
		lrp[i] = 1;
	}
	assert(k == nI1);
	for (; k < nI1 + nIs; k++) {
		// the next nIs elements form I*
		Is[k - nI1] = i = queue[k].i;
		lrp[i] = 2;
	}
	assert(k == nI1 + nIs);
	for (; k < nu; k++) {
		// the remaining nI0 elements form I0
		I0[k - (nI1 + nIs)] = i = queue[k].i;
		lrp[i] = 0;
	}


	// Begin tightening the bounds
label_begin:
	if (inst->fix) {
		if (node == NULL)
			return NULL;
		if (node_eval->z > ub) {
			node_destroy(node);
			return NULL;
		}

		// Recompute max1, mins, maxs, max0 if necessary
		if (nI1 > 0) {
			if (max1 == -1) {
				max1 = I1[0], max1k = 0;
				for (k = 1; k < nI1; k++)
					if (rc[I1[k]] > rc[max1])
						max1 = I1[k], max1k = k;
			} else {
				assert(lrp[max1] == 1), assert(max1k >= 0), assert(max1k < nI1), assert(I1[max1k] == max1);
			}
		} else {
			assert(max1 == -1), assert(max1k == -1);
		}

		if (nIs > 0) {
			if (mins == -1) {
				mins = Is[0], minsk = 0;
				for (k = 1; k < nIs; k++)
					if (rc[Is[k]] < rc[mins])
						mins = Is[k], minsk = k;
			} else {
				assert(lrp[mins] == 2), assert(minsk >= 0), assert(minsk < nIs), assert(Is[minsk] == mins);
			}
			if (maxs == -1) {
				maxs = Is[0], maxsk = 0;
				for (k = 1; k < nIs; k++)
					if (rc[Is[k]] > rc[maxs])
						maxs = Is[k], maxsk = k;
			} else {
				assert(lrp[maxs] == 2), assert(maxsk >= 0), assert(maxsk < nIs), assert(Is[maxsk] == maxs);
			}
		} else {
			assert(mins == -1), assert(minsk == -1), assert(maxs == -1), assert(maxsk == -1);
		}

		if (nI0 > 0) {
			if (min0 == -1) {
				min0 = I0[0], min0k = 0;
				for (k = 1; k < nI0; k++)
					if (rc[I0[k]] < rc[min0])
						min0 = I0[k], min0k = k;
			} else {
				assert(lrp[min0] == 0), assert(min0k >= 0), assert(min0k < nI0), assert(I0[min0k] == min0);
			}
		} else {
			assert(min0 == -1), assert(min0k == -1);
		}


		// Compute lower bounds on potential subnodes where max_open is set to min_open and vice versa
		zmaxmax = zminmin = node_eval->z;
		assert(node_eval->z <= ub);
		for (k = 0; k < nIs; k++) {
			assert(lrp[Is[k]] == 2);
			rci = rc[Is[k]];
			if (rci < 0)
				zminmin -= rci; 
			else
				zmaxmax += rci;
		}

		if (zminmin > ub) { 
			// Setting max_open to min_open would cause the lower bound to increase too much,
			// therefore min_open can be incremented by 1.
			assert(nIs > 0), assert(nI1 == 0 || rc[mins] >= rc[max1]);
			// Move mins from I* to I1
			lrp[mins] = 1;
			max1k = nI1;
			I1[nI1++] = max1 = mins;
			Is[minsk] = Is[--nIs];
			mins = maxs = minsk = maxsk = -1;
			// Reoptimize the lagrangian relaxation
			node_eval->x[max1] = 1;
			node_eval->z += (rc[max1] > 0) ? rc[max1] : 0;
			// Try again?
			goto label_begin;
		}
		if (zmaxmax > ub) {
			// Setting min_open to max_open would cause the lower bound to increase too much,
			// therefore max_open can be decremented by 1.
			assert(nIs > 0), assert(nI0 == 0 || rc[maxs] <= rc[min0]);
			// Move maxs from I* to I0
			lrp[maxs] = 0;
			min0k = nI0;
			I0[nI0++] = min0 = maxs;
			Is[maxsk] = Is[--nIs];
			mins = maxs = minsk = maxsk = -1;
			// Reoptimize the lagrangian relaxation
			node_eval->x[min0] = 0;
			node_eval->z -= (rc[min0] < 0) ? rc[min0] : 0;
			// Try again?
			goto label_begin;
		}

		// Actually update min_open and max_open
		k = node_eval->n_fixed_1 + nI1 - node->min_open;
		if (k > 0) {
			node = dual_node_increase_min_open(dual, node, k);
			if (node == NULL)
				return NULL;
		}
		k = node_eval->n_fixed_0 + nI0 - (inst->n - node->max_open);
		if (k > 0) {
			node = dual_node_decrease_max_open(dual, node, k);
			if (node == NULL)
				return NULL;
		}
		assert(node->min_open >= node_eval->n_fixed_1);
		assert(node->max_open <= inst->n - node_eval->n_fixed_0);

		if (dual->phase == 2) { 
			// Not in preprocessing

			// Attempt to fix locations in I1 to remain open
			for (k = 0; k < nI1; k++) {
				i = I1[k];
				assert(lrp[i] == 1);
				if (dual_bound_fix0I1(node_eval, i, min0, mins) > ub) {
					// We can fix location i as opened
					lrp[i] = -1; 
					I1[k] = I1[--nI1]; --nu;
					max1 = max1k = -1;
					node = dual_node_fix_open(dual, node, i);
					// Try again?
					goto label_begin;
				}
			}

			// Attempt to fix locations in I0 to remain closed
			for (k = 0; k < nI0; k++) {
				i = I0[k];
				assert(lrp[i] == 0);
				if (dual_bound_fix1I0(node_eval, i, max1, maxs) > ub) {
					// We can fix location i as closed
					lrp[i] = -1; 
					I0[k] = I0[--nI0]; --nu;
					min0 = min0k = -1;
					node = dual_node_fix_closed(dual, node, i);
					// Try again?
					goto label_begin;
				}
			}

			// Attempt to fix locations in I*
			for (k = 0; k < nIs; k++) {
				i = Is[k];
				assert(lrp[i] == 2);
				z0 = dual_bound_fix0Is(node_eval, i, min0);
				z1 = dual_bound_fix1Is(node_eval, i, max1);
				if (z0 > ub) {
					// We can fix location i as opened
					lrp[i] = -1; 
					node_eval->z = node->lb = z1; --nu;
					if (nI1 == 0) {
						Is[k] = Is[--nIs]; 
					} else {
						Is[k] = max1;
						lrp[max1] = 2;
						max1 = max1k = -1;
					}
					mins = maxs = minsk = maxsk = -1;
					node_eval->x[i] = 1;
					node = dual_node_fix_open(dual, node, i);
					// Try again?
					goto label_begin;
				} 
				if (z1 > ub) {
					// We can fix location i as closed
					lrp[i] = -1; 
					node_eval->z = node->lb = z0; --nu;
					if (nI0 == 0) {
						Is[k] = Is[--nIs];
					} else {
						Is[k] = min0;
						lrp[min0] = 2;
						min0 = min0k = -1;
					}
					mins = maxs = minsk = maxsk = -1;
					node_eval->x[i] = 0;
					node = dual_node_fix_closed(dual, node, i);
					// Try again?
					goto label_begin;
				} 
			}
		}
	}

	// Assertion tests
	/*
	   assert((max1 == -1 && max1k == -1 && nI1 == 0) || max1 == I1[max1k]);
	   for (k = 0; k < nI1; k++) {
	   assert(lrp[I1[k]] == 1);
	   assert(rc[max1] >= rc[I1[k]]);
	   }
	   assert((maxs == -1 && maxsk == -1 && nIs == 0) || maxs == Is[maxsk]);
	   assert((mins == -1 && minsk == -1 && nIs == 0) || mins == Is[minsk]);
	   for (k = 0; k < nIs; k++) {
	   assert(lrp[Is[k]] == 2);
	   assert(rc[mins] <= rc[Is[k]]);
	   assert(rc[maxs] >= rc[Is[k]]);
	   }
	   assert((min0 == -1 && min0k == -1 && nI0 == 0) || min0 == I0[min0k]);
	   for (k = 0; k < nI0; k++) {
	   assert(lrp[I0[k]] == 0);
	   assert(rc[min0] <= rc[I0[k]]);
	   }
	   for (k = inst->n, i = 0; i < inst->n; i++) {
	   if (lrp[i] == 0) {
	   for (l = 0; l < nI0; l++)
	   if (I0[l] == i)
	   break;
	   assert(l < nI0);
	   } else if (lrp[i] == 1) {
	   for (l = 0; l < nI1; l++)
	   if (I1[l] == i)
	   break;
	   assert(l < nI1);
	   } else if (lrp[i] == 2) {
	   for (l = 0; l < nIs; l++)
	   if (Is[l] == i)
	   break;
	   assert(l < nIs);
	   } else {
	   --k;
	   }
	   }
	   assert(k == nu);
	   assert(nu == nI0 + nI1 + nIs);
	   */


	// Reoptimize the lagrangian relaxation to update node_eval->x1 and node_eval->n_x1
	z = 0;
	for (k = 0; k < inst->m; k++)
		z += node->mul[k];
	node_eval->n_x1 = 0;
	for (i = 0; i < inst->n; i++) 
		if (node_eval->x[i]) {
			node_eval->x1[node_eval->n_x1++] = i;
			z += rc[i];
		}
	assert(fabs(z - node_eval->z) < bound_eps);

	node->lb = node_eval->z;
	if (node->lb > ub) {
		node_destroy(node);
		return NULL;
	}

	if (dual->phase == 2) {
		// Not preprocessing.
		// Select branching location for ExBranchP,
		// i.e. the unfixed location minimizing the absolute value of its reduced cost
		i = -1;
		for (k = 0; k < nIs; k++) 
			if (i == -1 || fabs(rc[Is[k]]) < fabs(rc[i]))
				i = Is[k];
		if (i == -1 && (max1 != -1 || min0 != -1)) {
			if (max1 == -1)
				i = min0;
			else if (min0 == -1)
				i = max1;
			else 
				i = (fabs(rc[max1]) < fabs(rc[min0])) ? max1 : min0;
		}
		node->branch_i = i;

		if (i != -1) {
			// Compute lower bounds for each subnode
			if (lrp[i] == 1) {
				node_eval->z_b1 = node_eval->z;
				node_eval->z_b0 = dual_bound_fix0I1(node_eval, i, min0, mins);
				assert(!inst->fix || node_eval->z_b0 <= ub);
			} else if (lrp[i] == 0) {
				node_eval->z_b0 = node_eval->z;
				node_eval->z_b1 = dual_bound_fix1I0(node_eval, i, max1, maxs); 
				assert(!inst->fix || node_eval->z_b1 <= ub);
			} else if (lrp[i] == 2) {
				node_eval->z_b0 = dual_bound_fix0Is(node_eval, i, min0);
				assert(!inst->fix || node_eval->z_b0 <= ub);
				node_eval->z_b1 = dual_bound_fix1Is(node_eval, i, max1);
				assert(!inst->fix || node_eval->z_b1 <= ub);
			} else 
				assert(0);
			assert(node_eval->z_b0 >= node_eval->z - bound_eps);
			assert(node_eval->z_b1 >= node_eval->z - bound_eps);
		}
	} else {
		node_eval->z_b0 = node_eval->z_b1 = node_eval->z;
	}

	return node;
}




// Implementation of the preprocessing rules of Goldengorin et al. 
// its performance is terrible and the code is ugly, but it works as specified
int dual_node_goldengorin (dual_t *dual, node_t *node)
{
	int nfix, i, j, k, imin1, imin2, ik, loop;
	instance_t *inst = dual->inst;
	node_eval_t *node_eval = dual->node_eval;
	double a, z;

	if (dual->phase == 1)
		return 0;

	do {
		loop = nfix = 0;
		for (i = 0; i < inst->n; i++) {
			if (get_bit(node->fix_mask, i)) {
				++nfix;
			} else {
				a = -inst->f[i];
				for (j = 0; j < inst->m; j++) {
					for (imin1 = imin2 = -1, k = 0; k < inst->n; k++) {
						ik = inst->inc[j][k];
						if (!get_bit(node->fix_mask, ik)) {
							if (imin1 == -1) {
								imin1 = ik;
								if (imin1 != i)
									break;
							} else {
								imin2 = ik;
								break;
							}
						}
					}
					if (imin1 == i && imin2 != -1) 
						a += inst->c[imin2][j] - inst->c[imin1][j];
				}
				if (a >= 0.0) {
					fprintf(stderr, "fix open %i %f\n", i, a);
					dual_node_fix_open(dual, node, i);
					loop = 1;
					++nfix;
				} else {
					for (j = 0; j < inst->m; j++) {
						for (imin2 = imin1 = -1, k = 0; k < inst->n; k++) {
							ik = inst->inc[j][k];
							if (!get_bit(node->fix_mask, ik)) {
								if (ik == i) {
									imin1 = ik;
									break;
								}
							}
						}
						if (imin1 != -1)
							a += inst->c[inst->inc[j][inst->n - 1]][j] - inst->c[imin1][j];
					}
					if (a <= 0.0) {
						fprintf(stderr, "fix closed %i\n", i);
						dual_node_fix_closed(dual, node, i);
						loop = 1;
						++nfix;
					}
				}
			}
		} 
	} while (loop && nfix < inst->n);

	if (nfix < inst->n)
		return 0;

	// compute primal cost
	for (j = 0; j < inst->m; j++)
		node_eval->min_c[j] = INFINITY;

	z = 0;
	for (i = 0; i < inst->n; i++) {
		node_eval->x[i] = 0;
		if (get_bit(node->fix_val, i)) {
			node_eval->x[i] = 1;
			z += inst->f[i];
			for (j = 0; j < inst->m; j++) 
				if (inst->c[i][j] < node_eval->min_c[j]) 
					node_eval->min_c[j] = inst->c[i][j]; 
		}
	}

	for (j = 0; j < inst->m; j++)
		z += node_eval->min_c[j];

	node->ub = node->lb = z;
	return 1;
}
