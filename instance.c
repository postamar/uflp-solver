// (C) 2011, Marius Posta (postamar@iro.umontreal.ca)
// Check LICENSE.txt for the legal blah-blah.

#include "instance.h"
#include <getopt.h>



// ** BLAS **

#ifdef __cplusplus
extern "C" double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy);
#else
extern double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy);
#endif

double ddot (int n, const double *x, const double *y)
{
	int one = 1;
	return ddot_(&n, x, &one, y, &one);
}


// ** Quickselect **

// Comparison function used by quickselect(...)
int ic_pair_cmp (const void *a, const void *b) 
{
	ic_pair_t *pa = (ic_pair_t*) a;
	ic_pair_t *pb = (ic_pair_t*) b;
	if (pa->c == pb->c && pa->i == pb->i)
		return 0;
	return (pa->c < pb->c || (pa->c == pb->c && pa->i < pb->i)) ? -1 : 1;
}

// Auxiliary function used by quickselect(...)
// same partition function as in quicksort
int partition (ic_pair_t *queue, int left, int right, int pivot_index)
{
	int i, store_index = left;
	ic_pair_t pivot = queue[pivot_index];

	ic_pair_t temp = queue[pivot_index];
	queue[pivot_index] = queue[right];
	queue[right] = temp;

	for (i = left; i < right; i++)
		if (-1 == ic_pair_cmp(&queue[i], &pivot)) {
			temp = queue[store_index];
			queue[store_index++] = queue[i];
			queue[i] = temp;
		}

	temp = queue[store_index];
	queue[store_index] = queue[right];
	queue[right] = temp;

	return store_index;
}

// Selects the k smallest elements among n in the array 'queue'.
// Pivots are selected randomly using the stream 'rng'
void quickselect (RngStream rng, ic_pair_t *queue, int n, int k)
{
	int pivot_index, pivot_new_index;
	int left = 0;
	int right = n - 1;
	if (n == 0 || k <= 0 || k >= n)
		return;
	--k;

	while (1) {
		pivot_index = RngStream_RandInt(rng, left, right);
		pivot_new_index = partition(queue, left, right, pivot_index);
		if (k < pivot_new_index)
			right = pivot_new_index - 1;
		else if (k > pivot_new_index)
			left = pivot_new_index + 1;
		else
			break;
	}
}


// ** instance_t **

void instance_alloc (instance_t *inst, int n, int m)
{
	int i, j;

	inst->n = n;
	inst->m = m;

	inst->n_words_in_bit_vector = 1 + n / ((int) sizeof(word_t) * 8);

	inst->f = (double*) calloc(n, sizeof(double));
	inst->c = (double**) calloc(n, sizeof(double*));
	for (i = 0; i < n; i++)
		inst->c[i] = (double*) calloc(m, sizeof(double));
	inst->inc = (int**) calloc(m, sizeof(int*));
	for (j = 0; j < m; j++)
		inst->inc[j] = (int*) calloc(n, sizeof(int));
}



instance_t *instance_copy (instance_t *inst)
{
	int i, j;
	instance_t *copy = (instance_t*) malloc(sizeof(instance_t));

	*copy = *inst;

	copy->f = (double*) calloc(copy->n, sizeof(double));
	copy->c = (double**) calloc(copy->n, sizeof(double*));
	for (i = 0; i < copy->n; i++)
		copy->c[i] = (double*) calloc(copy->m, sizeof(double));

	copy->inc = (int**) calloc(copy->m, sizeof(int*));
	for (j = 0; j < copy->m; j++)
		copy->inc[j] = (int*) calloc(copy->n, sizeof(int));

	return copy;
}


void instance_destroy (instance_t *inst)
{
	int i, j;

	for (j = 0; j < inst->m; j++)
		free(inst->inc[j]);
	free(inst->inc);

	for (i = 0; i < inst->n; i++)
		free(inst->c[i]);
	free(inst->c);
	free(inst->f);

	free(inst);
}

// Reads costs from stdin according to the appropriate format,
// returns 0 on failure, 1 on success.
int instance_read (instance_t *inst)
{
	int i, j, k;
	double fcost;
	char buf[256];

	switch (inst->fmt) {
		case inst_fmt_cap:
			for (i = 0; i < inst->n; i++) 
				if (2 != fscanf(stdin, "%s %lf", buf, &inst->f[i])) return 0;
			for (j = 0; j < inst->m; j++) {
				if (1 != fscanf(stdin, "%s", buf)) return 0;
				for (i = 0; i < inst->n; i++) 
					if (1 != fscanf(stdin, "%lf", &inst->c[i][j])) return 0;
			}
			break;
		case inst_fmt_orlib:
			for (k = 0; k < inst->n; k++) {
				if (1 != fscanf(stdin, "%i", &i)) return 0;
				if (k == 0 && i == 0) {
					--k;
					continue;
				}
				if (1 != fscanf(stdin, "%lf", &fcost)) return 0;
				inst->f[--i] = fcost;
				for (j = 0; j < inst->m; j++) 
					if (1 != fscanf(stdin, "%lf", &inst->c[i][j])) return 0;
			}
			break;
		default:
			assert(0);
	}


	return 1;
}


// Sorts locations according to service cost, for each client.
// Also computes an initial upper bound.
double instance_reset (instance_t *inst)
{
	int i, j;
	double maxc, ub;
	ic_pair_t *queue;

	queue = (ic_pair_t*) calloc(inst->n, sizeof(ic_pair_t));
	// get ordinals 
	for (j = 0; j < inst->m; j++) {
		for (i = 0; i < inst->n; i++) {
			queue[i].i = i;
			queue[i].c = inst->c[i][j];
		}
		qsort(queue, inst->n, sizeof(ic_pair_t), ic_pair_cmp);
		for (i = 0; i < inst->n; i++)
			inst->inc[j][i] = queue[i].i;
	}
	free(queue);

	ub = 0.0;
	for (i = 0; i < inst->n; i++)
		ub += inst->f[i];

	for (j = 0; j < inst->m; j++) {
		for (maxc = -INFINITY, i = 0; i < inst->n; i++)
			if (inst->c[i][j] > maxc && inst->c[i][j] < 1e19)
				maxc = inst->c[i][j];
		ub += maxc;
	}

	return (inst->ub < ub) ? inst->ub : ub;
}


// Copies search parameters and instance costs
void instance_sync (instance_t *master, instance_t *slave)
{
	int i, j;

	strncpy(slave->name, master->name, 64);
	slave->request_period = master->request_period;
	slave->bundle_max_size = master->bundle_max_size;
	slave->n_root_node_max_iter = master->n_root_node_max_iter;
	slave->n_node_max_iter_phase1 = master->n_node_max_iter_phase1;
	slave->n_node_max_iter_phase2 = master->n_node_max_iter_phase2;
	slave->max_utime = master->max_utime;
	slave->ub = master->ub;
	slave->fix = master->fix;
	slave->log = master->log;
	slave->printsol = master->printsol;
	slave->strategy = master->strategy;

	memcpy(slave->f, master->f, master->n * sizeof(double));
	for (i = 0; i < master->n; i++)
		memcpy(slave->c[i], master->c[i], master->m * sizeof(double));
	for (j = 0; j < master->m; j++)
		memcpy(slave->inc[j], master->inc[j], master->n * sizeof(int));
}



// ** shared_t **

shared_t *shared_create (instance_t *inst)
{
	int i;
	shared_t *shared = (shared_t*) malloc(sizeof(shared_t));

	shared->inst = inst;
	shared->n_solved = 0;
	shared->search_state = state_waiting;
	
	shared->improving_partial_x = (int*) calloc(inst->n, sizeof(int));
	shared->guiding_x = (int*) calloc(inst->n, sizeof(int));
	shared->best_x = (int*) calloc(inst->n, sizeof(int));
	shared->best_z = INFINITY;

	for (i = 0; i < inst->n; i++) {
		shared->best_x[i] = shared->guiding_x[i] = 0;
		shared->improving_partial_x[i] = -1;
	}

	return shared;
}

void shared_destroy (shared_t *shared)
{
	free(shared->best_x);
	free(shared->guiding_x);
	free(shared->improving_partial_x);

	free(shared);
}

void shared_reset (shared_t *shared, instance_t *inst, double ub)
{
	int i;

	instance_sync(inst, shared->inst);

	shared->best_z = ub;
	for (i = 0; i < shared->inst->n; i++) {
		shared->best_x[i] = shared->guiding_x[i] = 0;
		shared->improving_partial_x[i] = -1;
	}

	shared->global_lb = 0.0;
	shared->n_node_eval = 0;
	shared->n_lag_eval = 0;
	shared->n_moves = 0;

	shared->search_state = state_running;
}

void shared_output_results (shared_t *shared)
{
	instance_t *inst;
	double z;
	int i, j, best_i;

	inst = shared->inst;
	switch (inst->printsol) {
		case 2:
			// Output best known solution and its cost
			for (z = 0.0, i = 0; i < inst->n; i++)
				if (shared->best_x[i] == 1 || inst->f[i] < 0.0) 
					z += inst->f[i];

			for (j = 0; j < inst->m; j++) {
				for (best_i = -1, i = 0; i < inst->n; i++)
					if (shared->best_x[i] == 1 || inst->f[i] < 0.0) {
						if (best_i == -1 || inst->c[i][j] < inst->c[best_i][j])
							best_i = i;
					}
				z += inst->c[best_i][j];
				fprintf(stdout, "%i ", best_i);
			}
			fprintf(stdout, "%.2f\n", z);
			assert(fabs(z - shared->best_z) < bound_eps);
		case 1:
			// Output statistics
			fprintf(stdout, "%s\t", inst->name);
			fprintf(stdout, "%.2f\t", shared->global_lb);
			fprintf(stdout, "%.2f\t", shared->best_z);
			if (shared->search_state == state_solved)
				fprintf(stdout, "%.2f\t", shared->exec_time);
			else 
				fprintf(stdout, "timeout\t");
			fprintf(stdout, "%i\t%i\t%i\n", shared->n_node_eval, shared->n_lag_eval, shared->n_moves);
			fflush(stdout);
	}
}



// Print help string to stderr, using default parameter values in 'defarg'
void print_help (const search_param_t defarg) {
#define HELP(x) fprintf(stderr, "  %30s  ", (x));

	fprintf(stderr, "\nThis program performs the following loop:\n 1) reads solver parameters and instance data in stdin,\n 2) solves the instance,\n 3) prints the output in stdout.\n\nUsage:  solver [options]\nOptions:\n");
	HELP("-n NAME, --name=NAME"); fprintf(stderr, "Instance name,\n");
	HELP(""); fprintf(stderr, " string, default: %s.\n", defarg.name);
	HELP("-p BIAS, --primalbias=BIAS"); fprintf(stderr, "Bias toward primal iterations,\n");
	HELP(""); fprintf(stderr, " non-negative double; default: %.2f.\n", defarg.primalbias);
	HELP("-s STRAT, --strategy=TRAT"); fprintf(stderr, "Node selection strategy,\n");
	HELP(""); fprintf(stderr, " {bfs, dfs}; default: %s.\n", defarg.strategy);
	HELP("-c LIMIT, --cachesize=LIMIT"); fprintf(stderr, "Maximal subgradient cache size,\n");
	HELP(""); fprintf(stderr, " unsigned short; default: %d.\n", defarg.cachesize);
	HELP("-g PERIOD, --guiding=PERIOD"); fprintf(stderr, "Guiding solution update period,\n");
	HELP(""); fprintf(stderr, " positive int; default: %d.\n", defarg.guiding);
	HELP("-b LIMIT, --bundlesize=LIMIT"); fprintf(stderr, "Maximal bundle size,\n");
	HELP(""); fprintf(stderr, " positive int; default: %d.\n", defarg.bundlesize);
	HELP("-r LIMIT, --rootnodeiter=LIMIT"); fprintf(stderr, "Maximal bundle iterations at root node,\n");
	HELP(""); fprintf(stderr, " positive int; default: %d.\n", defarg.rootnodeiter);
	HELP("-i LIMIT, --subnodeiter=LIMIT"); fprintf(stderr, "Maximal bundle iterations at non-root node,\n");
	HELP(""); fprintf(stderr, " positive int; default: %d.\n", defarg.subnodeiter);
	HELP("-t LIMIT, --time=LIMIT"); fprintf(stderr, "Maximal execution time (utime in seconds),\n");
	HELP(""); fprintf(stderr, " non-negative double; default: %.2f.\n", defarg.time);
	HELP("-u LIMIT, --upperbound=LIMIT"); fprintf(stderr, "Objective value upper bound,\n");
	HELP(""); fprintf(stderr, " double; default: %.2f.\n", defarg.upperbound);
	HELP("-f STRAT, --implbranch=STRAT"); fprintf(stderr, "Implicit branching setting.\n");
	HELP(""); fprintf(stderr, " {on,off}; default: %s.\n", defarg.implbranch);
	HELP("-l, --log"); fprintf(stderr, "Logging setting (stderr),\n");
	HELP(""); fprintf(stderr, " {on,off}; default: %s.\n", defarg.log);
	HELP("-o, --output"); fprintf(stderr, "Best solution output setting (stdout),\n");
	HELP(""); fprintf(stderr, " {none,short,full}; default: %s.\n", defarg.output);
	HELP("-?, --help"); fprintf(stderr, "Prints this help to stderr.\n");
	HELP("-v, --version"); fprintf(stderr, "Prints version info to stderr.\n");
	fprintf(stderr, "Input:\n");
	HELP("quit"); fprintf(stderr, "Exits the solver.\n");
	HELP("set PARAM VALUE"); fprintf(stderr, "Sets parameter PARAM to new value VALUE.\n");
	HELP("INSTANCE_DATA"); fprintf(stderr, "Sets prices to INSTANCE_DATA and solves.\n");
	fprintf(stderr, "Output:\n  <NAME> <lb> <ub> <duration> <nnodes> <nlageval> <nmoves> \n\n");
}


// Sets search parameter if possible, complains to stderr otherwise
void instance_set_parameter (instance_t *inst, const char name[256], const char val[256]) 
{
	int o;
	double d;

	if (!strcmp(name, "name")) 
		strcpy(inst->name, val);
	else if (!strcmp(name, "strategy")) {
		if (!strcmp(val, "bfs"))
			inst->strategy = search_bf;
		else if (!strcmp(val, "dfs"))
			inst->strategy = search_df;
		else 
			fprintf(stderr, "WARNING:\tunknown node selection strategy.\n");
	} else if (!strcmp(name, "implbranch")) {
		if (!strcmp(val, "on"))
			inst->fix = 1;
		else if (!strcmp(val, "off"))
			inst->fix = 0;
		else 
			fprintf(stderr, "WARNING:\tunknown implicit branching strategy.\n");
	} else if (!strcmp(name, "log")) {
		if (!strcmp(val, "on"))
			inst->log= 1;
		else if (!strcmp(val, "off"))
			inst->log= 0;
		else 
			fprintf(stderr, "WARNING:\tunknown log setting.\n");
	} else if (!strcmp(name, "output")) {
		if (!strcmp(val, "short"))
			inst->printsol = 1;
		else if (!strcmp(val, "none"))
			inst->printsol = 0;
		else if (!strcmp(val, "full"))
			inst->printsol = 2;
		else 
			fprintf(stderr, "WARNING:\tunknown output setting.\n");
	} else if (!strcmp(name, "cachesize")) {
		if (1 != sscanf(val, "%i", &o) || o < 0 || o > 65535)
			fprintf(stderr, "WARNING:\tvalue not an integer in range [0, 65535].\n");
		else 
			inst->subg_cache_size = (o == 1) ? 0 : o;
	} else if (!strcmp(name, "guiding")) {
		if (1 != sscanf(val, "%i", &o) || o < 1 || o > 65535)
			fprintf(stderr, "WARNING:\tvalue not an integer in range [1, 65535].\n");
		else 
			inst->request_period = o;
	} else if (!strcmp(name, "bundlesize")) {
		if (1 != sscanf(val, "%i", &o) || o < 1 || o > 65535)
			fprintf(stderr, "WARNING:\tvalue not an integer in range [1, 65535].\n");
		else 
			inst->bundle_max_size = o;
	} else if (!strcmp(name, "rootnodeiter")) {
		if (1 != sscanf(val, "%i", &o) || o < 1 || o > 65535)
			fprintf(stderr, "WARNING:\tvalue not an integer in range [1, 65535].\n");
		else 
			inst->n_root_node_max_iter = o;
	} else if (!strcmp(name, "subnodeiter")) {
		if (1 != sscanf(val, "%i", &o) || o < 1 || o > 65535)
			fprintf(stderr, "WARNING:\tvalue not an integer in range [1, 65535].\n");
		else 
			inst->n_node_max_iter_phase1 = inst->n_node_max_iter_phase2 = o;
	} else if (!strcmp(name, "time")) {
		if (1 != sscanf(val, "%lf", &d) || d <= 0.0)
			fprintf(stderr, "WARNING:\tvalue not a positive double.\n");
		else 
			inst->max_utime = d;
	} else if (!strcmp(name, "upperbound")) {
		if (1 != sscanf(val, "%lf", &d))
			fprintf(stderr, "WARNING:\tvalue not a double.\n");
		else 
			inst->ub = d;
	} else if (!strcmp(name, "primalbias")) {
		if (1 != sscanf(val, "%lf", &d) || d < 0.0)
			fprintf(stderr, "WARNING:\tvalue not a non-negative double.\n");
		else 
			inst->primal_bias = d;
	} else 
	
		fprintf(stderr, "WARNING:\tparameter name not recognized.\n");
}


// Parse command line arguments for search parameters
int instance_parse_args (instance_t *inst, int argc, char **argv, const search_param_t defarg) 
{
	static struct option longopts[] = {
		{ "name",			required_argument,	NULL,	'n'},
		{ "primalbias",		required_argument,	NULL,	'p'},
		{ "strategy",		required_argument,	NULL,	's'},
		{ "cachesize",		required_argument,	NULL,	'c'},
		{ "guiding",		required_argument,	NULL,	'g'},
		{ "bundlesize",		required_argument,	NULL,	'b'},
		{ "rootnodeiter",	required_argument,	NULL,	'r'},
		{ "subnodeiter",	required_argument,  NULL,	'i'},
		{ "time",			required_argument,	NULL,	't'},
		{ "upperbound",		required_argument,	NULL,	'u'},
		{ "implbranch",		required_argument,	NULL,	'f'},
		{ "log",			required_argument,	NULL,	'l'},
		{ "output",			required_argument,	NULL,	'o'},
		{ "help",			no_argument,		NULL,	'?'},
		{ "version",		no_argument,		NULL,	'v'},
		{ NULL,				0,					NULL,	0}};
	int ch;
	char name[256];
	char val[256];
	instance_set_parameter(inst, "name", defarg.name);
	instance_set_parameter(inst, "strategy", defarg.strategy);
	instance_set_parameter(inst, "implbranch", defarg.implbranch);
	instance_set_parameter(inst, "log", defarg.log);
	instance_set_parameter(inst, "output", defarg.output);
	sprintf(val, "%f", defarg.primalbias); 
	instance_set_parameter(inst, "primalbias", val);
	sprintf(val, "%f", defarg.upperbound); 
	instance_set_parameter(inst, "upperbound", val);
	sprintf(val, "%f", defarg.time); 
	instance_set_parameter(inst, "time", val);
	sprintf(val, "%d", defarg.cachesize);
	instance_set_parameter(inst, "cachesize", val);
	sprintf(val, "%d", defarg.guiding);
	instance_set_parameter(inst, "guiding", val);
	sprintf(val, "%d", defarg.bundlesize);
	instance_set_parameter(inst, "bundlesize", val);
	sprintf(val, "%d", defarg.rootnodeiter);
	instance_set_parameter(inst, "rootnodeiter", val);
	sprintf(val, "%d", defarg.subnodeiter);
	instance_set_parameter(inst, "subnodeiter", val);

	while ((ch = getopt_long(argc, argv, "n:p:s:c:g:b:r:i:t:u:f:l:o:v?", longopts, NULL)) != -1) {
		switch (ch) {
			case 'n':
				strcpy(name, "name"); break;
			case 'p':
				strcpy(name, "primalbias"); break;
			case 's':
				strcpy(name, "strategy"); break;
			case 'c':
				strcpy(name, "cachesize"); break;
			case 'g':
				strcpy(name, "guiding"); break;
			case 'b':
				strcpy(name, "bundlesize"); break;
			case 'r':
				strcpy(name, "rootnodeiter"); break;
			case 'i':
				strcpy(name, "subnodeiter"); break;
			case 't':
				strcpy(name, "time"); break;
			case 'u':
				strcpy(name, "upperbound"); break;
			case 'f':
				strcpy(name, "implbranch"); break;
			case 'l':
				strcpy(name, "log"); break;
			case 'o':
				strcpy(name, "output"); break;
			case 'v':
				fprintf(stderr, "SPLP solver, built %s.\n", __DATE__);
				return 1;
			case '?':
			default:
				print_help(defarg);
				return 1;
		}	
		strcpy(val, optarg);
		fprintf(stderr, "set %s %s\n", name, val);
		instance_set_parameter(inst, name, val);
	}
	return 0;
}



// Reads instance data header,
// determines instance size if first_read is true,
// determintes instance cost format.
int instance_read_header (instance_t *inst, int first_read)
{
	int c, rval, n, m;
	char ws;
	char buf[256];
	char name[256];
	char val[256];

label_read_header:
	while ((c = getchar()) == EOF) {};
	ws = c;
	if (ws == ' ' || ws == '\t' || ws == '\r' || ws == '\n')
		goto label_read_header;
	rval = ungetc(c, stdin);
	assert(rval == c);

	n = fscanf(stdin, "%s", buf);
	if (n == 0)
		goto label_read_header;
	if (n != 1) {
		fprintf(stderr, "ERROR: failed to read stdin.\n");
		return 0;
	}

	if (!strncmp(buf, "quit", 4)) {
		return 0;
	}

	if (!strncmp(buf, "set", 3)) {
		val[0] = name[0] = '\0';
		if (2 != fscanf(stdin, "%s %s", name, val)) {
			fprintf(stderr, "WARNING\tinvalid syntax in command `set %s %s`\n", name, val);
			goto label_read_header;
		}
		instance_set_parameter(inst, name, val);
		goto label_read_header;
	}


	if (!strncmp(buf, "FILE:", 5)) {
		if (1 != fscanf(stdin, "%s", name)) {
			fprintf(stderr, "ERROR: failed to read line 1 of ORLIB-formatted instance data.\n");
			return 0;
		}
		if (2 != fscanf(stdin, "%i %i", &n, &m)) {
			fprintf(stderr, "ERROR: failed to read line 2 of ORLIB-formatted instance data.\n");
			return 0;
		}
		inst->fmt = inst_fmt_orlib;
	} else {
		if (1 != sscanf(buf, "%i", &n)) {
			fprintf(stderr, "ERROR: failed to read line 1 of CFL-formatted instance data.\n");
			return 0;
		}
		if (1 != fscanf(stdin, "%i", &m)) {
			fprintf(stderr, "ERROR: failed to read line 1 of CFL-formatted instance data.\n");
			return 0;
		}
		inst->fmt = inst_fmt_cap;
	}

	if (n < 1 || m < 1) {
		fprintf(stderr, "ERROR: n = %d or m = %d out of bounds.\n", n, m);
		return 0;
	}
	if (first_read) {
		inst->n = n, inst->m = m;
	}
	else if (n != inst->n || m != inst->m) {
		fprintf(stderr, "ERROR: n = %d or m = %d are not equal to %d and %d, respectively, as they should.\n", n, m, inst->n, inst->m);
		return 0;
	}

	return 1;
}


