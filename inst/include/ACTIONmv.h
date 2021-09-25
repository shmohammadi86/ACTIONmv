#ifndef STARTER_H
#define STARTER_H

#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>
#include <ctime>

#define ARMA_DONT_USE_WRAPPER
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE
#define ARMA_64BIT_WORD
#define ARMA_BLAS_LONG_LONG

#include <arma/armadillo>
#include <my_cblas.h>

using namespace std;
using namespace arma;

// min_{X} (|| AX - B ||) s.t. simplex constraint using ACTIVE Set Method
// mat run_simplex_regression(mat &A, mat &B);
mat run_simplex_regression(mat &A, mat &B, bool computeXtX);

struct ACTION_results {
	vector<uvec> selected_cols;
	vector<mat> H;
	vector<mat> C;
};


struct mvtrace_obj {
	vector<uvec> selected_cols;

	vector<mat> H_primary;
	vector<mat> C_primary;

	vector<mat> H_secondary;
	vector<mat> C_secondary;

	vector<mat> C_consensus;
};

struct full_trace {
	vector<mvtrace_obj> indiv_trace;
	vector<mat> H_consensus;
	
};

full_trace runACTION_muV(vector<mat> cell_signatures, int k_min, int k_max, vec alpha, double lambda = 1, int AA_iters = 50, int Opt_iters = 0);
field<mat> run_AA(mat &A, mat &W0, int max_it = 100, double min_delta = 1e-6);


#endif
