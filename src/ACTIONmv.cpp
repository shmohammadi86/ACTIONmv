#include <ACTIONmv.h>
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

#define ARMA_USE_CXX11_RNG

mat MWM(mat G);

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List run_ACTION_muV(const List& S, int k_min, int k_max, vec alpha, double lambda = 1, int AA_iters = 50, int Opt_iters= 100, int numThreads = 8) {
		
	int n_list = S.size();	
	vector<mat> cell_signatures(n_list);
	for(int i = 0; i < n_list; i++) {
		cell_signatures[i] = (as<mat>(S[i]));
	}
	
	
	full_trace run_trace = runACTION_muV(cell_signatures, k_min, k_max, alpha, lambda, AA_iters, Opt_iters, numThreads);
	
	
	List res;

	List H_consensus(k_max);
	for (int kk = k_min; kk <= k_max; kk++) {
		H_consensus[kk-1] = run_trace.H_consensus[kk];
	}
	res["H_consensus"] = H_consensus;	
		
	char ds_name[128];
	for(int i = 0; i < n_list; i++) {
		List individual_trace;

		List H_primary(k_max);
		for (int kk = k_min; kk <= k_max; kk++) {
			H_primary[kk-1] = run_trace.indiv_trace[kk].H_primary[i];
		}
		individual_trace["H_primary"] = H_primary;
		
		List H_secondary(k_max);
		for (int kk = k_min; kk <= k_max; kk++) {
			H_secondary[kk-1] = run_trace.indiv_trace[kk].H_secondary[i];
		}
		individual_trace["H_secondary"] = H_secondary;

		
		List C_primary(k_max);
		for (int kk = k_min; kk <= k_max; kk++) {
			C_primary[kk-1] = run_trace.indiv_trace[kk].C_primary[i];
		}
		individual_trace["C_primary"] = C_primary;	
	
	
		List C_consensus(k_max);
		for (int kk = k_min; kk <= k_max; kk++) {
			C_consensus[kk-1] = run_trace.indiv_trace[kk].C_consensus[i];
		}
		individual_trace["C_consensus"] = C_consensus;		
	
		sprintf(ds_name, "View%d_trace", i+1);
		res[ds_name] = individual_trace;
	}
	
		
	return res;	
}
