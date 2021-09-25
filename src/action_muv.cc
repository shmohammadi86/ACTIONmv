#include "ACTIONmv.h"



double *l1, *l2, *w;	
int *match1, *match2, *v1, *v2;	
int *s, *tt, *deg, *offset, *my_list;


/**
 * n the number of nodes
 * m the number of nodes
 * nedges the number of edges
 * vv1 is the source for each of the nedges 
 * vv2 is the target for each of the nedges
 * weight is the weight of each of the nedges
 * out1 is a vector of length at most min(n,m),
 * out2 is a vector of length at most min(n,m),
 * noutedges is the number of out edges
 */
double MWM_driver(int n, int m, int nedges, double *vv1, double *vv2, double *weight, double *out1, double *out2, int *noutedges) {	
	double ret, al;
	int i, j, k, p, q, r, t1, t2;
		

	for (i = 0; i < nedges; i++) {
		v1[i] = (int)(vv1[i] + .5);
		v2[i] = (int)(vv2[i] + .5);
	}
	for (i = 0; i < n; i++) {
		offset[i] = 0;
		deg[i] = 1;
	}
	for (i = 0; i < nedges; i++) deg[v1[i]]++;
	for (i = 1; i < n; i++) offset[i] = offset[i-1] + deg[i-1];
	for (i = 0; i < n; i++) deg[i] = 0;
	for (i = 0; i < nedges; i++) {
		my_list[offset[v1[i]] + deg[v1[i]]] = v2[i];
		w[offset[v1[i]] + deg[v1[i]]] = weight[i];
		deg[(int)v1[i]]++;
	}
	for (i = 0; i < n; i++) {
		my_list[offset[i] + deg[i]] = m + i;
		w[offset[i] + deg[i]] = 0;
		deg[i]++;
	}
	for (i = 0; i < n; i++) {
		l1[i] = 0;
		for (j = 0; j < deg[i]; j++) {
			if (w[offset[i]+j] > l1[i]) l1[i] = w[offset[i] + j];
		}
	}
	for (i = 0; i < n; i++) {
		match1[i] = -1;
	}
	for (i = 0; i < n + m; i++) {
		l2[i] = 0;
		match2[i] = -1;
	}
	for (i = 0; i < n; i++) {
		for(j = 0; j < n + m; j++) tt[j] = -1;
		s[p = q = 0] = i;
		for(; p <= q; p++) {
			if (match1[i] >= 0) break;
			k = s[p];
			for (r = 0; r < deg[k]; r++) {
				if (match1[i] >= 0) break;
				j = my_list[offset[k] + r];
				if (w[offset[k] + r] < l1[k] + l2[j] - 1e-8) continue;
				if (tt[j] < 0) {
					s[++q] = match2[j];
					tt[j] = k;
					if (match2[j] < 0) {
						for(; j>=0 ;) {
							k = match2[j] = tt[j];
							p = match1[k];
							match1[k] = j;
							j = p;
						}
					}
				}
			}
		}
		if (match1[i] < 0) {
			al = 1e20;
			for (j = 0; j < p; j++) {
				t1 = s[j];
				for (k = 0; k < deg[t1]; k++) {
					t2 = my_list[offset[t1] + k];
					if (tt[t2] < 0 && l1[t1] + l2[t2] - w[offset[t1] + k] < al) {
						al = l1[t1] + l2[t2] - w[offset[t1] + k];
					}
				}
			}
			for (j = 0; j < p; j++) l1[s[j]] -= al;
			for (j = 0; j < n + m; j++) if (tt[j] >= 0) l2[j] += al;
			i--;
			continue;
		}
	}
	ret = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < deg[i]; j++) {
			if (my_list[offset[i] + j] == match1[i]) {
				ret += w[offset[i] + j];
			}
		}
	}
    *noutedges = 0;
    for (i = 0; i < n; i++) {
        if (match1[i] < m) (*noutedges)++;
    }
    *noutedges = 0;
    for (i = 0; i < n; i++) {
        if (match1[i] < m) {
            out1[*noutedges] = i;
            out2[*noutedges] = match1[i];
            (*noutedges)++;
        }
    }
    
	return ret;
}


mat MWM(mat G) {
	// Memory allocation for bipartite matching function
	l1 = new double[(G.n_rows)];
	l2 = new double[(G.n_rows)+(G.n_cols)];
	v1 = new int[(G.n_rows)*(G.n_cols)];
	v2 = new int[(G.n_rows)*(G.n_cols)];
	s = new int[(G.n_rows)];
	tt = new int[(G.n_rows)+(G.n_cols)];
	match1 = new int[(G.n_rows)];
	match2 = new int[(G.n_rows)+(G.n_cols)];
	offset = new int[(G.n_rows)];
	deg = new int[(G.n_rows)];
	my_list = new int[(G.n_rows)*(G.n_cols) + (G.n_rows)];
	w = new double[(G.n_rows)*(G.n_cols) + (G.n_rows)];	
	
	
	mat G_matched = zeros(size(G));

	int n = G.n_rows;
	int m = G.n_cols;
	
	uvec idx = find(G);
	if(idx.n_elem == 0)
		return G_matched;
	
	int nedges = idx.n_elem;
		
	double *vv1 = new double[nedges];
	double *vv2 = new double[nedges];
	double *weight = new double[nedges];
	
	umat subs = ind2sub( size(G), idx );
	for (int i = 0; i < nedges; i++) {
		weight[i] = G(idx(i));
		vv1[i] = subs(0, i);
		vv2[i] = subs(1, i);		
	}
	
	int match_no = std::min(m, n);
	double *ii = new double[match_no];
	double *jj = new double[match_no];
	
	int matched_edge_no;
	
	MWM_driver(n, m, nedges, vv1, vv2, weight, ii, jj, &matched_edge_no);	
	
	for (int k = 0; k < matched_edge_no; k++) {
		G_matched(ii[k], jj[k]) = G(ii[k], jj[k]);
	}
	
	
	delete [] vv1;
	delete [] vv2;
	delete [] weight;
	delete [] ii;
	delete [] jj;

    delete [] l1;
    delete [] l2;
    delete [] v1;
    delete [] v2;
    delete [] s;
    delete [] tt;
	delete [] match1;
	delete [] match2;
    delete [] offset;
    delete [] deg;
    delete [] my_list;
    delete [] w;
    	
	return G_matched;
}

  field<mat> run_AA(mat &A, mat &W0, int max_it, double min_delta)
  {
    int sample_no = A.n_cols;
    int d = A.n_rows;  // input dimension
    int k = W0.n_cols; // AA components

    mat C = zeros(sample_no, k);
    mat H = zeros(k, sample_no);

    mat W = W0;
    vec c(sample_no);

    double old_RSS = 0;
    // printf("(New) %d- %d\n", k, max_it);

    for (int it = 0; it < max_it; it++)
    {
      mat C_old = C;
      mat H_old = H;

      double A_norm = norm(A, "fro");
      H = run_simplex_regression(W, A, true);
      mat R = A - W * H;
      mat Ht = trans(H);
      for (int i = 0; i < k; i++)
      {
        vec w = W.col(i);
        vec h = Ht.col(i);

        double norm_sq = arma::dot(h, h);
        if (norm_sq < double(10e-8))
        {
          // singular
          int max_res_idx = index_max(rowvec(sum(square(R), 0)));
          W.col(i) = A.col(max_res_idx);
          c.zeros();
          c(max_res_idx) = 1;
          C.col(i) = c;
        }
        else
        {
          // b = (1.0 / norm_sq) *R*ht + w;
          vec b = w;
          cblas_dgemv(CblasColMajor, CblasNoTrans, R.n_rows, R.n_cols,
                      (1.0 / norm_sq), R.memptr(), R.n_rows, h.memptr(), 1, 1,
                      b.memptr(), 1);

          C.col(i) = run_simplex_regression(A, b, false);

          vec w_new = A * C.col(i);
          vec delta = (w - w_new);

          // Rank-1 update: R += delta*h
          cblas_dger(CblasColMajor, R.n_rows, R.n_cols, 1.0, delta.memptr(), 1,
                     h.memptr(), 1, R.memptr(), R.n_rows);

          W.col(i) = w_new;
        }
      }
      double RSS = norm(R, "fro");
      double delta_RSS = abs(RSS - old_RSS) / A_norm;
      old_RSS = RSS;
      /*
    double delta_C = norm(C - C_old, "fro") / norm(C, "fro");
    double delta_H = norm(H - H_old, "fro") / norm(H, "fro");	
    printf("\t<%d, %d>- norm_RSS = %e, delta_RSS = %e, delta_C = %.3e, delta_H = %.3e\n", k, it, RSS/A_norm, delta_RSS, delta_C, delta_H);
    printf("\t<%d, %d>- norm_RSS = %e\n", k, it, delta_RSS);
    */

      if (delta_RSS < min_delta)
        break;
    }

    C = clamp(C, 0, 1);
    C = normalise(C, 1);
    H = clamp(H, 0, 1);
    H = normalise(H, 1);

    field<mat> decomposition(2, 1);
    decomposition(0) = C;
    decomposition(1) = H;

    return decomposition;
  }
  
uvec SPA(mat M, int k) {	

	int n = M.n_cols;
	uvec K(k); // selected columns from M
		
 
	rowvec normM = sum(M % M, 0); 
	rowvec normM1 = normM;
	
	mat U(M.n_rows, k);
	
	double eps = 1e-6;
	for (int i = 0; i < k; i++) {
		// Find the column with maximum norm. In case of having more than one column with almost very small diff in norm, pick the one that originally had the largest norm
		double a = max(normM); 
		
		uvec b = find((a*ones(1, n)-normM)/a <= eps); 
		
		if(b.n_elem > 1) {
			uword idx = index_max(normM1(b)); 
			K(i) = b(idx);
		}
		else {
			K(i) = b(0);			
		}			
		
		// Pick column
		U.col(i) = M.col(K(i));

		// Orthogonalize with respect to current basis
		for (int j = 0; j < i-1; j++) {
			U.col(i) = U.col(i) - dot(U.col(j), U.col(i)) * U.col(j);
		}
		U.col(i) = U.col(i)/ norm(U.col(i), 2); 
		
		// Update column norms
		vec u = U.col(i);            
		for (int j = i-1; 0 <= j; j--) {
			u = u - dot(U.col(j), u)*U.col(j); 
		}
		rowvec r = u.t()*M;
		normM = normM - (r % r);
	}
		
	return K;
}


/* alpha: sums to one and indicates relative importance of each view
 * 
 */
void findConsensus(vector<mat> S, full_trace &run_trace, int arch_no, vec alpha, double lambda, int max_it, double lambda2 = 1e-5, double epsilon = 1e-5) {
	printf("Find shared subspace\n");
	
	register int i;
	int ds_no = run_trace.indiv_trace[arch_no].H_primary.size(); // number of datasets ( ~ 2)	
	int cell_no = S[0].n_cols;
	vec c(cell_no);	
	

	// make sure it's a convex vector
	alpha.transform( [](double val) { return (max(0.0, val)); } );
	alpha = normalise(alpha, 1); 
	
	
	run_trace.indiv_trace[arch_no].C_secondary[0] = run_trace.indiv_trace[arch_no].C_primary[0]; 
	run_trace.indiv_trace[arch_no].H_secondary[0] = run_trace.indiv_trace[arch_no].H_primary[0]; 
	for(int ds = 1; ds < ds_no; ds++) {		
		mat G = 1 + cor(trans(run_trace.indiv_trace[arch_no].H_primary[0]), trans(run_trace.indiv_trace[arch_no].H_primary[ds]));
		mat G_matched = MWM(G);
		uvec perm = index_max(G_matched, 1);
		
		run_trace.indiv_trace[arch_no].C_secondary[ds] = run_trace.indiv_trace[arch_no].C_primary[ds].cols(perm); 
		run_trace.indiv_trace[arch_no].H_secondary[ds] = run_trace.indiv_trace[arch_no].H_primary[ds].rows(perm); 		
	}

	// Compute initial H_consensus
	run_trace.H_consensus[arch_no] = zeros(size(run_trace.indiv_trace[arch_no].H_primary[0]));	
	for(int ds = 0; ds < ds_no; ds++) {
		run_trace.H_consensus[arch_no] += alpha(ds)*run_trace.indiv_trace[arch_no].H_secondary[ds];
	}


	// Estimate relative ratio of error terms
	mat H_hat = run_trace.H_consensus[arch_no];
	double a = 0.0, b = 0.0, x, y;
	for(int ds = 0; ds < ds_no; ds++) {
		mat W = run_trace.indiv_trace[arch_no].C_secondary[ds];
		mat H = run_trace.indiv_trace[arch_no].H_secondary[ds];

		x = norm(S[ds] - S[ds]*W*H, "fro");
		y = norm(H - H_hat, "fro");
		a += (x*x);
		b += (alpha(ds)*y*y);
	}	
	double ratio = a / b;
	lambda *= ratio;
	
	
	// Main loop
	for(int it = 0; it < max_it; it++) {

		// Permute rows
		for(int ds = 1; ds < ds_no; ds++) {
			mat G = 1 + cor(trans(run_trace.indiv_trace[arch_no].H_secondary[0]), trans(run_trace.indiv_trace[arch_no].H_secondary[ds]));
			mat G_matched = MWM(G);
			uvec perm = index_max(G_matched, 1);
		 			
			run_trace.indiv_trace[arch_no].C_secondary[ds] = run_trace.indiv_trace[arch_no].C_secondary[ds].cols(perm); 
			run_trace.indiv_trace[arch_no].H_secondary[ds] = run_trace.indiv_trace[arch_no].H_secondary[ds].rows(perm); 		
		}
		
		
		// Compute shared subspace 
		run_trace.H_consensus[arch_no] = zeros(size(run_trace.indiv_trace[arch_no].H_primary[0]));	
		for(int ds = 0; ds < ds_no; ds++) {
			run_trace.H_consensus[arch_no] += alpha(ds)*run_trace.indiv_trace[arch_no].H_secondary[ds];
		}		
	
		
		// Recompute H_i
		for(int ds = 0; ds < ds_no; ds++) {
			mat I = eye(arch_no, arch_no);
			mat Z = S[ds]*run_trace.indiv_trace[arch_no].C_secondary[ds]; // Archetype matrix
			double weight = lambda*alpha[ds];
			
			mat A = join_vert(trans(Z)*Z, weight*I);
			mat B = join_vert(trans(Z)*S[ds], weight*run_trace.H_consensus[arch_no]);			
			
			
		    run_trace.indiv_trace[arch_no].H_secondary[ds] = run_simplex_regression(A, B, false);
		}		
		
		
		// Recompute C_i
		for(int ds = 0; ds < ds_no; ds++) {		
			mat W = S[ds]*run_trace.indiv_trace[arch_no].C_secondary[ds];
			mat H = run_trace.indiv_trace[arch_no].H_secondary[ds];
			mat R = S[ds] - W*H;
			for(int j = 0; j < arch_no; j++) {
				double norm_sq = sum(square(H.row(j)));
				vec h = trans(H.row(j))/norm_sq;
				vec b = R*h + W.col(j);

				c = run_simplex_regression(S[ds], b, false);

				R += (W.col(j) - S[ds]*c)*H.row(j);
				W.col(j) = S[ds]*c;
				run_trace.indiv_trace[arch_no].C_secondary[ds].col(j) = c;			
			}			
		}				
	}
}

full_trace runACTION_muV(vector<mat> S_r, int k_min, int k_max, vec alpha, double lambda, int AA_iters = 50, int Opt_iters= 0, int numThreads = 8) {
	register int i, kk;
	printf("Running ACTION\n");

	double lambda2 = 1e-5, epsilon = 1e-5;
	
	k_min = std::max(k_min, 2);
	k_max = std::min(k_max, (int)(S_r[0].n_cols));	

	int cell_no = S_r[0].n_cols;
	vec c(cell_no);
	
	
	full_trace run_trace;		
	run_trace.H_consensus.resize(k_max+1);
	run_trace.indiv_trace.resize(k_max+1);
	for (kk = 0; kk <= k_max; kk++) {
		run_trace.indiv_trace[kk].selected_cols.resize(S_r.size());		
		run_trace.indiv_trace[kk].H_primary.resize(S_r.size());
		run_trace.indiv_trace[kk].C_primary.resize(S_r.size());
		run_trace.indiv_trace[kk].H_secondary.resize(S_r.size());
		run_trace.indiv_trace[kk].C_secondary.resize(S_r.size());
		run_trace.indiv_trace[kk].C_consensus.resize(S_r.size());
	}
	

		
	// Normalize signature profiles
	for(i = 0; i < S_r.size(); i++) {
		S_r[i] = normalise(S_r[i], 1, 0); // norm-1 normalize across columns -- particularly important for SPA
	}	
	
	field<mat> AA_res(2,1);	
	for(int kk = k_min; kk <= k_max; kk++) {
		printf("K = %d\n", kk);
		// Solve ACTION for a fixed-k to "jump-start" the joint optimization problem.
		for(i = 0; i < S_r.size(); i++) {
			printf("\tRun ACTION for dataset %d/%d\n", i+1, S_r.size());
			
			run_trace.indiv_trace[kk].selected_cols[i] = SPA(S_r[i], kk);
			
			mat W = S_r[i].cols(run_trace.indiv_trace[kk].selected_cols[i]);
			
			//AA_res = AA(X_r, W);
			AA_res = run_AA(S_r[i], W);		
			
			mat C0 = AA_res(0);
			C0.transform( [](double val) { return (min(1.0, max(0.0, val))); } );
			C0 = normalise(C0, 1);					
			run_trace.indiv_trace[kk].C_primary[i] = C0;
			
			
			mat H0 = AA_res(1);
			H0.transform( [](double val) { return (min(1.0, max(0.0, val))); } );
			H0 = normalise(H0, 1);	
			run_trace.indiv_trace[kk].H_primary[i] = H0;
		}

		
		// Compute consensus latent subspace, H^*
		findConsensus(S_r, run_trace, kk, alpha, lambda, Opt_iters); // sets secondary and consensus objects

		
		// decouple to find individual consensus C matrices
		for(i = 0; i < S_r.size(); i++) {	
			mat S = S_r[i];
			mat C = run_trace.indiv_trace[kk].C_secondary[i];
			mat H = run_trace.indiv_trace[kk].H_secondary[i];
			mat W = S*C;
			
			mat R = S - W*H;			

			run_trace.indiv_trace[kk].C_consensus[i] = zeros(cell_no, kk);
			for(int j = 0; j < kk; j++) {
				vec w = W.col(j);
				vec h = trans(H.row(j));

				double norm_sq = arma::dot(h, h);
				if (norm_sq < double(10e-8))
				{
				  // singular
				  int max_res_idx = index_max(rowvec(sum(square(R), 0)));
				  W.col(j) = S.col(max_res_idx);
				  c.zeros();
				  c(max_res_idx) = 1;
				  C.col(j) = c;
				}
				else
				{
				  // b = (1.0 / norm_sq) *R*ht + w;
				  vec b = w;
				  cblas_dgemv(CblasColMajor, CblasNoTrans, R.n_rows, R.n_cols,
							  (1.0 / norm_sq), R.memptr(), R.n_rows, h.memptr(), 1, 1,
							  b.memptr(), 1);

				  C.col(j) = run_simplex_regression(S, b, false);

				  vec w_new = S * C.col(j);
				  vec delta = (w - w_new);

				  // Rank-1 update: R += delta*h
				  cblas_dger(CblasColMajor, R.n_rows, R.n_cols, 1.0, delta.memptr(), 1,
							 h.memptr(), 1, R.memptr(), R.n_rows);

				  W.col(j) = w_new;
				}				

				run_trace.indiv_trace[kk].C_consensus[i].col(j) = C.col(j);
			}		
		}			
		
	}	

	return run_trace;
}
