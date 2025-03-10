functions {
          // Function to do partial pooling on correlation matrices
          // See https://discourse.mc-stan.org/t/hierarchical-prior-for-partial-pooling-on-correlation-matrices/4852/27
          matrix convex_combine_cholesky(matrix global_chol_cor, matrix local_chol_cor, real alpha){
            int dim = rows(local_chol_cor);
            int global_dim = rows(global_chol_cor);
            matrix[global_dim,global_dim] global_cor = multiply_lower_tri_self_transpose(global_chol_cor);
            matrix[dim,dim] local_cor = multiply_lower_tri_self_transpose(local_chol_cor);
            matrix[dim,dim] global_cor_use;
            matrix[dim,dim] combined_chol_cor;
            
            if(dim < rows(global_cor)){
              global_cor_use = global_cor[1:dim,1:dim];
            } else {
              global_cor_use = global_cor;
            }
            combined_chol_cor = cholesky_decompose((1 - alpha)*global_cor_use + (alpha)*local_cor);
            return(combined_chol_cor);
  }
          real corr(int num, int[] vec1, int[] vec2) {
            real pearson_corr =      (num*sum(to_vector(vec1) .* to_vector(vec2)) -  sum(to_vector(vec1))*sum(to_vector(vec2)) )  /  
                              sqrt( (num*sum(to_vector(vec1) .* to_vector(vec1))   -  sum(to_vector(vec1))^2) * (num*sum(to_vector(vec2) .*  to_vector(vec2))   - sum(to_vector(vec2))^2 ) ); 
            return(pearson_corr);
          }
  real phi_logit_approx(real b) {
          real a;
          a = inv_logit( 1.702*b );
            return(a);
          }
  real inv_phi_logit_approx(real b) {
          real a;
          a = (1/1.702)*logit(b);    
            return(a);
          }

 real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}

data {
          int total_n;
          int n_binary_tests;
          int<lower=1> nt; // number of tests
          int<lower=1> n_studies;	// total # of studies
          int<lower=1> ns[n_studies]; // number of individuals in each study 
          int<lower=0> y[max(ns), nt, n_studies]; // N individuals and nt tests, n_studies studies 
          int n_patterns;
          int<lower=0> pa[max(ns), n_studies]; 
          int numg;
          int ns_cumsum[n_studies];
          int ind[n_studies];

          int n_thresholds;
          int n_thresholds_study[n_studies];
          int max_cutpoint;
          int length_cutpoint_vec;
          int cutpoints[n_studies, length_cutpoint_vec];

          int twoway_thr[n_studies];
          int<lower=0, upper=1> primary_care_ind[n_studies]; //indicator for whether the study comes from primary care or not
          int r[n_studies, 4, choose(nt,2), n_thresholds];
          int num_refs; // number of reference tests
          int ref[n_studies]; // ref. test used in each of the n_studies studies
          int re[n_studies];
          vector[500] thresh;
}

transformed data {
              int r2[choose(nt,2), n_studies, 4, n_thresholds];
      
       for (a in 1:4)
          for (s in 1:n_studies)
            for (thr in 1:n_thresholds)
              r2[1,s,a,thr] = r[s,a,1,thr];
}

parameters {
              vector[2] a1_m_raw[num_refs]; // Se > Fp for each reference test
              vector[2] a2_m_raw;   // don't order this due to potential thr. interction
              vector<lower=0>[2] sd1[num_refs];  
              vector<lower=0>[2] sd2; 
              real<lower=0,upper=1> u_d[total_n,nt]; // nuisance that absorbs inequality constraints
              real<lower=0,upper=1> u_nd[total_n,nt]; // nuisance that absorbs inequality constraints
              cholesky_factor_corr[2] L_Omega_bs1[num_refs];
              cholesky_factor_corr[2] L_Omega_bs2;
              vector[2] z1[n_studies];
              vector[2] z2[n_studies];
              vector[2] b_primary_care; // coefficients for primary care (1 for MMSE Se and 1 for MMSE Sp)
              real<lower=0, upper=1> p[n_studies]; 
              vector<lower=0>[n_thresholds+1] alpha;
              ordered[n_thresholds] C_d[n_studies];

              cholesky_factor_corr[nt] L_Omega_nd[n_studies];
              cholesky_factor_corr[nt] L_Omega_d[n_studies];

         //     real<lower=0,upper=1> alpha_d;
         //     real<lower=0,upper=1> alpha_nd;
         //     cholesky_factor_corr[nt] L_Omega_global_nd;
        //      cholesky_factor_corr[nt] L_Omega_global_d;
        //      cholesky_factor_corr[nt] R_diff_L_nd[n_studies];
        //      cholesky_factor_corr[nt] R_diff_L_d[n_studies];

}

transformed parameters {
              vector[2] mu[num_refs + 1];  
              matrix[2, nt] nu[n_studies];
              vector[total_n] log_lik; 
              vector[n_studies] lp1[n_patterns];
              vector[n_studies] lp0[n_patterns];
              vector[nt] z_d[n_patterns];
              vector[nt] z_nd[n_patterns]; 
              vector[nt] y1d[n_patterns];
              vector[nt] y1nd[n_patterns];
           //   cholesky_factor_corr[nt] L_Omega_nd[n_studies];
           //   cholesky_factor_corr[nt] L_Omega_d[n_studies];

                     //     for (s in 1:n_studies) {
                     //       L_Omega_d[s]  = convex_combine_cholesky(L_Omega_global_d,  R_diff_L_d[s], alpha_d);
                     //       L_Omega_nd[s] = convex_combine_cholesky(L_Omega_global_nd, R_diff_L_nd[s], alpha_nd); 
                     //     }

                          // index test (MMSE)
                          mu[num_refs + 1, 1] =  a2_m_raw[1];
                          mu[num_refs + 1, 2] =  a2_m_raw[2];
                          
          for (s in 1:n_studies) { 
                          // ref test (Se > Fp)
                          mu[ref[s], 1] =   a1_m_raw[ref[s], 2];       // Se  -- 26 - 99% CrI, median = 76
                          mu[ref[s], 2] =   a1_m_raw[ref[s], 1];       // Fp -- 88-99% CrI, median - 98%

                 // between-study model
                  if (ref[s] == 1 || ref[s] == 2) {
                         nu[s, 1:2 , 1] =  mu[ref[s], 1:2]   + diag_pre_multiply(sd1[ref[s],1:2], L_Omega_bs1[ref[s],]) * z1[s,1:2]; // (refmod 1 and 2)
                  }
                  else {
                        nu[s, 1:2 , 1] =  mu[3, 1:2]   + re[s]*diag_pre_multiply(sd1[3,1:2], L_Omega_bs1[3,]) * z1[s,1:2]; // refmod2
                  }
                      // index test random effect - w/ covariate for primary care
                   nu[s, 1:2 , 2] =  mu[num_refs + 1 , 1:2 ]   + 
                                     primary_care_ind[s]*b_primary_care[1:2]  + 
                                     diag_pre_multiply(sd2, L_Omega_bs2) * z2[s,1:2]; 
           }
                      // Parameters for likelihood function. Based on code upoaded by Ben Goodrich which uses the 
                      // GHK algorithm for generating TruncMVN. See:
                      // https://github.com/stan-dev/example-models/blob/master/misc/multivariate-probit/probit-multi-good.stan#L11
                      {
                 for (s in 1:n_studies) { 
                    for (n in 1:ns[s]) {
                     real prev_d = 0;
                     real prev_nd = 0;
                      for (i in 1:nt) { // loop over each test
                         real u_d1 =   u_d[n + ns_cumsum[s-1]*ind[s] ,i];
                         real u_nd1 = u_nd[n + ns_cumsum[s-1]*ind[s] ,i];
                           real bound_d_bin     = Phi_approx(    -  ( nu[s,1,1]   + prev_d )  / L_Omega_d[s,i,i]   ); // se
                           real bound_nd_bin    = Phi_approx(    -  ( nu[s,2,1]  + prev_nd )  / L_Omega_nd[s,i,i]  );  // fp
                           vector[ n_thresholds_study[s] + 1 ] bound_d;
                           vector[ n_thresholds_study[s] + 1 ] bound_nd;
                           bound_d[ n_thresholds_study[s] + 1 ] = 1;
                           bound_nd[ n_thresholds_study[s] + 1 ] = 1;

                           for (j in 1:(n_thresholds_study[s])) {
                             bound_d[j] =   Phi_approx(    (  C_d[s,cutpoints[s,j]]  -  ( nu[s,1,2]  +  prev_d ) ) /  L_Omega_d[s,i,i] ); // se
                             bound_nd[j] =  Phi_approx(    (  C_d[s,cutpoints[s,j]]  -  ( nu[s,2,2]  + prev_nd ) ) / L_Omega_nd[s,i,i] ); // fp
                          }
         if (i == 1) {
                                if (y[n,1,s] == 1) {
                                  z_d[pa[n,s],i]   = inv_Phi(bound_d_bin + (1 - bound_d_bin)*u_d1);      
                                  z_nd[pa[n,s],i]  = inv_Phi(bound_nd_bin + (1 - bound_nd_bin)*u_nd1);    
                                  y1d[pa[n,s],i]   = log1m(bound_d_bin);  // Jacobian adjustment
                                  y1nd[pa[n,s],i]  = log1m(bound_nd_bin); // Jacobian adjustment
                                }
                                if (y[n,1,s] == 0) { 
                                  z_d[pa[n,s],i]   = inv_Phi(bound_d_bin*u_d1);
                                  z_nd[pa[n,s],i]  = inv_Phi(bound_nd_bin*u_nd1);
                                  y1d[pa[n,s],i]   = log(bound_d_bin);  // Jacobian adjustment
                                  y1nd[pa[n,s],i]  = log(bound_nd_bin); // Jacobian adjustment
                                }
               }
          else {
              if ( n_thresholds_study[s] ==  1 )  {
                                if (  y[n,i,s] == 2 ) { 
                                  z_d[pa[n,s],i]   = inv_Phi(bound_d[1] + (1 - bound_d[1])*u_d1);  
                                  z_nd[pa[n,s],i]  = inv_Phi(bound_nd[1] + (1 - bound_nd[1])*u_nd1);
                                  y1d[pa[n,s],i]   = log1m(bound_d[1]); // Jacobian adjustment
                                  y1nd[pa[n,s],i]  = log1m(bound_nd[1]); // Jacobian adjustment
                                }
                                if ( y[n,i,s] == 1 ) {  
                                  z_d[pa[n,s],i]   = inv_Phi(bound_d[1]*u_d1);  
                                  z_nd[pa[n,s],i]  = inv_Phi(bound_nd[1]*u_nd1);
                                  y1d[pa[n,s],i]   = log(bound_d[1]); // Jacobian adjustment
                                  y1nd[pa[n,s],i]  = log(bound_nd[1]); // Jacobian adjustment
                                  }
                     }
              else { 
                if (y[n,i,s] == (n_thresholds_study[s] + 1) ) {
                      int K =   (n_thresholds_study[s] + 1) ; 
                        if ( y[n,i,s] == K ) { // Y from (num_thresh_params[s] + 1)  to 2
                                      z_d[pa[n,s],i]   = inv_Phi(bound_d[K-1] + (bound_d[K] - bound_d[K-1])*u_d1); 
                                      z_nd[pa[n,s],i]  = inv_Phi(bound_nd[K-1] + (bound_nd[K] - bound_nd[K-1])*u_nd1); 
                                      y1d[pa[n,s],i]   = log1m(bound_d[K-1]);  // Jacobian adjustment
                                      y1nd[pa[n,s],i]  = log1m(bound_nd[K-1]);  // Jacobian adjustment
                        }   
                }
                  for (J in 2:n_thresholds_study[s]) {
                    if ( y[n,i,s] == J ) { //  Y from 2 to num_thresh_params[s]  
                              z_d[pa[n,s],i]    = inv_Phi(bound_d[J-1] + (bound_d[J] - bound_d[J-1])*u_d1); 
                              z_nd[pa[n,s],i]   = inv_Phi(bound_nd[J-1] + (bound_nd[J] - bound_nd[J-1])*u_nd1); 
                              y1d[pa[n,s],i]    = log(bound_d[J] - bound_d[J-1]);  // Jacobian adjustment
                              y1nd[pa[n,s],i]   = log(bound_nd[J] - bound_nd[J-1]);  // Jacobian adjustment
                    }
                   }
                    if (y[n,i,s] == 1 ) {  
                      z_d[pa[n,s],i]   = inv_Phi(bound_d[1]*u_d1); 
                      z_nd[pa[n,s],i]  = inv_Phi(bound_nd[1]*u_nd1); 
                      y1d[pa[n,s],i]   = log(bound_d[1]); // Jacobian adjustment
                      y1nd[pa[n,s],i]  = log(bound_nd[1]); // Jacobian adjustment
                    }
                 }
             } 
                             if (i < nt)   prev_d    = L_Omega_d[s,i+1,1:i] * head(z_d[pa[n,s],] ,i);  
                             if (i < nt)   prev_nd  = L_Omega_nd[s,i+1,1:i] * head(z_nd[pa[n,s],] ,i);
          } // end i block
                            lp1[pa[n,s],s]  = sum(y1d[pa[n,s],])  +   bernoulli_lpmf(1 | p[s]); 
                            lp0[pa[n,s],s]  = sum(y1nd[pa[n,s],]) +   bernoulli_lpmf(0 | p[s]); 
                           log_lik[n + ns_cumsum[s-1]*ind[s]] =  log_sum_exp(lp1[pa[n,s],s] , lp0[pa[n,s],s] ); 
                        }
                       }
                     }
}

model {
        for (s in 1:n_studies) 
             if (primary_care_ind[s] == 0)  p[s] ~ beta(1,5); 

      b_primary_care ~ normal(0, 1);
  {
        real a = 0.75; 
        real b = 0.65; 
        a1_m_raw[1, 1]     ~ normal(-a, b); // ~ std_normal(); 
        a1_m_raw[1, 2]     ~ normal(a, b);  // ~ std_normal(); 

        a1_m_raw[2, 1]  ~ normal(-a, b); // ~ std_normal(); 
        a1_m_raw[2, 2]  ~ normal(a, b);  // ~ std_normal(); 

      for (t in 3:7) {
        a1_m_raw[t, 1]  ~ normal(-a, b); // ~ std_normal(); 
        a1_m_raw[t, 2]  ~ normal(a, b);  // ~ std_normal(); 
      }
  }
           alpha  ~ normal(1, 5);
       //    alpha  ~ normal(1, 10);

            for (s in 1:n_studies) 
                   C_d[s, ] ~   induced_dirichlet(alpha, 0);
          
       //   alpha_d  ~ beta(2,2);
       //   alpha_nd ~ beta(2,2);

       //   L_Omega_global_d  ~ lkj_corr_cholesky(4);
       //   L_Omega_global_nd ~ lkj_corr_cholesky(4);

         for (s in 1:n_studies) {
        //    R_diff_L_d[s,] ~ lkj_corr_cholesky(exp(alpha_d));
        //    R_diff_L_nd[s,] ~ lkj_corr_cholesky(exp(alpha_nd));
          L_Omega_d[s,]  ~ lkj_corr_cholesky(5);
          L_Omega_nd[s,] ~ lkj_corr_cholesky(5);
          }

            a2_m_raw  ~ normal(0, 1); 
            sd2   ~ normal(0, 0.60);  
            L_Omega_bs2 ~ lkj_corr_cholesky(2); 

        for (s in 1:n_studies) {
            sd1[ref[s],] ~ normal(0, 0.60);
            L_Omega_bs1[ref[s],] ~ lkj_corr_cholesky(2); 
            z1[s, ] ~ std_normal(); 
            z2[s, ] ~ std_normal(); 
        }
          
          for (s in 1:n_studies)  // likelihood
             for (n in 1:ns[s]) 
               target +=  log_lik[ n + ns_cumsum[s-1]*ind[s] ];
}

generated quantities {
  
      vector[num_refs] SeR; 
      vector[num_refs] SpR; 
      vector[n_thresholds] SeI; 
      vector[n_thresholds] SpI; 
      vector[nt] se[n_studies]; 
      vector[nt] sp[n_studies]; 
      vector[nt] fp[n_studies]; 
      vector[n_thresholds] sp2[n_studies];
      vector[n_thresholds] se2[n_studies];
      vector[n_thresholds] fp2[n_studies];
      // vector[n_studies] log_lik_study;
      vector[n_thresholds+1] p_dm_sim[numg];
      vector[n_thresholds+1] p_dm;
      vector[n_thresholds] C_dm;
      vector[n_thresholds] C_dm2;
      vector[500] fpr;
      vector[500] tpr;

     for (i in 1:500) { 
        fpr[i]  =        Phi_approx(    thresh[i] - mu[num_refs+1,2] ) ; 
        tpr[i]  =        Phi_approx(    thresh[i] - mu[num_refs+1,1] ) ; 
     }

    // log_lik_study[1] = sum(log_lik[1:ns[1]]);
    //  for (s in 2:n_studies)  {
    //      log_lik_study[s] =  sum(log_lik[ ns[s-1]:ns[s] ]); 
    //  }
                     
     for (i in 1:numg) { 
          p_dm_sim[i,]  =  dirichlet_rng(alpha);
     }
     for (i in 1:(n_thresholds+1)) {
         p_dm[i] = mean(p_dm_sim[,i]);
     }

       C_dm[1] =   0 - logit(1 - p_dm[1]); // logit
         
     for (i in 2:n_thresholds)  {
         C_dm[i] =    0 - logit(inv_logit(0 -  C_dm[i-1]) -  p_dm[i]); // logit
     }

       C_dm2[1] =   0 - logit(1 - p_dm_sim[1,1]); // logit 
           
     for (i in 2:n_thresholds)  {
           C_dm2[i] =    0 - logit(inv_logit(0 -  C_dm2[i-1]) -  p_dm_sim[1,i]); // logit
     }

     // Global Estimates
     for (s in 1:n_studies)  {
        SeR[ref[s]]  =  Phi_approx(    mu[ref[s],1]  );
        SpR[ref[s]]  =  Phi_approx(  - mu[ref[s],2]  );
     }
   
     for (i in 1:n_thresholds) { 
        SpI[i]  =   1 -  Phi_approx(   ( C_dm[i]  - mu[num_refs+1,2] ) ) ;
        SeI[i]  =        Phi_approx(   ( C_dm[i]  - mu[num_refs+1,1] ) ) ;
     }
 
 
     // study-specific accuracy estimates 
     // for index test
     for (s in  1:n_studies) {
         
            sp[s,2]   =    1 -   Phi_approx(  ( C_d[s, twoway_thr[s]]  - nu[s,2,2] ) );     
            se[s,2]   =          Phi_approx(  ( C_d[s, twoway_thr[s]]  - nu[s,1,2] ) );   
    
            for (thr in 1:n_thresholds) {
              sp2[s,thr] =  1 -     Phi_approx(  ( C_d[s, thr]  - nu[s,2,2] ) );  
              se2[s,thr] =          Phi_approx(  ( C_d[s, thr]  - nu[s,1,2] ) );
              fp2[s, ] = 1-sp2[s, ];
              }
     }
    // for ref test
    for (s in 1:n_studies) { 
        se[s,1]   =       Phi_approx(  nu[s,1,1] ); 
        sp[s,1]   =       Phi_approx( -nu[s,2,1] ); 
        fp[s, ] = 1-sp[s, ];
    
      
      }    
      
}

