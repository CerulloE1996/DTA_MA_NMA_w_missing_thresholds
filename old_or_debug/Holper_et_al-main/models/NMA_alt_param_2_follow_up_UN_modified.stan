// Model is from Nyaga et al - arm-based NMA model, "traditional" (non-copula) normal-binomial

// Code based on that provided from Nyaga et al. supp. material - with the some small modifications:
// using more efficient non-centered param for between-study model, uses LKJ corr priors, and calculation
// of parameters of interest in the generated quantities block, allows calculation of K-fold CV which
// is more reliable than WAIC and can be used when LOO fails (i.e. many pareto-k's > 0.7)
 

functions {  // make function to count # of unique elements in a vector - equiv to unique() in R. 
            int unique(int[] x) {
                int res[num_elements(x)] = sort_desc(x);   // Sort the elements of x in descending orde
                int count = 1;
                for (i in 2:num_elements(res))  {
                   if ( res[i] != res[i-1])  {
                     count = count + 1;
                    } else { 
                      count = count;
                    }
                }
                return(count); 
            }
}
                

data {
  
              int N; // number of obs. / rows 
              int n_tests; // total number of tests
              int n_test_thresholds;
              int n_thresholds[n_tests]; 
              int fu_summary;
              int TP[N];
              int FN[N];
              int FP[N];
              int TN[N];
              int Study[N];
              int Test[N];
              int Threshold[N];
              int Test_Threshold[N];
              vector[N] fu;
              vector[100] cts_cov_points; 
              real prior_intercept_thr_1_sd;
              real prior_intercept_thr_1_diseased_mean;
              real prior_intercept_thr_1_non_diseased_mean;
              real prior_intercept_thr_2_sd;
              real prior_intercept_thr_2_diseased_mean;
              real prior_intercept_thr_2_non_diseased_mean;
              real prior_intercept_other_thr_sd;
              real prior_coeff_sd;
              real prior_sigma_sd;
              real prior_tau_sd;
              vector<lower=0,upper=1>[5] prev; // prev for each test, to calculate population prevelance-dependent statistics e.g. PPV's and NPV's
              vector<lower=-1,upper=1>[1] youden_weight; // weight to calculate weighted youden index, 0.5 -> same as standard youden (equal weight on Se and Sp)
              int<lower=0,upper=1> holdout[N]; //index whether the observation should be used (for K-fold CV for model comparison)
              
}


transformed data {
  
              int max_Threshold;
              int n_studies; // total number of studies 
              int pos[N];
              int neg[N];
              
                for (n in 1:N) { 
                  pos[n] = TP[n] + FN[n];
                  neg[n] = TN[n] + FP[n];
                }
              
              max_Threshold = max(n_thresholds);
              n_studies = unique(Study);
              
}


parameters {
  
              vector[2] intercept[n_test_thresholds]; 
              cholesky_factor_corr[2] L_Omega; 
              vector<lower=0>[2] sigma; 
              matrix<lower=0>[n_tests, 2] tau; 
              matrix[n_studies, 2] z_nu;
              matrix[N, 2] z_tau;
              
}


transformed parameters {
  
        vector[2] logit_pi[N];
        vector[2] log_lik[N];
        matrix[N, 2] delta;
        vector[2] eta[n_studies];  
                    
          //// eta[s,] ~ normal(0, sigma).
          //// captures between-study var (aka the "study effect") which is SHARED BETWEEN TESTS, even if using UN (i.e. test-specific) structure 
          //// - hence NMA as using indirect evidence (borrowing info on other tests from other studies)
          for (s in 1:n_studies) {
             eta[s, ] =   diag_pre_multiply(sigma, L_Omega) * to_vector(z_nu[s, ]);
          }
                                                                                       
                                                                        
          for (n in 1:N) { //// likelihood is only implemented for the observed data points (hence loop through the N observations only)
    
                  for (d in 1:2) {
                         delta[n, d] = 0.0 + z_tau[n, d] * tau[Test[n], 1];
                  }

                  logit_pi[n, 1] = intercept[Test_Threshold[n], 1] +  eta[Study[n], 1] + delta[n, 1] ;
                  logit_pi[n, 2] = intercept[Test_Threshold[n], 2] +  eta[Study[n], 2] + delta[n, 2] ;
                       
                  //// Likelihood Model -  Pointwise (i.e. observation level) log-Likelihood - evaluated only at the observed data (missing data is imputed) hence looping through the N observations only 
                  log_lik[n, 1]  =  binomial_logit_lpmf(TP[n]  | pos[n] , logit_pi[n, 1] );
                  log_lik[n, 2]  =  binomial_logit_lpmf(TN[n]  | neg[n] , logit_pi[n, 2] );
          }
          
}


model {
    
        for (n in 1:N) {
           if ( Threshold[n] == 1) { 
             intercept[Test_Threshold[n],1] ~ normal(prior_intercept_thr_1_diseased_mean, prior_intercept_thr_1_sd); 
             intercept[Test_Threshold[n],2] ~ normal(prior_intercept_thr_1_non_diseased_mean, prior_intercept_thr_1_sd);
           }
           else if ( Threshold[n] == 2) { 
             intercept[Test_Threshold[n],1] ~ normal(prior_intercept_thr_2_diseased_mean, prior_intercept_thr_2_sd); 
             intercept[Test_Threshold[n],2] ~ normal(prior_intercept_thr_2_non_diseased_mean, prior_intercept_thr_2_sd);
          }
           else {                 
             intercept[Test_Threshold[n],1] ~ normal(0, prior_intercept_other_thr_sd); 
             intercept[Test_Threshold[n],2] ~ normal(0, prior_intercept_other_thr_sd); 
           }
        }
    
             to_vector(fu_coeff) ~ normal(0, prior_coeff_sd); 
           
        
    
             to_vector(tau) ~ normal(0, prior_tau_sd); 
             sigma ~ normal(0, prior_sigma_sd); 
             L_Omega ~ lkj_corr_cholesky(2);
             
             to_vector(z_nu)  ~ std_normal();
             to_vector(z_tau) ~ std_normal();
               
    
    
            // Likelihood Model
              for (n in 1:N) {
                if(holdout[n] == 0) {
                  target += log_lik[n,  1];
                  target += log_lik[n,  2];
                 }
              }
          
}

// 
//         // //// Declare Study-level random effects (eta in Nyaga notation) -  eta[s, 1:2] ~ multi_normal({0, 0}, Sigma):
//         // ////
//         // vector[n_studies] beta_eta;  // first eta corresponds to beta and 2nd eta corresponds to raw_scale.
//         // vector[n_studies] raw_scale_eta;  // first eta corresponds to beta and 2nd eta corresponds to raw_scale.
//         // ////
//         // //// Compute Study-level random effects (eta in Nyaga notation) -  eta[s, 1:2] ~ multi_normal({0, 0}, Sigma).
//         // ////
//         // for (s in 1:n_studies) {
//         //     beta_eta[s]      = beta_sigma      * beta_eta_z[s];      // NOTE: 1st eta's ("beta_eta") correspond to shared (between tests) component of "beta"
//         //     raw_scale_eta[s] = raw_scale_sigma * raw_scale_eta_z[s]; // NOTE: 2nd eta's ("raw_scale_eta") correspond to shared (between tests) component of "raw_scale"
//         // }
//         ////
//         //// Declare test-specific deviations ("delta" in Nyaga notation) - delta_{s, c, t} ~ normal(0, tau_{c, t}):
//         ////
//         // matrix[n_studies, n_tests] beta_delta;
//         // matrix[n_studies, n_tests] raw_scale_delta;
//         ////
//         //// Compute test-specific deviations ("delta" in Nyaga notation) - delta_{s, c, t} ~ normal(0, tau_{c, t}):
//         ////
//         for (t in 1:n_tests) {
//           for (s in 1:n_studies) {
//             if (test_in_study[s, t] == 1) {
//                   beta_delta[s, t]      = 0.0 + beta_delta_z[s, t]      * beta_tau[t];       ////  NOTE: 1st delta's ("beta_delta") correspond to shared (between tests) component of "beta"   -   delta_{s, c, t} ~ normal(0, tau_{c, t}):
//                   raw_scale_delta[s, t] = 0.0 + raw_scale_delta_z[s, t] * raw_scale_tau[t];  ////  NOTE: 2nd delta's ("raw_scale_delta") correspond to shared (between tests) component of "raw_scale" - delta_{s, c, t} ~ normal(0, tau_{c, t}):            }
//           }
//         }
//         
//         


generated quantities { 
  
                 matrix[n_tests, max_Threshold] Se; 
                 matrix[n_tests, max_Threshold] Sp;
                 vector[100] Se_range[n_tests, max_Threshold]; 
                 vector[100] Sp_range[n_tests, max_Threshold]; 
                 vector[N] se; 
                 vector[N] sp; 
                 matrix[n_tests, max_Threshold] lSe; 
                 matrix[n_tests, max_Threshold] lSp; 
                 vector[2] pred_nu;
                 matrix[n_tests,max_Threshold] pred_delta[2];
                 matrix[n_tests,max_Threshold] pred[2];
                 matrix[n_tests, max_Threshold] Se_pred; 
                 matrix[n_tests, max_Threshold] Sp_pred;  
                 matrix[n_tests, max_Threshold] lSe_pred; 
                 matrix[n_tests, max_Threshold] lSp_pred; 
                 corr_matrix[2] Omega; 
                 matrix[2,2] Sigma; 
                 matrix[n_tests, max_Threshold] Fp; 
                 matrix[n_tests, max_Threshold] LRpos;  
                 matrix[n_tests, max_Threshold] LRneg; 
                 matrix[n_tests, max_Threshold] DOR; 
                 matrix[n_tests, max_Threshold] PPV; 
                 matrix[n_tests, max_Threshold] NPV; 
                 matrix[n_tests, max_Threshold] Youden_index; 
                 matrix[n_tests, max_Threshold] Youden_index_weighted; 
                 matrix[n_tests,2] tausq;
                 vector[2] sigmabsq;
                 matrix[n_tests, n_tests] sigmasq[2];
                 matrix[n_tests, n_tests] rho[2];
                 matrix[n_tests, n_tests] rho12[2];

        // Between-study variance/cor parameters          
                  Omega = multiply_lower_tri_self_transpose(L_Omega); // correlation matrix (same across tests)
                  Sigma = quad_form_diag(Omega, sigma);               // var-cov matrix (same across tests)
                
                sigmabsq[1] = Sigma[1,1];
                sigmabsq[2] = Sigma[2,2];
                
        
            for (d in 1:2) 
               for (k in 1:n_tests)
                  tausq[k,d] = tau[k,d] * tau[k,d]; // for model w/ diff SD's across tests
          
                    
            for (d in 1:2) {
                for (k in 1:n_tests) {
                    for (l in 1:n_tests) {
                        sigmasq[d,k,l] = (sigmabsq[d] + tausq[k,d]) * ((sigmabsq[d] + tausq[l,d])); 
                        rho[d,k,l] = sigmabsq[d] / sqrt(sigmasq[d,k,l]);
                        // rho12 is the correlation between the t1'th and t2'th test (t1=t2 and t1 =/=t2 both possible) 
                       rho12[d,k,l] =      Omega[1,1]*sqrt(Sigma[1,1])*sqrt(Sigma[2,2]) /
                                           sqrt( (Sigma[1,1] + tausq[k,d]) * (Sigma[2,2] + tausq[l,d]) );
                    }
                }
            }
            
            
          pred_nu[1:2] = multi_normal_rng(rep_vector(0,2), Sigma);
         
          
         for (d in 1:2) {
            for (n in 1:N) { 
                     pred_delta[d,Test[n],Threshold[n]] =   normal_rng(0, tausq[Test[n],d]);
                     pred[d,Test[n],Threshold[n]] = intercept[Test_Threshold[n],d] + pred_nu[d] + pred_delta[d,Test[n],Threshold[n]];
             }
           }
        
        
         for (n in 1:N) { 
                Se[Test[n],Threshold[n]] = inv_logit(intercept[Test_Threshold[n],1] + fu_summary*fu_coeff[Test_Threshold[n],1]) ; 
                Sp[Test[n],Threshold[n]] = inv_logit(intercept[Test_Threshold[n],2] + fu_summary*fu_coeff[Test_Threshold[n],2]) ; 
                
                for (i in 1:100) { 
                  Se_range[Test[n],Threshold[n],i] = inv_logit(intercept[Test_Threshold[n],1] + cts_cov_points[i]*fu_coeff[Test_Threshold[n],1]) ; 
                  Sp_range[Test[n],Threshold[n],i] = inv_logit(intercept[Test_Threshold[n],2] + cts_cov_points[i]*fu_coeff[Test_Threshold[n],2]) ; 
                }
                
                
                
                lSe[Test[n],Threshold[n]] = intercept[Test_Threshold[n],1];
                lSp[Test[n],Threshold[n]] = intercept[Test_Threshold[n],2];
                lSe_pred[Test[n],Threshold[n]] = pred[1,Test[n],Threshold[n]];
                lSp_pred[Test[n],Threshold[n]] = pred[2,Test[n],Threshold[n]];
                Se_pred[Test[n],Threshold[n]] = inv_logit(lSe_pred[Test[n],Threshold[n]]);
                Sp_pred[Test[n],Threshold[n]] = inv_logit(lSp_pred[Test[n],Threshold[n]]);
                se[n] =     inv_logit(logit_pi[n, 1]) ; 
                sp[n] =     inv_logit(logit_pi[n, 2]) ; 
           }
           


        // Other summary non-comparative estimates 
                for (n in 1:N) { 
                          Fp[Test[n],Threshold[n]] = 1 - Sp[Test[n],Threshold[n]];
                          
                        // LR's
                         	LRpos[Test[n],Threshold[n]] = Se[Test[n],Threshold[n]]/(1-Sp[Test[n],Threshold[n]]);
          	              LRneg[Test[n],Threshold[n]] = (1-Se[Test[n],Threshold[n]])/Sp[Test[n],Threshold[n]];
          	              
          	            // PPV and NPV (using inputted prev data vector)
          	            PPV[Test[n],Threshold[n]] = Se[Test[n],Threshold[n]]*prev[Test[n]]/( Se[Test[n],Threshold[n]]*prev[Test[n]] + (1-Sp[Test[n],Threshold[n]])*(1-prev[Test[n]]) );
          	            NPV[Test[n],Threshold[n]] = Sp[Test[n],Threshold[n]]*prev[Test[n]]/( (1-Se[Test[n],Threshold[n]])*prev[Test[n]] + Sp[Test[n],Threshold[n]]*(1-prev[Test[n]]) );
                       }
	              
                 
        // Rank statistics
              for (n in 1:N) { 
                       
                        // Youden index
                        Youden_index[Test[n],Threshold[n]] =  Se[Test[n],Threshold[n]] + Sp[Test[n],Threshold[n]] - 1;
                        
                        // Weighted Youden index (weight = 0.5 -> same as youden, weight > 0.5 -> more importance on Se)
                        Youden_index_weighted[Test[n],Threshold[n]] = 2 * ( youden_weight[1]*Se[Test[n],Threshold[n]] + (1 - youden_weight[1])*Sp[Test[n],Threshold[n]] ) ;
                        
                        // DOR
                      	DOR[Test[n],Threshold[n]] = (Se[Test[n],Threshold[n]]*Sp[Test[n],Threshold[n]])/((1-Se[Test[n],Threshold[n]])*(1-Sp[Test[n],Threshold[n]]));
              }

        
}

