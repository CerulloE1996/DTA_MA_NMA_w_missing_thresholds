
//// NOTE: This verion is more suitable for "real world" use than the other "FIXED_thr" (fixed thresholds)
//// version because it allows the Sp's to vary between studies eventhough the thresholds are still fixed 
//// between studies. 
//// This is because we have introduced mean parameters for the within-study ordinal-probit models in
//// BOTH groups (i.e. not only in the diseased group) - we can do this without identifiability issues 
//// - especially in our case where the within-study models are:
//// (i) univariate, and:
//// (ii) NOT latent class (i.e. assuming a perfect gold standard) 
//// - this makes identification of parameters (in general) much easier. 
////
//// Also, please see the Michael Betancourt case study to understand why - particulaly in the Bayesian case 
//// and especially when we use an Induced Dirichlet model (or prior model) - that identifiability is not as
//// much of an issue and that we can freely estimate means AND cutpoint parameters in BOTH groups, link is here: 
//// https://betanalpha.github.io/assets/case_studies/ordinal_regression.html
//// Specifically, in this case study he he demonstrates that whilst there is indeed a theoretical 
//// location/mean "non-identifiability"" (since shifting means and cutpoints by a constant preserves the 
//// ordinal probabilities), in practice having an informative prior on one set of parameters
//// (e.g. the mean) sufficiently anchors the model - especially one based on the Induced Dirichlet prior model
//// (which allows us to set priors * directly * on the induced ordinal probabilities). 
////
//// Note that this cannot be said about frequentist versions of our model nor can it be said about models
//// which naively use so-called "vague" priors. 



functions {
  
        //// Induced-Dirichlet ("ind_dir") log-density function:
        //// NOTE: You can use this for both ind_dir PRIORS and ind_dir MODELS:
        real induced_dirichlet_v2_lpdf(  vector p_ord,
                                         vector p_cumul,
                                         vector alpha) {
                  
                int n_cat = num_elements(p_ord);
                matrix[n_cat, n_cat] J = rep_matrix(0.0, n_cat, n_cat);
                
                //// Jacobian computation:
                for (k in 1:n_cat) {
                        J[k, 1] = 1.0;
                }
                for (k in 2:n_cat) {
                        real rho = exp(std_normal_lpdf(inv_Phi( p_cumul[k - 1]))); //  p_cumul[k - 1] * (1.0 - p_cumul[k - 1]); // original
                        J[k, k] = +rho;
                        J[k - 1, k] = -rho;
                }
                
                return dirichlet_lpdf(p_ord | alpha) + log_determinant(J);
          
        }
        
        real median(vector x) {
      
                int n = num_elements(x);
                vector[n] sorted_x = sort_asc(x);
                
                if (n % 2 == 0) {   // For even number of elements, average the middle two
                    int left_element = to_int(n/2.0);
                    int right_element = to_int(n/2.0 + 1.0);
                    return (sorted_x[left_element] + sorted_x[right_element]) / 2.0;
                } else {            // For odd number of elements, return the middle one
                    int middle_element = to_int((n+1.0)/2.0);
                    return sorted_x[middle_element];
                }
       }
       
       real median(row_vector x) {
      
                int n = num_elements(x);
                row_vector[n] sorted_x = sort_asc(x);
                
                if (n % 2 == 0) {   // For even number of elements, average the middle two
                    int left_element = to_int(n/2.0);
                    int right_element = to_int(n/2.0 + 1.0);
                    return (sorted_x[left_element] + sorted_x[right_element]) / 2.0;
                } else {    // For odd number of elements, return the middle one
                    int middle_element = to_int((n+1.0)/2.0);
                    return sorted_x[middle_element];
                }
       }

}




data {
    
        int<lower=1> n_studies;
        int<lower=1> n_thr;
        array[2, n_studies] int n_cutpoints;
        array[2] matrix[n_studies, n_thr] x_with_missings;
        array[2] matrix[n_studies, n_thr] n;
        array[2] matrix[n_studies, n_thr] x;
        array[2] matrix[n_studies, n_thr] cutpoint_index;
        ////
        array[2] matrix[n_studies, n_thr + 1] x_cat;
        //// Priors:
        vector<lower=0>[n_thr + 1] prior_alpha; // Induced-Dirichlet prior vector 
        vector<lower=0.0>[2] prior_kappa_mean;
        vector<lower=0.0>[2] prior_kappa_SD;
        vector[2] prior_beta_mu_mean;
        vector[2] prior_beta_mu_SD;
        vector[2] prior_beta_SD_mean;
        vector[2] prior_beta_SD_SD;
        vector[2] prior_log_scale_mu_mean;
        vector[2] prior_log_scale_mu_SD;
        vector[2] prior_log_scale_SD_mean;
        vector[2] prior_log_scale_SD_SD;
        //// Other:
        int prior_only;
        int estimate_scales;
        real log_alpha_lb;
        int use_empirical_cutpoint_means;
  
}

transformed data { 
    
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        int n_sets_of_cutpoints = 1;
        real alpha_lb = exp(log_alpha_lb);
}
  


parameters {
  
        vector[2] beta_mu;
        vector<lower=0.0>[2] beta_SD;
        matrix[2, n_studies] beta_z;    // Study-specific random effects for beta (off-centered parameterisation)
        cholesky_factor_corr[2] beta_L_Omega;
        ////
        array[n_studies] ordered[n_thr] C_nd;  // Global cutpoints
        ////
        matrix<lower=alpha_lb>[n_sets_of_cutpoints, n_cat] alpha; // For induced-Dirichlet hierarchical model
        
}


transformed parameters { 
    
        matrix[2, n_studies] beta;
        cholesky_factor_cov[2] beta_L_Sigma;
        ////
        array[2] matrix[n_studies, n_thr] logit_cumul_prob; // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cumul_prob;  // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cond_prob;   // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] log_lik;     // log_lik storage
        ////
        array[n_sets_of_cutpoints] matrix[n_studies, n_thr] Ind_Dir_cumul_prob;
        array[n_sets_of_cutpoints] matrix[n_studies, n_cat] Ind_Dir_ord_prob;
        ////
        array[n_sets_of_cutpoints] matrix[n_studies, n_thr] C;
        ////
        ////
        //// Initialise matrices:
        ////
        for (c in 1:2) {
                 logit_cumul_prob[c]   = rep_matrix(positive_infinity(), n_studies, n_thr);
                 cumul_prob[c]         = rep_matrix(1.0, n_studies, n_thr);
                 cond_prob[c]          = rep_matrix(1.0, n_studies, n_thr);
                 log_lik[c]            = rep_matrix(0.0, n_studies, n_thr);
        }
        for (c in 1:n_sets_of_cutpoints) {
                 Ind_Dir_cumul_prob[c] = rep_matrix(0.0, n_studies, n_thr);
                 Ind_Dir_ord_prob[c]   = rep_matrix(0.0, n_studies, n_cat);
        }
        ////
        //// Between-study model for the location parameters ("beta") - models between-study correlation:
        ////
        beta_L_Sigma = diag_pre_multiply(beta_SD, beta_L_Omega);
        for (s in 1:n_studies) {
             beta[, s] = beta_mu + beta_L_Sigma * beta_z[, s];
             // beta[2, s] = 0.0; //////// FOR DEBUGGING
        }
        ////
        //// Cutpoints:
        ////
        for (s in 1:n_studies) {
               for (k in 1:n_thr) {
                    C[1][s, k] = C_nd[s][k];
               }
        }
        ////
        //// Induced-Dirichlet model cumulative probs:
        ////
        for (c in 1:n_sets_of_cutpoints) {
                real ID_anchor = 0.0;
                Ind_Dir_cumul_prob[c] = Phi(C[c] - ID_anchor);
                for (s in 1:n_studies) {
                        //// Induced-Dirichlet ordinal probs:
                        Ind_Dir_ord_prob[c][s, 1] = Ind_Dir_cumul_prob[c][s, 1] - 0.0;
                        for (k in 2:n_thr) {
                           Ind_Dir_ord_prob[c][s, k] = Ind_Dir_cumul_prob[c][s, k] - Ind_Dir_cumul_prob[c][s, k - 1]; // since cutpoints are increasing with k
                        }
                        Ind_Dir_ord_prob[c][s, n_cat] = 1.0 - Ind_Dir_cumul_prob[c][s, n_cat - 1];
                }
        }
        ////
        //// Likelihood using binomial factorization:
        ////
        for (s in 1:n_studies) {
                  for (c in 1:2) {
                      for (cut_i in 1:to_int(n_cutpoints[c, s])) {
                             int k = to_int(cutpoint_index[c][s, cut_i]);
                             real cutpoint = C[1][s, k];
                             if (n_sets_of_cutpoints == 2) cutpoint = C[c][s, k];
                             logit_cumul_prob[c][s, cut_i] = (cutpoint - beta[c, s]);
                      }
                  }
        }
        ////
        //// Calculate CUMULATIVE probabilities (vectorised):
        ////
        for (c in 1:2) {
            cumul_prob[c] = Phi(logit_cumul_prob[c]); //// INCREASING sequence (as C_k > C_{k - 1})
        }
        ////
        //// ------- Binomial likelihood:
        ////
        for (s in 1:n_studies) {
                for (c in 1:2) {
                      for (cut_i in 1:n_cutpoints[c, s]) {
                              // Current and next cumulative counts
                              int x_current = to_int(x[c][s, cut_i]);
                              int x_next    = to_int(n[c][s, cut_i]);
                              
                              // Skip if the current count is zero (no observations to classify)
                              if (x_current != 0)  {
                              
                                    // Conditional probability of being at or below the current cutpoint - given being at or below the next cutpoint
                                    if (cut_i == n_cutpoints[c, s]) { 
                                             cond_prob[c][s, cut_i]  = cumul_prob[c][s, cut_i] / 1.0;
                                    } else {
                                          if (x_next > 0) { 
                                             cond_prob[c][s, cut_i] = cumul_prob[c][s, cut_i] / cumul_prob[c][s, cut_i + 1];
                                          } else { 
                                             cond_prob[c][s, cut_i] = 1.0;
                                          }
                                    }
                                    
                                    // Binomial for observations at or below current cutpoint out of those at or below next cutpoint
                                    log_lik[c][s, cut_i] = binomial_lpmf(x_current | x_next, cond_prob[c][s, cut_i]);
                                    
                              }
                      }
                }
        }
}


model {
      
        // Priors:
        {
            //// locations:
            beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
            beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
            beta_L_Omega ~ lkj_corr_cholesky(2.0);
        }
        //// Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        for (c in 1:n_sets_of_cutpoints) {
            target += normal_lpdf(alpha[c, 1:n_cat] | prior_kappa_mean[c], prior_kappa_SD[c]); // putting prior DIRECTLY on parameter so NO JACOBIAN needed here !!
        }
        for (c in 1:n_sets_of_cutpoints) {
          for (s in 1:n_studies) {
            target += induced_dirichlet_v2_lpdf(to_vector(Ind_Dir_ord_prob[c][s, 1:n_cat]) | to_vector(Ind_Dir_cumul_prob[c][s, 1:n_thr]), to_vector(alpha[c, 1:n_cat]));
          }
        }
        ////
        // Likelihood / Model:
        {
            to_vector(beta_z) ~ std_normal();   // (part of between-study model, NOT prior)
        }
        //// Likelihood using binomial factorization:
        if (prior_only == 0) {
            for (c in 1:2) {
              target +=  sum(log_lik[c]);
            }
        }
  
}




 

generated quantities {

          //// Summary accuracy parameters (at each threshold):
          vector[n_thr] Se;
          vector[n_thr] Sp;
          vector[n_thr] Fp;
          vector[n_thr] Se_pred;
          vector[n_thr] Sp_pred;
          vector[n_thr] Fp_pred;
          matrix[n_studies, n_thr] se;
          matrix[n_studies, n_thr] sp;
          matrix[n_studies, n_thr] fp;
          matrix[n_studies, n_thr] x_hat_nd = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] x_hat_d  = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_nd   = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_d    = rep_matrix(-1, n_studies, n_thr);
          array[2] matrix[n_studies, n_thr] x_hat;
          array[2] matrix[n_studies, n_thr] dev;
          matrix[n_sets_of_cutpoints, n_thr] C_mu;       
          matrix[n_sets_of_cutpoints, n_thr] C_mu_empirical; 
          matrix[n_sets_of_cutpoints, n_thr] prob_cumul_mu;  
          matrix[n_sets_of_cutpoints, n_cat] prob_ord_mu;  
          corr_matrix[2] beta_Omega;
          int n_sims = 1000;
             
          //// Compute between-study correlation matrix for location parameters:
          beta_Omega  =  multiply_lower_tri_self_transpose(beta_L_Omega);
  
          //// Initialise containers:
          for (c in 1:2) {
              x_hat[c]  = rep_matrix(-1, n_studies, n_thr);
              dev[c]    = rep_matrix(-1, n_studies, n_thr); 
          }
          
          for (c in 1:n_sets_of_cutpoints) {

                // Expected value of Dirichlet is proportional to alpha parameters:
                vector[n_cat] expected_p_ord;
                for (k in 1:n_cat) {
                    expected_p_ord[k] = alpha[c][k] / sum(alpha[c]);
                }

                // Calculate cumulative probabilities:
                vector[n_thr] expected_p_cumul;
                expected_p_cumul[1] = expected_p_ord[1];
                for (k in 2:n_thr) {
                    expected_p_cumul[k] = expected_p_cumul[k - 1] + expected_p_ord[k];
                }

                // Transform to cutpoints:
                for (k in 1:n_thr) {
                       real prob = expected_p_cumul[k];
                      if (prob < 0.0001) {
                            C_mu[c, k] = -10.0;
                      } else if (prob > 0.99999) {
                            C_mu[c, k] = +10.0;
                      } else {
                            C_mu[c, k] = 0.0 + inv_Phi(prob);
                      }
          
                }
          }
          
          // {
          // 
          //      array[n_sets_of_cutpoints] matrix[n_cat, n_sims] prob_ord_mu_sim;
          //      array[n_sets_of_cutpoints] matrix[n_thr, n_sims] prob_cumul_mu_sim;
          //      array[n_sets_of_cutpoints] matrix[n_thr, n_sims] C_mu_sim;
          // 
          //      for (c in 1:n_sets_of_cutpoints) {
          //              for (i in 1:n_sims) {
          // 
          //                             // // //// Induced-Dirichlet ordinal probs are structured like this:
          //                             //  Ind_Dir_cumul_prob[c] = Phi(C[c] - Ind_Dir_anchor[c]);
          //                             //  for (s in 1:n_studies) {
          //                             //          //// Induced-Dirichlet ordinal probs:
          //                             //          Ind_Dir_ord_prob[c][s, 1] = Ind_Dir_cumul_prob[c][s, 1] - 0.0;
          //                             //          for (k in 2:n_thr) {
          //                             //             Ind_Dir_ord_prob[c][s, k] = Ind_Dir_cumul_prob[c][s, k] - Ind_Dir_cumul_prob[c][s, k - 1]; // since cutpoints are increasing with k
          //                             //          }
          //                             //          Ind_Dir_ord_prob[c][s, n_cat] = 1.0 - Ind_Dir_cumul_prob[c][s, n_cat - 1];
          //                             //  }
          //                             // Simulate from Dirichlet by using the summary "alpha" parameters:
          //                             prob_ord_mu_sim[c][1:n_cat, i]   =  dirichlet_rng(to_vector(alpha[c, ]));
          //                             ////
          //                             //// Calculate cumulative probabilities:
          //                             prob_cumul_mu_sim[c][1, i] = prob_ord_mu_sim[c][1, i];
          //                             for (k in 2:n_thr) {
          //                                 prob_cumul_mu_sim[c][k, i] = prob_cumul_mu_sim[c][k - 1, i] + prob_ord_mu_sim[c][k, i];
          //                             }
          //                             ////
          //                             //// Transform to cutpoints:
          //                             real anchor_for_summaries = 0.0;
          //                             for (k in 1:n_thr) {
          //                                   real prob_1;
          //                                   if (prob_cumul_mu_sim[c][k, i] < 1e-38) { 
          //                                         prob_1 = 1e-38;
          //                                         C_mu_sim[c][k, i] = -10.0;
          //                                   } else if (prob_cumul_mu_sim[c][k, i] > 0.9999999999999) { 
          //                                         prob_1 = 0.9999999999999;
          //                                         C_mu_sim[c][k, i] = +10.0;
          //                                   } else {
          //                                         prob_1 = inv_Phi(prob_cumul_mu_sim[c][k, i]);
          //                                         C_mu_sim[c][k, i] = anchor_for_summaries + prob_1;
          //                                   }
          //                              
          //                             }
          // 
          //              }
          //      }
          // 
          // 
          //      for (c in 1:n_sets_of_cutpoints) {
          // 
          //            for (k in 1:n_thr) {
          //                prob_ord_mu[c, k] = median(prob_ord_mu_sim[c][k, 1:n_sims]);
          //                C_mu[c, k]  = median(C_mu_sim[c][k, 1:n_sims]);
          //                prob_cumul_mu[c, k] = median(prob_cumul_mu_sim[c][k, 1:n_sims]);
          //            }
          //            prob_ord_mu[c, n_cat] = median(prob_ord_mu_sim[c][n_cat, 1:n_sims]);
          //            // prob_cumul_mu[c, n_cat] = median(prob_cumul_mu_sim[c][n_cat, ]);
          // 
          //      }
          // 
          // }
          ////
          //// Empirical-mean cutpoints:     
          for (c in 1:n_sets_of_cutpoints) {
              for (k in 1:n_thr) {
                    C_mu_empirical[c, k] = median(C[c][, k]);
              }
          }
          //// Calculate study-specific accuracy:
          for (s in 1:n_studies) {
                for (k in 1:n_thr) { 
                    fp[s, k] =   1.0 - cumul_prob[1][s, k];
                    sp[s, k] =   1.0 - fp[s, k];
                    se[s, k] =   1.0 - cumul_prob[2][s, k];
                }
          }
          //// Calculate summary accuracy (using mean parameters):
          for (k in 1:n_thr) {
                Fp[k] =   1.0 - Phi(C_mu[1, k] - beta_mu[1]);
                Sp[k] =   1.0 - Fp[k];
                Se[k] =   1.0 - Phi(C_mu[1, k] - beta_mu[2]);
          }
          //// Calculate predictive accuracy:
          {
                vector[2] beta_pred =  to_vector(multi_normal_cholesky_rng(beta_mu, beta_L_Sigma));
                for (k in 1:n_thr) {
                      Fp_pred[k] =   1.0 - Phi(C_mu[1, k] - beta_pred[1]);
                      Sp_pred[k] =   1.0 - Fp_pred[k];
                      Se_pred[k] =   1.0 - Phi(C_mu[1, k] - beta_pred[2]);
                }
          }
          //// Model-predicted ("re-constructed") data:
          for (s in 1:n_studies) {
              for (c in 1:2) {
                 for (cut_i in 1:to_int(n_cutpoints[c, s])) {
                  
                      //// Model-estimated data:
                      x_hat[c][s, cut_i] = cond_prob[c][s, cut_i] * n[c][s, cut_i];  	 // Fitted values
                    
                      //// Compute residual deviance contribution:
                      real n_i =  (n[c][s, cut_i]);
                      real x_i =  (x[c][s, cut_i]);
                      real x_hat_i =  (x_hat[c][s, cut_i]);
                      real log_x_minus_log_x_hat = log(x_i) - log(x_hat_i);
                      real log_diff_n_minus_x = log(n_i - x_i);
                      real log_diff_n_minus_x_hat = log(abs(n_i - x_hat_i));
                       
                      dev[c][s, cut_i] = 2.0 * ( x_i * log_x_minus_log_x_hat + (n_i - x_i) * (log_diff_n_minus_x - log_diff_n_minus_x_hat) ); 
                      
                 }
              }
          }
          
         x_hat_nd = x_hat[1];
         dev_nd = dev[1];
         x_hat_d = x_hat[2];
         dev_d = dev[2];
       
}

















