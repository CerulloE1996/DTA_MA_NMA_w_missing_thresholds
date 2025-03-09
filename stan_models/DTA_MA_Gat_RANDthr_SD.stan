


functions {
  
          // Normal PDF for a real input
          real normal_pdf( real x, 
                           real mu, 
                           real sigma) {
                             
               real sqrt_2_pi = sqrt(2 * pi());
               return (1.0 / (sigma * sqrt_2_pi)) * exp(-0.5 * ((x - mu) / sigma)^2);
            
          }
        
        
          // Normal PDF for a vector input
          vector normal_pdf( vector x, 
                             real mu, 
                             real sigma) {

                real sqrt_2_pi = sqrt(2 * pi());
                return (1.0 / (sigma .* sqrt_2_pi)) * exp(-0.5 * ((x - mu) ./ sigma)^2);

          }
          
          // Overload for vector mean, vector sigma
          vector normal_pdf( vector x, 
                             vector mu, 
                             vector sigma) {
            
                real sqrt_2_pi = sqrt(2 * pi());
                return (1.0 / (sigma * sqrt_2_pi)) .* exp(-0.5 * ((x - mu) ./ sigma)^2);
          
        }
              //// Induced-Dirichlet ("ind_dir") log-density function:
        //// NOTE: You can use this for both ind_dir PRIORS and ind_dir MODELS:
        real induced_dirichlet_v2_lpdf(  vector p_ord,
                                         vector rho,
                                         vector alpha) {
                  
                int n_cat = num_elements(p_ord);
                matrix[n_cat, n_cat] J = rep_matrix(0.0, n_cat, n_cat);
                
                //// Jacobian computation:
                for (k in 1:n_cat) {
                        J[k, 1] = 1.0;
                }
                for (k in 2:n_cat) {
                        // real rho =  normal_pdf(inv_Phi( p_cumul[k - 1])); //  p_cumul[k - 1] * (1.0 - p_cumul[k - 1]); // original
                        J[k, k] = +rho[k - 1];
                        J[k - 1, k] = -rho[k - 1];
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
                    int middle_element = to_int((n + 1.0)/2.0);
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
        //// Priors:
        real prior_beta_mu_mean;
        real prior_beta_mu_SD;
        real prior_beta_SD_mean;
        real prior_beta_SD_SD;
        //// Priors for raw_scale:
        real prior_raw_scale_mu_mean;
        real prior_raw_scale_mu_SD;
        real prior_raw_scale_SD_mean;
        real prior_raw_scale_SD_SD;
        //// Induced-Dirichlet priors:
        vector[n_thr + 1] prior_dirichlet_cat_means_alpha;
        vector[n_thr + 1] prior_dirichlet_cat_SDs_mean;
        vector<lower=0.0>[n_thr + 1] prior_dirichlet_cat_SDs_SD;
        //// Other:
        int prior_only;
        real kappa_lb;
  
}

transformed data { 
    
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        
}
  


parameters {
  
        real beta_mu;    
        real<lower=0.0> beta_SD;   
        vector[n_studies] beta_z;
        ////
        real raw_scale_mu;    
        real<lower=0.0> raw_scale_SD;   
        vector[n_studies] raw_scale_z;
        ////
        array[n_studies] ordered[n_thr] C_nd; // study-specific cutpoints
        simplex[n_cat] dirichlet_cat_means_phi;
        real<lower=kappa_lb> kappa;  
        
}


transformed parameters { 
  
        matrix[2, n_studies] locations;
        matrix[2, n_studies] scales;
        ////
        array[2] matrix[n_studies, n_thr] logit_cumul_prob; // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cumul_prob;  // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cond_prob;   // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] log_lik;     // log_lik storage
        ////
        matrix[n_studies, n_thr] Ind_Dir_anchor;
        matrix[n_studies, n_thr] Ind_Dir_cumul_prob;
        matrix[n_studies, n_cat] Ind_Dir_ord_prob;
        real anchor = 0.0;
        ////
        matrix[n_studies, n_thr] C;
        real log_kappa = log(kappa);
        vector[n_cat] log_dirichlet_cat_means_phi = log(dirichlet_cat_means_phi);
        ////
        //// Compute alpha:
        ////
        vector[n_cat] log_alpha = log_kappa + log_dirichlet_cat_means_phi;
        vector[n_cat] alpha = exp(log_alpha);
        real alpha_0 = sum(alpha);
        ////
        //// Other quantities needed (for Jacobian for the (alpha_k, alpha_0) -> dirichlet_cat_SDs_sigma transformation):
        ////
        real Jacobian_for_alpha_k_to_category_SDs = 0.0;
        real log_half = log(0.5);
        real alpha_0_sq = square(alpha_0);
        real alpha_0_cube = alpha_0_sq * alpha_0;
        ////
        //// NOTE: SD for each category probability in a Dirichlet is √(α_k(α_0-α_k)/(α_0^2(α_0+1))) - where α_0 is the sum of all alphas:
        ////
        vector[n_cat] dirichlet_cat_SDs_sigma = sqrt((alpha .* (alpha_0 - alpha) ) ./ (alpha_0_sq * (alpha_0 + 1.0))); // We are putting a prior on this!!
        // vector[n_cat] Phi_dirichlet_cat_SDs_sigma = Phi(dirichlet_cat_SDs_sigma); // We are putting a prior on this!!
        ////
        //// Then compute the Jacobian adjustment:
        ////
        {
            vector[n_cat] log_sigma = log(dirichlet_cat_SDs_sigma);
            vector[n_cat] alpha_sq = square(alpha);
            vector[n_cat] deriv_A = alpha_sq + alpha_0 - 2.0*alpha;
            vector[n_cat] deriv_B = (1.0 ./ alpha) .* (3.0*alpha_0_sq*alpha + 2.0*alpha*alpha_0);
            vector[n_cat] deriv_var  = deriv_A .* (alpha_0_cube + alpha_0_sq) + deriv_B .* (alpha*alpha_0 - alpha_sq);
            vector[n_cat] log_abs_deriv_var_wrt_alpha = log(abs(deriv_var));
            vector[n_cat] log_abs_deriv_SD_wrt_alpha = log_half - log_sigma + log_abs_deriv_var_wrt_alpha;
            Jacobian_for_alpha_k_to_category_SDs += sum(log_abs_deriv_SD_wrt_alpha);
        }
        ////
        //// Initialise matrices:
        ////
        for (c in 1:2) {
                 logit_cumul_prob[c]   = rep_matrix(positive_infinity(), n_studies, n_thr);
                 cumul_prob[c]         = rep_matrix(1.0, n_studies, n_thr);
                 cond_prob[c]          = rep_matrix(1.0, n_studies, n_thr);
                 log_lik[c]            = rep_matrix(0.0, n_studies, n_thr);
        }
        {
           Ind_Dir_cumul_prob = rep_matrix(0.0, n_studies, n_thr);
           Ind_Dir_ord_prob   = rep_matrix(0.0, n_studies, n_cat);
           Ind_Dir_anchor     = rep_matrix(anchor, n_studies, n_thr);
        }
        ////
        //// Between-study model for location and scale:
        ////
        {
            ////
            for (s in 1:n_studies) {
                locations[1, s] = -0.5*(beta_mu      + beta_SD      * beta_z[s]); // = -0.5*(beta_mu + beta_SD * beta_z[s])
                locations[2, s] = +0.5*(beta_mu      + beta_SD      * beta_z[s]);
                scales[1, s] = log1p_exp(-0.5*(raw_scale_mu + raw_scale_SD * raw_scale_z[s]));
                scales[2, s] = log1p_exp(+0.5*(raw_scale_mu + raw_scale_SD * raw_scale_z[s]));
            }
        }
        // ////
        // //// Between-study model for the location parameters ("beta") - models between-study correlation:
        // ////
        // {
        //     matrix[2, n_studies] bs_mat;
        //     for (s in 1:n_studies) {
        //          bs_mat[1, s] = beta_mu      + beta_SD      * beta_z[s];       // beta_{s}     ~ normal(beta_mu,      beta_SD)
        //          bs_mat[2, s] = raw_scale_mu + raw_scale_SD * raw_scale_z[s];  // raw_scale{s} ~ normal(raw_scale_mu, raw_scale_SD)
        //     }
        //     ////
        //     for (s in 1:n_studies) {
        //         locations[1, s] = -0.5*bs_mat[1, s];
        //         locations[2, s] = +0.5*bs_mat[1, s];
        //         scales[1, s] = log1p_exp(-0.5*bs_mat[2, s]);
        //         scales[2, s] = log1p_exp(+0.5*bs_mat[2, s]);
        //     }
        // }
        ////
        //// Cutpoints:
        ////
        for (s in 1:n_studies) {
               for (k in 1:n_thr) {
                    C[s, k] = C_nd[s][k];
               }
        }
        ////
        //// Induced-Dirichlet model cumulative probs:
        ////
        {
                Ind_Dir_cumul_prob = Phi(C - Ind_Dir_anchor);
                for (s in 1:n_studies) {
                        //// Induced-Dirichlet ordinal probs:
                        Ind_Dir_ord_prob[s, 1] = Ind_Dir_cumul_prob[s, 1] - 0.0;
                        for (k in 2:n_thr) {
                           Ind_Dir_ord_prob[s, k] = Ind_Dir_cumul_prob[s, k] - Ind_Dir_cumul_prob[s, k - 1]; // since cutpoints are increasing with k
                        }
                        Ind_Dir_ord_prob[s, n_cat] = 1.0 - Ind_Dir_cumul_prob[s, n_cat - 1];
                }
        }
        ////
        //// Likelihood using binomial factorization:
        ////
        for (s in 1:n_studies) {
                  for (c in 1:2) {
                      for (cut_i in 1:to_int(n_cutpoints[c, s])) {
                             int k = to_int(cutpoint_index[c][s, cut_i]);
                             logit_cumul_prob[c][s, cut_i] = (C[s, k] - locations[c, s])/scales[c, s];
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
      
        ////
        //// Priors:
        ////
        {
            //// locations:
            beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
            beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
            //// scales:
            raw_scale_mu ~ normal(prior_raw_scale_mu_mean, prior_raw_scale_mu_SD);
            raw_scale_SD ~ normal(prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
        }
        ////
        //// Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        ////
        {
            //// target += normal_lpdf(kappa | prior_kappa_mean, prior_kappa_SD);
            target += dirichlet_lpdf( dirichlet_cat_means_phi | prior_dirichlet_cat_means_alpha ); // "flat" prior on the simplex dirichlet_cat_means_phi. 
            target += normal_lpdf(dirichlet_cat_SDs_sigma | prior_dirichlet_cat_SDs_mean, prior_dirichlet_cat_SDs_SD );
        }
        ////
        {
            for (s in 1:n_studies) {
                vector[n_thr] rho =  normal_pdf(to_vector(Ind_Dir_cumul_prob[s, 1:n_thr]), 0.0, 1.0);   //  p_cumul[k - 1] * (1.0 - p_cumul[k - 1]); // original
                target += induced_dirichlet_v2_lpdf(to_vector(Ind_Dir_ord_prob[s, 1:n_cat]) | rho, alpha[1:n_cat]);
            }
            if (prior_only == 0) {
                  for (k in 1:n_cat) {
                      target += log_kappa;
                      target += log_dirichlet_cat_means_phi[k];
                  }
            }
            target += Jacobian_for_alpha_k_to_category_SDs;
        }
        ////
        //// Likelihood / Model:
        ////
        {
            to_vector(beta_z) ~ std_normal();   // (part of between-study model, NOT prior)
            to_vector(raw_scale_z) ~ std_normal();  // (part of between-study model, NOT prior)
        }
        ////
        //// Increment the log-likelihood:
        ////
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
          vector[n_thr] C_mu_pred;
          matrix[n_studies, n_thr] se;
          matrix[n_studies, n_thr] sp;
          matrix[n_studies, n_thr] fp;
          matrix[n_studies, n_thr] x_hat_nd = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] x_hat_d  = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_nd   = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_d    = rep_matrix(-1, n_studies, n_thr);
          array[2] matrix[n_studies, n_thr] x_hat;
          array[2] matrix[n_studies, n_thr] dev;
          vector[n_thr] C_MU;       
          vector[n_thr] C_MU_empirical; 
          vector[n_thr] prob_cumul_mu;  
          vector[n_cat] prob_ord_mu;  
          real<lower=kappa_lb> precision = kappa; 
          real dispersion = 1.0 / precision;
          vector[n_thr] expected_p_cumul;
  
          //// Initialise containers:
          for (c in 1:2) {
              x_hat[c]  = rep_matrix(-1, n_studies, n_thr);
              dev[c]    = rep_matrix(-1, n_studies, n_thr); 
          }
          
          {
                // Calculate cumulative probabilities:
                expected_p_cumul[1] = dirichlet_cat_means_phi[1];
                for (k in 2:n_thr) {
                    expected_p_cumul[k] = expected_p_cumul[k - 1] + dirichlet_cat_means_phi[k];
                }

                // Transform to cutpoints:
                for (k in 1:n_thr) {
                       real prob = expected_p_cumul[k];
                      if (prob < 0.0001) {
                            C_MU[k] = -10.0;
                      } else if (prob > 0.99999) {
                            C_MU[k] = +10.0;
                      } else {
                            C_MU[k] = 0.0 + inv_Phi(prob);
                      }
          
                }
          }
          ////
          //// Generate a prediction for each cutpoint:
          ////
          {

                  vector[n_cat] prob_ord_mu_pred;
                  vector[n_thr] prob_cumul_mu_pred;

                  // // //// Induced-Dirichlet ordinal probs are structured like this:
                  //  Ind_Dir_cumul_prob = Phi(C - Ind_Dir_anchor);
                  //  for (s in 1:n_studies) {
                  //          //// Induced-Dirichlet ordinal probs:
                  //          Ind_Dir_ord_prob[s, 1] = Ind_Dir_cumul_prob[s, 1] - 0.0;
                  //          for (k in 2:n_thr) {
                  //             Ind_Dir_ord_prob[s, k] = Ind_Dir_cumul_prob[s, k] - Ind_Dir_cumul_prob[s, k - 1]; // since cutpoints are increasing with k
                  //          }
                  //          Ind_Dir_ord_prob[s, n_cat] = 1.0 - Ind_Dir_cumul_prob[s, n_cat - 1];
                  //  }
                  // Simulate from Dirichlet by using the summary "alpha" parameters:
                  prob_ord_mu_pred[1:n_cat]   =  dirichlet_rng(alpha);
                  ////
                  //// Calculate cumulative probabilities:
                  prob_cumul_mu_pred[1] = prob_ord_mu_pred[1];
                  for (k in 2:n_thr) {
                      prob_cumul_mu_pred[k] = prob_cumul_mu_pred[k - 1] + prob_ord_mu_pred[k];
                  }
                  ////
                  //// Transform to cutpoints:
                  real anchor_for_summaries = 0.0;
                  for (k in 1:n_thr) {
                        real prob_1;
                        if (prob_cumul_mu_pred[k] < 1e-38) {
                              prob_1 = 1e-38;
                              C_mu_pred[k] = -10.0;
                        } else if (prob_cumul_mu_pred[k] > 0.9999999999999) {
                              prob_1 = 0.9999999999999;
                              C_mu_pred[k] = +10.0;
                        } else {
                              prob_1 = inv_Phi(prob_cumul_mu_pred[k]);
                              C_mu_pred[k] = anchor_for_summaries + prob_1;
                        }

                  }


          }
          ////
          //// Empirical-mean cutpoints:
          ////
          {
              for (k in 1:n_thr) {
                    C_MU_empirical[k] = median(C[, k]);
              }
          }
          ////
          //// Calculate study-specific accuracy:
          ////
          for (s in 1:n_studies) {
                for (k in 1:n_thr) { 
                    fp[s, k] =   1.0 - cumul_prob[1][s, k];
                    sp[s, k] =   1.0 - fp[s, k];
                    se[s, k] =   1.0 - cumul_prob[2][s, k];
                }
          }
          ////
          //// Calculate summary accuracy (using mean parameters):
          ////
          for (k in 1:n_thr) {
                Fp[k] =   1.0 - Phi((C_MU[k] - (-0.5)*beta_mu)/log1p_exp((-0.5)*raw_scale_mu));
                Sp[k] =   1.0 - Fp[k];
                Se[k] =   1.0 - Phi((C_MU[k] - (+0.5)*beta_mu)/log1p_exp((+0.5)*raw_scale_mu));
          }
          ////
          //// Calculate predictive accuracy:
          ////
          {
                real beta_pred      = normal_rng(beta_mu,       beta_SD);
                real raw_scale_pred = normal_rng(raw_scale_mu, raw_scale_SD);
                for (k in 1:n_thr) {
                      Fp_pred[k] =   1.0 - Phi((C_mu_pred[k] - (-0.5)*beta_pred)/log1p_exp((-0.5)*raw_scale_pred));
                      Sp_pred[k] =   1.0 - Fp_pred[k];
                      Se_pred[k] =   1.0 - Phi((C_mu_pred[k] - (+0.5)*beta_pred)/log1p_exp((+0.5)*raw_scale_pred));
                }
          }
          ////
          //// Model-predicted ("re-constructed") data:
          ////
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














