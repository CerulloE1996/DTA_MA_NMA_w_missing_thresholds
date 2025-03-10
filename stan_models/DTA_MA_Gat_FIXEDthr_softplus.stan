

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
        array[n_studies] int n_obs_cutpoints;  //// OBSERVED cutpoints for test in study s
        ////
        array[2] matrix[n_studies, n_thr] x_with_missings;
        array[2] matrix[n_studies, n_thr] x;
        array[2] matrix[n_studies, n_thr] n;
        array[2] matrix[n_studies, n_thr] cutpoint_index;
        //// Priors for beta:
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
        vector<lower=0.0>[n_thr + 1] prior_alpha; // Induced-Dirichlet prior vector 
        //// Other:
        int prior_only;
  
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
        ordered[n_thr] C;  // Global cutpoints
        
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
        vector[n_thr] Ind_Dir_anchor;
        vector[n_thr] Ind_Dir_cumul;
        vector[n_thr] Ind_Dir_cumul_prob;
        vector[n_cat] Ind_Dir_ord_prob;
        ////
        //// Initialise matrices:
        ////
        for (c in 1:2) {
                 logit_cumul_prob[c]   = rep_matrix(positive_infinity(), n_studies, n_thr);
                 cumul_prob[c]         = rep_matrix(1.0, n_studies, n_thr);
                 cond_prob[c]          = rep_matrix(1.0, n_studies, n_thr);
                 log_lik[c]            = rep_matrix(0.0, n_studies, n_thr);
        }
        ////
        //// Induced-Dirichlet ** prior model ** stuff:
        ////
        {
                 //// Ind-Dir cumulative probs:
                 Ind_Dir_anchor     = rep_vector(0.0, n_thr);
                 Ind_Dir_cumul[1:n_thr] =  (C - Ind_Dir_anchor[1:n_thr]);
                 Ind_Dir_cumul_prob[1:n_thr] = Phi(Ind_Dir_cumul[1:n_thr]);
                 //// Ind-Dir ordinal probs:
                 Ind_Dir_ord_prob[1] = Ind_Dir_cumul_prob[1] - 0.0;
                 for (k in 2:n_thr) {
                     Ind_Dir_ord_prob[k] = Ind_Dir_cumul_prob[k] - Ind_Dir_cumul_prob[k - 1]; // since probs are INCREASING with k
                 }
                 Ind_Dir_ord_prob[n_cat] =  1.0 - Ind_Dir_cumul_prob[n_cat - 1];
        }
        ////
        //// Between-study model for location and scale:
        ////
        //// real Jacobian_raw_scale_to_scale = 0.0;
        {
            ////
            for (s in 1:n_studies) {
                real raw_beta_baseline = beta_mu     + beta_SD      * beta_z[s];
                real raw_scale_baseline = raw_scale_mu + raw_scale_SD * raw_scale_z[s];
                ////
                locations[1, s] = (-1)*(-0.5)*raw_beta_baseline;
                locations[2, s] = (-1)*(+0.5)*raw_beta_baseline;
                ////
                scales[1, s] = log1p_exp(-0.5*raw_scale_baseline);
                scales[2, s] = log1p_exp(+0.5*raw_scale_baseline);
                //// Jacobian for raw_scale -> scale:
                //// Jacobian_raw_scale_to_scale += +log(2) - 0.5*raw_scale_baseline; //// deriv of exp(-0.5*raw_scale_baseline) = -0.5*exp(-0.5*raw_scale_baseline) -> log_abs_deriv = +log(2) -0.5*raw_scale_baseline; //// log_inv_logit(raw_scale_nd);
                //// Jacobian_raw_scale_to_scale += +log(2) + 0.5*raw_scale_baseline; //// deriv of exp(+0.5*raw_scale_baseline) = +0.5*exp(+0.5*raw_scale_baseline) -> log_abs_deriv = -log(2) +0.5*raw_scale_baseline; //// log_inv_logit(raw_scale_nd);
            }
        }
        ////
        //// Likelihood using binomial factorization:
        ////
        for (s in 1:n_studies) {
                  for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                             int k = to_int(cutpoint_index[c][s, cut_i]);
                             logit_cumul_prob[c][s, cut_i] = (C[k] - locations[c, s])/scales[c, s];
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
                      for (cut_i in 1:n_obs_cutpoints[s]) {
                              // Current and next cumulative counts
                              int x_current = to_int(x[c][s, cut_i]);
                              int x_next    = to_int(n[c][s, cut_i]);
                              
                              // Skip if the current count is zero (no observations to classify)
                              if (x_current != 0)  {
                              
                                    // Conditional probability of being at or below the current cutpoint - given being at or below the next cutpoint
                                    if (cut_i == n_obs_cutpoints[s]) { 
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
      
        //// Priors:
        {
            //// locations:
            beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
            beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
            //// scales:
            raw_scale_mu ~ normal(prior_raw_scale_mu_mean, prior_raw_scale_mu_SD);
            raw_scale_SD ~ normal(prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
        }
        ////
        //// Induced-dirichlet ** Prior ** model:
        ////
        {
               vector[n_thr] rho =  normal_pdf(Ind_Dir_cumul[1:n_thr], 0.0, 1.0);   //  p_cumul[k - 1] * (1.0 - p_cumul[k - 1]); // original
               target += induced_dirichlet_v2_lpdf( Ind_Dir_ord_prob[1:n_cat] | rho, prior_alpha);
        }
        ////
        //// Likelihood / Model:
        ////
        //// target += Jacobian_raw_scale_to_scale;
        {
            to_vector(beta_z) ~ std_normal();       // (part of between-study model, NOT prior)
            to_vector(raw_scale_z) ~ std_normal();  // (part of between-study model, NOT prior)
        }
        //// Increment the log-likelihood:
        if (prior_only == 0) {
            for (c in 1:2) {
              target +=  sum(log_lik[c]);
            }
        }
  
}




// 
//       real raw_scale_baseline = raw_scale_mu + raw_scale_SD * raw_scale_z[s];
//                 ////
//                 locations[1, s] = (-1)*(-0.5)*raw_beta_baseline;
//                 locations[2, s] = (-1)*(+0.5)*raw_beta_baseline;
//                 ////
//                 scales[1, s] = exp(-0.5*raw_scale_baseline);
//                 scales[2, s] = exp(+0.5*raw_scale_baseline);
//                 //// Jacobian for raw_scale -> scale:
//                 Jacobian_raw_scale_to_scale += +log(2) - 0.5*raw_scale_baseline; //// deriv of exp(-0.5*raw_scale_baseline) = -0.5*exp(-0.5*raw_scale_baseline) -> log_abs_deriv = +log(2) -0.5*raw_scale_baseline; //// log_inv_logit(raw_scale_nd);
//                 Jacobian_raw_scale_to_scale += +log(2) + 0.5*raw_scale_baseline; //// deriv of exp(+0.5*raw_scale_baseline) = +0.5*exp(+0.5*raw_scale_baseline) -> log_abs_deriv = -log(2) +0.5*raw_scale_baseline; //// log_inv_logit(raw_scale_nd);
//             }
//         }
//  
//                              logit_cumul_prob[c][s, cut_i] = (C[k] - locations[c, s])/scales[c, s];
                             
                             
 

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
          int n_sims = 1000;
  
          //// Initialise containers:
          for (c in 1:2) {
              x_hat[c]  = rep_matrix(-1, n_studies, n_thr);
              dev[c]    = rep_matrix(-1, n_studies, n_thr); 
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
          {
              real location_nd = (-1)*(-0.5)*beta_mu;
              real location_d  = (-1)*(+0.5)*beta_mu;
              real scale_nd = log1p_exp(-0.5*raw_scale_mu);
              real scale_d  = log1p_exp(+0.5*raw_scale_mu); // This corresponds to the MEDIAN of log-normal (could also use mean which is e.g.,: "scale_nd = exp(-0.5*(raw_scale_mu + 0.5*raw_scale_SD"))"
              for (k in 1:n_thr) {
                    Fp[k] =   1.0 - Phi((C[k] - location_nd)/scale_nd);
                    Sp[k] =   1.0 - Fp[k];
                    Se[k] =   1.0 - Phi((C[k] - location_d)/scale_d);
              }
          }
          ////
          //// Calculate predictive accuracy:
          ////
          {
              real beta_pred      = normal_rng(beta_mu,      beta_SD);
              real raw_scale_pred = normal_rng(raw_scale_mu, raw_scale_SD);
              ////
              real location_nd = (-1)*(-0.5)*beta_pred;
              real location_d  = (-1)*(+0.5)*beta_pred;
              real scale_nd = log1p_exp(-0.5*raw_scale_pred);
              real scale_d  = log1p_exp(+0.5*raw_scale_pred); // This corresponds to the MEDIAN of log-normal (could also use mean which is e.g.,: "scale_nd = exp(-0.5*(raw_scale_pred + 0.5*raw_scale_SD"))"
              {
                 
                    for (k in 1:n_thr) {
                          Fp_pred[k] =   1.0 - Phi((C[k] - location_nd)/scale_nd);
                          Sp_pred[k] =   1.0 - Fp_pred[k];
                          Se_pred[k] =   1.0 - Phi((C[k] - location_d)/scale_d);
                    }
          }
          }
          ////
          //// Model-predicted ("re-constructed") data:
          ////
          for (s in 1:n_studies) {
              for (c in 1:2) {
                 for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                  
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














