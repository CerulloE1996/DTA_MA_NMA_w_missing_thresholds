

functions {
        
      //// Box-Cox transform function:
      real box_cox( real x, 
                    real lambda) {
            
            if (lambda == 0.0) {
                return log(x);
            } else {
                return (pow(x, lambda) - 1.0) / lambda;
            }
          
      }
      
      //// Vectorized Box-Cox transform function:
      vector box_cox( vector x, 
                      real lambda) {
                               
            int N = num_elements(x);
            vector[N] result;
            
            for (n in 1:N) {
              
                if (lambda == 0.0) {
                     result[n] = log(x[n]); 
                } else {
                     result[n] = (pow(x[n], lambda) - 1.0) / lambda;
                }
                
            }
            
            return result;
          
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
      vector[n_thr] cts_thr_values;
      //// Priors:
      vector[2] prior_beta_mu_mean;
      vector[2] prior_beta_mu_SD;
      vector[2] prior_beta_SD_mean;
      vector[2] prior_beta_SD_SD;
      vector[2] prior_raw_scale_mu_mean;
      vector[2] prior_raw_scale_mu_SD;
      vector[2] prior_raw_scale_SD_mean;
      vector[2] prior_raw_scale_SD_SD;
      //// Other:
      int prior_only;
      int use_box_cox;
      vector[2] prior_boxcox_lambda_mean;
      vector<lower=0.0>[2] prior_boxcox_lambda_SD;
  
}


transformed data { 
    
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
   
}


parameters {
  
      vector[2] beta_mu;    
      vector<lower=0.0>[2] beta_SD;   
      matrix[2, n_studies] beta_z;    // Study-specific random effects for beta (off-centered parameterisation)
      cholesky_factor_corr[2] beta_L_Omega;
      ////
      vector[2] raw_scale_mu;    
      vector<lower=0.0>[2] raw_scale_SD;   
      matrix[2, n_studies] raw_scale_z;    // Study-specific random effects for beta (off-centered parameterisation)
      cholesky_factor_corr[2] raw_scale_L_Omega;
      ////
      real<lower=-5.0, upper = 5.0> lambda;   // box-cox params
  
}


transformed parameters { 
  
        //// Construct the study-specific random effects (off-centered param.):
        matrix[2, n_studies] beta;
        matrix[2, n_studies] raw_scale;
        matrix[2, n_studies] scale;
        ////
        array[2] matrix[n_studies, n_thr] logit_cumul_prob; // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cumul_prob; // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cond_prob;  // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] log_lik;  // log_lik storage
        ////
        vector[n_thr] C;
        //// 
        real Jacobian_raw_scale_to_scale = 0.0;
        ////
        //// Initialise matrices:
        ////
        for (c in 1:2) {
                 logit_cumul_prob[c] = rep_matrix(1.0, n_studies, n_thr);
                 cumul_prob[c]       = rep_matrix(1.0, n_studies, n_thr);
                 cond_prob[c]        = rep_matrix(1.0, n_studies, n_thr);
                 log_lik[c]          = rep_matrix(0.0, n_studies, n_thr);
        }
        ////
        if (use_box_cox == 0) {
              C  = log(cts_thr_values);
        } else { 
              C  = box_cox(cts_thr_values, lambda);
        }
        ////
        //// Between-study model for the location parameters ("beta") and the raw_scale parameters - models between-study correlation:
        ////
        for (s in 1:n_studies) {
             beta[, s]      =  beta_mu      + diag_pre_multiply(beta_SD, beta_L_Omega) * beta_z[, s];
             raw_scale[, s] =  raw_scale_mu + diag_pre_multiply(raw_scale_SD, raw_scale_L_Omega) * raw_scale_z[, s];
        }
        
        for (c in 1:2) {
            scale[c, ] = log1p_exp(raw_scale[c, ]);
            Jacobian_raw_scale_to_scale += sum(log_inv_logit(raw_scale[c, ]));
        }
        ////
        //// Likelihood using binomial factorization:
        ////
        for (s in 1:n_studies) {
                  for (c in 1:2) {
                      for (cut_i in 1:to_int(n_cutpoints[c, s])) {
                             int k = to_int(cutpoint_index[c][s, cut_i]);
                             logit_cumul_prob[c][s, cut_i] = (C[k] - beta[c, s])/scale[c, s];
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
                                             cond_prob[c][s, cut_i] = cumul_prob[c][s, cut_i] / 1.0;
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
        beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
        beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
        beta_L_Omega ~ lkj_corr_cholesky(2.0);
        ////
        raw_scale_mu ~ normal(prior_raw_scale_mu_mean, prior_raw_scale_mu_SD);
        raw_scale_SD ~ normal(prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
        raw_scale_L_Omega ~ lkj_corr_cholesky(2.0);
        
        //// Likelihood / Model:
        for (c in 1:2) {
             to_vector(beta_z[c, ]) ~ std_normal(); // (part of between-study model, NOT prior)
             to_vector(raw_scale_z[c, ]) ~ std_normal(); // (part of between-study model, NOT prior)
        }
        
        //// For box-cox parameters:
        target +=  normal_lpdf(lambda | prior_boxcox_lambda_mean, prior_boxcox_lambda_SD);
        
        //// Jacobian adjustments needed:
        { 
            target +=  Jacobian_raw_scale_to_scale;
            //// Also need the Jacobian to go from (raw_scale_z, raw_scale_mu, raw_scale_SD) -> raw_scale! (using chain rule):
            if (abs(sum(raw_scale_SD)) != 0.0) target += log(abs(sum(raw_scale_SD)));      // double-checked the log-derivative of this by hand (correct)
            if (abs(sum(raw_scale_z)) != 0.0)  target += log(abs(sum(raw_scale_z)));  // double-checked the log-derivative of this by hand (correct)
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
          array[n_studies] vector[n_thr] se;
          array[n_studies] vector[n_thr] sp;
          array[n_studies] vector[n_thr] fp;
          matrix[n_studies, n_thr] x_hat_nd = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] x_hat_d  = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_nd   = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_d    = rep_matrix(-1, n_studies, n_thr);
          array[2] matrix[n_studies, n_thr] x_hat;
          array[2] matrix[n_studies, n_thr] dev;
          corr_matrix[2] beta_Omega;
          corr_matrix[2] raw_scale_Omega;
          vector[2] scale_mu = log1p_exp(raw_scale_mu);
          
          //// Compute between-study correlation matrices:
          beta_Omega       =  multiply_lower_tri_self_transpose(beta_L_Omega);
          raw_scale_Omega  =  multiply_lower_tri_self_transpose(raw_scale_L_Omega);
          
          //// Initialise containers:
          for (c in 1:2) {
              x_hat[c]  = rep_matrix(-1, n_studies, n_thr);
              dev[c]    = rep_matrix(-1, n_studies, n_thr);
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
                Fp[k] =   1.0 - Phi((C[k] - beta_mu[1])/scale_mu[1]);
                Sp[k] =   1.0 - Fp[k];
                Se[k] =   1.0 - Phi((C[k] - beta_mu[2])/scale_mu[2]);
          }
          ////
          //// Calculate predictive accuracy:
          ////
          {
                vector[2] beta_pred      =  to_vector(multi_normal_cholesky_rng(beta_mu, diag_pre_multiply(beta_SD, beta_L_Omega)));
                vector[2] raw_scale_pred =  to_vector(multi_normal_cholesky_rng(raw_scale_mu, diag_pre_multiply(raw_scale_SD, raw_scale_L_Omega)));
                vector[2] scale_pred = log1p_exp(raw_scale_pred);
                ////
                for (k in 1:n_thr) {
                      Fp_pred[k] =   1.0 - Phi((C[k] - beta_pred[1])/scale_pred[1]);
                      Sp_pred[k] =   1.0 - Fp_pred[k];
                      Se_pred[k] =   1.0 - Phi((C[k] - beta_pred[2])/scale_pred[2]);
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























