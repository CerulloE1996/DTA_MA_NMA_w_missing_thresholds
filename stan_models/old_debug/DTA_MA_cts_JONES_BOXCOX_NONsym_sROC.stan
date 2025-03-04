
//// NOTE: This verion is more flexible than the other "RANDOM_thr" (fixed thresholds) version of the
//// model because we are not assuming that the mean/location parameters in the non-diseasded group 
//// are all equal to zero and also allows for between-study heterogenity in the study-specific 
//// mean parameters in BOTH groups. 
//// We can typically do this without identifiability issues, especially in our case since the 
/// within-study models are:
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
  
      // //// Basic bounding function
      // real fn_bound_10( real x ) { 
      //                         
      //        vector[2] container_bound_below;
      //        container_bound_below[1] = -10.0;
      //        container_bound_below[2] = x;
      //        
      //        vector[2] container_bound_above;
      //        container_bound_above[1] = max(container_bound_below);
      //        container_bound_above[2] = +10.0;
      //        
      //        return min(container_bound_above);
      //  
      // }
        
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
      vector[n_thr] cts_thr_values_nd; // = logit(cts_thr_values_nd);   // Global cutpoints 
      vector[n_thr] cts_thr_values_d; //  = logit(cts_thr_values_d);    // Global cutpoints
      //// Priors:
      vector[2] prior_beta_mu_mean;
      vector[2] prior_beta_mu_SD;
      vector[2] prior_beta_SD_mean;
      vector[2] prior_beta_SD_SD;
      vector[2] prior_log_scale_mu_mean;
      vector[2] prior_log_scale_mu_SD;
      vector[2] prior_log_scale_SD_mean;
      vector[2] prior_log_scale_SD_SD;
      //// Other:
      int use_box_cox;
  
}


transformed data { 
    
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
   
}


parameters {
  
      vector[2] beta_mu;    
      vector<lower=0.0>[2] beta_SD;   
      matrix[2, n_studies] beta_z;    // Study-specific random effects for beta (off-centered parameterisation)
      ////
      vector[2] log_scale_mu;    
      vector<lower=0.0>[2] log_scale_SD;   
      matrix[2, n_studies] log_scale_z;    // Study-specific random effects for beta (off-centered parameterisation)
      ////
      vector<lower=-5.0, upper = 5.0>[2] lambda;   // box-cox params
  
}


transformed parameters { 
  
        //// Construct the study-specific random effects (off-centered param.):
        matrix[2, n_studies] beta;
        matrix[2, n_studies] log_scale;
        matrix[2, n_studies] scale;
        matrix[2, n_thr] cutpoints;  // Global cutpoints
        ////
        array[2] matrix[n_studies, n_thr] logit_cumul_prob; // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cumul_prob; // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cond_prob;  // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] log_lik;  // log_lik storage
      
        //// Initialise matrices:
        for (c in 1:2) {
                 logit_cumul_prob[c] = rep_matrix(1.0, n_studies, n_thr);
                 cumul_prob[c]       = rep_matrix(1.0, n_studies, n_thr);
                 cond_prob[c]        = rep_matrix(1.0, n_studies, n_thr);
                 log_lik[c]          = rep_matrix(0.0, n_studies, n_thr);
        }
      
        if (use_box_cox == 0) {
              cutpoints[1, ]  = to_row_vector(log(cts_thr_values_nd));
              cutpoints[2, ]  = to_row_vector(log(cts_thr_values_d));
        } else { 
              cutpoints[1, ]  = to_row_vector(box_cox(cts_thr_values_nd, lambda[1]));
              cutpoints[2, ]  = to_row_vector(box_cox(cts_thr_values_d,  lambda[2]));
        }
        
        for (c in 1:2) {
            beta[c, ]  = to_row_vector(beta_mu[c] + beta_z[c, ] * beta_SD[c]);
            log_scale[c, ] = to_row_vector(log_scale_mu[c] + log_scale_z[c, ] * log_scale_SD[c]);
            scale[c, ] = exp(log_scale[c, ]);
        }
       
        //// Likelihood using binomial factorization:
        for (s in 1:n_studies) {
                    for (c in 1:2) {
                            for (cut_i in 1:to_int(n_cutpoints[c, s])) {
                                    int k = to_int(cutpoint_index[c][s, cut_i]);
                                    logit_cumul_prob[c][s, cut_i] = (beta[c, s] - cutpoints[c, k])/scale[c, s];
                            }
                    }
        }
        //// Calculate CUMULATIVE probabilities (vectorised):
        for (c in 1:2) {
            cumul_prob[c] = Phi_approx(logit_cumul_prob[c]);
        }
        for (s in 1:n_studies) {
                    for (c in 1:2) {
                          real previous_cumul_prob_temp = 1.0;
                          for (cut_i in 1:to_int(n_cutpoints[c, s])) {
                            
                                  real cumul_prob_temp = cumul_prob[c][s, cut_i];
                                  cond_prob[c][s, cut_i] = cumul_prob_temp / previous_cumul_prob_temp;
                                  if (cut_i == 1) { //// First non-missing cutpoint
                                       int success_prob = to_int(n[c][s, 1]);
                                       int x_i = to_int(x[c][s, 1]);
                                       log_lik[c][s, cut_i] = binomial_lpmf( x_i | success_prob, cond_prob[c][s, 1] );
                                  } else if (cut_i > 1) { 
                                       int success_prob = to_int(n[c][s, cut_i]);
                                       int x_i = to_int(x[c][s, cut_i]);
                                       log_lik[c][s, cut_i] = binomial_lpmf( x_i | success_prob,  cond_prob[c][s, cut_i] );
                                  }
                                  previous_cumul_prob_temp = cumul_prob_temp;
                        
                          }
                    }
         }
      
}


model {
      
        //// Priors:
        beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
        beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
        log_scale_mu ~ normal(prior_log_scale_mu_mean, prior_log_scale_mu_SD);
        log_scale_SD ~ normal(prior_log_scale_SD_mean, prior_log_scale_SD_SD);
        
        //// Likelihood / Model:
        for (c in 1:2) {
             to_vector(beta_z[c, ]) ~ std_normal(); // (part of between-study model, NOT prior)
             to_vector(log_scale_z[c, ]) ~ std_normal(); // (part of between-study model, NOT prior)
        }
        
        //// For box-cox parameters:
        target +=  normal_lpdf(lambda | 0.00, 0.50);
        
        //// Jacobian adjustments needed:
        { 
            target += sum(log_scale);                // double-checked the log-derivative of this by hand (correct)
            if (abs(sum(log_scale_SD)) != 0.0) {  // just in case inits are set to exactly zero (sampler will otherwise fail)
                 target += log(abs(sum(log_scale_SD)));      // double-checked the log-derivative of this by hand (correct)
            } 
            if (abs(sum(log_scale_z)) != 0.0) {  // just in case inits are set to exactly zero (sampler will otherwise fail)
                 target += log(abs(sum(log_scale_z)));  // double-checked the log-derivative of this by hand (correct)
            }
        }
                    
        //// Likelihood using binomial factorization:
        for (c in 1:2) {
          target +=  sum(log_lik[c]);
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
          
          //// Initialise containers:
          for (c in 1:2) {
              x_hat[c]  = rep_matrix(-1, n_studies, n_thr);
              dev[c]    = rep_matrix(-1, n_studies, n_thr);
          }
              
          //// Calculate study-specific accuracy:
          for (s in 1:n_studies) {
             for (k in 1:n_thr) {
                    fp[s, k] =   Phi_approx((beta[1, s] - cutpoints[1][k])/scale[1, s]);
                    sp[s, k] =   1.0 - fp[s, k];
                    se[s, k] =   Phi_approx((beta[2, s] - cutpoints[2][k])/scale[2, s]);
                }
          }
          
          //// Calculate summary accuracy (using mean parameters):
          for (k in 1:n_thr) {
                Fp[k] =   Phi_approx((beta_mu[1] - cutpoints[1][k])/exp(log_scale_mu[1]));
                Sp[k] =   1.0 - Fp[k];
                Se[k] =   Phi_approx((beta_mu[2] - cutpoints[2][k])/exp(log_scale_mu[2]));
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























