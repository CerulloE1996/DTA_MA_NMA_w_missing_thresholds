

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
          int<lower=1> n_tests;
          ////
          array[n_tests] int<lower=1> n_thr;
          int<lower=1> n_thr_max; // corresponding to the test with the most thresholds/cutpoints. 
          array[n_studies, n_tests] int n_obs_cutpoints; //// OBSERVED cutpoints for test t in study s
          ////
          array[n_studies] int<lower=0, upper=n_tests> n_tests_per_study;  // Tests per study
          array[n_studies, n_tests] int<lower=0, upper=1> indicator_test_in_study;   // Binary indicator if test t is in study s
          ////
          //// Data:
          ////
          array[n_tests, 2] matrix[n_studies, n_thr_max] x_with_missings;
          array[n_tests, 2] matrix[n_studies, n_thr_max] n;
          array[n_tests, 2] matrix[n_studies, n_thr_max] x;
          array[n_tests, 2] matrix[n_studies, n_thr_max] cutpoint_index;
          ////
          // matrix[n_studies, n_tests] obs_indicator; // BINARY indicator for whether test t is observed in study s
          ////
          //// Priors for beta:
          ////
          vector[n_tests] prior_beta_mu_mean;
          vector[n_tests] prior_beta_mu_SD;
          vector[n_tests] prior_beta_tau_SD;
          real prior_beta_sigma_SD; //// "sigma's (Nyaga et al. notation) are shared between tests.
          ////
          //// Priors for raw_scale:
          ////
          vector[n_tests] prior_raw_scale_mu_mean;
          vector[n_tests] prior_raw_scale_mu_SD;
          vector[n_tests] prior_raw_scale_tau_SD;
          real prior_raw_scale_sigma_SD; //// "sigma's (Nyaga et al. notation) are shared between tests.
          ////
          //// Induced-Dirichlet priors:
          ////
          array[n_tests] vector<lower=0.0>[n_thr_max + 1] prior_alpha; // Induced-Dirichlet prior vector ("staggered" matrix using "n_thr[t]" to index correctly).
          //// Other:
          int prior_only;
  
}


transformed data { 
      
          array[n_tests] int<lower=1> n_cat;
          int n_cat_max = n_thr_max + 1; //// Number of ordinal categories for index test
          int n_total_cutpoints = 0;     //// total C across ALL tests
          int N_obs = 0;                 //// The TOTAL number of observed data points across ALL studies and ALL tests and ALL thresholds
          
          for (t in 1:n_tests) {
              n_cat[t] = n_thr[t] + 1;
          }
          
          for (c in 1:2) {
              for (s in 1:n_studies) {
                  for (t in 1:n_tests) {
                      for (k in 1:n_thr[t]) {
                            int obs_ind = indicator_test_in_study[s, t];
                            if (obs_ind == 1) { //// Only increment if data point is actually observed
                                N_obs += 1;
                            }
                      }
                  }
              }
          }
          
          for (t in 1:n_tests) {
              n_total_cutpoints += n_thr[t];
          }
          
}


parameters {
  
          vector[n_tests] beta_mu;    
          vector[n_tests] raw_scale_mu;    
          ////
          vector[n_total_cutpoints] C_raw_vec; //// RAW LOG-DIFFERENCES - Global cutpoints for each test (staggered array/matrix using "n_thr[t]" to index correctly)
          ////
          //// "NMA" params:
          ////
          vector[n_studies] beta_eta_z;       //// Standard normal RVs for study-level effects - eta[s, 1:2] ~ multi_normal({0, 0}, Sigma).
          vector[n_studies] raw_scale_eta_z;  //// Standard normal RVs for study-level effects - eta[s, 1:2] ~ multi_normal({0, 0}, Sigma).
          real beta_sigma;                    //// Variance component - Between-study SD (Nyaga's σ)
          real raw_scale_sigma;               //// Variance component - Between-study SD (Nyaga's σ)
          ////
          vector<lower=0>[n_tests] beta_tau;             //// Variance component - Test-specific SD (Nyaga's τ) - delta_{s, c, t} ~ normal(0, tau_{c, t}).
          vector<lower=0>[n_tests] raw_scale_tau;        //// Variance component - Test-specific SD (Nyaga's τ) - delta_{s, c, t} ~ normal(0, tau_{c, t}).
          matrix[n_studies, n_tests] beta_delta_z;       //// Standard normal RVs for test-specific effects
          matrix[n_studies, n_tests] raw_scale_delta_z;  //// Standard normal RVs for test-specific effects
        
}


transformed parameters {
    
          ////
          //// -------- Cutpoint params:
          ////
          array[n_tests] vector[n_thr_max] C;  // Global cutpoints for each test (staggered array/matrix using "n_thr[t]" to index correctly)
          ////
          //// Initialise matrices:
          ////
          for (t in 1:n_tests) {
              C[t]         = rep_vector(0.0, n_thr_max);
          }
          ////
          //// Construct cutpoints for each test:
          ////
          {
            int counter = 0;
            for (t in 1:n_tests) {
              for (k in 1:n_thr[t]) {
                 counter += 1;
                 C[t][k] = C_raw_vec[counter];
              }
            }
          }
          //// 1st cutpoint (for each test) is unconstrained (no Jacobian needed):
          { 
            int counter = 0;
            for (t in 1:n_tests) {
               counter += 1;
               C[t][1] = C_raw_vec[counter];
            }
          }
          //// Rest of cutpoints are made using LOG-differences:
          real Jacobian_C = 0;
          {
            int counter = 0;
            for (t in 1:n_tests) {
              for (k in 2:n_thr[t]) {
                   counter += 1;
                   //// C[t][k] = C[t][k - 1] + exp(C_raw_vec[counter]);
                   //// Jacobian_C += C_raw_vec[counter]; // Jacobian for transformation C_raw -> C. 
                   C[t][k] = C[t][k - 1] + log1p_exp(C_raw_vec[counter]);
                   Jacobian_C += log_inv_logit(C_raw_vec[counter]); //// Jacobian for transformation C_raw -> C. (derivative of softplus(x) = log(1.0 + exp(x)) is inv_logit(x) so log(deriv) = log_inv_logit(x)). 
              }
            }
          }
          ////
          array[n_tests] matrix[2, n_studies] locations;
          array[n_tests] matrix[2, n_studies] scales;
          ////
          array[n_tests, 2] matrix[n_studies, n_thr_max] logit_cumul_prob; // Ordinal probs for the likelihood (staggered array/matrix using "n_thr[t]" to index correctly)
          array[n_tests, 2] matrix[n_studies, n_thr_max] cumul_prob;       // Ordinal probs for the likelihood (staggered array/matrix using "n_thr[t]" to index correctly)
          array[n_tests, 2] matrix[n_studies, n_thr_max] cond_prob;        // Ordinal probs for the likelihood (staggered array/matrix using "n_thr[t]" to index correctly)
          array[n_tests, 2] matrix[n_studies, n_thr_max] log_lik;          // log_lik storage (staggered array/matrix using "n_thr[t]" to index correctly)
          ////
          array[n_tests] vector[n_thr_max] Ind_Dir_anchor;
          array[n_tests] vector[n_thr_max] Ind_Dir_cumul;
          array[n_tests] vector[n_thr_max] Ind_Dir_cumul_prob;
          array[n_tests] vector[n_cat_max] Ind_Dir_ord_prob;
          ////
          //// Initialise matrices:
          ////
          for (t in 1:n_tests) { 
              for (c in 1:2) {
                       logit_cumul_prob[t, c]   = rep_matrix(positive_infinity(), n_studies, n_thr_max);
                       cumul_prob[t, c]         = rep_matrix(1.0, n_studies, n_thr_max);
                       cond_prob[t, c]          = rep_matrix(1.0, n_studies, n_thr_max);
                       log_lik[t, c]            = rep_matrix(0.0, n_studies, n_thr_max);
              }
          }
          ////
          //// Induced-Dirichlet ** prior model ** stuff:
          ////
          for (t in 1:n_tests) { 
                   //// Ind-Dir cumulative probs:
                   Ind_Dir_anchor[t]    = rep_vector(0.0, n_thr_max);
                   Ind_Dir_cumul[t][1:n_thr_max] =  (C[t] - Ind_Dir_anchor[t][1:n_thr_max]);
                   Ind_Dir_cumul_prob[t][1:n_thr_max] = Phi(Ind_Dir_cumul[t][1:n_thr_max]);
                   //// Ind-Dir ordinal probs:
                   Ind_Dir_ord_prob[t][1] = Ind_Dir_cumul_prob[t][1] - 0.0;
                   for (k in 2:n_thr[t]) {
                       Ind_Dir_ord_prob[t][k] = Ind_Dir_cumul_prob[t][k] - Ind_Dir_cumul_prob[t][k - 1]; // since probs are INCREASING with k
                   }
                   Ind_Dir_ord_prob[t][n_cat_max] =  1.0 - Ind_Dir_cumul_prob[t][n_cat[t] - 1];
          }
          ////
          //// -------- "NMA" params:
          ////
          //// Declare Study-level random effects (eta in Nyaga notation) -  eta[s, 1:2] ~ multi_normal({0, 0}, Sigma):
          ////
          vector[n_studies] beta_eta;  // first eta corresponds to beta and 2nd eta corresponds to raw_scale.
          vector[n_studies] raw_scale_eta;  // first eta corresponds to beta and 2nd eta corresponds to raw_scale.
          ////
          //// Compute Study-level random effects (eta in Nyaga notation)
          ////
          for (s in 1:n_studies) {
              beta_eta[s]      = 0.0 + beta_sigma      * beta_eta_z[s];      // NOTE: 1st eta's ("beta_eta") correspond to shared (between tests) component of "beta"           -  eta_{s, i} ~ normal(0, sigma_{i}).
              raw_scale_eta[s] = 0.0 + raw_scale_sigma * raw_scale_eta_z[s]; // NOTE: 2nd eta's ("raw_scale_eta") correspond to shared (between tests) component of "raw_scale" -  eta_{s, i} ~ normal(0, sigma_{i}).
          }
          ////
          //// Declare test-specific deviations ("delta" in Nyaga notation) - delta_{s, c, t} ~ normal(0, tau_{c, t}):
          ////
          matrix[n_studies, n_tests] beta_delta;
          matrix[n_studies, n_tests] raw_scale_delta;
          ////
          //// Compute test-specific deviations ("delta" in Nyaga notation) - delta_{s, c, t} ~ normal(0, tau_{c, t}):
          ////
          for (t in 1:n_tests) {
              for (s in 1:n_studies) {
                  if (indicator_test_in_study[s, t] == 1) {
                        beta_delta[s, t]      = 0.0 + beta_delta_z[s, t]      * beta_tau[t];       ////  NOTE: 1st delta's ("beta_delta") correspond to shared (between tests) component of "beta"   -   delta_{s, c, t} ~ normal(0, tau_{c, t}):
                        raw_scale_delta[s, t] = 0.0 + raw_scale_delta_z[s, t] * raw_scale_tau[t];  ////  NOTE: 2nd delta's ("raw_scale_delta") correspond to shared (between tests) component of "raw_scale" - delta_{s, c, t} ~ normal(0, tau_{c, t}):            }
                  }
              }
          }
          ////  
          //// Compute probit CUMULATIVE probabilities:
          ////
          for (t in 1:n_tests) {
                for (s in 1:n_studies) {
                      if (indicator_test_in_study[s, t] == 1) {
                        
                                  for (cut_i in 1:n_obs_cutpoints[s, t]) { //// only loop through the OBSERVED cutpoints for test t in study s
                                        
                                        int k = to_int(cutpoint_index[t, 1][s, cut_i]);
                                        
                                        real link_pi_for_beta      = (beta_mu[t]      + beta_eta[s]      + beta_delta[s, t]);      //// pi_{s, c, t} = Phi(mu_{c, t} + eta_{s, c} + delta_{s, c, t}).
                                        real link_pi_for_raw_scale = (raw_scale_mu[t] + raw_scale_eta[s] + raw_scale_delta[s, t]); //// pi_{s, c, t} = Phi(mu_{c, t} + eta_{s, c} + delta_{s, c, t}).
                                        
                                        //// Use opposite signs for non-diseased vs diseased:
                                        locations[t][1, s] = -0.5*link_pi_for_beta; //// D- group ////  -0.5*(beta_mu + beta_SD * beta_z[s]);
                                        locations[t][2, s] = +0.5*link_pi_for_beta; //// D+ group ////  +0.5*(beta_mu + beta_SD * beta_z[s]);
                                        scales[t][1, s] = log1p_exp(-0.5*link_pi_for_raw_scale); //// D- group //// log1p_exp(-0.5*(raw_scale_mu + raw_scale_SD * raw_scale_z[s]));
                                        scales[t][2, s] = log1p_exp(+0.5*link_pi_for_raw_scale); //// D+ group //// log1p_exp(+0.5*(raw_scale_mu + raw_scale_SD * raw_scale_z[s]));
                                        
                                        //// logit_cumul_prob[c, t, s, cut_i] = (C[t][k] - location)/scale;
                                        for (c in 1:2) {
                                            logit_cumul_prob[t, c][s, cut_i] = (C[t][k] - locations[t][c, s])/scales[t][c, s];
                                        }
                                    
                                  }
                      }
                      
               } 
          }
          ////
          //// Calculate CUMULATIVE probabilities (vectorised):
          ////
          for (t in 1:n_tests) {
              for (c in 1:2) {
                  cumul_prob[t, c] = Phi(logit_cumul_prob[t, c]); //// INCREASING sequence (as C_k > C_{k - 1})
              }
          }
          ////
          //// ------- Likelihood using binomial factorization:
          ////
          for (t in 1:n_tests) {
              for (s in 1:n_studies) {
                    
                     //// Log-Likelihood is evaluated only at the OBSERVED data:
                     if (indicator_test_in_study[s, t] == 1) {
                          for (c in 1:2) {
                              for (cut_i in 1:n_obs_cutpoints[s, t]) {
                              
                                      // //// Log-Likelihood is evaluated only at the OBSERVED data:
                                      // if (observed[t, c][s, cut_i] == 1) {
                                        
                                            //// Current and next cumulative counts:
                                            int x_current = to_int(x[t, c][s, cut_i]);
                                            int x_next    = to_int(n[t, c][s, cut_i]);
                                            
                                            //// Skip if the current count is zero (no observations to classify):
                                            if (x_current != 0)  {
                                            
                                                  //// Conditional probability of being at or below the current cutpoint - given being at or below the next cutpoint:
                                                  if (cut_i == n_obs_cutpoints[s, t]) { 
                                                           cond_prob[t, c][s, cut_i] = cumul_prob[t, c][s, cut_i] / 1.0;
                                                  } else {
                                                        if (x_next > 0) { 
                                                           cond_prob[t, c][s, cut_i] = cumul_prob[t, c][s, cut_i] / cumul_prob[t, c][s, cut_i + 1];
                                                        } else { 
                                                           cond_prob[t, c][s, cut_i] = 1.0;
                                                        }
                                                  }
                                                  
                                                  //// Binomial for observations at or below current cutpoint out of those at or below next cutpoint:
                                                  log_lik[t, c][s, cut_i] = binomial_lpmf(x_current | x_next, cond_prob[t, c][s, cut_i]);
                                                  
                                            }
                                      // }
                              
                              }
                          }
                   }
               
            }
            
          }
      
}


model {
          ////
          //// Priors:
          ////
          for (t in 1:n_tests) {
              beta_mu[t]      ~ normal(prior_beta_mu_mean[t], prior_beta_mu_SD[t]);
              raw_scale_mu[t] ~ normal(prior_raw_scale_mu_mean[t], prior_raw_scale_mu_SD[t]);
          }
          ////
          for (t in 1:n_tests) {
               beta_tau[t]        ~ normal(0.0, prior_beta_tau_SD[t]);       //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
               raw_scale_tau[t]   ~ normal(0.0, prior_raw_scale_tau_SD[t]);  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
          }
          ////
          beta_sigma      ~ normal(0.0, prior_beta_sigma_SD);      //// eta[s, i] ~ normal(0, sigma[i]):
          raw_scale_sigma ~ normal(0.0, prior_raw_scale_sigma_SD); //// eta[s, i] ~ normal(0, sigma[i]):
          ////
          //// Induced-dirichlet ** Prior ** model (NOT actual between-study model because the cutpoints are fixed between studies for this model variation). 
          ////
          for (t in 1:n_tests) {
                 vector[n_thr[t]] rho =  normal_pdf(Ind_Dir_cumul[t][1:n_thr[t]], 0.0, 1.0);   //  p_cumul[k - 1] * (1.0 - p_cumul[k - 1]); // original
                 target += induced_dirichlet_v2_lpdf(Ind_Dir_ord_prob[t][1:n_cat[t]] | rho, prior_alpha[t][1:n_thr[t]]);
                 target += Jacobian_C; // Jacobian for transfomation C_raw -> C.
          }
          ////
          //// Likelihood / Model:
          //// 
          {
              beta_eta_z      ~ std_normal(); // (part of between-test / between-study model, NOT prior) - eta[s, i] ~ normal(0, sigma[i]):
              raw_scale_eta_z ~ std_normal(); // (part of between-test / between-study model, NOT prior) - eta[s, i] ~ normal(0, sigma[i]):
          }
          ////
          to_vector(beta_delta_z)      ~ std_normal(); // (part of between-test / between-study model, NOT prior)  -   delta_{s, i, t} ~ normal(0, tau_{i, t}):
          to_vector(raw_scale_delta_z) ~ std_normal(); // (part of between-test / between-study model, NOT prior)  -   delta_{s, i, t} ~ normal(0, tau_{i, t}):
          ////
          //// Increment the log-likelihood:
          ////
          if (prior_only == 0) {
            for (t in 1:n_tests) {
              for (c in 1:2) {
                target +=  sum(log_lik[t][c]);
              }
            }
          }
    
}




 

generated quantities {

          //// Summary accuracy parameters (at each threshold):
          array[n_tests] vector[n_thr_max] Se;
          array[n_tests] vector[n_thr_max] Sp;
          array[n_tests] vector[n_thr_max] Fp;
          array[n_tests] vector[n_thr_max] Se_pred;
          array[n_tests] vector[n_thr_max] Sp_pred;
          array[n_tests] vector[n_thr_max] Fp_pred;
          array[n_tests] matrix[n_studies, n_thr_max] se;
          array[n_tests] matrix[n_studies, n_thr_max] sp;
          array[n_tests] matrix[n_studies, n_thr_max] fp;
          array[n_tests] matrix[n_studies, n_thr_max] x_hat_nd;// = rep_matrix(-1, n_studies, n_thr_max);
          array[n_tests] matrix[n_studies, n_thr_max] x_hat_d;//  = rep_matrix(-1, n_studies, n_thr_max);
          array[n_tests] matrix[n_studies, n_thr_max] dev_nd;//    = rep_matrix(-1, n_studies, n_thr_max);
          array[n_tests] matrix[n_studies, n_thr_max] dev_d;//     = rep_matrix(-1, n_studies, n_thr_max);
          array[n_tests, 2] matrix[n_studies, n_thr_max] x_hat;
          array[n_tests, 2] matrix[n_studies, n_thr_max] dev;
          int n_sims = 1000;
          ////
          //// Initialise containers:
          ////
          for (t in 1:n_tests) {
            for (c in 1:2) {
                x_hat[t, c]  = rep_matrix(-1, n_studies, n_thr_max);
                dev[t, c]    = rep_matrix(-1, n_studies, n_thr_max); 
            }
          }
          ////
          //// Calculate study-specific accuracy:
          ////
          for (s in 1:n_studies) {
              for (t in 1:n_tests) {
                  for (k in 1:n_thr[t]) { 
                      fp[t][s, k] =   1.0 - cumul_prob[t, 1][s, k];
                      sp[t][s, k] =   1.0 - fp[t][s, k];
                      se[t][s, k] =   1.0 - cumul_prob[t, 2][s, k];
                  }
             }
          }
          ////
          //// Calculate summary accuracy (using mean parameters):
          ////
          for (t in 1:n_tests) {
            for (k in 1:n_thr[t]) {
                  Fp[t][k] =   1.0 - Phi((C[t][k] - (-0.5)*beta_mu[t])/log1p_exp((-0.5)*raw_scale_mu[t]));
                  Sp[t][k] =   1.0 - Fp[t][k];
                  Se[t][k] =   1.0 - Phi((C[t][k] - (+0.5)*beta_mu[t])/log1p_exp((+0.5)*raw_scale_mu[t]));
            }
          }
          ////
          //// Calculate predictive accuracy:
          ////
          real beta_eta_pred      = normal_rng(0.0, beta_sigma);      //// shared between tests
          real raw_scale_eta_pred = normal_rng(0.0, raw_scale_sigma); //// shared between tests
          ////
          for (t in 1:n_tests) {
                ////
                real beta_delta_t_pred      = normal_rng(0.0, beta_tau[t]);
                real raw_scale_delta_t_pred = normal_rng(0.0, raw_scale_tau[t]);
                ////
                real beta_pred      = beta_eta_pred       + beta_delta_t_pred;
                real raw_scale_pred = raw_scale_eta_pred  + raw_scale_delta_t_pred;
                ////
                for (k in 1:n_thr[t]) {
                      Fp_pred[k] =   1.0 - Phi((C[k] - (-0.5)*beta_pred)/log1p_exp((-0.5)*raw_scale_pred));
                      Sp_pred[k] =   1.0 - Fp_pred[k];
                      Se_pred[k] =   1.0 - Phi((C[k] - (+0.5)*beta_pred)/log1p_exp((+0.5)*raw_scale_pred));
                }
          }
          ////
          //// Model-predicted ("re-constructed") data:
          ////
          for (t in 1:n_tests) {
              for (s in 1:n_studies) {
                    if (indicator_test_in_study[s, t] == 1) {
                          
                              for (c in 1:2) {
                                 for (cut_i in 1:to_int(n_obs_cutpoints[s, t])) {
                                  
                                          //// Model-estimated data:
                                          x_hat[t, c][s, cut_i] = cond_prob[t, c][s, cut_i] * n[t, c][s, cut_i];  	 //// Fitted values
                                        
                                          //// Compute residual deviance contribution:
                                          real n_i =  (n[t, c][s, cut_i]);
                                          real x_i =  (x[t, c][s, cut_i]);
                                          real x_hat_i =  (x_hat[t, c][s, cut_i]);
                                          real log_x_minus_log_x_hat = log(x_i) - log(x_hat_i);
                                          real log_diff_n_minus_x = log(n_i - x_i);
                                          real log_diff_n_minus_x_hat = log(abs(n_i - x_hat_i));
                                           
                                          dev[t, c][s, cut_i] = 2.0 * ( x_i * log_x_minus_log_x_hat + (n_i - x_i) * (log_diff_n_minus_x - log_diff_n_minus_x_hat) ); 
                                      
                                 }
                              }
                      
                    }
              }
          }
          ////
          //// Store deviance and x_hat split by disease status:
          ////
          for (t in 1:n_tests) {
               x_hat_nd[t] = x_hat[t][1];
               dev_nd[t]   = dev[t][1];
               x_hat_d[t]  = x_hat[t][2];
               dev_d[t]    = dev[t][2];
          }
       
}














