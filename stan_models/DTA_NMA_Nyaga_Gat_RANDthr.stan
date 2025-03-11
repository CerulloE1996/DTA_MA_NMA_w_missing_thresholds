

////
//// Include files to compile the necessary Stan functions:
////
#include "Stan_fns_basic.stan"
#include "Stan_fns_included_Dirichlet.stan"
#include "Stan_fns_simplex.Stan"



data {
      
          int<lower=1> n_studies;
          int<lower=1> n_index_tests;
          ////
          array[n_index_tests] int<lower=1> n_thr;
          array[n_index_tests] int<lower=2> n_cat;
          int<lower=1> n_thr_max; // corresponding to the test with the most thresholds/cutpoints. 
          int<lower=1> n_cat_max; // corresponding to the test with the most thresholds/cutpoints. 
          array[n_studies, n_index_tests] int n_obs_cutpoints; //// OBSERVED cutpoints for test t in study s
          int n_total_cutpoints;
          int n_total_raw_simplex_elements;
          ////
          // array[n_studies] int<lower=0, upper=n_index_tests> n_index_tests_per_study;  // Tests per study
          array[n_studies, n_index_tests] int<lower=0, upper=1> indicator_index_test_in_study;   // Binary indicator if test t is in study s
          ////
          //// Data:
          ////
          array[n_index_tests, 2, n_studies, n_thr_max] int x_with_missings;
          array[n_index_tests, 2, n_studies, n_thr_max] int x;
          array[n_index_tests, 2, n_studies, n_thr_max] int n;
          array[n_index_tests, 2, n_studies, n_thr_max] int cutpoint_index;
          ////
          // matrix[n_studies, n_index_tests] obs_indicator; // BINARY indicator for whether test t is observed in study s
          ////
          //// Priors for beta:
          ////
          vector[n_index_tests] prior_beta_mu_mean;
          vector<lower=0.0>[n_index_tests] prior_beta_mu_SD;
          vector<lower=0.0>[n_index_tests] prior_beta_tau_SD;
          real<lower=0.0> prior_beta_sigma_SD; //// "sigma's (Nyaga et al. notation) are shared between tests.
          ////
          //// Priors for raw_scale:
          ////
          vector[n_index_tests] prior_raw_scale_mu_mean;
          vector<lower=0.0>[n_index_tests] prior_raw_scale_mu_SD;
          vector<lower=0.0>[n_index_tests] prior_raw_scale_tau_SD;
          real<lower=0.0> prior_raw_scale_sigma_SD; //// "sigma's (Nyaga et al. notation) are shared between tests.
          ////
          //// Induced-Dirichlet priors:
          ////
          array[n_index_tests] vector[n_thr_max + 1] prior_dirichlet_cat_means_alpha;
          array[n_index_tests] vector[n_thr_max + 1] prior_dirichlet_cat_SDs_mean;
          array[n_index_tests] vector<lower=0.0>[n_thr_max + 1] prior_dirichlet_cat_SDs_SD;
          //// Other:
          int prior_only;
          vector[n_index_tests] kappa_lb;
          int use_softplus_for_scales;
  
}



parameters {
  
          vector[n_index_tests] beta_mu;    
          vector[n_index_tests] raw_scale_mu;    
          ////
          vector<lower=-15.0, upper=15.0>[n_total_cutpoints] C_raw_vec;  //// RAW LOG-DIFFERENCES - Global cutpoints for each test (staggered array/matrix using "n_thr[t]" to index correctly)
          // array[n_index_tests] vector[n_cat_max - 1] dirichlet_cat_means_phi_raw_mat;//// Need to construct the simplex'es MANUALLY since (potentially) using staggered-arrays
          vector[n_total_raw_simplex_elements] dirichlet_cat_means_phi_raw_vec;
          vector<lower=kappa_lb>[n_index_tests] kappa;  
          ////
          //// "NMA" params:
          ////
          vector[n_studies] beta_eta_z;       //// Standard normal RVs for study-level effects - eta[s, 1:2] ~ multi_normal({0, 0}, Sigma).
          real<lower=0.0>   beta_sigma;                    //// Variance component - Between-study SD (Nyaga's σ)
          vector[n_studies] raw_scale_eta_z;  //// Standard normal RVs for study-level effects - eta[s, 1:2] ~ multi_normal({0, 0}, Sigma).
          real<lower=0.0>   raw_scale_sigma;               //// Variance component - Between-study SD (Nyaga's σ)
          ////
          vector<lower=0>[n_index_tests]   beta_tau;             //// Variance component - Test-specific SD (Nyaga's τ) - delta_{s, c, t} ~ normal(0, tau_{c, t}).
          matrix[n_studies, n_index_tests] beta_delta_z;       //// Standard normal RVs for test-specific effects
          vector<lower=0>[n_index_tests]   raw_scale_tau;        //// Variance component - Test-specific SD (Nyaga's τ) - delta_{s, c, t} ~ normal(0, tau_{c, t}).
          matrix[n_studies, n_index_tests] raw_scale_delta_z;  //// Standard normal RVs for test-specific effects
        
}


transformed parameters {
    
          ////
          //// -------- Cutpoint params:
          ////
          array[n_index_tests] matrix[n_studies, n_thr_max] C;  //// Global cutpoints for each test ("staggered" array/matrix using "n_thr[t]" to index correctly) - "staggered" array
          array[n_index_tests] vector[n_cat_max] dirichlet_cat_means_phi; //// "staggered" array
          ////
          //// Initialise matrices:
          ////
          for (t in 1:n_index_tests) {
                C[t] = rep_matrix(0.0, n_studies, n_thr_max);
                dirichlet_cat_means_phi[t] = rep_vector(0.0, n_cat_max);
          }
          ////
          //// ---- Construct cutpoints for each test:
          ////
          real Jacobian_C = 0.0;
          {
                int counter = 1;
                ////
                //// 1st cutpoint (for each test) is unconstrained (no Jacobian needed):
                ////
                for (t in 1:n_index_tests) {
                      for (s in 1:n_studies) {
                         real C_raw = C_raw_vec[counter]; //// C_raw_mat[t][s, 1];
                         C[t][s, 1] = C_raw; //// No Jacobian needd here. 
                         counter += 1;
                      }
                }
                ////
                //// Rest of cutpoints are made using LOG-differences:
                ////
                {
                      for (t in 1:n_index_tests) {
                        
                            int n_thr_t = n_thr[t];
                            for (s in 1:n_studies) {
                              
                                    for (k in 2:n_thr_t) {
                                           if (counter <= num_elements(C_raw_vec)) {
                                               real C_raw = C_raw_vec[counter]; //// C_raw_mat[t][s, k];
                                               C[t][s, k] = C[t][s, k - 1] + exp(C_raw);
                                               //// Jacobian_C += log_inv_logit(C_raw); //// Jacobian for trans. C_raw -> C. (deriv of SP(x) = log(1.0 + exp(x)) is inv_logit(x) so log(deriv) = log_inv_logit(x)). 
                                               Jacobian_C += C_raw;
                                               counter += 1;
                                           }
                                    }
                                    
                             }
                        
                      }
                }
          }
          ////
          //// ---- Construct the Dirichlet simplex for each test (manually as using potentially "staggered" arrays):
          ////
          real Jacobian_simplex = 0.0;
          { 
                //// Now create the simplex for the t-th test using the "stickbreaking_logistic_simplex_constrain" fn:
                int counter = 1;
                for (t in 1:n_index_tests) {
                       int n_cat_t = n_cat[t];
                       vector[n_cat_t - 1] simplex_raw_vec; /// = dirichlet_cat_means_phi_raw_vec; //// dirichlet_cat_means_phi_raw_mat[t][1:(n_cat[t] - 1)];
                       for (k in 1:(n_cat_t - 1)) {
                           simplex_raw_vec[k] = dirichlet_cat_means_phi_raw_vec[counter];
                           counter += 1;
                       }
                       vector[n_cat_t + 1] out_vec = stickbreaking_logistic_simplex_constrain(simplex_raw_vec);
                       // vector[n_cat_t + 1] out_vec = stickbreaking_logistic_simplex_constrain( dirichlet_cat_means_phi_raw_mat[t][1:(n_cat[t] - 1)] );
                       Jacobian_simplex += out_vec[1];
                       vector[n_cat_t] simplex_out = out_vec[2:(n_cat_t + 1)];
                       dirichlet_cat_means_phi[t][1:n_cat_t] = simplex_out;
                } 
          }
          ////
          ////
          ////
          vector[n_index_tests] log_kappa = log(kappa);
          array[n_index_tests] vector[n_cat_max] log_dirichlet_cat_means_phi;
          for (t in 1:n_index_tests) {
              int n_cat_t = n_cat[t];
              log_dirichlet_cat_means_phi[t][1:n_cat_t] = log(dirichlet_cat_means_phi[t][1:n_cat_t]);
          }
          ////
          //// Compute alpha:
          ////
          array[n_index_tests] vector[n_cat_max] log_alpha;
          array[n_index_tests] vector[n_cat_max] alpha;
          vector[n_index_tests] alpha_0;
          for (t in 1:n_index_tests) {
              int n_cat_t = n_cat[t];
              log_alpha[t][1:n_cat_t] = log_kappa[t] + log_dirichlet_cat_means_phi[t][1:n_cat_t];
              alpha[t][1:n_cat_t]     = exp(log_alpha[t][1:n_cat_t]);
              alpha_0[t]   = sum(alpha[t][1:n_cat_t]);
          }
          ////
          //// Other quantities needed (for Jacobian for the (alpha_k, alpha_0) -> dirichlet_cat_SDs_sigma transformation):
          ////
          real log_half = log(0.5);
          ////
          //// NOTE: SD for each category probability in a Dirichlet is √(α_k(α_0-α_k)/(α_0^2(α_0+1))) - where α_0 is the sum of all alphas:
          ////
          array[n_index_tests] vector[n_cat_max] dirichlet_cat_SDs_sigma; ////  = sqrt((alpha[t] .* (alpha_0 - alpha[t]) ) ./ (alpha_0_sq * (alpha_0 + 1.0))); // We are putting a prior on this!!
          for (t in 1:n_index_tests) {
               int n_cat_t = n_cat[t];
               real alpha_0_sq = square(alpha_0[t]);
               dirichlet_cat_SDs_sigma[t][1:n_cat_t] = sqrt((alpha[t][1:n_cat_t] .* (alpha_0[t] - alpha[t][1:n_cat_t]) ) ./ (alpha_0_sq * (alpha_0[t] + 1.0))); // We are putting a prior on this!!
          }
          ////
          //// Then compute the Jacobian adjustment:
          ////
          real Jacobian_for_alpha_k_to_category_SDs = 0.0;
          for (t in 1:n_index_tests) {
              int n_cat_t = n_cat[t];
              real alpha_0_t = alpha_0[t];
              real alpha_0_sq = square(alpha_0_t);
              real alpha_0_cube = alpha_0_sq * alpha_0_t;
              vector[n_cat_t] alpha_test_t = alpha[t][1:n_cat_t];
              vector[n_cat_t] alpha_t_sq = square(alpha_test_t);
              ////
              vector[n_cat_t] log_sigma = log(dirichlet_cat_SDs_sigma[t][1:n_cat_t]);
              vector[n_cat_t] deriv_A = alpha_t_sq + alpha_0_t - 2.0*alpha_test_t;
              vector[n_cat_t] deriv_B = (1.0 ./ alpha_test_t) .* (3.0*alpha_0_sq*alpha_test_t + 2.0*alpha_test_t*alpha_0_t);
              vector[n_cat_t] deriv_var  = deriv_A .* (alpha_0_cube + alpha_0_sq) + deriv_B .* (alpha_test_t*alpha_0_t - alpha_t_sq);
              vector[n_cat_t] log_abs_deriv_var_wrt_alpha = log(abs(deriv_var));
              vector[n_cat_t] log_abs_deriv_SD_wrt_alpha = log_half - log_sigma + log_abs_deriv_var_wrt_alpha;
              Jacobian_for_alpha_k_to_category_SDs += sum(log_abs_deriv_SD_wrt_alpha);
          }
          ////
          ////
          ////
          array[n_index_tests] matrix[2, n_studies] locations;
          array[n_index_tests] matrix[2, n_studies] scales;
          ////
          array[n_index_tests, 2] matrix[n_studies, n_thr_max] logit_cumul_prob; // Ordinal probs for the likelihood (staggered array/matrix using "n_thr[t]" to index correctly)
          array[n_index_tests, 2] matrix[n_studies, n_thr_max] cumul_prob;       // Ordinal probs for the likelihood (staggered array/matrix using "n_thr[t]" to index correctly)
          array[n_index_tests, 2] matrix[n_studies, n_thr_max] cond_prob;        // Ordinal probs for the likelihood (staggered array/matrix using "n_thr[t]" to index correctly)
          array[n_index_tests, 2] matrix[n_studies, n_thr_max] log_lik;          // log_lik storage (staggered array/matrix using "n_thr[t]" to index correctly)
          ////
          array[n_index_tests] matrix[n_studies, n_thr_max] Ind_Dir_anchor;
          array[n_index_tests] matrix[n_studies, n_thr_max] Ind_Dir_cumul_prob;
          array[n_index_tests] matrix[n_studies, n_cat_max] Ind_Dir_ord_prob;
          ////
          //// Initialise matrices:
          ////
          for (t in 1:n_index_tests) { 
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
          for (t in 1:n_index_tests) { 
                   //// Ind-Dir cumulative probs:
                   Ind_Dir_anchor[t]    = rep_matrix(0.0, n_studies, n_thr_max);
                   for (s in 1:n_studies) {
                         Ind_Dir_cumul_prob[t][s, 1:n_thr_max] = Phi((C[t][s, 1:n_thr_max] - Ind_Dir_anchor[t][s, 1:n_thr_max]));
                         //// Ind-Dir ordinal probs:
                         Ind_Dir_ord_prob[t][s, 1] = Ind_Dir_cumul_prob[t][s, 1] - 0.0;
                         for (k in 2:n_thr[t]) {
                             Ind_Dir_ord_prob[t][s, k] = Ind_Dir_cumul_prob[t][s, k] - Ind_Dir_cumul_prob[t][s, k - 1]; // since probs are INCREASING with k
                         }
                         Ind_Dir_ord_prob[t][s, n_cat[t]] =  1.0 - Ind_Dir_cumul_prob[t][s, n_cat[t] - 1];
                   }
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
          matrix[n_studies, n_index_tests] beta_delta;
          matrix[n_studies, n_index_tests] raw_scale_delta;
          ////
          //// Compute test-specific deviations ("delta" in Nyaga notation) - delta_{s, c, t} ~ normal(0, tau_{c, t}):
          ////
          for (t in 1:n_index_tests) {
              for (s in 1:n_studies) {
                  if (indicator_index_test_in_study[s, t] == 1) {
                        beta_delta[s, t]      = 0.0 + beta_delta_z[s, t]      * beta_tau[t];       ////  NOTE: 1st delta's ("beta_delta") correspond to shared (between tests) component of "beta"   -   delta_{s, c, t} ~ normal(0, tau_{c, t}):
                        raw_scale_delta[s, t] = 0.0 + raw_scale_delta_z[s, t] * raw_scale_tau[t];  ////  NOTE: 2nd delta's ("raw_scale_delta") correspond to shared (between tests) component of "raw_scale" - delta_{s, c, t} ~ normal(0, tau_{c, t}):            }
                  }
              }
          }
          ////  
          //// Compute probit CUMULATIVE probabilities:
          ////
          real Jacobian_raw_scale_to_scale = 0.0;
          for (t in 1:n_index_tests) {
                for (s in 1:n_studies) {
                      if (indicator_index_test_in_study[s, t] == 1) {
                        
                                  for (cut_i in 1:n_obs_cutpoints[s, t]) { //// only loop through the OBSERVED cutpoints for test t in study s
                                        
                                        int k = to_int(cutpoint_index[t, 1][s, cut_i]);
                                        
                                        real raw_beta_baseline  = beta_mu[t]      + beta_eta[s]      + beta_delta[s, t];      //// pi_{s, c, t} = Phi(mu_{c, t} + eta_{s, c} + delta_{s, c, t}).
                                        real raw_scale_baseline = raw_scale_mu[t] + raw_scale_eta[s] + raw_scale_delta[s, t]; //// pi_{s, c, t} = Phi(mu_{c, t} + eta_{s, c} + delta_{s, c, t}).
                                        
                                        //// Use opposite signs for non-diseased vs diseased:
                                        locations[t][1, s] =  (-1)*(-0.5)*raw_beta_baseline; //// D- group ////  -0.5*(beta_mu + beta_SD * beta_z[s]);
                                        locations[t][2, s] =  (-1)*(+0.5)*raw_beta_baseline; //// D+ group ////  +0.5*(beta_mu + beta_SD * beta_z[s]);
                                        scales[t][1, s] = exp(-0.5*raw_scale_baseline); //// D- group //// log1p_exp(-0.5*(raw_scale_mu + raw_scale_SD * raw_scale_z[s]));
                                        scales[t][2, s] = exp(+0.5*raw_scale_baseline); //// D+ group //// log1p_exp(+0.5*(raw_scale_mu + raw_scale_SD * raw_scale_z[s]));
                                        //// Jacobian for raw_scale -> scale:
                                        Jacobian_raw_scale_to_scale += +log(2) - 0.5*raw_scale_baseline; //// deriv of exp(-0.5*raw_scale_baseline) = -0.5*exp(-0.5*raw_scale_baseline) -> log_abs_deriv = +log(2) -0.5*raw_scale_baseline; //// log_inv_logit(raw_scale_nd);
                                        Jacobian_raw_scale_to_scale += +log(2) + 0.5*raw_scale_baseline; //// deriv of exp(+0.5*raw_scale_baseline) = +0.5*exp(+0.5*raw_scale_baseline) -> log_abs_deriv = -log(2) +0.5*raw_scale_baseline; //// log_inv_logit(raw_scale_nd);
                                        
                                        //// logit_cumul_prob[c, t, s, cut_i] = (C[t][k] - location)/scale;
                                        for (c in 1:2) {
                                            logit_cumul_prob[t, c][s, cut_i] = (C[t][s, k] - locations[t][c, s])/scales[t][c, s];
                                        }
                                    
                                  }
                      }
                      
               } 
          }
          
            // if (abs(sum(raw_scale_SD)) != 0.0) target += log(abs(sum(raw_scale_SD)));      // double-checked the log-derivative of this by hand (correct)
            // if (abs(sum(raw_scale_z)) != 0.0)  target += log(abs(sum(raw_scale_z)));  // double-checked the log-derivative of this by hand (correct)
            
          ////
          //// Calculate CUMULATIVE probabilities (vectorised):
          ////
          for (t in 1:n_index_tests) {
              for (c in 1:2) {
                  cumul_prob[t, c] = Phi(logit_cumul_prob[t, c]); //// INCREASING sequence (as C_k > C_{k - 1})
              }
          }
          ////
          //// ------- Likelihood using binomial factorization:
          ////
          for (t in 1:n_index_tests) {
                    for (s in 1:n_studies) {
                           for (c in 1:2) {
                          
                                 //// Log-Likelihood is evaluated only at the OBSERVED data:
                                 if (indicator_index_test_in_study[s, t] == 1) {
                                              
                                                    for (cut_i in 1:n_obs_cutpoints[s, t]) {
                                                              
                                                                  //// Current and next cumulative counts:
                                                                  int x_current = x[t, c, s, cut_i];
                                                                  int x_next    = n[t, c, s, cut_i];
                                                                  
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
                                                                        
                                                                        if (cond_prob[t, c][s, cut_i] < 1.00001) {
                                                                        
                                                                              //// Binomial for observations at or below current cutpoint out of those at or below next cutpoint:
                                                                              log_lik[t, c][s, cut_i] = binomial_lpmf(x_current | x_next, cond_prob[t, c][s, cut_i]);
                                                                        
                                                                        }
                                                                        
                                                                  }
                                                    
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
          for (t in 1:n_index_tests) {
              beta_mu[t]      ~ normal(prior_beta_mu_mean[t], prior_beta_mu_SD[t]);
              raw_scale_mu[t] ~ normal(prior_raw_scale_mu_mean[t], prior_raw_scale_mu_SD[t]);
          }
          ////
          for (t in 1:n_index_tests) {
               beta_tau[t]        ~ normal(0.0, prior_beta_tau_SD[t]);       //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
               raw_scale_tau[t]   ~ normal(0.0, prior_raw_scale_tau_SD[t]);  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
          }
          ////
          beta_sigma      ~ normal(0.0, prior_beta_sigma_SD);      //// eta[s, i] ~ normal(0, sigma[i]):
          raw_scale_sigma ~ normal(0.0, prior_raw_scale_sigma_SD); //// eta[s, i] ~ normal(0, sigma[i]):
          ////
          //// Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
          ////
          for (t in 1:n_index_tests) {
                int n_cat_t = n_cat[t];
                target += dirichlet_lpdf( dirichlet_cat_means_phi[t][1:n_cat_t] | prior_dirichlet_cat_means_alpha[t][1:n_cat_t] ); // "flat" prior on the simplex dirichlet_cat_means_phi. 
                target += normal_lpdf(dirichlet_cat_SDs_sigma[t][1:n_cat_t] | prior_dirichlet_cat_SDs_mean[t][1:n_cat_t], prior_dirichlet_cat_SDs_SD[t][1:n_cat_t] );
          }
          //// Increment Jacobian adjustments:
          target += Jacobian_C; //// For C_raw -> C
          target += Jacobian_simplex; //// For simplex_raw -> simplex
          target += Jacobian_for_alpha_k_to_category_SDs; //// For (alpha, kappa) -> Dirichlet_SD's
          ////
          for (t in 1:n_index_tests) {
                for (s in 1:n_studies) {
                    if (indicator_index_test_in_study[s, t] == 1) {
                        int n_thr_t = n_thr[t];
                        int n_cat_t = n_cat[t];
                        vector[n_thr_t] rho =  normal_pdf(to_vector(Ind_Dir_cumul_prob[t][s, 1:n_thr_t]), 0.0, 1.0);
                        target += induced_dirichlet_v2_lpdf(to_vector(Ind_Dir_ord_prob[t][s, 1:n_cat_t]) | rho, alpha[t][1:n_cat_t]);
                    }
                }
                ////
                for (k in 1:n_cat[t]) {
                    target += log_kappa[t];
                    target += log_dirichlet_cat_means_phi[t][k];
                }
          }
          ////
          //// Likelihood / Model:
          //// 
          target += Jacobian_raw_scale_to_scale;
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
            for (t in 1:n_index_tests) {
              for (c in 1:2) {
                target +=  sum(log_lik[t, c]);
              }
            }
          }
          dirichlet_cat_means_phi_raw_vec ~ normal(0, 10); //// TEMP
          C_raw_vec ~ normal(0, 10); //// TEMP
    
}




 


generated quantities {

          //// Summary accuracy parameters (at each threshold):
          array[n_index_tests] vector[n_thr_max] Se;
          array[n_index_tests] vector[n_thr_max] Sp;
          array[n_index_tests] vector[n_thr_max] Fp;
          ////
          array[n_index_tests] vector[n_thr_max] Se_pred;
          array[n_index_tests] vector[n_thr_max] Sp_pred;
          array[n_index_tests] vector[n_thr_max] Fp_pred;
          array[n_index_tests] vector[n_thr_max] C_mu_pred;
          ////
          array[n_index_tests] matrix[n_studies, n_thr_max] se;
          array[n_index_tests] matrix[n_studies, n_thr_max] sp;
          array[n_index_tests] matrix[n_studies, n_thr_max] fp;
          ////
          array[n_index_tests] matrix[n_studies, n_thr_max] x_hat_nd;// = rep_matrix(-1, n_studies, n_thr_max);
          array[n_index_tests] matrix[n_studies, n_thr_max] x_hat_d;//  = rep_matrix(-1, n_studies, n_thr_max);
          array[n_index_tests] matrix[n_studies, n_thr_max] dev_nd;//    = rep_matrix(-1, n_studies, n_thr_max);
          array[n_index_tests] matrix[n_studies, n_thr_max] dev_d;//     = rep_matrix(-1, n_studies, n_thr_max);
          array[2, n_index_tests] matrix[n_studies, n_thr_max] x_hat;
          array[2, n_index_tests] matrix[n_studies, n_thr_max] dev;
          ////
          array[n_index_tests] vector[n_thr_max] C_MU;
          array[n_index_tests] vector[n_thr_max] C_MU_empirical;
          // array[n_index_tests] vector[n_thr_max] prob_cumul_mu;
          // real<lower=kappa_lb> precision = kappa;
          // real dispersion = 1.0 / precision;
          array[n_index_tests] vector[n_thr_max] expected_p_cumul;
          ////
          //// Initialise containers:
          ////
          for (t in 1:n_index_tests) {
            for (c in 1:2) {
                x_hat[t, c]  = rep_matrix(-1, n_studies, n_thr_max);
                dev[t, c]    = rep_matrix(-1, n_studies, n_thr_max);
            }
          }
          ////
          ////
          ////
          for (t in 1:n_index_tests) {
                //// Calculate cumulative probabilities:
                expected_p_cumul[t][1] = dirichlet_cat_means_phi[t][1];
                for (k in 2:n_thr[t]) {
                    expected_p_cumul[t][k] = expected_p_cumul[t][k - 1] + dirichlet_cat_means_phi[t][k];
                }

                // Transform to cutpoints:
                for (k in 1:n_thr[t]) {
                       real prob = expected_p_cumul[t][k];
                      if (prob < 0.0001) {
                            C_MU[t][k] = -10.0;
                      } else if (prob > 0.99999) {
                            C_MU[t][k] = +10.0;
                      } else {
                            C_MU[t][k] = 0.0 + inv_Phi(prob);
                      }

                }
          }
          ////
          //// Generate a prediction for each cutpoint:
          ////
          for (t in 1:n_index_tests) {

                  vector[n_cat[t]] prob_ord_mu_t_pred;
                  vector[n_thr[t]] prob_cumul_mu_t_pred;
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
                  prob_ord_mu_t_pred[1:n_cat[t]]   =  dirichlet_rng(alpha[t]);
                  ////
                  //// Calculate cumulative probabilities:
                  prob_cumul_mu_t_pred[1] = prob_ord_mu_t_pred[1];
                  for (k in 2:n_thr[t]) {
                      prob_cumul_mu_t_pred[k] = prob_cumul_mu_t_pred[k - 1] + prob_ord_mu_t_pred[k];
                  }
                  ////
                  //// Transform to cutpoints:
                  real anchor_for_summaries = 0.0;
                  for (k in 1:n_thr[t]) {
                        real prob_1;
                        if (prob_cumul_mu_t_pred[k] < 1e-38) {
                              prob_1 = 1e-38;
                              C_mu_pred[t][k] = -10.0;
                        } else if (prob_cumul_mu_t_pred[k] > 0.9999999999999) {
                              prob_1 = 0.9999999999999;
                              C_mu_pred[t][k] = +10.0;
                        } else {
                              prob_1 = inv_Phi(prob_cumul_mu_t_pred[k]);
                              C_mu_pred[t][k] = anchor_for_summaries + prob_1;
                        }

                  }

          }
          ////
          //// Empirical-mean cutpoints:
          ////
          for (t in 1:n_index_tests) {
              for (k in 1:n_thr[t]) {
                    C_MU_empirical[t][k] = median(C[t][, k]);
              }
          }

          ////
          //// Calculate study-specific accuracy:
          ////
          for (s in 1:n_studies) {
              for (t in 1:n_index_tests) {
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
          for (t in 1:n_index_tests) {
            for (k in 1:n_thr[t]) {
                  Fp[t][k] =   1.0 - Phi((C_MU[t][k] - (-0.5)*beta_mu[t])/log1p_exp((-0.5)*raw_scale_mu[t]));
                  Sp[t][k] =   1.0 - Fp[t][k];
                  Se[t][k] =   1.0 - Phi((C_MU[t][k] - (+0.5)*beta_mu[t])/log1p_exp((+0.5)*raw_scale_mu[t]));
            }
          }
          ////
          //// Calculate predictive accuracy:
          ////
          real beta_eta_pred      = normal_rng(0.0, beta_sigma);      //// shared between tests
          real raw_scale_eta_pred = normal_rng(0.0, raw_scale_sigma); //// shared between tests
          ////
          for (t in 1:n_index_tests) {
                ////
                real beta_delta_t_pred      = normal_rng(0.0, beta_tau[t]);
                real raw_scale_delta_t_pred = normal_rng(0.0, raw_scale_tau[t]);
                ////
                real beta_pred      = beta_eta_pred       + beta_delta_t_pred;
                real raw_scale_pred = raw_scale_eta_pred  + raw_scale_delta_t_pred;
                ////
                for (k in 1:n_thr[t]) {
                      Fp_pred[t][k] =   1.0 - Phi((C_mu_pred[t][k] - (-0.5)*beta_pred)/log1p_exp((-0.5)*raw_scale_pred));
                      Sp_pred[t][k] =   1.0 - Fp_pred[t][k];
                      Se_pred[t][k] =   1.0 - Phi((C_mu_pred[t][k] - (+0.5)*beta_pred)/log1p_exp((+0.5)*raw_scale_pred));
                }
          }
          ////
          //// Model-predicted ("re-constructed") data:
          ////
          for (t in 1:n_index_tests) {
              for (s in 1:n_studies) {
                    if (indicator_index_test_in_study[s, t] == 1) {

                              for (c in 1:2) {
                                 for (cut_i in 1:to_int(n_obs_cutpoints[s, t])) {

                                          //// Model-estimated data:
                                          x_hat[c][t][s, cut_i] = cond_prob[c][t][s, cut_i] * n[c][t][s, cut_i];  	 //// Fitted values

                                          //// Compute residual deviance contribution:
                                          real n_i =  (n[t, c, s, cut_i]);
                                          real x_i =  (x[t, c, s, cut_i]);
                                          real x_hat_i =  (x_hat[c][t][s, cut_i]);
                                          real log_x_minus_log_x_hat = log(x_i) - log(x_hat_i);
                                          real log_diff_n_minus_x = log(n_i - x_i);
                                          real log_diff_n_minus_x_hat = log(abs(n_i - x_hat_i));

                                          dev[c][t][s, cut_i] = 2.0 * ( x_i * log_x_minus_log_x_hat + (n_i - x_i) * (log_diff_n_minus_x - log_diff_n_minus_x_hat) );

                                 }
                              }

                    }
              }
          }
          ////
          //// Store deviance and x_hat split by disease status:
          ////
          for (t in 1:n_index_tests) {
               x_hat_nd[t] = x_hat[t][1];
               dev_nd[t]   = dev[t][1];
               x_hat_d[t]  = x_hat[t][2];
               dev_d[t]    = dev[t][2];
          }

}














