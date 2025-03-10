// 
// 
// 
// 

data {

    int<lower=1> n_studies;
    int<lower=1> n_tests;          // Total number of tests across all studies
    int<lower=1> n_thr;            // Maximum number of thresholds
    array[n_studies] int<lower=1> n_tests_per_study;  // Number of tests in each study
    array[n_studies, n_tests] int test_in_study;      // Binary indicator if test t is in study s
    array[2, n_studies, n_tests] int n_cutpoints;     // Cutpoints for each test in each study
    array[2] matrix[n_studies, n_tests, n_thr] x_with_missings;
    array[2] matrix[n_studies, n_tests, n_thr] n;
    array[2] matrix[n_studies, n_tests, n_thr] x;
    array[2] matrix[n_studies, n_tests, n_thr] cutpoint_index;
    // Priors...

}





parameters {

    // Test-specific parameters
    vector[n_tests] beta_mu;      // Mean effect for each test
    vector<lower=0>[n_tests] beta_SD;  // Between-study SD for each test
    matrix[n_studies, n_tests] beta_z;  // Study-specific z-scores for each test

    // Scale parameters
    vector[n_tests] raw_scale_mu;
    vector<lower=0>[n_tests] raw_scale_SD;
    matrix[n_studies, n_tests] raw_scale_z;

    // Between-study correlation
    cholesky_factor_corr[2] bs_L_Omega;

    // Global cutpoints for each test
    array[n_tests] ordered[n_thr] C;  // Test-specific global cutpoints

}






transformed parameters {

    array[n_tests] matrix[2, n_studies] locations;
    array[n_tests] matrix[2, n_studies] scales;

    // Between-study model for each test
    for (t in 1:n_tests) {
        for (s in 1:n_studies) {
            if (test_in_study[s, t] == 1) {
                locations[t, 1, s] = -0.5 * (beta_mu[t] + beta_SD[t] * beta_z[s, t]);
                locations[t, 2, s] = +0.5 * (beta_mu[t] + beta_SD[t] * beta_z[s, t]);

                scales[t, 1, s] = log1p_exp(-0.5 * (raw_scale_mu[t] + raw_scale_SD[t] * raw_scale_z[s, t]));
                scales[t, 2, s] = log1p_exp(+0.5 * (raw_scale_mu[t] + raw_scale_SD[t] * raw_scale_z[s, t]));
            }
        }
    }

    // Calculate ordinal probabilities for each test in each study
    array[2, n_tests] matrix[n_studies, n_thr] logit_cumul_prob;
    array[2, n_tests] matrix[n_studies, n_thr] cumul_prob;
    array[2, n_tests] matrix[n_studies, n_thr] cond_prob;
    array[2, n_tests] matrix[n_studies, n_thr] log_lik;

    // Initialize matrices
    for (t in 1:n_tests) {
        for (c in 1:2) {
            logit_cumul_prob[c, t] = rep_matrix(positive_infinity(), n_studies, n_thr);
            cumul_prob[c, t] = rep_matrix(1.0, n_studies, n_thr);
            cond_prob[c, t] = rep_matrix(1.0, n_studies, n_thr);
            log_lik[c, t] = rep_matrix(0.0, n_studies, n_thr);
        }
    }

    // Calculate probabilities for each test in each study
    for (s in 1:n_studies) {
        for (t in 1:n_tests) {
            if (test_in_study[s, t] == 1) {
                for (c in 1:2) {
                    for (cut_i in 1:n_cutpoints[c, s, t]) {
                        int k = cutpoint_index[c, s, t, cut_i];
                        logit_cumul_prob[c, t, s, cut_i] = (C[t, k] - locations[t, c, s]) / scales[t, c, s];
                    }
                }
            }
        }
    }

    // Calculate cumulative probabilities
    for (t in 1:n_tests) {
        for (c in 1:2) {
            cumul_prob[c, t] = Phi(logit_cumul_prob[c, t]);
        }
    }

    // Calculate likelihood using binomial factorization
    for (s in 1:n_studies) {
        for (t in 1:n_tests) {
            if (test_in_study[s, t] == 1) {
                for (c in 1:2) {
                    for (cut_i in 1:n_cutpoints[c, s, t]) {
                        // Calculate conditional probabilities and likelihood as in your original model
                        // ...
                    }
                }
            }
        }
    }

}






model {

    // Priors
    for (t in 1:n_tests) {
        beta_mu[t] ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
        beta_SD[t] ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
        raw_scale_mu[t] ~ normal(prior_raw_scale_mu_mean, prior_raw_scale_mu_SD);
        raw_scale_SD[t] ~ normal(prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
    }

    // Between-study correlation
    bs_L_Omega ~ lkj_corr_cholesky(2.0);

    // Induced-dirichlet prior for each test's cutpoints
    for (t in 1:n_tests) {
        vector[n_thr] rho = normal_pdf(C[t], 0.0, 1.0);
        target += induced_dirichlet_v2_lpdf(ord_prob_from_C(C[t]) | rho, prior_alpha);
    }

    // Likelihood
    for (t in 1:n_tests) {
        for (s in 1:n_studies) {
            if (test_in_study[s, t] == 1) {
                beta_z[s, t] ~ std_normal();
                raw_scale_z[s, t] ~ std_normal();
            }
        }
    }

    // Increment log-likelihood
    if (prior_only == 0) {
        for (t in 1:n_tests) {
            for (c in 1:2) {
                target += sum(log_lik[c, t]);
            }
        }
    }

}




// 
// 
// 
// 
// 
// generated quantities {
//   
//     // Summary accuracy parameters for each test at each threshold
//     array[n_tests] vector[n_thr] Se;
//     array[n_tests] vector[n_thr] Sp;
//     array[n_tests] vector[n_thr] Fp;
// 
//     // Study-specific accuracy for each test
//     array[n_tests] matrix[n_studies, n_thr] se;
//     array[n_tests] matrix[n_studies, n_thr] sp;
//     array[n_tests] matrix[n_studies, n_thr] fp;
// 
//     // Calculate test-specific summary measures
//     for (t in 1:n_tests) {
//         for (k in 1:n_thr) {
//             Fp[t, k] = 1.0 - Phi((C[t, k] - 0.5*beta_mu[t])/log1p_exp(-0.5*raw_scale_mu[t]));
//             Sp[t, k] = 1.0 - Fp[t, k];
//             Se[t, k] = 1.0 - Phi((C[t, k] + 0.5*beta_mu[t])/log1p_exp(+0.5*raw_scale_mu[t]));
//         }
//     }
// 
//     // Calculate study-specific accuracy for each test
//     for (t in 1:n_tests) {
//         for (s in 1:n_studies) {
//             if (test_in_study[s, t] == 1) {
//                 for (k in 1:n_thr) {
//                     fp[t, s, k] = 1.0 - cumul_prob[1, t, s, k];
//                     sp[t, s, k] = 1.0 - fp[t, s, k];
//                     se[t, s, k] = 1.0 - cumul_prob[2, t, s, k];
//                 }
//             }
//         }
//     }
// 
//     // Calculate relative effects between tests
//     matrix[n_tests, n_tests] test_differences;
//     for (t1 in 1:n_tests) {
//         for (t2 in 1:n_tests) {
//             test_differences[t1, t2] = beta_mu[t1] - beta_mu[t2];
//         }
//     }
//     
// }
// 
// 
// 











data {
  
    int<lower=1> n_studies;
    int<lower=1> n_tests;
    int<lower=1> n_thr;
    array[n_studies] int<lower=0, upper=n_tests> n_tests_per_study;
    array[n_studies, n_tests] int<lower=0, upper=1> test_in_study;
    array[2, n_studies, n_tests] int<lower=0> n_cutpoints;
    array[2, n_studies, n_tests, n_thr] int x_with_missings;
    array[2, n_studies, n_tests, n_thr] int n;
    array[2, n_studies, n_tests, n_thr] int x;
    array[2, n_studies, n_tests, n_thr] int cutpoint_index;
    
    // Priors
    real prior_beta_mu_mean;
    real prior_beta_mu_SD;
    real prior_beta_SD_mean;
    real prior_beta_SD_SD;
    real prior_raw_scale_mu_mean;
    real prior_raw_scale_mu_SD;
    real prior_raw_scale_SD_mean;
    real prior_raw_scale_SD_SD;
    vector<lower=0>[n_thr + 1] prior_alpha;
    int prior_only;
    
    // Reference test (optional - for identifiability)
    int<lower=1, upper=n_tests> reference_test;
    
}



parameters {
  
    // Test-specific mean parameters
    vector[n_tests-1] beta_mu_raw;  // Using n_tests-1 to enable constraints
    
    // Between-study variation parameters (test-specific)
    vector<lower=0>[n_tests] beta_SD;
    
    // Study random effects (shared across all tests for each disease status)
    matrix[2, n_studies] eta;  // This is η in Nyaga's notation
    
    // Test-specific deviations within studies
    array[n_tests] matrix[2, n_studies] delta;  // This is δ in Nyaga's notation
    
    // Correlation structure for study random effects
    cholesky_factor_corr[2] L_Omega;
    
    // Parameters for scale
    vector[n_tests] raw_scale_mu;
    vector<lower=0>[n_tests] raw_scale_SD;
    
    // Cutpoints for each test
    array[n_tests] ordered[n_thr] C;
    
}




transformed parameters {
  
    // Full test-specific effects
    vector[n_tests] beta_mu;
    
    // Incorporate the constraint on reference test
    beta_mu[reference_test] = 0;
    for (t in 1:(reference_test-1))
        beta_mu[t] = beta_mu_raw[t];
    for (t in (reference_test+1):n_tests)
        beta_mu[t] = beta_mu_raw[t-1];
    
    // Construct covariance matrix from cholesky factor
    matrix[2, 2] Sigma = multiply_lower_tri_self_transpose(L_Omega);
    
    // Storage for test-specific parameters in each study
    array[n_tests] matrix[2, n_studies] locations;
    array[n_tests] matrix[2, n_studies] scales;
    
    // Likelihood components
    array[2, n_tests, n_studies, n_thr] real logit_cumul_prob;
    array[2, n_tests, n_studies, n_thr] real cumul_prob;
    array[2, n_tests, n_studies, n_thr] real cond_prob;
    array[2, n_tests, n_studies, n_thr] real log_lik;
    
    // Initialize arrays
    for (c in 1:2)
        for (t in 1:n_tests)
            for (s in 1:n_studies)
                for (k in 1:n_thr) {
                    logit_cumul_prob[c,t,s,k] = positive_infinity();
                    cumul_prob[c,t,s,k] = 1.0;
                    cond_prob[c,t,s,k] = 1.0;
                    log_lik[c,t,s,k] = 0.0;
                }
    
    // Calculate total effects for each test, study and disease status
    // Following Nyaga: logit(π) = μₖ + ηᵢⱼ + δᵢⱼₖ
    for (t in 1:n_tests) {
        for (s in 1:n_studies) {
            if (test_in_study[s,t] == 1) {
                // Construct locations following Nyaga's model structure
                // For non-diseased (specificity) - c=1
                locations[t][1,s] = -0.5*beta_mu[t] + eta[1,s] + delta[t][1,s];
                
                // For diseased (sensitivity) - c=2
                locations[t][2,s] = 0.5*beta_mu[t] + eta[2,s] + delta[t][2,s];
                
                // Scale parameters - not directly in Nyaga's model but needed for your probit approach
                scales[t][1,s] = log1p_exp(-0.5*(raw_scale_mu[t] + raw_scale_SD[t]));
                scales[t][2,s] = log1p_exp(0.5*(raw_scale_mu[t] + raw_scale_SD[t]));
            }
        }
    }
    
    // Calculate ordinal probabilities for likelihood
    for (t in 1:n_tests) {
        for (s in 1:n_studies) {
            if (test_in_study[s,t] == 1) {
                for (c in 1:2) {
                    for (cut_i in 1:n_cutpoints[c,s,t]) {
                        int k = cutpoint_index[c,s,t,cut_i];
                        logit_cumul_prob[c,t,s,cut_i] = (C[t][k] - locations[t][c,s])/scales[t][c,s];
                        cumul_prob[c,t,s,cut_i] = Phi(logit_cumul_prob[c,t,s,cut_i]);
                    }
                }
            }
        }
    }
    
    // Calculate conditional probabilities and log-likelihood
    for (t in 1:n_tests) {
        for (s in 1:n_studies) {
            if (test_in_study[s,t] == 1) {
                for (c in 1:2) {
                    for (cut_i in 1:n_cutpoints[c,s,t]) {
                        // Current and next cumulative counts
                        int x_current = x[c,s,t,cut_i];
                        int x_next = n[c,s,t,cut_i];
                        
                        // Skip if the current count is zero
                        if (x_current != 0) {
                            // Conditional probability for binomial model
                            if (cut_i == n_cutpoints[c,s,t]) {
                                cond_prob[c,t,s,cut_i] = cumul_prob[c,t,s,cut_i] / 1.0;
                            } else {
                                if (x_next > 0) {
                                    cond_prob[c,t,s,cut_i] = cumul_prob[c,t,s,cut_i] / cumul_prob[c,t,s,cut_i+1];
                                } else {
                                    cond_prob[c,t,s,cut_i] = 1.0;
                                }
                            }
                            
                            // Binomial likelihood
                            log_lik[c,t,s,cut_i] = binomial_lpmf(x_current | x_next, cond_prob[c,t,s,cut_i]);
                        }
                    }
                }
            }
        }
    }
    
}



model {
  
    // Priors for test-specific parameters
    beta_mu_raw ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
    beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
    raw_scale_mu ~ normal(prior_raw_scale_mu_mean, prior_raw_scale_mu_SD);
    raw_scale_SD ~ normal(prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
    
    // Prior for correlation matrix
    L_Omega ~ lkj_corr_cholesky(2.0);
    
    // Priors for cutpoints using induced Dirichlet
    for (t in 1:n_tests) {
        vector[n_thr] rho = normal_pdf(to_vector(C[t]), 0.0, 1.0);
        vector[n_thr+1] ord_probs = ord_prob_from_C(to_vector(C[t]));
        target += induced_dirichlet_v2_lpdf(ord_probs | rho, prior_alpha);
    }
    
    // Study random effects model - following Nyaga
    for (s in 1:n_studies) {
        eta[,s] ~ multi_normal_cholesky(rep_vector(0, 2), L_Omega);
    }
    
    // Test-specific effects - following Nyaga
    for (t in 1:n_tests) {
        for (s in 1:n_studies) {
            if (test_in_study[s,t] == 1) {
                for (c in 1:2) {
                    delta[t][c,s] ~ normal(0, beta_SD[t]);
                }
            }
        }
    }
    
    // Likelihood
    if (prior_only == 0) {
        for (t in 1:n_tests) {
            for (s in 1:n_studies) {
                if (test_in_study[s,t] == 1) {
                    for (c in 1:2) {
                        for (cut_i in 1:n_cutpoints[c,s,t]) {
                            if (x[c,s,t,cut_i] > 0) {
                                target += log_lik[c,t,s,cut_i];
                            }
                        }
                    }
                }
            }
        }
    }
}

generated quantities {
    // Summary accuracy parameters for each test
    array[n_tests] vector[n_thr] Se;
    array[n_tests] vector[n_thr] Sp;
    
    // Relative accuracy between tests
    matrix[n_tests, n_tests] relative_accuracy;
    
    // Calculate summary accuracy
    for (t in 1:n_tests) {
        for (k in 1:n_thr) {
            Sp[t][k] = Phi((C[t][k] - (-0.5*beta_mu[t]))/log1p_exp(-0.5*raw_scale_mu[t]));
            Se[t][k] = 1 - Phi((C[t][k] - (0.5*beta_mu[t]))/log1p_exp(0.5*raw_scale_mu[t]));
        }
    }
    
    // Compute relative effects between all pairs of tests
    for (t1 in 1:n_tests) {
        for (t2 in 1:n_tests) {
            relative_accuracy[t1, t2] = beta_mu[t1] - beta_mu[t2];
        }
    }
    
    // Correlation matrix
    matrix[2, 2] Omega = multiply_lower_tri_self_transpose(L_Omega);
}


















functions {
  // Your existing functions
  real normal_pdf(real x, real mu, real sigma) {
    real sqrt_2_pi = sqrt(2 * pi());
    return (1.0 / (sigma * sqrt_2_pi)) * exp(-0.5 * ((x - mu) / sigma)^2);
  }
  
  vector normal_pdf(vector x, real mu, real sigma) {
    real sqrt_2_pi = sqrt(2 * pi());
    return (1.0 / (sigma * sqrt_2_pi)) * exp(-0.5 * ((x - mu) ./ sigma)^2);
  }
  
  real induced_dirichlet_v2_lpdf(vector p_ord, vector rho, vector alpha) {
    int n_cat = num_elements(p_ord);
    matrix[n_cat, n_cat] J = rep_matrix(0.0, n_cat, n_cat);
    
    // Jacobian computation
    for (k in 1:n_cat) {
      J[k, 1] = 1.0;
    }
    for (k in 2:n_cat) {
      J[k, k] = +rho[k - 1];
      J[k - 1, k] = -rho[k - 1];
    }
    
    return dirichlet_lpdf(p_ord | alpha) + log_determinant(J);
  }
  
  // Function to get ordinal probabilities from cutpoints
  vector ord_prob_from_C(vector C) {
    int n_thr = num_elements(C);
    vector[n_thr + 1] p;
    vector[n_thr] cumul_p = Phi(C);
    
    p[1] = cumul_p[1];
    for (k in 2:n_thr)
      p[k] = cumul_p[k] - cumul_p[k-1];
    p[n_thr + 1] = 1 - cumul_p[n_thr];
    
    return p;
  }
}




data {
  int<lower=1> n_studies;                 // Number of studies
  int<lower=1> n_tests;                   // Number of tests
  int<lower=1> n_thr;                     // Number of thresholds
  array[n_studies] int<lower=0, upper=n_tests> n_tests_per_study;  // Tests per study
  array[n_studies, n_tests] int<lower=0, upper=1> test_in_study;   // Binary indicator if test t is in study s
  array[2, n_studies, n_tests] int<lower=0> n_cutpoints;           // Cutpoints for each test in each study
  array[2, n_studies, n_tests, n_thr] int x_with_missings;         // Missing indicator
  array[2, n_studies, n_tests, n_thr] int n;                       // Denominators
  array[2, n_studies, n_tests, n_thr] int x;                       // Numerators
  array[2, n_studies, n_tests, n_thr] int cutpoint_index;          // Cutpoint indices
  
  // Priors
  real prior_beta_mu_mean;
  real prior_beta_mu_SD;
  real prior_tau_SD;                      // Prior for test-specific SD (Nyaga's τ)
  real prior_sigma_SD;                    // Prior for study-level SD (Nyaga's σ)
  real prior_raw_scale_mu_mean;
  real prior_raw_scale_mu_SD;
  vector<lower=0>[n_thr + 1] prior_alpha; // Induced-Dirichlet prior
  int prior_only;
  
  // Reference test (optional)
  int<lower=1, upper=n_tests> reference_test;
}




parameters {
  // Test-specific effects
  vector[n_tests-1] beta_mu_raw;          // Test-specific means
  
  // Study random effects (shared across tests)
  matrix[n_studies, 2] z_nu;              // Standard normal RVs for study effects
  
  // Observation-level random effects
  array[n_tests] matrix[n_studies, 2] z_tau;  // Standard normal RVs for test-specific effects
  
  // Between-study correlation
  cholesky_factor_corr[2] L_Omega;
  
  // Variance components
  vector<lower=0>[2] sigma;               // Between-study SD (Nyaga's σ)
  matrix<lower=0>[n_tests, 2] tau;        // Test-specific SD (Nyaga's τ)
  
  // Scale parameters for probit model
  vector[n_tests] raw_scale_mu;
  
  // Global cutpoints for each test
  array[n_tests] ordered[n_thr] C;
}




transformed parameters {
  // Full test-specific effects vector
  vector[n_tests] beta_mu;
  
  // Study-level random effects (eta in Nyaga notation)
  array[n_studies] vector[2] eta;
  
  // Test-specific deviations (delta in Nyaga notation)
  array[n_tests, n_studies] matrix[2, n_thr] delta;
  
  // Storage for logit probabilities and likelihood
  array[2, n_tests, n_studies, n_thr] real logit_cumul_prob;
  array[2, n_tests, n_studies, n_thr] real cumul_prob;
  array[2, n_tests, n_studies, n_thr] real cond_prob;
  array[2, n_tests, n_studies, n_thr] real log_lik;
  
  // Initialize arrays
  for (t in 1:n_tests)
    for (s in 1:n_studies)
      for (c in 1:2)
        for (k in 1:n_thr) {
          delta[t, s][c, k] = 0.0;
          logit_cumul_prob[c, t, s, k] = positive_infinity();
          cumul_prob[c, t, s, k] = 1.0;
          cond_prob[c, t, s, k] = 1.0;
          log_lik[c, t, s, k] = 0.0;
        }
  
  // Apply constraint on reference test
  beta_mu[reference_test] = 0;
  for (t in 1:(reference_test-1))
    beta_mu[t] = beta_mu_raw[t];
  for (t in (reference_test+1):n_tests)
    beta_mu[t] = beta_mu_raw[t-1];
  
  // Compute study random effects (Nyaga's η) - shared across tests within study
  for (s in 1:n_studies) {
    eta[s] = diag_pre_multiply(sigma, L_Omega) * to_vector(z_nu[s, ]);
  }
  
  // Compute test-specific effects (Nyaga's δ)
  for (t in 1:n_tests) {
    for (s in 1:n_studies) {
      if (test_in_study[s, t] == 1) {
        for (c in 1:2) {
          for (k in 1:n_cutpoints[c, s, t]) {
            delta[t, s][c, k] = z_tau[t][s, c] * tau[t, c];
          }
        }
      }
    }
  }
  
  // Calculate logit cumulative probabilities
  for (t in 1:n_tests) {
    for (s in 1:n_studies) {
      if (test_in_study[s, t] == 1) {
        for (c in 1:2) {
          for (cut_i in 1:n_cutpoints[c, s, t]) {
            int k = cutpoint_index[c, s, t, cut_i];
            
            // Use opposite signs for non-diseased vs diseased
            real location;
            real scale;
            
            if (c == 1) {  // Non-diseased (specificity)
              location = -0.5*beta_mu[t] + eta[s, 1] + delta[t, s][1, cut_i];
              scale = log1p_exp(-0.5*raw_scale_mu[t]);
            } else {       // Diseased (sensitivity)
              location = 0.5*beta_mu[t] + eta[s, 2] + delta[t, s][2, cut_i];
              scale = log1p_exp(0.5*raw_scale_mu[t]);
            }
            
            logit_cumul_prob[c, t, s, cut_i] = (C[t][k] - location)/scale;
            cumul_prob[c, t, s, cut_i] = Phi(logit_cumul_prob[c, t, s, cut_i]);
          }
        }
      }
    }
  }
  
  // Calculate conditional probabilities and log-likelihood
  for (t in 1:n_tests) {
    for (s in 1:n_studies) {
      if (test_in_study[s, t] == 1) {
        for (c in 1:2) {
          for (cut_i in 1:n_cutpoints[c, s, t]) {
            // Current and next cumulative counts
            int x_current = x[c, s, t, cut_i];
            int x_next = n[c, s, t, cut_i];
            
            // Skip if the current count is zero
            if (x_current != 0) {
              // Conditional probability for binomial model
              if (cut_i == n_cutpoints[c, s, t]) {
                cond_prob[c, t, s, cut_i] = cumul_prob[c, t, s, cut_i] / 1.0;
              } else {
                if (x_next > 0) {
                  cond_prob[c, t, s, cut_i] = cumul_prob[c, t, s, cut_i] / cumul_prob[c, t, s, cut_i+1];
                } else {
                  cond_prob[c, t, s, cut_i] = 1.0;
                }
              }
              
              // Binomial likelihood
              log_lik[c, t, s, cut_i] = binomial_lpmf(x_current | x_next, cond_prob[c, t, s, cut_i]);
            }
          }
        }
      }
    }
  }
}




model {
  // Priors for test-specific parameters
  beta_mu_raw ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
  raw_scale_mu ~ normal(prior_raw_scale_mu_mean, prior_raw_scale_mu_SD);
  
  // Priors for variance components (following your Nyaga implementation)
  to_vector(tau) ~ normal(0, prior_tau_SD);
  sigma ~ normal(0, prior_sigma_SD);
  
  // Prior for correlation matrix
  L_Omega ~ lkj_corr_cholesky(2);
  
  // Standard normal priors for z parameters
  to_vector(z_nu) ~ std_normal();
  for (t in 1:n_tests) {
    to_vector(z_tau[t]) ~ std_normal();
  }
  
  // Priors for cutpoints using induced Dirichlet
  for (t in 1:n_tests) {
    vector[n_thr] rho = normal_pdf(to_vector(C[t]), 0.0, 1.0);
    vector[n_thr+1] ord_probs = ord_prob_from_C(to_vector(C[t]));
    target += induced_dirichlet_v2_lpdf(ord_probs | rho, prior_alpha);
  }
  
  // Likelihood
  if (prior_only == 0) {
    for (t in 1:n_tests) {
      for (s in 1:n_studies) {
        if (test_in_study[s, t] == 1) {
          for (c in 1:2) {
            for (cut_i in 1:n_cutpoints[c, s, t]) {
              if (x[c, s, t, cut_i] > 0) {
                target += log_lik[c, t, s, cut_i];
              }
            }
          }
        }
      }
    }
  }
}




generated quantities {
  // Summary accuracy parameters for each test
  array[n_tests] vector[n_thr] Se;
  array[n_tests] vector[n_thr] Sp;
  
  // Relative accuracy between tests
  matrix[n_tests, n_tests] relative_accuracy;
  
  // Correlation matrix
  matrix[2, 2] Omega = multiply_lower_tri_self_transpose(L_Omega);
  
  // Calculate summary accuracy
  for (t in 1:n_tests) {
    for (k in 1:n_thr) {
      Sp[t][k] = Phi((C[t][k] - (-0.5*beta_mu[t]))/log1p_exp(-0.5*raw_scale_mu[t]));
      Se[t][k] = 1 - Phi((C[t][k] - (0.5*beta_mu[t]))/log1p_exp(0.5*raw_scale_mu[t]));
    }
  }
  
  // Compute relative effects between all pairs of tests
  for (t1 in 1:n_tests) {
    for (t2 in 1:n_tests) {
      relative_accuracy[t1, t2] = beta_mu[t1] - beta_mu[t2];
    }
  }
}



















