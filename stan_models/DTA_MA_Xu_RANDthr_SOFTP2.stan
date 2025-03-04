
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
        array[2, n_studies] int n_cutpoints;
        array[2] matrix[n_studies, n_thr] x_with_missings;
        array[2] matrix[n_studies, n_thr] n;
        array[2] matrix[n_studies, n_thr] x;
        array[2] matrix[n_studies, n_thr] cutpoint_index;
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
        array[2] vector[n_thr + 1] prior_dirichlet_phi;
        //// Other:
        int prior_only;
        int estimate_scales;
        real log_alpha_lb;
        real log_alpha_ub;
        ////
  
}


transformed data { 
    
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        int n_sets_of_cutpoints = 1;
        
}
  

parameters {
  
        vector[2] beta_mu;    
        vector<lower=0.0>[2] beta_SD;   
        matrix[2, n_studies] beta_z;    // Study-specific random effects for beta (off-centered parameterisation)
        cholesky_factor_corr[2] beta_L_Omega;
        ////
        ordered[n_thr] C_MU;
        vector<lower=0.0>[n_thr] C_SD;
        matrix[n_studies, n_thr - 1] raw_inc;
        ////
        vector[n_studies] first_C;
        
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
        array[n_sets_of_cutpoints] vector[n_thr] Ind_Dir_cumul;
        array[n_sets_of_cutpoints] vector[n_thr] Ind_Dir_cumul_prob;
        array[n_sets_of_cutpoints] vector[n_cat] Ind_Dir_ord_prob;
        ////
        array[n_sets_of_cutpoints] matrix[n_studies, n_thr] C;
        // matrix[n_sets_of_cutpoints, n_thr] C_MED;
        real Jacobian_for_unc_C_to_trans_C_MU = 0.0;
        real Jacobian_for_unc_C_to_trans_C = 0.0;
        ////
        vector<lower=0.0>[n_thr - 1] gap;
        matrix<lower=0.0>[n_studies, n_thr - 1] inc;
        //// 
        //// study-specific cutpoints (construct study-specific trans_C and then study-specific cutpoints)
        //// 
        C[1][, 1] = first_C;
        for (k in 2:n_thr) {
                gap[k - 1] = C_MU[k] - C_MU[k - 1];
                for (s in 1:n_studies) {
                    //// Ensure the study-specific cutpoint maintains ordering:
                    inc[s, k - 1] =  log1p_exp(raw_inc[s, k - 1]);
                    C[1][s, k] = C[1][s, k - 1] + inc[s, k - 1];
                    //// Jacobian for   raw_inc -> inc transformation:
                    Jacobian_for_unc_C_to_trans_C_MU += log_inv_logit(raw_inc[s, k - 1]);
                    //// Jacobian for   inc -> C transformation:
                    Jacobian_for_unc_C_to_trans_C_MU += 0.0; //// No Jacobian as inc -> C is linear transformation. 
                }
        }
        // // Subsequent cutpoints:
        // {
        //       //// ---- For constructing each study-specific trans_C just in terms of unc_C_normal[s, k]:
        //       ////
        //       trans_C[, 2:n_thr] = log1p_exp(unc_C_normal[, 2:n_thr]);
        //       //// deriv of trans_C w.r.t unc_C_normal:
        //       Jacobian_for_unc_C_to_trans_C +=  sum(log_inv_logit(unc_C_normal[, 2:n_thr]));
        //       // // //// deriv of unc_C_normal w.r.t unc_C_normal_MU, unc_C_normal_SD and unc_C_normal_z:
        //       // // Jacobian_for_unc_C_to_trans_C += 0.0;
        //       // // Jacobian_for_unc_C_to_trans_C += log(sum(abs(unc_C_normal_z[, 2:n_thr])));
        //       // // Jacobian_for_unc_C_to_trans_C += log(sum(abs(unc_C_normal_SD[2:n_thr])));
        //       ////
        //       // //// ---- For constructing each study-specific trans_C in terms of unc_C_normal[s, k] AND ALSO in terms of unc_C_normal_SD[k]:
        //       // matrix[n_studies, n_thr] unc_C_normal_SD_sq_mat;
        //       // for (s in 1:n_studies) {
        //       //     unc_C_normal_SD_sq_mat[s, ] = to_row_vector(unc_C_normal_SD_sq);
        //       // }
        //       // trans_C[, 2:n_thr]  = log1p_exp(unc_C_normal[, 2:n_thr] + 0.5 * unc_C_normal_SD_sq_mat[, 2:n_thr]);
        //       // //// deriv of trans_C w.r.t unc_C_normal:
        //       // Jacobian_for_unc_C_to_trans_C += sum(log_inv_logit(unc_C_normal[, 2:n_thr]));
        //       // // //// deriv of unc_C_normal w.r.t unc_C_normal_MU, unc_C_normal_SD and unc_C_normal_z:
        //       // // Jacobian_for_unc_C_to_trans_C += 0.0;
        //       // // Jacobian_for_unc_C_to_trans_C += log(sum(abs(unc_C_normal_z[, 2:n_thr])));
        //       // // Jacobian_for_unc_C_to_trans_C += log(sum(abs(unc_C_normal_SD[2:n_thr])));
        //       // //// deriv of trans_C w.r.t unc_C_normal_SD:
        //       // Jacobian_for_unc_C_to_trans_C += n_studies*log(sum(abs(unc_C_normal_SD[2:n_thr])));
        // }
        // ////
        // //// construct positive-constrained trans_C_MU (SUMMARY) params:
        // ////
        // trans_C_MED[1] = unc_C_normal_MU[1]; // first cutpoint is unconstrained
        // trans_C_MU[1]  = unc_C_normal_MU[1]; // first cutpoint is unconstrained
        // //// 
        // for (k in 2:n_thr) {
        //     // real stuff_to_softplus =   unc_C_normal_MU[k] + 0.5 * unc_C_normal_SD_sq[k];
        //     trans_C_MU[k]  = log1p_exp(unc_C_normal_MU[k] + 0.5 * unc_C_normal_SD_sq[k]);
        //     real stuff_to_softplus = unc_C_normal_MU[k];
        //     trans_C_MED[k] =  log1p_exp(unc_C_normal_MU[k]);
        //     //// deriv of trans_C_MU w.r.t unc_C_normal_MU:
        //     //// deriv of trans_C_MU[k]  = exp(stuff_to_softplus) w.r.t  unc_C_normal_MU is:
        //     //// exp(stuff_to_softplus) * deriv_of_correction_wrt_unc_C_normal_MU = exp(stuff_to_softplus) * 1;
        //     Jacobian_for_unc_C_to_trans_C_MU += log_inv_logit(stuff_to_softplus);  //// deriv of trans_C_MU w.r.t unc_C_normal_MU
        //     // //// deriv of trans_C_MU w.r.t unc_C_normal_SD:
        //     // real log_abs_deriv_of_correction_wrt_SD =  log(abs(unc_C_normal_SD[k]));
        //     // Jacobian_for_unc_C_to_trans_C_MU += log_inv_logit(stuff_to_softplus) + log_abs_deriv_of_correction_wrt_SD;  
        //     // //// Jacobian_for_MU_to_MED += log_inv_logit(unc_C_normal_MU[k]); //// since derivative of log1p_exp is inv_logit!
        // }
        // ////
        // //// Summary cutpoints:
        // ////
        // // First cutpoint is directly parameterized:
        // C_MED[1, 1] = trans_C_MED[1];
        // C_MU[1, 1]  = trans_C_MU[1];
        // for (k in 2:n_thr) {
        //         C_MED[1, k] = C_MED[1, k - 1] + trans_C_MED[k]; // No Jacobian needed for transformation trans_C_MU -> C_MED
        //         C_MU[1, k]  = C_MU[1, k - 1]  + trans_C_MU[k]; // No Jacobian needed for transformation trans_C_MU -> C_MED
        // }
        //// Initialise matrices:
        ////
        for (c in 1:2) {
                 logit_cumul_prob[c]   = rep_matrix(positive_infinity(), n_studies, n_thr);
                 cumul_prob[c]         = rep_matrix(1.0, n_studies, n_thr);
                 cond_prob[c]          = rep_matrix(1.0, n_studies, n_thr);
                 log_lik[c]            = rep_matrix(0.0, n_studies, n_thr);
        }
        for (c in 1:n_sets_of_cutpoints) {
                 Ind_Dir_cumul[c]      = rep_vector(0.0, n_thr);
                 Ind_Dir_cumul_prob[c] = rep_vector(0.0, n_thr);
                 Ind_Dir_ord_prob[c]   = rep_vector(0.0, n_cat);
        }
        ////
        //// Between-study model for the location parameters ("beta") - models between-study correlation:
        ////
        beta_L_Sigma = diag_pre_multiply(beta_SD, beta_L_Omega);
        for (s in 1:n_studies) {
             beta[, s] = beta_mu + beta_L_Sigma * beta_z[, s];
        }
        ////
        //// Induced-Dirichlet model cumulative probs:
        ////
        for (c in 1:n_sets_of_cutpoints) {
                for (k in 1:n_thr) {
                         Ind_Dir_cumul[c][k] = (C_MU[k] - 0.0);
                }
                Ind_Dir_cumul_prob[c] = Phi(Ind_Dir_cumul[c]); /// NaN here ???
                {
                    //// Induced-Dirichlet ordinal probs:
                    Ind_Dir_ord_prob[c][1] = Ind_Dir_cumul_prob[c][1] - 0.0;
                    for (k in 2:n_thr) {
                       Ind_Dir_ord_prob[c][k] = Ind_Dir_cumul_prob[c][k] - Ind_Dir_cumul_prob[c][k - 1]; // since cutpoints are increasing with k
                    }
                    Ind_Dir_ord_prob[c][n_cat] = 1.0 - Ind_Dir_cumul_prob[c][n_cat - 1];
                }
        }
        ////
        //// Likelihood using binomial factorization:
        ////
        for (s in 1:n_studies) {
                  for (c in 1:2) {
                      for (cut_i in 1:to_int(n_cutpoints[c, s])) {
                             int k = to_int(cutpoint_index[c][s, cut_i]);
                             logit_cumul_prob[c][s, cut_i] = (C[1][s, k] - beta[c, s]);
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
      
        //// location priors:
        {
            beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
            beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
            beta_L_Omega ~ lkj_corr_cholesky(2.0);
        }
        //// Cutpoint priors:
        C_MU[1] ~ normal(0, 5.0);
        C_SD    ~ normal(0, 1.0);
        // unc_C_normal_SD[2:n_thr] ~ normal(0, 2.5); //// need this prior if parameterising summary cutpoints in terms of MEDIANS of unc_C (NOT means). 
        ////
        ////
        for (s in 1:n_studies) {
            first_C[s] ~ normal(C_MU[1], C_SD[1]);
        }
        //// Rest of cutpoints:
        for (k in 2:n_thr) {
              for (s in 1:n_studies) {
                  raw_inc[s, k - 1] ~ normal(gap[k - 1], C_SD[k - 1]);
              }
        }
        ////
        //// Jacobian adjustment for transforming from unc_C -> C
        ////
        target += Jacobian_for_unc_C_to_trans_C;
        target += Jacobian_for_unc_C_to_trans_C_MU;
        ////
        //// Cutpoint between study model:
        ////
        for (c in 1:n_sets_of_cutpoints) {
            vector[n_thr] rho =  normal_pdf(to_vector(Ind_Dir_cumul[c][1:n_thr]), 0.0, 1.0);   //  p_cumul[k - 1] * (1.0 - p_cumul[k - 1]); // original
            target += induced_dirichlet_v2_lpdf(to_vector(Ind_Dir_ord_prob[c][1:n_cat]) | rho, prior_alpha);
        }
        ////
        //// Likelihood / Model:
        {
            to_vector(beta_z) ~ std_normal();   // (part of between-study model, NOT prior)
            // to_vector(unc_C_normal_z) ~ std_normal();
        }
        //// Increment the log-likelihood:
        if (prior_only == 0) {
            for (c in 1:2) {
              target +=  sum(log_lik[c]);
            }
        }
  
}




 

generated quantities {

          vector[n_thr] Fp;
          vector[n_thr] Sp;
          vector[n_thr] Se;
          ////
          vector[n_thr] Fp_MED;
          vector[n_thr] Sp_MED;
          vector[n_thr] Se_MED;
          ////
          vector[n_thr] Fp_MU;
          vector[n_thr] Sp_MU;
          vector[n_thr] Se_MU;
          ////
          vector[n_thr] Fp_EMP;
          vector[n_thr] Sp_EMP;
          vector[n_thr] Se_EMP;
          ////
          vector[n_thr] Fp_SIM_MED;
          vector[n_thr] Sp_SIM_MED;
          vector[n_thr] Se_SIM_MED;
          ////
          vector[n_thr] Fp_SIM_MU;
          vector[n_thr] Sp_SIM_MU;
          vector[n_thr] Se_SIM_MU;
          ////
          vector[n_thr] Se_pred;
          vector[n_thr] Sp_pred;
          vector[n_thr] Fp_pred;
          ////
          matrix[n_studies, n_thr] se;
          matrix[n_studies, n_thr] sp;
          matrix[n_studies, n_thr] fp;
          matrix[n_studies, n_thr] x_hat_nd = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] x_hat_d  = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_nd   = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_d    = rep_matrix(-1, n_studies, n_thr);
          array[2] matrix[n_studies, n_thr] x_hat;
          array[2] matrix[n_studies, n_thr] dev;
          matrix[n_sets_of_cutpoints, n_thr] C_mu_empirical; 
          corr_matrix[2] beta_Omega;
          ////
          //// Compute between-study correlation matrix for location parameters:
          beta_Omega  =  multiply_lower_tri_self_transpose(beta_L_Omega);
          ////
          matrix[n_sets_of_cutpoints, n_thr] C_SIM_medians;  
          matrix[n_sets_of_cutpoints, n_thr] C_SIM_means;
          matrix[n_sets_of_cutpoints, n_thr] trans_C_medians; 
          matrix[n_sets_of_cutpoints, n_thr] trans_C_means; 
                     
          // //// Summary cutpoints:
          // {
          //       int n_sims = 1000;
          //       matrix[n_sims, n_thr] C_sim;
          //       matrix[n_sims, n_thr] trans_C_sim;
          //       for (i in 1:n_sims) {
          //                C_sim[i, 1] = normal_rng(unc_C_normal_MU[1], unc_C_normal_SD[1]);
          //                trans_C_sim[i, 1] = C_sim[i, 1];
          //                // Subsequent cutpoints:
          //                for (k in 2:n_thr) {
          //                       real unc_C_normal_sim = normal_rng(unc_C_normal_MU[k], unc_C_normal_SD[k]);
          //                       trans_C_sim[i, k] = log1p_exp(unc_C_normal_sim);
          //                       C_sim[i, k] = C_sim[i, k - 1] + trans_C_sim[i, k];
          //                }
          //       }
          // 
          //       // Compute summary cutpoints as medians:
          //       for (k in 1:n_thr) {
          //           trans_C_medians[1, k] = median(trans_C_sim[, k]);
          //           trans_C_means[1, k]   = mean(trans_C_sim[, k]);
          //           C_SIM_medians[1, k] = median(C_sim[, k]);
          //           C_SIM_means[1, k]   = mean(C_sim[, k]);
          //       }
          // }
          //// Initialise containers:
          for (c in 1:2) {
              x_hat[c]  = rep_matrix(-1, n_studies, n_thr);
              dev[c]    = rep_matrix(-1, n_studies, n_thr); 
          }
          ////
          //// Empirical-mean cutpoints:
          for (c in 1:n_sets_of_cutpoints) {
              for (k in 1:n_thr) {
                    C_mu_empirical[c, k] = median(C[c][, k]);
              }
          }
          ////
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
                // Fp_MED[k] =   1.0 - Phi(C_MED[1, k] - beta_mu[1]);
                // Sp_MED[k] =   1.0 - Fp_MED[k];
                // Se_MED[k] =   1.0 - Phi(C_MED[1, k] - beta_mu[2]);
                ////
                Fp_MU[k] =   1.0 - Phi(C_MU[k] - beta_mu[1]);
                Sp_MU[k] =   1.0 - Fp_MU[k];
                Se_MU[k] =   1.0 - Phi(C_MU[k] - beta_mu[2]);
                ////
                Fp_EMP[k] =   1.0 - Phi(C_mu_empirical[1, k] - beta_mu[1]);
                Sp_EMP[k] =   1.0 - Fp_EMP[k];
                Se_EMP[k] =   1.0 - Phi(C_mu_empirical[1, k] - beta_mu[2]);
                // ////
                // Fp_SIM_MED[k] =   1.0 - Phi(C_SIM_medians[1, k] - beta_mu[1]);
                // Sp_SIM_MED[k] =   1.0 - Fp_SIM_MED[k];
                // Se_SIM_MED[k] =   1.0 - Phi(C_SIM_medians[1, k] - beta_mu[2]);
                // ////
                // Fp_SIM_MU[k] =   1.0 - Phi(C_SIM_means[1, k] - beta_mu[1]);
                // Sp_SIM_MU[k] =   1.0 - Fp_SIM_MU[k];
                // Se_SIM_MU[k] =   1.0 - Phi(C_SIM_means[1, k] - beta_mu[2]);
                ////
                Fp[k] = Fp_MED[k];
                Sp[k] = Sp_MED[k];
                Se[k] = Se_MED[k];
          }
          //// Calculate predictive accuracy:
          {
                vector[2] beta_pred =  to_vector(multi_normal_cholesky_rng(beta_mu, beta_L_Sigma));
                for (k in 1:n_thr) {
                      Fp_pred[k] =   1.0 - Phi(C_MU[k] - beta_pred[1]);
                      Sp_pred[k] =   1.0 - Fp_pred[k];
                      Se_pred[k] =   1.0 - Phi(C_MU[k] - beta_pred[2]);
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






















