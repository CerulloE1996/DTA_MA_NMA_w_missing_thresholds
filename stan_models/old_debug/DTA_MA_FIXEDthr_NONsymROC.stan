
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
                            // real rho =  p_cumul[k - 1] * (1.0 - p_cumul[k - 1]); // original
                            // J[k, k] = -rho;     // original
                            // J[k - 1, k] = +rho; // original
                            real rho =  p_cumul[k - 1] * (1.0 - p_cumul[k - 1]); // original
                            J[k, k] = +rho;
                            J[k - 1, k] = -rho;
                    }
                    
                    return dirichlet_lpdf(p_ord | alpha) + log_determinant(J);
          
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
        //// Other:
        int prior_only;
        int estimate_scales;
    
  
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
        real log_scale_d_mu;    
        real<lower=0.0> log_scale_d_SD;   
        vector[n_studies] log_scale_d_z;    // Study-specific random effects for beta (off-centered parameterisation)
        ////
        ordered[n_thr] C_nd;  // Global cutpoints
        ordered[n_thr] C_d;   // Global cutpoints
      
}


transformed parameters { 
    
        matrix[2, n_studies] beta;
        vector[2] log_scale_mu;    
        vector<lower=0.0>[2] log_scale_SD;   
        matrix[2, n_studies] log_scale_z;    // Study-specific random effects for beta (off-centered parameterisation)
        matrix[2, n_studies] log_scale;
        matrix[2, n_studies] scale;
        ////
        array[2] matrix[n_studies, n_thr] logit_cumul_prob; // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cumul_prob; // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cond_prob;  // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] log_lik;  // log_lik storage
        ////
        matrix[2, n_thr] Ind_Dir_anchor;
        matrix[2, n_thr] Ind_Dir_cumul_prob;
        matrix[2, n_cat] Ind_Dir_ord_prob;
        ////
        matrix[2, n_thr] C;
        
        //// Initialise matrices:
        for (c in 1:2) {
                   logit_cumul_prob[c] = rep_matrix(1.0, n_studies, n_thr);
                   cumul_prob[c]       = rep_matrix(1.0, n_studies, n_thr);
                   cond_prob[c]        = rep_matrix(1.0, n_studies, n_thr);
                   log_lik[c]          = rep_matrix(0.0, n_studies, n_thr);
        }
        //// Cutpoints:
        C[1, ] = to_row_vector(C_nd);
        C[2, ] = to_row_vector(C_d);
        //// Induced-Dirichlet ** prior model ** stuff:
        for (c in 1:2) {
                
                   //// Ind-Dir cumulative probs:
                   Ind_Dir_anchor[c, ]       = rep_row_vector(0.0, n_thr);
                   Ind_Dir_cumul_prob[c, ] = inv_logit(Ind_Dir_anchor[c, ] - C[c, ]);
                   //// Ind-Dir ordinal probs:
                   Ind_Dir_ord_prob[c, 1] = 1.0 - Ind_Dir_cumul_prob[c, 1];
                   for (k in 2:n_thr) {
                       Ind_Dir_ord_prob[c, k] = Ind_Dir_cumul_prob[c, k - 1] - Ind_Dir_cumul_prob[c, k]; // since probs are decreasing with k
                   }
                   Ind_Dir_ord_prob[c, n_cat] =  Ind_Dir_cumul_prob[c, n_cat - 1];
            
        }
         
        //// Between-study model:
        log_scale_mu[1] = 0.0;
        log_scale_SD[1] = 0.0;
        log_scale_z[1, ]  = rep_row_vector(0.0, n_studies);
        
        if (estimate_scales == 1) { 
              log_scale_mu[2]   = log_scale_d_mu;
              log_scale_SD[2]   = log_scale_d_SD;
              log_scale_z[2, ]  = to_row_vector(log_scale_d_z);
        } else { 
              log_scale_mu[2] = 0.0;
              log_scale_SD[2] = 0.0;
              log_scale_z[2, ]  = rep_row_vector(0.0, n_studies);
        }
        
        //// Between-study model for the location parameters ("beta") - models between-study correlation:
        for (s in 1:n_studies) {
             beta[, s] = beta_mu + diag_pre_multiply(beta_SD, beta_L_Omega) * beta_z[, s];
        }
        
        for (c in 1:2) {
            log_scale[c, ] = to_row_vector(log_scale_mu[c] + log_scale_z[c, ] * log_scale_SD[c]);
            scale[c, ] = exp(log_scale[c, ]);
        }
        //// Likelihood using binomial factorization:
        for (s in 1:n_studies) {
                for (c in 1:2) {
                      for (cut_i in 1:to_int(n_cutpoints[c, s])) {
                              int k = to_int(cutpoint_index[c][s, cut_i]);
                              logit_cumul_prob[c][s, cut_i] = ( beta[c, s] - C[c, k] )/scale[c, s];
                      } 
                }
        }
        //// Calculate CUMULATIVE probabilities (vectorised):
        for (c in 1:2) {
            cumul_prob[c] = inv_logit(logit_cumul_prob[c]);
        }
        //// ------- Binomial likelihood:
        for (s in 1:n_studies) {
                    for (c in 1:2) {
                            real previous_cumul_prob = 1.0;
                            for (cut_i in 1:to_int(n_cutpoints[c, s])) {
                              
                                        cond_prob[c][s, cut_i] = cumul_prob[c][s, cut_i] / previous_cumul_prob;
                                        if (cut_i == 1) { //// First non-missing cutpoint
                                                  int success_prob = to_int(n[c][s, 1]);
                                                  int x_i = to_int(x[c][s, 1]);
                                                  log_lik[c][s, cut_i] = binomial_lpmf( x_i | success_prob, cond_prob[c][s, 1] );
                                        } else if (cut_i > 1) { 
                                                  int success_prob = to_int(n[c][s, cut_i]);
                                                  int x_i = to_int(x[c][s, cut_i]);
                                                  log_lik[c][s, cut_i] = binomial_lpmf( x_i | success_prob,  cond_prob[c][s, cut_i] );
                                        }
                                        previous_cumul_prob = cond_prob[c][s, cut_i];
                          
                            }
                    }
         }
      
}


model {
      
        //// Priors:
        beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
        beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
        beta_L_Omega ~ lkj_corr_cholesky(2.0);
        log_scale_d_mu ~ normal(prior_log_scale_mu_mean[2], prior_log_scale_mu_SD[2]);
        log_scale_d_SD ~ normal(prior_log_scale_SD_mean[2], prior_log_scale_SD_SD[2]);
        
        //// Likelihood / Model:
        to_vector(log_scale_d_z) ~ std_normal(); // part of between-study model, NOT prior
        to_vector(beta_z) ~ std_normal();        // part of between-study model, NOT prior
     
        //// Induced-dirichlet ** Prior ** model:
        for (c in 1:2) {
               target += induced_dirichlet_v2_lpdf( to_vector(Ind_Dir_ord_prob[c, 1:n_cat]) | to_vector(Ind_Dir_cumul_prob[c, 1:n_thr]), prior_alpha);
        }
        
        //// Jacobian adjustments needed for scales (if estimating them):
        if (estimate_scales == 1) { 
            target += sum(log_scale);                // double-checked the log-derivative of this by hand (correct)
            if (abs(log_scale_d_SD) != 0.0)     target += log(abs(log_scale_d_SD));      // double-checked the log-derivative of this by hand (correct)
            if (abs(sum(log_scale_d_z)) != 0.0) target += log(abs(sum(log_scale_d_z)));  // double-checked the log-derivative of this by hand (correct)
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
          corr_matrix[2] beta_Omega;
          vector[2] scale_mu = exp(log_scale_mu);
          
          //// Compute between-study correlation matrix:
          beta_Omega = multiply_lower_tri_self_transpose(beta_L_Omega);

          //// Initialise containers:
          for (c in 1:2) {
              x_hat[c]  = rep_matrix(-1, n_studies, n_thr);
              dev[c]    = rep_matrix(-1, n_studies, n_thr);
          }
          //// Calculate study-specific accuracy:
          for (s in 1:n_studies) {
                for (k in 1:n_thr) {
                    fp[s, k] =   cumul_prob[1][s, k];
                    sp[s, k] =   1.0 - fp[s, k];
                    se[s, k] =   cumul_prob[2][s, k];
                }
          }
          //// Calculate summary accuracy (using mean parameters):
          for (k in 1:n_thr) {
                Fp[k] =   inv_logit((beta_mu[1] - C[1, k])/scale_mu[1]);
                Sp[k] =   1.0 - Fp[k];
                Se[k] =   inv_logit((beta_mu[2] - C[2, k])/scale_mu[2]);
          }
          //// Calculate predictive accuracy:
          {
                vector[2] beta_pred =  to_vector(normal_rng(beta_mu[1:2], beta_SD[1:2]));
                vector[2] log_cale_pred;
                log_cale_pred[1]  = 0.0;
                if (estimate_scales == 1) { 
                   log_cale_pred[2]  = normal_rng(log_scale_mu[2], log_scale_SD[2]);
                } else { 
                   log_cale_pred[2]  = 0.0;
                }
    
                vector[2] scale_pred = exp(log_cale_pred);
       
                for (k in 1:n_thr) {
                      Fp_pred[k] =   inv_logit((beta_pred[1] - C[1, k])/scale_pred[1]);
                      Sp_pred[k] =   1.0 - Fp_pred[k];
                      Se_pred[k] =   inv_logit((beta_pred[2] - C[2, k])/scale_pred[2]);
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






















