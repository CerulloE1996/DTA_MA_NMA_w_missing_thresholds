
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
        
}


        

data {
    
        int<lower=1> n_studies;
        int<lower=1> n_thr;
        array[2, n_studies] int n_cutpoints;
        array[2] matrix[n_studies, n_thr] x_with_missings;
        array[2] matrix[n_studies, n_thr] n;
        array[2] matrix[n_studies, n_thr] x;
        array[2] matrix[n_studies, n_thr] cutpoint_index;
        //// ////
        //// array[2] matrix[n_studies, n_thr + 1] x_cat;
        //// Priors:
        vector<lower=0>[n_thr + 1] prior_alpha; // Induced-Dirichlet prior vector 
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
        int estimate_scales;
    
  
}

transformed data { 
    
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        // int n_sets_of_cutpoints = 1;
}
  


parameters {
  
        vector[2] beta_mu;    
        vector<lower=0.0>[2] beta_SD;   
        matrix[2, n_studies] beta_z;    // Study-specific random effects for beta (off-centered parameterisation)
        cholesky_factor_corr[2] beta_L_Omega;
        ////
        real raw_scale_d_mu;    
        real<lower=0.0> raw_scale_d_SD;   
        vector[n_studies] raw_scale_d_z;
        // ////
        ordered[n_thr] C;  // Global cutpoints
      
}


transformed parameters { 
    
        matrix[2, n_studies] beta;
        cholesky_factor_cov[2] beta_L_Sigma;
        ////
        vector[2] raw_scale_mu;    
        vector<lower=0.0>[2] raw_scale_SD;   
        matrix[2, n_studies] raw_scale_z;
        matrix[2, n_studies] raw_scale;
        matrix[2, n_studies] scale;
        ////
        array[2] matrix[n_studies, n_thr] logit_cumul_prob; // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cumul_prob; // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cond_prob;  // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] log_lik;  // log_lik storage
        ////
        vector[n_thr] Ind_Dir_anchor;
        vector[n_thr] Ind_Dir_cumul;
        vector[n_thr] Ind_Dir_cumul_prob;
        vector[n_cat] Ind_Dir_ord_prob;
        ////
        real Jacobian_raw_scale_to_scale = 0.0;
        ////
        //// Initialise matrices:
        ////
        for (c in 1:2) {
                   logit_cumul_prob[c] = rep_matrix(positive_infinity(), n_studies, n_thr);
                   cumul_prob[c]       = rep_matrix(1.0, n_studies, n_thr);
                   cond_prob[c]        = rep_matrix(0.0, n_studies, n_thr);
                   log_lik[c]          = rep_matrix(0.0, n_studies, n_thr);
        }
        ////
        //// Between-study model for the location parameters ("beta") - models between-study correlation:
        ////
        beta_L_Sigma = diag_pre_multiply(beta_SD, beta_L_Omega);
        for (s in 1:n_studies) {
             beta[, s] = beta_mu + beta_L_Sigma * beta_z[, s];
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
        //// Between-study model for raw_sscale parameters:
        ////
        raw_scale_mu[1] = 0.0;
        raw_scale_SD[1] = 0.0;
        raw_scale_z[1, ]  = rep_row_vector(0.0, n_studies);
        ////
        if (estimate_scales == 1) { 
              raw_scale_mu[2]   = raw_scale_d_mu;
              raw_scale_SD[2]   = raw_scale_d_SD;
              raw_scale_z[2, ]  = to_row_vector(raw_scale_d_z);
        } else { 
              raw_scale_mu[2] = 0.0;
              raw_scale_SD[2] = 0.0;
              raw_scale_z[2, ]  = rep_row_vector(0.0, n_studies);
        }
        for (c in 1:2) {
            raw_scale[c, ] = to_row_vector(raw_scale_mu[c] + raw_scale_z[c, ] * raw_scale_SD[c]);
            scale[c, ] = log1p_exp(raw_scale[c, ]);
        }
        if (estimate_scales == 1) { 
            Jacobian_raw_scale_to_scale += sum(log_inv_logit(raw_scale[2, ])); // Only in diseased class!
        }
        //// Likelihood using binomial factorization:
        for (s in 1:n_studies) {
                for (c in 1:2) {
                      for (cut_i in 1:to_int(n_cutpoints[c, s])) {
                              int k = to_int(cutpoint_index[c][s, cut_i]);
                              logit_cumul_prob[c][s, cut_i] = (C[k] - beta[c, s])/scale[c, s];
                      } 
                }
        }
        //// Calculate CUMULATIVE probabilities (vectorised):
        for (c in 1:2) {
            cumul_prob[c] = Phi(logit_cumul_prob[c]);
        }
        //// ------- Binomial likelihood:
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
      
        //// Priors:
        beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
        beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
        beta_L_Omega ~ lkj_corr_cholesky(2.0);
        //// scales:
        raw_scale_d_mu ~ normal(prior_raw_scale_mu_mean[2], prior_raw_scale_mu_SD[2]);
        raw_scale_d_SD ~ normal(prior_raw_scale_SD_mean[2], prior_raw_scale_SD_SD[2]);
        
        //// Likelihood / Model:
        to_vector(raw_scale_d_z) ~ std_normal(); // part of between-study model, NOT prior
        to_vector(beta_z) ~ std_normal();        // part of between-study model, NOT prior
     
        //// Induced-dirichlet ** Prior ** model:
        {
               vector[n_thr] rho =  normal_pdf(Ind_Dir_cumul[1:n_thr], 0.0, 1.0);   //  p_cumul[k - 1] * (1.0 - p_cumul[k - 1]); // original
               target += induced_dirichlet_v2_lpdf( Ind_Dir_ord_prob[1:n_cat] | rho, prior_alpha);
        }
        
        //// Jacobian adjustments needed for scales (if estimating them):
        if (estimate_scales == 1) { 
            target +=  Jacobian_raw_scale_to_scale;
            //// Also need the Jacobian to go from (raw_scale_z, raw_scale_mu, raw_scale_SD) -> raw_scale! (using chain rule):
            if (abs(raw_scale_d_SD) != 0.0)     target += log(abs(raw_scale_d_SD));      // double-checked the log-derivative of this by hand (correct)
            if (abs(sum(raw_scale_d_z)) != 0.0) target += log(abs(sum(raw_scale_d_z)));  // double-checked the log-derivative of this by hand (correct)
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
          vector[2] scale_mu = log1p_exp(raw_scale_mu);
          
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
                    fp[s, k] =   1.0 - cumul_prob[1][s, k];
                    sp[s, k] =   1.0 - fp[s, k];
                    se[s, k] =   1.0 - cumul_prob[2][s, k];
                }
          }
          //// Calculate summary accuracy (using mean parameters):
          for (k in 1:n_thr) {
                Fp[k] =   1.0 - Phi((beta_mu[1] - C[k])/scale_mu[1]);
                Sp[k] =   1.0 - Fp[k];
                Se[k] =   1.0 - Phi((beta_mu[2] - C[k])/scale_mu[2]);
          }
          //// Calculate predictive accuracy:
          {
                vector[2] beta_pred =  to_vector(multi_normal_cholesky_rng(beta_mu, beta_L_Sigma));
                vector[2] log_cale_pred;
                log_cale_pred[1]  = 0.0;
                if (estimate_scales == 1) { 
                   log_cale_pred[2]  = normal_rng(raw_scale_mu[2], raw_scale_SD[2]);
                } else { 
                   log_cale_pred[2]  = 0.0;
                }
    
                vector[2] scale_pred = exp(log_cale_pred);
       
                for (k in 1:n_thr) {
                      Fp_pred[k] =   1.0 - Phi((beta_pred[1] - C[k])/scale_pred[1]);
                      Sp_pred[k] =   1.0 - Fp_pred[k];
                      Se_pred[k] =   1.0 - Phi((beta_pred[2] - C[k])/scale_pred[2]);
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






















