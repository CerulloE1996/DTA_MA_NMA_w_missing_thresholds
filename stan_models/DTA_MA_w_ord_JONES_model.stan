
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
  
        
      //// Box-Cox transform function:
      real box_cox( real x, 
                    real lambda) {
            
            if (lambda == 0) {
                return log(x);
            } else {
                return (pow(x, lambda) - 1) / lambda;
            }
          
      }
      
      //// Vectorized version:
      vector box_cox( vector x, 
                      real lambda) {
                               
            int N = num_elements(x);
            vector[N] result;
            
            for (n in 1:N) {
              
                if (lambda == 0) {
                     result[n] = log(x[n]); 
                } else {
                     result[n] = (pow(x[n], lambda) - 1) / lambda;
                }
                
            }
            
            return result;
          
      }
      
      // //// Log probability for Box-Cox normal model (useful for estimating optimal lambda):
      // real box_cox_normal_lpdf( vector y, 
      //                           real mu, 
      //                           real sigma, 
      //                           real lambda) {
      //                             
      //       int N = num_elements(y);
      //       vector[N] transformed_y = box_cox(y, lambda);
      //       
      //       // Add Jacobian adjustment for transformation:
      //       real log_jacobian = (lambda - 1) * sum(log(y));
      //       
      //       // Return log probability:
      //       return normal_lpdf(transformed_y | mu, sigma) + log_jacobian;
      //   
      // }
        
  
}


data {
  
      int<lower=1> n_studies;
      int<lower=1> n_thr;
      vector<lower=0>[n_thr + 1] alpha_non_diseased; // Induced-Dirichlet prior vector 
      vector<lower=0>[n_thr + 1] alpha_diseased; // Induced-Dirichlet prior vector 
      array[n_studies] int<lower=1> n_non_diseased;  // number in non-diseased group in each study
      array[n_studies] int<lower=1> n_diseased;  // number in diseased group in each study
      array[n_studies, n_thr] int x_non_diseased;  // counts for non-diseased ("missing" threshold data recorded as 999)
      array[n_studies, n_thr] int x_diseased;  // counts for diseased     ("missing" threshold data recorded as 999)
      vector[n_thr] cts_thr_values_nd; // = logit(cts_thr_values_nd);   // Global cutpoints 
      vector[n_thr] cts_thr_values_d; //  = logit(cts_thr_values_d);    // Global cutpoints
  
}

transformed data { 
  
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        matrix[n_studies, n_cat] x_inc_n_non_diseased; 
        matrix[n_studies, n_cat] x_inc_n_diseased; 
        
        for (s in 1:n_studies) {
        
              x_inc_n_non_diseased[s, 1] = n_non_diseased[s];
              x_inc_n_diseased[s, 1]     = n_diseased[s];
              
               for (k in 2:n_cat) {
                  x_inc_n_non_diseased[s, k] = x_non_diseased[s, k - 1];
                  x_inc_n_diseased[s, k]     = x_diseased[s, k - 1];
               }
             
        }
  
}


parameters {
  
      real beta_d_mu;
      real beta_nd_mu; 
      real log_scale_d_mu;
      real log_scale_nd_mu;
      real<lower=0> beta_d_SD;  
      real<lower=0> beta_nd_SD;   
      real<lower=0> log_scale_d_SD;
      real<lower=0> log_scale_nd_SD;
      vector[n_studies] beta_d_z;    // Study-specific random effects for beta (off-centered parameterisation)
      vector[n_studies] beta_nd_z;    // Study-specific random effects for beta (off-centered parameterisation)
      vector[n_studies] log_scale_d_z;    // Study-specific random effects for scale (off-centered parameterisation)
      vector[n_studies] log_scale_nd_z;    // Study-specific random effects for scale (off-centered parameterisation)
      //// box-cox params:
      real<lower=-3.0,upper=3.0> lambda_nd;
      real<lower=-3.0,upper=3.0> lambda_d;
  
}


transformed parameters { 
  
      //// Construct the study-specific random effects (off-centered param.):
      vector[n_studies] beta_nd =     beta_nd_mu      + beta_nd_z      .* beta_nd_SD; 
      vector[n_studies] beta_d =      beta_d_mu       + beta_d_z       .* beta_d_SD;  
      vector[n_studies] log_scale_nd = log_scale_nd_mu + log_scale_nd_z .* log_scale_nd_SD;  
      vector[n_studies] log_scale_d  = log_scale_d_mu  + log_scale_d_z  .* log_scale_d_SD;  
      vector[n_studies] scale_d  = exp(log_scale_d);
      vector[n_studies] scale_nd = exp(log_scale_nd);
      //// Cutpoints:
      vector[n_thr] cutpoints_nd = box_cox(cts_thr_values_nd, lambda_nd);   // Global cutpoints 
      vector[n_thr] cutpoints_d  = box_cox(cts_thr_values_d,  lambda_d);    // Global cutpoints
      //// Cumulative ordinal probs for the likelihood:
      matrix<lower=0.0>[n_studies, n_thr] cumul_cond_prob_nd = rep_matrix(0.0, n_studies, n_thr);
      matrix<lower=0.0>[n_studies, n_thr] cumul_cond_prob_d  = rep_matrix(0.0, n_studies, n_thr);
      
      
              //// Likelihood using binomial factorization:
        for (s in 1:n_studies) {
                
                //// Non-diseased group (fixed parameters)
                if (x_non_diseased[s, 1] != 999) { //// If non-missing (i.e. only if study reports @ this ordinal cutpoint!)
                         cumul_cond_prob_nd[s, 1] = inv_logit((beta_nd[s] - cutpoints_nd[1])/scale_nd[s]);
                }
                
                for (k in 2:n_thr) {
                     if (x_non_diseased[s, k] != 999) { //// If non-missing (i.e. only if study reports @ this ordinal cutpoint!)
                         cumul_cond_prob_nd[s, k] = inv_logit((beta_nd[s] - cutpoints_nd[k])/scale_nd[s]) / inv_logit((beta_nd[s] - cutpoints_nd[k - 1])/scale_nd[s]);
                     }
                }
                
                //// Diseased group (D+):
                if (x_diseased[s, 1] != 999) { //// If non-missing (i.e. only if study reports @ this ordinal cutpoint!)
                         cumul_cond_prob_d[s, 1] = inv_logit((beta_d[s] - cutpoints_d[1])/scale_d[s]);
                }
                
                for (k in 2:n_thr) {
                   if (x_diseased[s, k] != 999) { //// If non-missing (i.e. only if study reports @ this ordinal cutpoint!)
                         cumul_cond_prob_d[s, k] = inv_logit((beta_d[s] - cutpoints_d[k])/scale_d[s]) /  inv_logit((beta_d[s] - cutpoints_d[k - 1])/scale_d[s]);
                   }
                }
          
        }
      
}


model {
      
      //// Between-study heterogenity priors:
      target += normal_lpdf(beta_d_mu      | 0.0, 2.00);
      target += normal_lpdf(beta_d_SD      | 0.0, 1.00);
      target += std_normal_lpdf(beta_d_z);
      target += normal_lpdf(log_scale_d_mu | 0.0, 1.00);
      target += normal_lpdf(log_scale_d_SD | 0.0, 0.50);
      target += std_normal_lpdf(log_scale_d_z);
      
      //// For non-diseased group only: 
      target += normal_lpdf( beta_nd_mu     | 0.0, 2.00);
      target += normal_lpdf( beta_nd_SD     | 0.0, 1.00);
      target += std_normal_lpdf( beta_nd_z);
      target += normal_lpdf(log_scale_nd_mu | 0.0, 1.00);
      target += normal_lpdf(log_scale_nd_SD | 0.0, 0.50);
      target += std_normal_lpdf(log_scale_nd_z);
      
      target +=  normal_lpdf(lambda_nd | 0.0, 1.0);
      target +=  normal_lpdf(lambda_d  | 0.0, 1.0);
                  
      //// Likelihood using binomial factorization:
        //// Likelihood using binomial factorization:
        for (s in 1:n_studies) {
                
                for (k in 2:n_cat) {
                  
                       //// Non-diseased group (D+):
                       int x_nd = to_int(x_inc_n_non_diseased[s, k]);
                       int x_nd_previous = to_int(x_inc_n_non_diseased[s, k - 1]);
                       if (x_nd != 999) { //// If non-missing (i.e. only if study reports @ this ordinal cutpoint!)
                             target += binomial_lpmf( x_nd | x_nd_previous, cumul_cond_prob_nd[s, k - 1] );
                       }
                       //// Diseased group (D+):
                       int x_d = to_int(x_inc_n_diseased[s, k]);
                       int x_d_previous = to_int(x_inc_n_diseased[s, k - 1]);
                       if (x_d != 999) { //// If non-missing (i.e. only if study reports @ this ordinal cutpoint!)
                             target += binomial_lpmf( x_d | x_d_previous, cumul_cond_prob_d[s, k - 1] );
                       }
                   
                }
          
        }
  
}


generated quantities {
  
        //// Summary accuracy parameters (at each threshold):
        vector[n_thr] Se;
        vector[n_thr] Sp;
        vector[n_thr] Se_pred;
        vector[n_thr] Sp_pred;
        array[n_studies] vector[n_thr] se;
        array[n_studies] vector[n_thr] sp;
        matrix[n_studies, n_thr] x_non_diseased_est =  rep_matrix(0.0, n_studies, n_thr);
        matrix[n_studies, n_thr] x_diseased_est =      rep_matrix(0.0, n_studies, n_thr);
        matrix[n_studies, n_thr] x_non_diseased_diff = rep_matrix(0.0, n_studies, n_thr);
        matrix[n_studies, n_thr] x_diseased_diff =     rep_matrix(0.0, n_studies, n_thr);
        
        //// Calculate study-specific accuracy:
        for (s in 1:n_studies) {
              for (k in 1:n_thr) {
                  se[s, k] =        inv_logit((cutpoints_d[k]      - beta_d[s]) /scale_d[s]);
                  sp[s, k] =  1.0 - inv_logit((cutpoints_nd[k]     - beta_nd[s])/scale_nd[s]);
              }
        }
        
        //// Calculate summary accuracy (using mean parameters):
        for (k in 1:n_thr) {
            Se[k] =         inv_logit((cutpoints_d[k]     - beta_d_mu) /exp(log_scale_d_mu));
            Sp[k] =  1.0 -  inv_logit((cutpoints_nd[k]    - beta_nd_mu)/exp(log_scale_d_mu));
        }
        
        //// Calculate predictive accuracy:
        {
            real beta_d_pred  = normal_rng(beta_d_mu, beta_d_SD);
            real scale_d_pred = exp(normal_rng(log_scale_d_mu, log_scale_d_SD));
            //// For non-diseased group only: 
            real beta_nd_pred  = normal_rng(beta_nd_mu, beta_nd_SD);
            real scale_nd_pred = exp(normal_rng(log_scale_nd_mu, log_scale_nd_SD));
            
            for (k in 1:n_thr) {
                Se_pred[k] =        inv_logit((cutpoints_d[k]     - beta_d_pred) /scale_d_pred);
                Sp_pred[k] =  1.0 - inv_logit((cutpoints_nd[k]    - beta_nd_pred)/scale_nd_pred);
            }
        }
        
                  //// Model-predicted ("re-constructed") data:
          for (s in 1:n_studies) {
  
                    for (k in 2:n_cat) {
  
                           //// Non-diseased group (D+):
                           int x_nd =          to_int(x_inc_n_non_diseased[s, k]);
                           int x_nd_previous = to_int(x_inc_n_non_diseased[s, k - 1]);
                           if (x_nd != 999) { //// If non-missing (i.e. only if study reports @ this ordinal cutpoint!)
                                   x_non_diseased_est[s, k - 1]  = binomial_rng( x_nd_previous, cumul_cond_prob_nd[s, k - 1] );
                                   x_non_diseased_diff[s, k - 1] = abs( x_non_diseased_est[s, k - 1] - x_nd );
                           }
                           //// Diseased group (D+):
                           int x_d =          to_int(x_inc_n_diseased[s, k]);
                           int x_d_previous = to_int(x_inc_n_diseased[s, k - 1]);
                           if (x_d != 999) { //// If non-missing (i.e. only if study reports @ this ordinal cutpoint!)
                                  x_diseased_est[s, k - 1]  = binomial_rng( x_d_previous, cumul_cond_prob_d[s, k - 1] );
                                  x_diseased_diff[s, k - 1] = abs( x_diseased_est[s, k - 1] - x_d );
                           }
  
                    }
  
        }
    
}























