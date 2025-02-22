
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
  
        real induced_dirichlet_lpdf(  vector c, 
                                      vector alpha, 
                                      real anchor) {
          
              int K = num_elements(c) + 1;
              vector[K - 1] anchoredcutoffs = c - anchor;
              vector[K] sigma;
              vector[K] p;
              matrix[K, K] J = rep_matrix(0, K, K);
              
              sigma[1:(K - 1)] = inv_logit(anchoredcutoffs);
              sigma[K] = 1;
              
              p[1] = sigma[1];
              for (k in 2:(K - 1)) {
                 p[k] = sigma[k] - sigma[k - 1];
              }
              p[K] = 1 - sigma[K - 1];
              
              // Jacobian computation
              for (k in 1:K) {
                  J[k, 1] = 1;
              }
              for (k in 2:K) {
                  real rho = sigma[k - 1] * (1.0 - sigma[k - 1]);
                  J[k, k] = -rho;
                  J[k - 1, k] = rho;
              }
              
              return dirichlet_lpdf(p | alpha) + log_determinant(J);
          
        }
  
}


data {
  
        int<lower=1> n_studies;
        int<lower=1> n_thr;
        // vector<lower=0>[n_thr + 1] alpha; // Induced-Dirichlet prior vector 
        array[n_studies] int<lower=1> n_non_diseased;  // number in D- group in each study
        array[n_studies] int<lower=1> n_diseased;  // number in diseased group in each study
        array[n_studies, n_thr] int x_non_diseased;  // counts for D- ("missing" threshold data recorded as 999)
        array[n_studies, n_thr] int x_diseased;  // counts for diseased     ("missing" threshold data recorded as 999)
  
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
  
        array[n_studies] ordered[n_thr] cutpoints_nd;  // Global cutpoints
        array[n_studies] ordered[n_thr] cutpoints_d;   // Global cutpoints
        real beta_d_mu;      // Mean in diseased group (0 in D- for identifiability)
        real log_scale_d_mu;  // Scale in diseased group (1 in D- for identifiability))
        real<lower=0> beta_d_SD;   
        vector[n_studies] beta_d_z;    // Study-specific random effects for beta (off-centered parameterisation)
        //// For D- group only: 
        real beta_nd_mu; 
        real<lower=0> beta_nd_SD;   
        vector[n_studies] beta_nd_z;    // Study-specific random effects for beta (off-centered parameterisation)
        //// for induced-dirichlet between-study cutpoint model:
        real<lower=0> kappa_d;
        real<lower=0> kappa_nd;
        simplex[n_cat] phi_d;
        simplex[n_cat] phi_nd;
  
}


transformed parameters { 
    
        //// Construct the study-specific random effects (off-centered param.):
        vector[n_studies] beta_d  =     beta_d_mu       + beta_d_z       .* beta_d_SD;  
        //// Estimate means in D- (but consider fixing the scales in this group for identifiability - especially if you also want cutpoints to vary between studies):
        vector[n_studies] beta_nd =     beta_nd_mu      + beta_nd_z      .* beta_nd_SD;
        //// for induced-dirichlet between-study cutpoint model:
        vector<lower=0>[n_thr + 1] alpha_d  =  phi_d*kappa_d;
        vector<lower=0>[n_thr + 1] alpha_nd =  phi_nd*kappa_nd;
        //// Cumulative ordinal probs for the likelihood:
        matrix<lower=0.0>[n_studies, n_thr] cumul_cond_prob_nd = rep_matrix(0.0, n_studies, n_thr);
        matrix<lower=0.0>[n_studies, n_thr] cumul_cond_prob_d  = rep_matrix(0.0, n_studies, n_thr);
      
        
        //// Likelihood using binomial factorization:
        for (s in 1:n_studies) {
                
                //// Non-diseased group (fixed parameters)
                if (x_non_diseased[s, 1] != 999) { //// If non-missing (i.e. only if study reports @ this ordinal cutpoint!)
                         cumul_cond_prob_nd[s, 1] = inv_logit(beta_nd[s] - cutpoints_nd[s, 1]);
                }
                
                for (k in 2:n_thr) {
                     if (x_non_diseased[s, k] != 999) { //// If non-missing (i.e. only if study reports @ this ordinal cutpoint!)
                         cumul_cond_prob_nd[s, k] = inv_logit(beta_nd[s] - cutpoints_nd[s, k]) / inv_logit(beta_nd[s] - cutpoints_nd[s, k - 1]);
                     }
                }
                
                //// Diseased group (D+):
                if (x_diseased[s, 1] != 999) { //// If non-missing (i.e. only if study reports @ this ordinal cutpoint!)
                         cumul_cond_prob_d[s, 1] = inv_logit((beta_d[s] - cutpoints_d[s, 1]));
                }
                
                for (k in 2:n_thr) {
                   if (x_diseased[s, k] != 999) { //// If non-missing (i.e. only if study reports @ this ordinal cutpoint!)
                         cumul_cond_prob_d[s, k] = inv_logit((beta_d[s] - cutpoints_d[s, k])) /  inv_logit((beta_d[s] - cutpoints_d[s, k - 1]));
                   }
                }
          
        }
      
}


model {
      
        //// For D+ group:
        target += std_normal_lpdf(beta_d_z);
        target += normal_lpdf(beta_d_mu      | 0.0, 2.0);
        target += normal_lpdf(beta_d_SD      | 0.0, 1.0);
        
        //// For D- group:
        target += std_normal_lpdf(beta_nd_z);
        target += normal_lpdf( beta_nd_mu | 0.0, 2.0);
        target += normal_lpdf( beta_nd_SD | 0.0, 1.0);
        
        //// Induced-Dirichlet priors for cutpoints:
        target += normal_lpdf(kappa_d   | 0.0, 5.0);
        target += normal_lpdf(kappa_nd  | 0.0, 5.0);
           
        // //// Induced-Dirichlet priors for cutpoints (FIXED between studies):
        // target += induced_dirichlet_lpdf( cutpoints | alpha, 0.0);
        
        //// Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        for (s in 1:n_studies)  {
                   target += induced_dirichlet_lpdf(cutpoints_d[s, 1:n_thr]  | alpha_d, 0.0);
                   target += induced_dirichlet_lpdf(cutpoints_nd[s, 1:n_thr] | alpha_nd, 0.0);
        }
      
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
          vector<lower=0.0>[n_thr] Se;
          vector<lower=0.0>[n_thr] Sp;
          vector<lower=0.0>[n_thr] Se_pred;
          vector<lower=0.0>[n_thr] Sp_pred;
          matrix<lower=0.0>[n_studies, n_thr] se;
          matrix<lower=0.0>[n_studies, n_thr] sp;
          matrix[n_studies, n_thr] x_non_diseased_est =  rep_matrix(0.0, n_studies, n_thr);
          matrix[n_studies, n_thr] x_diseased_est =      rep_matrix(0.0, n_studies, n_thr);
          matrix[n_studies, n_thr] x_non_diseased_diff = rep_matrix(0.0, n_studies, n_thr);
          matrix[n_studies, n_thr] x_diseased_diff =     rep_matrix(0.0, n_studies, n_thr);
          //// For Mean thr. parameters:
          vector[n_thr] cutpoints_nd_mu;
          vector[n_thr] cutpoints_d_mu;
          vector[n_cat] prob_ord_nd_mu;
          vector[n_cat] prob_ord_d_mu;
          
          {
          
                 int n_sims = 100;
                 array[n_sims] vector[n_cat] prob_ord_nd_mu_sim;
                 array[n_sims] vector[n_cat] prob_ord_d_mu_sim;
                    
                 for (i in 1:n_sims) {
                     prob_ord_nd_mu_sim[i,]  =  dirichlet_rng(alpha_nd);
                     prob_ord_d_mu_sim[i,]   =  dirichlet_rng(alpha_d);
                 }
                  
                 for (k in 1:(n_thr + 1)) {
                     prob_ord_nd_mu[k]       = mean(prob_ord_nd_mu_sim[, k]); 
                     prob_ord_d_mu[k]        = mean(prob_ord_d_mu_sim[, k]); 
                 }
                
                 cutpoints_d_mu[1] =         logit(prob_ord_nd_mu[1]); 
                 cutpoints_nd_mu[1] =        logit(prob_ord_d_mu[1]); 
                 for (k in 2:n_thr) {
                       cutpoints_d_mu[k] =         logit( prob_ord_nd_mu[k]   + inv_logit(cutpoints_d_mu[k - 1]) ); 
                       cutpoints_nd_mu[k] =        logit( prob_ord_d_mu[k]    + inv_logit(cutpoints_nd_mu[k - 1]) ); 
                 }
    
          }
  
          //// Calculate study-specific accuracy:
          for (s in 1:n_studies) {
                for (k in 1:n_thr) {
                    se[s, k] =       inv_logit(cutpoints_d[s, k]   - beta_d[s]);
                    sp[s, k] = 1.0 - inv_logit(cutpoints_nd[s, k]  - beta_nd[s]);
                }
          }
  
          //// Calculate summary accuracy (using mean parameters):
          for (k in 1:n_thr) {
              real scale_d_mu = 1.0;
              Se[k] =        inv_logit((cutpoints_d_mu[k]  - beta_d_mu)/scale_d_mu);
              Sp[k] =  1.0 - inv_logit(cutpoints_nd_mu[k]  - beta_nd_mu);
          }
  
          //// Calculate predictive accuracy:
          {
              real beta_d_pred  = normal_rng(beta_d_mu, beta_d_SD);
              // real scale_d_pred = exp(normal_rng(log_scale_d_mu, log_scale_d_SD));
              real scale_d_pred = 1.0;
              
              //// For D- group only:
              real beta_nd_pred = normal_rng(beta_nd_mu, beta_nd_SD);
  
              for (k in 1:n_thr) {
                  Se_pred[k] =        inv_logit((cutpoints_d_mu[k]  - beta_d_pred)/scale_d_pred);
                  Sp_pred[k] =  1.0 - inv_logit(cutpoints_nd_mu[k]  - beta_nd_pred);
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






















