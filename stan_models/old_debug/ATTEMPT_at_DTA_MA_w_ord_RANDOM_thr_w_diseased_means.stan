functions {
  
      real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
        int K = num_elements(c) + 1;
        vector[K - 1] sigma = inv_logit(phi - c);
        vector[K] p;
        matrix[K, K] J = rep_matrix(0, K, K);
        
        // Induced ordinal probabilities
        p[1] = 1 - sigma[1];
        for (k in 2:(K - 1))
          p[k] = sigma[k - 1] - sigma[k];
        p[K] = sigma[K - 1];
        
        // Jacobian matrix construction (same as your current model)
        for (k in 1:K) J[k, 1] = 1;
        for (k in 2:K) {
          real rho = sigma[k - 1] * (1 - sigma[k - 1]);
          J[k, k] = - rho;
          J[k - 1, k] = rho;
        }
        
        return dirichlet_lpdf(p | alpha) + log_determinant(J);
      }
      
}


data {
  
      int<lower=1> n_studies;
      int<lower=1> n_thr;
      vector<lower=0>[n_thr + 1] alpha_nd; 
      vector<lower=0>[n_thr + 1] alpha_d;
      //// Ordinal test data:
      array[n_studies] int<lower=1> n_non_diseased;
      array[n_studies] int<lower=1> n_diseased; 
      array[n_studies, n_thr] int x_non_diseased;
      array[n_studies, n_thr] int x_diseased;
      //// Binary reference test data:
      array[n_studies] int<lower=1> n_total; // total patients per study
      array[n_studies] int<lower=0> n_ref_pos; // number of pts. testing +'ve on reference test
  
}


parameters {
  
      // Study-specific parameters
      vector[n_studies] prev_z; // prevalence
      vector[n_studies] beta_d_z;
      vector[n_studies] beta_nd_z;
      vector[n_studies] log_scale_d_z;
      
      // Population parameters
      real prev_mu;
      real<lower=0> prev_sigma;
      real beta_d_mu;
      real beta_nd_mu;
      real log_scale_d_mu;
      real<lower=0> beta_d_sigma;
      real<lower=0> beta_nd_sigma;
      real<lower=0> log_scale_d_sigma;
      
      // Reference test accuracy
      real<lower=0,upper=1> Se_ref; // sensitivity
      real<lower=0,upper=1> Sp_ref; // specificity
      
      // Cutpoint parameters  
      array[n_studies] ordered[n_thr] C_d;
      array[n_studies] ordered[n_thr] C_nd;
      
      // Dirichlet parameters
      real<lower=0> kappa_d;
      real<lower=0> kappa_nd;
  
}


transformed parameters {
      
      vector[n_studies] prev = inv_logit(prev_mu + prev_sigma * prev_z);
      vector[n_studies] beta_d = beta_d_mu + beta_d_sigma * beta_d_z;
      vector[n_studies] beta_nd = beta_nd_mu + beta_nd_sigma * beta_nd_z;
      vector[n_studies] log_scale_d = log_scale_d_mu + log_scale_d_sigma * log_scale_d_z;
      vector[n_studies] scale_d = exp(log_scale_d);
  
}


model {
  
      // Priors
      prev_z ~ normal(0, 1);
      beta_d_z ~ normal(0, 1);
      beta_nd_z ~ normal(0, 1);
      log_scale_d_z ~ normal(0, 1);
      
      prev_mu ~ normal(0, 1);
      prev_sigma ~ normal(0, 1);
      beta_d_mu ~ normal(0, 1);
      beta_nd_mu ~ normal(0, 1);
      log_scale_d_mu ~ normal(0, 1);
      beta_d_sigma ~ normal(0, 1);
      beta_nd_sigma ~ normal(0, 1);
      log_scale_d_sigma ~ normal(0, 1);
      
      Se_ref ~ beta(2, 2);
      Sp_ref ~ beta(2, 2);
      
      kappa_d ~ normal(0, 50);
      kappa_nd ~ normal(0, 50);
      
      
      for(s in 1:n_studies) {
        
            for(k in 1:n_thr) {
              
                    if(ref_pos_ord_above[s,k] != 999) {
                      
                            vector[2] lp;
                            
                            // D=1 component
                            real p_ref_pos_d = Se_ref; 
                            real p_ord_above_d = Phi((beta_d[s] - C_d[s,k])/scale_d[s]);
                            real joint_prob_d = p_ref_pos_d * p_ord_above_d + cov_d[s];
                            
                            // D=0 component  
                            real p_ref_pos_nd = 1 - Sp_ref;
                            real p_ord_above_nd = Phi(beta_nd[s] - C_nd[s,k]);
                            real joint_prob_nd = p_ref_pos_nd * p_ord_above_nd + cov_nd[s];
                            
                            // Mixture components
                            lp[1] = log(prev[s]) + multinomial_lpmf(...); // Need to complete joint distribution
                            lp[2] = log1m(prev[s]) + multinomial_lpmf(...);
                            
                            target += log_sum_exp(lp);
                            
                            // Add constraints on covariances:
                            real min_d = min(Se_ref, p_ord_above_d);
                            real min_nd = min(1-Sp_ref, p_ord_above_nd);
                            
                            target += (cov_d[s] <= min_d) ? 0 : negative_infinity();
                            target += (cov_nd[s] <= min_nd) ? 0 : negative_infinity();
                            // Add lower bound constraints too
                            
                    }
                    
            }
            
  }
  
  
}
























