
functions {
  
        real induced_dirichlet_lpdf(  vector c, 
                                      vector alpha, 
                                      real phi) {
          
              int K = num_elements(c) + 1;
              vector[K - 1] anchoredcutoffs = c - phi;
              vector[K] sigma;
              vector[K] p;
              matrix[K, K] J = rep_matrix(0, K, K);
              
              sigma[1:(K-1)] = Phi(anchoredcutoffs);
              sigma[K] = 1;
              
              p[1] = sigma[1];
              for (k in 2:(K - 1))
                p[k] = sigma[k] - sigma[k - 1];
              p[K] = 1 - sigma[K - 1];
              
              // Jacobian computation
              for (k in 1:K) J[k, 1] = 1;
              for (k in 2:K) {
                real rho = 1.702 * sigma[k - 1] * (1 - sigma[k - 1]);
                J[k, k] = -rho;
                J[k - 1, k] = rho;
              }
              
              return dirichlet_lpdf(p | alpha) + log_determinant(J);
          
        }
  
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
  
}


parameters {
  
      array[n_studies] ordered[n_thr] C_non_diseased;
      array[n_studies] ordered[n_thr] C_diseased;
      real beta_d_mu;      // Mean in diseased group (0 in non-diseased for identifiability)
      real log_scale_d_mu;  // Scale in diseased group (1 in non-diseased for identifiability))
      real<lower=0> beta_d_SD;   
      real<lower=0> log_scale_d_SD;
      vector[n_studies] beta_d_z;    // Study-specific random effects for beta (off-centered parameterisation)
      vector[n_studies] log_scale_d_z;    // Study-specific random effects for scale (off-centered parameterisation)
      simplex[n_thr + 1] phi_d;  // for induced-dirichlet between-study cutpoint model
      simplex[n_thr + 1] phi_nd; // for induced-dirichlet between-study cutpoint model
      real<lower=0> kappa_d;  // for induced-dirichlet between-study cutpoint model
      real<lower=0> kappa_nd; // for induced-dirichlet between-study cutpoint model
  
}


transformed parameters { 
  
      //// Construct the study-specific random effects (off-centered param.):
      vector[n_studies] beta_d =      beta_d_mu      + beta_d_z      .* beta_d_SD;  
      vector[n_studies] log_scale_d = log_scale_d_mu + log_scale_d_z .* log_scale_d_SD;  
      vector[n_studies] scale_d = exp(log_scale_d);
      //// for induced-dirichlet between-study cutpoint model:
      vector<lower=0>[n_thr + 1] alpha_d =  phi_d*kappa_d;
      vector<lower=0>[n_thr + 1] alpha_nd = phi_nd*kappa_nd;
      
}


model {
      
      //// Between-study heterogenity priors:
      target += normal_lpdf(beta_d_z       | 0.0, 1.0);
      target += normal_lpdf(log_scale_d_z  | 0.0, 1.0);
      target += normal_lpdf(beta_d_mu      | 0.0, 1.0);
      target += normal_lpdf(log_scale_d_mu | 0.0, 1.0);
      target += normal_lpdf(beta_d_SD      | 0.0, 1.0);
      target += normal_lpdf(log_scale_d_SD | 0.0, 1.0);
      //// Induced-Dirichlet priors for cutpoints (FIXED between studies):
      target += normal_lpdf(kappa_d  | 0.0, 50.0);
      target += normal_lpdf(kappa_nd | 0.0, 50.0);
                  
      //// Likelihood using binomial factorization:
      //// Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
      for (s in 1:n_studies)  {
                 target += induced_dirichlet_lpdf(C_diseased[s, 1:n_thr] | alpha_d, 0.0);
                 target += induced_dirichlet_lpdf(C_non_diseased[s, 1:n_thr] | alpha_nd, 0.0);
      }
      for (s in 1:n_studies) {
              
              //// Non-diseased group (fixed parameters)
              if (x_non_diseased[s, 1] != 999) {
                       target += binomial_lpmf(  x_non_diseased[s, 1] | n_non_diseased[s], Phi(C_non_diseased[s, 1]) );
              }
              
              for (k in 2:n_thr) {
                   if (x_non_diseased[s, k] != 999) {
                       target += binomial_lpmf(  x_non_diseased[s, k] | x_non_diseased[s, k - 1], Phi(C_non_diseased[s, k]) / Phi(C_non_diseased[s, k - 1]) );
                   }
              }
              
              //// Diseased group (D+):
              if (x_diseased[s, 1] != 999) {
                       target += binomial_lpmf( x_diseased[s, 1] | n_diseased[s], Phi((beta_d[s] - C_diseased[s, 1])/scale_d[s]) );
              }
              
              for (k in 2:n_thr) {
                 if (x_diseased[s, k] != 999) {
                       target += binomial_lpmf(  x_diseased[s, k] | x_diseased[s, k - 1], Phi((beta_d[s] - C_diseased[s, k])/scale_d[s]) /  Phi((beta_d[s] - C_diseased[s, k - 1])/scale_d[s]) );
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
        //// Maan threshold parameters:
        vector[n_thr] C_diseased_mu;
        vector[n_thr] C_non_diseased_mu;
        vector[n_thr + 1] p_diseased_mu;
        vector[n_thr + 1] p_non_diseased_mu;
        
        {
          
               int n_sims = 1000;
               vector[n_thr + 1] p_diseased_mu_sim[n_sims];
               vector[n_thr + 1] p_non_diseased_mu_sim[n_sims];
                
               for (i in 1:n_sims) {
                  p_non_diseased_mu_sim[i,]  =  dirichlet_rng(alpha_nd);
                  p_diseased_mu_sim[i,]      =  dirichlet_rng(alpha_d);
               }
                
               for (k in 1:(n_thr + 1)) {
                   p_non_diseased_mu[k] = mean(p_non_diseased_mu_sim[, k]); 
                   p_diseased_mu[k]     = mean(p_diseased_mu_sim[, k]); 
               }
              
               C_non_diseased_mu[1] =    inv_Phi(p_non_diseased_mu[1]); 
               C_diseased_mu[1] =        inv_Phi(p_diseased_mu[1]); 
               for (k in 2:n_thr) {
                     C_non_diseased_mu[k] =    inv_Phi(p_non_diseased_mu[k] + Phi(C_non_diseased_mu[k - 1])); 
                     C_diseased_mu[k] =        inv_Phi(p_diseased_mu[k] +     Phi(C_diseased_mu[k - 1])); 
               }
    
        }
        
        //// Calculate study-specific accuracy:
        for (s in 1:n_studies) {
              for (k in 1:n_thr) {
                  se[s][k] =       Phi((C_diseased[s, k] - beta_d[s])/scale_d[s]);
                  sp[s][k] = 1.0 - Phi(C_non_diseased[s, k]);
              }
        }
        
        //// Calculate summary accuracy (using mean parameters):
        for (k in 1:n_thr) {
            Se[k] =       Phi((C_diseased_mu[k] - beta_d_mu)/exp(log_scale_d_mu));
            Sp[k] = 1.0 - Phi(C_non_diseased_mu[k]);
        }
        
        //// Calculate predictive accuracy:
        {
            real beta_d_pred =  normal_rng(beta_d_mu, beta_d_SD);
            real scale_d_pred = exp(normal_rng(log_scale_d_mu, log_scale_d_SD));
            
            for (k in 1:n_thr) {
                Se_pred[k] =       Phi((C_diseased_mu[k] - beta_d_pred)/scale_d_pred);
                Sp_pred[k] = 1.0 - Phi(C_non_diseased_mu[k]);
            }
        }
    
}























