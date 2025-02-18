
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
  
      ordered[n_thr] C_non_diseased;  // Global cutpoints (FIXED between studies)
      ordered[n_thr] C_diseased;  // Global cutpoints (FIXED between studies)
      real beta_d_mu;      // Mean in diseased group (0 in non-diseased for identifiability)
      real log_scale_d_mu;  // Scale in diseased group (1 in non-diseased for identifiability))
      real<lower=0> beta_d_SD;   
      real<lower=0> log_scale_d_SD;
      vector[n_studies] beta_d_z;    // Study-specific random effects for beta (off-centered parameterisation)
      vector[n_studies] log_scale_d_z;    // Study-specific random effects for scale (off-centered parameterisation)
  
}

transformed parameters { 
  
      //// Construct the study-specific random effects (off-centered param.):
      vector[n_studies] beta_d =      beta_d_mu      + beta_d_z      .* beta_d_SD;  
      vector[n_studies] log_scale_d = log_scale_d_mu + log_scale_d_z .* log_scale_d_SD;  
      vector[n_studies] scale_d = exp(log_scale_d);
      
}
  


model {
      
      //// Between-study heterogenity priors:
      beta_d_z        ~ normal(0.0, 1.0);
      log_scale_d_z   ~ normal(0.0, 1.0);
      beta_d_mu       ~ normal(0.0, 1.0);
      log_scale_d_mu  ~ normal(0.0, 1.0);
      beta_d_SD       ~ normal(0.0, 1.0);
      log_scale_d_SD  ~ normal(0.0, 1.0);
      //// Induced-Dirichlet priors for cutpoints (FIXED between studies):
      C_non_diseased ~ induced_dirichlet(alpha_non_diseased, 0.0);
      C_diseased     ~ induced_dirichlet(alpha_diseased, 0.0);
       
      //// Likelihood using binomial factorization:
      for (s in 1:n_studies) {
              
              //// Non-diseased group (fixed parameters)
              if (x_non_diseased[s, 1] != 999) {
                  x_non_diseased[s, 1] ~ binomial( n_non_diseased[s], Phi(C_non_diseased[1]) );
              }
              
              for (k in 2:n_thr) {
                   if (x_non_diseased[s, k] != 999) {
                    x_non_diseased[s, k] ~ binomial( x_non_diseased[s, k-1], Phi(C_non_diseased[k]) / Phi(C_non_diseased[k-1]) );
                   }
              }
              
              //// Diseased group (D+):
              if (x_diseased[s, 1] != 999) {
                 x_diseased[s, 1] ~ binomial( n_diseased[s], Phi((beta_d[s] - C_diseased[1])/scale_d[s]) );
              }
              
              for (k in 2:n_thr) {
                 if (x_diseased[s, k] != 999) {
                      x_diseased[s, k] ~ binomial( x_diseased[s,k-1], Phi((beta_d[s] - C_diseased[k])/scale_d[s]) /  Phi((beta_d[s] - C_diseased[k-1])/scale_d[s]) );
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
        
        //// Calculate study-specific accuracy:
        for (s in 1:n_studies) {
              for (k in 1:n_thr) {
                  se[s][k] = 1.0 - Phi((C_diseased[k] - beta_d[s])/scale_d[s]);
                  sp[s][k] = Phi(C_non_diseased[k]);
              }
        }
        
        //// Calculate summary accuracy (using mean parameters):
        for (k in 1:n_thr) {
            Se[k] = 1.0 - Phi((C_diseased[k] - beta_d_mu)/exp(log_scale_d_mu));
            Sp[k] = Phi(C_non_diseased[k]);
        }
        
        //// Calculate predictive accuracy:
        {
            real beta_d_pred = normal_rng(beta_d_mu, beta_d_SD);
            real scale_d_pred = exp(normal_rng(log_scale_d_mu, log_scale_d_SD));
            
            for (k in 1:n_thr) {
                Se_pred[k] = 1 - Phi((C_diseased[k] - beta_d_pred)/scale_d_pred);
                Sp_pred[k] = Phi(C_non_diseased[k]);
            }
        }
    
}























