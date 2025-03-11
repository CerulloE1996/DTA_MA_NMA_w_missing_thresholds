



functions {
  
       /////
       //// Function to manually construct a simplex from raw unconstrained parameters.
       //// Adapted from: ////https://github.com/mjhajharia/transforms/blob/main/transforms/simplex/StickbreakingLogistic.stan
       ////
       vector stickbreaking_logistic_simplex_constrain(vector y) {
            
                    //// Initialise Jacobian:
                    real Jacobian = 0.0;
                    
                    //// Other variables:
                    int N = rows(y) + 1;
                    vector[N] x;
                    real log_zi;
                    real log_xi;
                    real log_cum_prod = 0;
                    
                    //// Construct simplex:
                    for (i in 1:(N - 1)) {
                        log_zi = log_inv_logit(y[i] - log(N - i)); // logistic_lcdf(y[i] | log(N - i), 1)
                        log_xi = log_cum_prod + log_zi;
                        x[i] = exp(log_xi);
                        log_cum_prod += log1m_exp(log_zi);
                        Jacobian += log_xi;
                    }   
                    
                    x[N] = exp(log_cum_prod);
                    
                    //// Increment Jacobian adjustment:
                    Jacobian += log_cum_prod;
                    
                    //// Output vector (1st element is the Jacobian adjustment):
                    vector[N + 1] out_vec;
                    out_vec[1] = Jacobian;
                    out_vec[2:(N + 1)] = x;
                    return out_vec;
        
       }
       
       
}
