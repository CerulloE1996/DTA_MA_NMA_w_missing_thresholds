


functions {
  
            
          //// Box-Cox transform function:
          real box_cox( real x, 
                        real lambda) {
                
                if (lambda == 0.0) {
                  
                  if (x != 0.0) {
                    return log(x);
                  } else { 
                    return -700.0;
                  }
                  
                } else {
                  
                    return (pow(x, lambda) - 1.0) / lambda;
                    
                }
              
          }
          
          //// Vectorized Box-Cox transform function:
          vector box_cox( vector x, 
                          real lambda) {
                                   
                int N = num_elements(x);
                vector[N] result;
                
                for (n in 1:N) {
                  
                    if (lambda == 0.0) {
                        if (x[n] != 0.0) {
                             result[n] = log(x[n]); 
                        } else { 
                             result[n] = -700.0;
                        }
                    } else {
                         result[n] = (pow(x[n], lambda) - 1.0) / lambda;
                    }
                    
                }
                
                return result;
              
          }
          
          
          
          
          
        
  
}


 