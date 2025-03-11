


functions {
  
          //// Normal PDF for a real input (overloaded fn):
          real normal_pdf( real x, 
                           real mu, 
                           real sigma) {
                             
               real sqrt_2_pi = sqrt(2 * pi());
               return (1.0 / (sigma * sqrt_2_pi)) * exp(-0.5 * ((x - mu) / sigma)^2);
            
          }
        
        
          //// Normal PDF for a vector input (overloaded fn):
          vector normal_pdf( vector x, 
                             real mu, 
                             real sigma) {

                real sqrt_2_pi = sqrt(2 * pi());
                return (1.0 / (sigma .* sqrt_2_pi)) * exp(-0.5 * ((x - mu) ./ sigma)^2);

          }
          
          //// Overload for vector mean, vector sigma (overloaded fn):
          vector normal_pdf( vector x, 
                             vector mu, 
                             vector sigma) {
            
                real sqrt_2_pi = sqrt(2 * pi());
                return (1.0 / (sigma * sqrt_2_pi)) .* exp(-0.5 * ((x - mu) ./ sigma)^2);
          
        }
        
        //// Median for vector (overloaded fn):
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
       
       //// Median for a row_vector (overloaded fn)
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
 