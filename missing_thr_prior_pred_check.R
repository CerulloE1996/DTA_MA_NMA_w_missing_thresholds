# K <- 3
# prior_kappa_sd <- 50
# 
# 
# #############################
# # prior predictive check plot
# 
# require(ggplot2)
# # require(hexbin)
# # require(latex2exp)
# 
# 
# # N <- 1e4
# # K <- 5
# # prior_kappa_sd <- 10

require(ggplot2)
require(dplyr)


# use_alpha_directly = TRUE
# use_log_alpha = FALSE
# use_log_kappa = FALSE
# log_alpha_lb = Stan_data_list$log_alpha_lb
# log_alpha_ub = Stan_data_list$log_alpha_ub
# prior_mean = Stan_data_list$prior_kappa_mean[1]
# prior_sd = Stan_data_list$prior_kappa_SD[1]
# n_cat = n_cat
# N = 10000

## method is either "kappa" or "alpha"

induced_Dirichlet_ppc_plot <- function(  method = "kappa", 
                                         use_log_alpha,
                                         log_alpha_lb,
                                         log_alpha_ub,
                                         use_log_kappa, 
                                         prior_mean,
                                         prior_sd, 
                                         prior_dirichlet_phi,
                                         n_cat, 
                                         N = 5000) {
  
            p_if_uniform <- 1/n_cat
              
            message(paste("p_if_uniform = ", p_if_uniform))
            
            kappa <- NULL
            res   <- array(NA, dim = c(n_cat, N))
            alpha <- array(NA, dim = c(n_cat, N))
            log_alpha <- array(NA, dim = c(n_cat, N))

            if (method == "kappa") {
             
                              if (use_log_kappa == TRUE) {
                                
                                    log_kappa <- truncnorm::rtruncnorm(n = N, mean = prior_mean, sd = prior_sd, a = log_alpha_lb, b = log_alpha_ub)
                                    kappa <- exp(log_kappa)
                                    #### log_kappa <- rnorm(n = N, mean = prior_mean, sd = prior_sd)
                                 
                              } else { 
                                
                                    kappa <- truncnorm::rtruncnorm(N, a = exp(log_alpha_lb), b = exp(log_alpha_ub), mean = prior_mean, sd = prior_sd)
                                    log_kappa <<- log(kappa)
                                    # kappa <- rgamma( n = N, 
                                    #                  shape = prior_mean, 
                                    #                  rate = prior_sd) ##  (N, a = 0.0, mean = prior_mean, sd = prior_sd)
                                    
                              }
                             
                              phi <- MCMCpack::rdirichlet(N, prior_dirichlet_phi) 
                              log_phi <- log(phi)
                                    
                              for (i in 1:N) {
                                 
                                    log_alpha[,i] <- log_kappa[i] + log_phi[i, 1:n_cat]
                                    ## Compute alpha:
                                    alpha[,i] <- exp(log_alpha[,i])
                                    ##
                                    res[,i] <- MCMCpack::rdirichlet(n = 1, alpha[,i])
                                    
                              }
                     
            } else if (method == "alpha") {
              
                            for (i in 1:N) {
                              
                                    # {
                                      
                                      if (use_log_alpha == TRUE) {
                                        log_alpha[,i] <- truncnorm::rtruncnorm(n = n_cat, mean = prior_mean, sd = prior_sd, a = log_alpha_lb, b = log_alpha_ub)
                                        alpha[,i] <- exp(log_alpha[,i])
                                      } else { 
                                        alpha[,i] <- truncnorm::rtruncnorm(n = n_cat, mean = prior_mean, sd = prior_sd, a = exp(log_alpha_lb), b = exp(log_alpha_ub))
                                      }
                                      
                                    #  else { 
                                    #   
                                    #   alpha[,i] <- log_kappa[i] + log(phi[i, 1:n_cat])
                                    #   
                                    # }
                                    
                                    res[,i] <- MCMCpack::rdirichlet(n = 1, alpha[,i])
                              
                            }
                            
            }
            

            
            df <- data.frame(p1 = res[1,], p2 = res[2,], dist = 1) %>% 
                  dplyr::filter(!is.na(p1), !is.na(p2))     
            
            g1 <- ggplot(df, aes(p1, p2)) + 
              geom_hex() + 
              scale_fill_continuous(trans = "log10") +
              theme_bw() + 
              # xlab(TeX("$P_{i}$")) + 
              # ylab(TeX("$P_{j}$")) + 
              xlim(0, 3.0*p_if_uniform) + 
              ylim(0, 3.0*p_if_uniform) + 
              geom_hline(yintercept = p_if_uniform) + 
              geom_vline(xintercept = p_if_uniform)
            
            print(g1)
            
            return(list(res = res, 
                        alpha = alpha,
                        kappa = kappa,
                        log_alpha = log_alpha))

}










