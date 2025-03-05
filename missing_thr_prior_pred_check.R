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


# method = "sigma"
# use_log_alpha = FALSE
# use_log_kappa = FALSE
# log_alpha_lb = Stan_data_list$log_alpha_lb
# log_alpha_ub = Stan_data_list$log_alpha_ub
# prior_mean = Stan_data_list$prior_kappa_mean[1]
# prior_sd = Stan_data_list$prior_kappa_SD[1]
# n_cat = n_cat
# N = 5000
# prior_dirichlet_cat_means_alpha <- rep(1, n_cat)
# prior_dirichlet_cat_SDs_mean <- rep(0.0, n_cat)
# prior_dirichlet_cat_SDs_SD <- rep(0.10, n_cat)





## method is either "sigma" or kappa" or "alpha":
induced_Dirichlet_ppc_plot <- function(  method = "sigma", 
                                         N = 5000,
                                         n_cat,
                                         other_args_list
) {
  
            p_if_uniform <- 1/n_cat
              
            message(paste("p_if_uniform = ", p_if_uniform))
            
            kappa <- NULL
            res   <- array(NA, dim = c(n_cat, N))
            alpha <- array(NA, dim = c(n_cat, N))
            log_alpha <- array(NA, dim = c(n_cat, N))

            
            if (method == "sigma") { 
              
                              ## Extract arguments from other_args_list:
                              prior_dirichlet_cat_means_alpha <- other_args_list$prior_dirichlet_cat_means_alpha
                              prior_dirichlet_cat_SDs_mean <- other_args_list$prior_dirichlet_cat_SDs_mean
                              prior_dirichlet_cat_SDs_SD <- other_args_list$prior_dirichlet_cat_SDs_SD
                      
                              # Initialize arrays:
                              phi <- matrix(NA, nrow = N, ncol = n_cat)
                              category_sds <- matrix(NA, nrow = n_cat, ncol = N)
                              kappa <- numeric(N)
                              
                              # Generate means from Dirichlet distribution
                              for (i in 1:N) {
                                phi[i,] <- gtools::rdirichlet(1, prior_dirichlet_cat_means_alpha)
                              }
                              
                              # Generate SD values from truncated normal (ensure they're positive)
                              for (k in 1:n_cat) {
                                category_sds[k,] <- truncnorm::rtruncnorm(N, 
                                                                          a = 0, 
                                                                          mean = prior_dirichlet_cat_SDs_mean[k], sd = prior_dirichlet_cat_SDs_SD[k])
                              }
                              
                              
                              ## Define rgw error function to optimise (to find alpha given phi and SD's):
                              error_function <- function(alpha_0, 
                                                         current_phi, 
                                                         SDs) {
                                
                                temp_alpha <- current_phi * alpha_0
                                calc_sds <- sqrt((temp_alpha * (alpha_0 - temp_alpha)) /  (alpha_0^2 * (alpha_0 + 1.0)))
                                
                                return(sum((calc_sds - SDs)^2))
                                
                              }
                              
                              ## Work backwards to find alpha GIVEN the SD and phi (i.e. category mean probabilities)
                              for (i in 1:N) {
                                
                                      ## Get the current phi and target_sds for this iteration:
                                      current_phi <- phi[i,]
                                      SDs <- category_sds[,i]
                                      
                                      # Optimize to find alpha_0:
                                      opt_result <- optimize(error_function, 
                                                             interval = c(exp(log_alpha_lb), 1000), 
                                                             tol = 1e-8,
                                                             current_phi = current_phi,
                                                             SDs = SDs)
                                      
                                      alpha_0 <- opt_result$minimum
                                      
                                      ## Store kappa (alpha_0):
                                      kappa[i] <- alpha_0
                                      
                                      ## Calculate alpha using alpha_0 and phi:
                                      alpha[,i] <- current_phi * alpha_0
                                      log_alpha[,i] <- log(alpha[,i])
                                      
                                      ## Sample from Dirichlet with these alphas:
                                      res[,i] <- gtools::rdirichlet(n = 1, alpha[,i])
                                
                              }
              
            } else if (method == "kappa") {
              
             
                              if (use_log_kappa == TRUE) {
                                
                                        ## Extract arguments from other_args_list:
                                        prior_log_kappa_mean <- other_args_list$prior_log_kappa_mean
                                        prior_log_kappa_sd <- other_args_list$prior_log_kappa_sd                             
                                        log_kappa_lb <- other_args_list$log_kappa_lb
                                        log_kappa_ub <- other_args_list$log_kappa_ub
                                        
                                        if (is.null(log_kappa_ub)) { 
                                          log_kappa_ub <- Inf
                                        }
                                        if (is.null(log_kappa_lb)) { 
                                          log_kappa_lb <- -Inf
                                        }
                                    
                                        log_kappa <- truncnorm::rtruncnorm(n = N, 
                                                                           a = log_kappa_lb, b = log_kappa_ub,
                                                                           mean = prior_log_kappa_mean, sd = prior_log_kappa_sd)
                                                                          
                                        kappa <- exp(log_kappa)
                                 
                              } else { 
                                
                                        ## Extract arguments from other_args_list:
                                        prior_kappa_mean <- other_args_list$prior_kappa_mean
                                        prior_kappa_sd <- other_args_list$prior_kappa_sd                             
                                        kappa_lb <- other_args_list$kappa_lb
                                        kappa_ub <- other_args_list$kappa_ub
                                        
                                        if (is.null(kappa_ub)) { 
                                          kappa_ub <- Inf
                                        }
                                        if (is.null(kappa_lb)) { 
                                          kappa_lb <- -Inf
                                        }
                                        
                                        kappa <- truncnorm::rtruncnorm(N, 
                                                                       a = kappa_lb, b = kappa_ub,
                                                                       mean = prior_kappa_mean, sd = prior_kappa_sd)
                                        log_kappa <<- log(kappa)
                                    
                              }
                             
                              phi <- gtools::rdirichlet(N, prior_dirichlet_cat_means_alpha) 
                         
                              log_phi <- log(phi)
                                    
                              for (i in 1:N) {
                                 
                                    log_alpha[,i] <- log_kappa[i] + log_phi[i, 1:n_cat]
                                    ## Compute alpha:
                                    alpha[,i] <- exp(log_alpha[,i])
                                    ##
                                    res[,i] <- gtools::rdirichlet(n = 1, alpha[,i])
                                    
                              }
                     
            } else if (method == "alpha") {
              
                            for (i in 1:N) {
                                      
                                      if (use_log_alpha == TRUE) {
                                        
                                                ## Extract arguments from other_args_list:
                                                prior_log_alpha_mean <- other_args_list$prior_log_alpha_mean
                                                prior_log_alpha_sd <- other_args_list$prior_log_alpha_sd                             
                                                log_alpha_lb <- other_args_list$log_alpha_lb
                                                log_alpha_ub <- other_args_list$log_alpha_ub
                                                
                                                if (is.null(log_alpha_lb)) { 
                                                  log_alpha_lb <- Inf
                                                }
                                                if (is.null(log_alpha_ub)) { 
                                                  log_alpha_ub <- -Inf
                                                }
                                                
                                                log_alpha[,i] <- truncnorm::rtruncnorm(n = n_cat, 
                                                                                       mean = prior_log_alpha_mean, sd = prior_log_alpha_sd, 
                                                                                       a = log_alpha_lb, b = log_alpha_ub)
                                                alpha[,i] <- exp(log_alpha[,i])
                                        
                                      } else { 
                                        
                                                ## Extract arguments from other_args_list:
                                                prior_alpha_mean <- other_args_list$prior_alpha_mean
                                                prior_alpha_sd <- other_args_list$prior_alpha_sd                             
                                                alpha_lb <- other_args_list$alpha_lb
                                                alpha_ub <- other_args_list$alpha_ub
                                                
                                                if (is.null(alpha_lb)) { 
                                                  alpha_lb <- Inf
                                                }
                                                if (is.null(alpha_ub)) { 
                                                  alpha_ub <- -Inf
                                                }
                                            
                                                alpha[,i] <- truncnorm::rtruncnorm(n = n_cat, 
                                                                                   mean = prior_alpha_mean, sd = prior_alpha_sd, 
                                                                                   a = alpha_lb, b = alpha_ub)
                                        
                                      }
                                      
                                    
                                    res[,i] <- gtools::rdirichlet(n = 1, alpha[,i])
                              
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










