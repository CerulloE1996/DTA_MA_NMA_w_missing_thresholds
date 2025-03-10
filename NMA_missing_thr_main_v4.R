


 

## Set wd:
setwd("/home/enzocerullo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")
setwd("/home/enzo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")




require(BayesMVP)
require(TruncatedNormal)

# install.packages("TruncatedNormal")
# remotes::install_github("tylermorganwall/spacefillr")
# 
# 
# install.packages("spacefillr", 
#                  configure.args = "--with-include-dir=include",
#                  type = "source")
# 
# Sys.setenv(CPPFLAGS = paste(Sys.getenv("CPPFLAGS"), "-I./include", sep = " "))
# Sys.setenv(CPPFLAGS = paste(Sys.getenv("CPPFLAGS"), "-I./inst/include", sep = " "))
# install.packages("spacefillr", type = "source")
# 
# 
# Sys.setenv(CXX="g++")
# Sys.setenv(CC="gcc")
# install.packages("MCMCpack")
# 
# 
# install.packages("/home/enzo/Downloads/spacefillr_download/spacefillr_fixed.tar.gz", repos = NULL, type = "source")


source("R_fn_load_data_ordinal_NMA_LC_MVP_sim.R")
source("NMA_missing_thr_prep_Stan_data.R")
source("missing_thr_prior_pred_check.R")
source("R_fn_compile_Stan_model.R")

options(scipen = 999999999999)

 


prior_vis_test_index <- 3
n_thr_t <- n_thr[prior_vis_test_index]
n_cat_t <- n_thr_t + 1

#### To overwrite some model options:
{
    
    ## Model options for JONES model:
    ##
    Stan_data_list$use_box_cox <- 1
    ##  
    ## Model options for "Cerullo" model:
    ##
    Stan_data_list$estimate_scales <- 0
    #### Stan_data_list$same_cutpoints_between_groups <- 1
    ##
    Stan_data_list$prior_only <- 0
    ##
    Stan_data_list$use_empirical_cutpoint_means <- 0
    ##
    
      method <- "sigma"
    # method <- "kappa"
    # method <- "alpha"
    
    
    if (method == "sigma") { 
      
            prior_dirichlet_cat_SDs_mean <- rep(0.00, n_cat_t)
            # prior_dirichlet_cat_SDs_SD   <- rep(0.50, n_cat) ## example of when "vague" prior is very NOT vague and bad !!!
            prior_dirichlet_cat_SDs_SD   <- rep(0.10, n_cat_t) ## example of when "vague" prior is very NOT vague and bad !!!
            prior_dirichlet_cat_means_alpha <- rep(1.00, n_cat_t)
            ##
            Stan_data_list$prior_dirichlet_cat_means_alpha <- prior_dirichlet_cat_means_alpha
            Stan_data_list$prior_dirichlet_cat_SDs_mean <- prior_dirichlet_cat_SDs_mean
            Stan_data_list$prior_dirichlet_cat_SDs_SD <- prior_dirichlet_cat_SDs_SD
            ##
            other_args_list <- list(prior_dirichlet_cat_means_alpha = prior_dirichlet_cat_means_alpha,
                                    prior_dirichlet_cat_SDs_mean = prior_dirichlet_cat_SDs_mean,
                                    prior_dirichlet_cat_SDs_SD = prior_dirichlet_cat_SDs_SD)
            ##
            inf_dir_samples <- induced_Dirichlet_ppc_plot(   method = method,
                                                             N = 5000,
                                                             n_cat = n_cat_t,
                                                             other_args_list = other_args_list)
      
      
    } else if (method == "kappa") {
  
      
            prior_kappa_mean <- 0
            prior_kappa_SD <-   200
            prior_dirichlet_cat_means_alpha <- rep(1, n_cat_t)
            ##
            Stan_data_list$prior_kappa_mean <-  rep(prior_kappa_mean, 1)
            Stan_data_list$prior_kappa_SD <-    rep(prior_kappa_SD, 1)
            Stan_data_list$prior_dirichlet_cat_means_alpha <- prior_dirichlet_cat_means_alpha
            ##
            kappa_lb <- 1
            Stan_data_list$log_alpha_lb <- log(kappa_lb)
            Stan_data_list$log_alpha_ub <- +Inf
            print(paste("alpha_lb = ")) ; print(exp( Stan_data_list$log_alpha_lb ))
            print(paste("alpha_ub = ")) ; print(exp( Stan_data_list$log_alpha_ub ))
            ##
            inf_dir_samples <- induced_Dirichlet_ppc_plot(   method = method,
                                                             use_log_alpha = FALSE,
                                                             use_log_kappa = FALSE,
                                                             log_alpha_lb = Stan_data_list$log_alpha_lb,
                                                             log_alpha_ub = Stan_data_list$log_alpha_ub,
                                                             prior_mean = Stan_data_list$prior_kappa_mean,
                                                             prior_sd = Stan_data_list$prior_kappa_SD,
                                                             prior_dirichlet_phi =   prior_dirichlet_cat_means_alpha,
                                                             n_cat = n_cat_t,
                                                             N = 5000)
    
    } else if (method == "alpha") {
            
            # ##
            # prior_alpha_mean <-  0
            # prior_alpha_SD   <-  10
            # ##
            # prior_dirichlet_cat_means_alpha <- rep(1, n_cat_t) #'# dummy
            # Stan_data_list$prior_kappa_mean <-  rep(prior_alpha_mean, 2)
            # Stan_data_list$prior_kappa_SD <-    rep(prior_alpha_SD, 2)
            # Stan_data_list$prior_dirichlet_cat_means_alpha <- list(prior_dirichlet_cat_means_alpha, prior_dirichlet_cat_means_alpha)
            # ##
            # Stan_data_list$log_alpha_lb <- log(0.1)
            # Stan_data_list$log_alpha_ub <- +Inf
            # print(paste("alpha_lb = ")) ; print(exp( Stan_data_list$log_alpha_lb ))
            # print(paste("alpha_ub = ")) ; print(exp( Stan_data_list$log_alpha_ub ))
            # 
            # 
            # inf_dir_samples <- induced_Dirichlet_ppc_plot(   method = method,
            #                                                  use_log_alpha = FALSE,
            #                                                  use_log_kappa = FALSE,
            #                                                  log_alpha_lb = Stan_data_list$log_alpha_lb,
            #                                                  log_alpha_ub = Stan_data_list$log_alpha_ub,
            #                                                  prior_mean = prior_alpha_mean,
            #                                                  prior_sd = prior_alpha_SD,
            #                                                  n_cat = n_cat_t,
            #                                                  N = 5000)
    
    }
    
    kappa <- inf_dir_samples$kappa
    alpha <- inf_dir_samples$alpha
    # 
    # plot(density((alpha)))
    # plot(xlim = c(0, 25), density((kappa)))
    
    # plot(density(inf_dir_samples$kappa))
    
    # plot(density(inf_dir_samples[1, ])) 
    # plot(density(inf_dir_samples[5, ]))
    # plot(density(inf_dir_samples[10, ]))
    # 
    # rowMeans(inf_dir_samples)
    # 
    
    print(mean(inf_dir_samples$alpha))
    
    
    # str(inf_dir_samples)
    
    # alpha
    # 

    # # print(mean(alpha))
    # print(quantile(alpha, c(0.025, 0.50, 0.975)))
    
    # dev.off()
    # plot(density(inf_dir_samples$alpha), xlim = c(0, 100))
    print(mean(inf_dir_samples$alpha))
    print(quantile(inf_dir_samples$alpha, c(0.025, 0.50, 0.975)))

    # # plot(density(exp(inf_dir_samples$log_alpha)))
    # print(mean(exp(inf_dir_samples$log_alpha)))
    # print(quantile(exp(inf_dir_samples$log_alpha), c(0.025, 0.50, 0.975)))
    
}





{
      Stan_data_list$prior_dirichlet_cat_means_alpha <- list()
      Stan_data_list$prior_dirichlet_cat_SDs_mean    <- list()
      Stan_data_list$prior_dirichlet_cat_SDs_SD      <- list()
      for (t in 1:n_index_tests) {
          Stan_data_list$prior_dirichlet_cat_means_alpha[[t]] <- rep(1.0, n_thr_max + 1)
          Stan_data_list$prior_dirichlet_cat_SDs_mean[[t]]    <- rep(0.0, n_thr_max + 1)
          Stan_data_list$prior_dirichlet_cat_SDs_SD[[t]]      <- rep(0.10, n_thr_max + 1)
      }
      ##
      ##
      Stan_data_list$kappa_lb <- c()
      for (t in 1:n_index_tests) {
         Stan_data_list$kappa_lb[t] <- 0
      }
}


n_total_raw_simplex_elements <- sum(Stan_data_list$n_thr)






simplex_to_unconstrained <- function(x) {

      N <- length(x)
      y <- numeric(N-1)

      # Calculate the stick-breaking proportions
      z <- numeric(N-1)
      remaining_mass <- 1.0

      for (i in 1:(N-1)) {
        z[i] <- x[i] / remaining_mass
        remaining_mass <- remaining_mass - x[i]
      }

      # Transform to unconstrained space using the inverse logit
      for (i in 1:(N-1)) {
        y[i] <- qlogis(z[i]) + log(N-i)  # inverse of log_inv_logit(y[i] - log(N - i))
      }

      return(y)

}

simplex_to_unconstrained <- function(x) {
      
      N <- length(x)
      y <- numeric(N-1)
      
      # Small constant to avoid exact 0 or 1 values
      epsilon <- 1e-10
      
      # Ensure x is a valid simplex and avoid exact 0s and 1s
      x <- pmax(x, epsilon)
      x <- pmin(x, 1-epsilon)
      x <- x / sum(x)  # Re-normalize to ensure sum is 1
      
      # Calculate the stick-breaking proportions
      remaining_mass <- 1.0
      for (i in 1:(N-1)) {
        # Calculate proportion with numerical safeguards
        z_i <- min(max(x[i] / remaining_mass, epsilon), 1-epsilon)
        
        # Transform to unconstrained space
        y[i] <- log(z_i/(1-z_i)) + log(N-i)
        
        # Update remaining mass
        remaining_mass <- remaining_mass - x[i]
        
        # Avoid numerical issues when remaining mass gets close to 0
        if (remaining_mass < epsilon) {
          # Fill remaining elements with reasonable defaults
          if (i < (N-1)) {
            y[(i+1):(N-1)] <- rnorm(N-1-i, 0, 1)
          }
          break
        }
      }
      
      return(y)
  
}




dirichlet_cat_means_phi_raw_mat <- array(0, dim = c(n_index_tests, n_cat_max - 1))
for (t in 1:n_index_tests) {
  
      n_cat_t <- Stan_data_list$n_cat[t]
      uniform_simplex <- rep(1/n_cat_t, n_cat_t)  # Start with uniform
      # Add some variation while ensuring it remains a valid simplex
      noise <- runif(n_cat_t, -0.01, 0.01)
      informative_simplex <- uniform_simplex + noise
      # Make sure it's still a valid simplex
      informative_simplex <- informative_simplex / sum(informative_simplex)
      
      dirichlet_cat_means_phi_raw_mat[t, 1:(n_cat_t - 1)] <- simplex_to_unconstrained(informative_simplex)
  
}

dirichlet_cat_means_phi_raw_mat

dirichlet_cat_means_phi_raw_vec <- c()
counter <- 1
for (t in 1:n_index_tests) {
  for (k in 1:Stan_data_list$n_thr[t]) {
   dirichlet_cat_means_phi_raw_vec[counter] <- dirichlet_cat_means_phi_raw_mat[t, k]
   counter <- counter + 1
  }
}

length(dirichlet_cat_means_phi_raw_vec)
 
{

    Stan_data_list$n_total_cutpoints <- n_studies*(sum(n_thr[2:length(n_thr)]))
    Stan_data_list$n_total_raw_simplex_elements <- sum(n_thr[2:length(n_thr)])
    
}

Stan_data_list$kappa_lb <- rep(1.0, n_index_tests)


##  | ------   Initial values:  -------------------------------------------------------------------------
{
        Stan_init_list <- list()
        ##
        Stan_init_list$beta_mu      <- rep(0.001, n_index_tests)
        Stan_init_list$raw_scale_mu <- rep(0.001, n_index_tests)
        ##
        Stan_init_list$C_raw_vec <- rep(-2.0, Stan_data_list$n_total_cutpoints)
        Stan_init_list$C_raw_mat <- list()
        for (t in 1:n_index_tests) {
           Stan_init_list$C_raw_mat[[t]] <- array(-2.0, dim = c(n_studies, n_thr_max))
        }
        Stan_init_list$dirichlet_cat_means_phi_raw_mat <- dirichlet_cat_means_phi_raw_mat
        Stan_init_list$dirichlet_cat_means_phi_raw_vec <- dirichlet_cat_means_phi_raw_vec
        Stan_init_list$kappa <- rep(10, n_index_tests)
        ##
        ## For "NMA" (Nyaga-like) params:
        ##
        Stan_init_list$beta_eta_z      <- rep(0.001, n_studies)
        Stan_init_list$beta_sigma      <- 0.001
        Stan_init_list$raw_scale_eta_z <- rep(0.001, n_studies)
        Stan_init_list$raw_scale_sigma <- 0.001
        ##
        Stan_init_list$beta_tau          <- rep(0.001, n_index_tests)
        Stan_init_list$beta_delta_z      <- array(0.001, dim = c(n_studies, n_index_tests))
        Stan_init_list$raw_scale_tau     <- rep(0.001, n_index_tests)
        Stan_init_list$raw_scale_delta_z <- array(0.001, dim = c(n_studies, n_index_tests))
}





#### --------- Select NMA model type:
# Model_type <- "Jones_Nyaga"
##
# Model_type <- "Cerullo_Nyaga_Xu_FIXED_cutpoints"
# Model_type <- "Cerullo_Nyaga_Xu_RANDOM_cutpoints"
##
# Model_type <- "Cerullo_Nyaga_Gat_FIXED_cutpoints"
Model_type <- "Cerullo_Nyaga_Gat_RANDOM_cutpoints"
 
 



##
#  cutpoint_param <- "kappa"
cutpoint_param <- "sigma"








## - | ----------  Compile Stan model --------------------------------------------------------------------
{   
  
           if (Model_type == "Jones_Nyaga") {
                             ##
                             file <- file.path(getwd(), "stan_models", "DTA_MA_JONES_BOXCOX.stan")
                         
           } else if (Model_type == "Cerullo_Nyaga_Xu_FIXED_cutpoints") {

                             # Stan_init_list$C_nd <-  c(array(dim = c(n_thr, 1), data = seq(from = -2.0, to = 2.0, length = n_thr)))
                             # Stan_init_list$C_d  <-  c(array(dim = c(n_thr, 1), data = seq(from = -2.0, to = 2.0, length = n_thr)))
                           
                             file <- file.path(getwd(), "stan_models", "DTA_NMA_Nyaga_Xu_FIXEDthr.stan")

           } else if (Model_type == "Cerullo_Nyaga_Xu_RANDOM_cutpoints") {

                         
                            if (cutpoint_param == "kappa") {

                                        file <- file.path(getwd(), "stan_models", "DTA_NMA_Nyaga_Xu_RANDthr_kappa.stan")

                            } else if (cutpoint_param == "sigma") {

                                        file <- file.path(getwd(), "stan_models", "DTA_NMA_Nyaga_Xu_RANDthr.stan")
                            }
                           
             } else if (Model_type == "Cerullo_Nyaga_Gat_FIXED_cutpoints") {
                           
                             # Stan_init_list$C <-   seq(from = -2.0, to = 2.0, length = n_thr)
                             # Stan_init_list$raw_scale_mu <-  0.55
                             # Stan_init_list$raw_scale_SD <-  0.001
                             # Stan_init_list$raw_scale_z <- rep(0.001, n_studies)
                             # Stan_init_list$bs_L_Omega <- t(chol(diag(2)))
                             # ##
                             # Stan_data_list$prior_beta_mu_mean <- 0.0
                             # Stan_data_list$prior_beta_mu_SD   <- 1.0
                             # Stan_data_list$prior_beta_SD_mean <- 0.0
                             # Stan_data_list$prior_beta_SD_SD   <- 0.50
                             # ##
                             # Stan_data_list$prior_raw_scale_mu_mean <- 0.50
                             # Stan_data_list$prior_raw_scale_mu_SD   <- 1.00
                             # Stan_data_list$prior_raw_scale_SD_mean <- 0.0
                             # Stan_data_list$prior_raw_scale_SD_SD   <- 0.50
                             ##
                             file <- file.path(getwd(), "stan_models", "DTA_NMA_Nyaga_Gat_FIXEDthr.stan")
                           
             } else if (Model_type == "Cerullo_Nyaga_Gat_RANDOM_cutpoints") {

                             file <- file.path(getwd(), "stan_models", "DTA_NMA_Nyaga_Gat_RANDthr.stan")
                           
             }
           
             ##
             ## Compile Stan model:
             ##
             outs <- R_fn_compile_Stan_model(Stan_model_file_path = file, 
                                           CCACHE_PATH = "/usr/bin/ccache")
             ##
             ## Extract Stan model object:
             ##
             mod <- outs$mod

}










#### | ------  Run model - using Stan  ---------------------------------------------------------------------------------------------------------
{
          
          seed <- 123
          
          n_chains <- 4
          init_lists_per_chain <- rep(list(Stan_init_list), n_chains) 
          
          n_burnin <- 500
          n_iter   <- 500
          
          tictoc::tic()
          
          Stan_mod_sample <- mod$sample( seed = seed,
                                         data = Stan_data_list,
                                         init =   init_lists_per_chain, 
                                         chains = n_chains,
                                         parallel_chains = n_chains, 
                                         refresh = 10,
                                         iter_sampling = n_iter,
                                         iter_warmup = n_burnin,
                                         max_treedepth = 10, 
                                         adapt_delta = 0.80,
                                         metric = "diag_e")
          
          try({
                {
                  print(tictoc::toc(log = TRUE))
                  log.txt <- tictoc::tic.log(format = TRUE)
                  tictoc::tic.clearlog()
                  total_time <- unlist(log.txt)
                }
          })
          
          try({
                total_time <- as.numeric( substr(start = 0, stop = 100,  strsplit(  strsplit(total_time, "[:]")[[1]], "[s]")  [[1]][1] ) )
          })
          
}










{
        ## Predicted data:
        #### dev    <- Stan_mod_sample$summary(c("dev"))  ####    %>% print(n = 25)
        dev_nd <- Stan_mod_sample$summary(c("dev_nd"))  #### %>% print(n = 25)
        dev_d  <- Stan_mod_sample$summary(c("dev_d"))   #### %>% print(n = 25)
        dev_nd_mat <- round( array(-999999, dim = c(n_studies, n_thr)), 3)
        dev_d_mat  <- round( array(-999999, dim = c(n_studies, n_thr)), 3)
        
        
        dev_nd_medians <- ifelse(dev_nd$median == -1, 0, dev_nd$median)
        dev_nd_means   <- ifelse(dev_nd$mean == -1, 0, dev_nd$mean)
        dev_d_medians <- ifelse(dev_d$median == -1, 0, dev_d$median)
        dev_d_means   <- ifelse(dev_d$mean == -1, 0, dev_d$mean)
        
        # sum(dev_nd_medians, na.rm = TRUE)
        # sum(dev_nd_means,   na.rm = TRUE)
        # sum(dev_d_medians, na.rm = TRUE)
        # sum(dev_d_means,   na.rm = TRUE)
  
        x_hat_nd <- Stan_mod_sample$summary(c("x_hat_nd")) ##  %>% print(n = 10)
        x_hat_d  <- Stan_mod_sample$summary(c("x_hat_d")) ##   %>% print(n = 10)
        x_hat_nd_mat <- round( array(-1, dim = c(n_studies, n_thr)), 3)
        x_hat_d_mat  <- round( array(-1, dim = c(n_studies, n_thr)), 3)
        
        x_nd <-  Stan_data_list$x[[1]]
        x_d  <-  Stan_data_list$x[[2]]
        
        x_nd <- ifelse(x_nd == -1, 0, x_nd)
        x_d  <- ifelse(x_d ==  -1, 0, x_d)
        
        counter <- 0 
        for (k in 1:n_thr) {
          for (s in 1:n_studies) {
                  counter <- counter + 1
                  x_hat_nd_mat[s, k] <-  (x_hat_nd$median[counter])
                  x_hat_d_mat[s, k]  <-  (x_hat_d$median[counter])
                  dev_nd_mat[s, k]   <-  (dev_nd$median[counter])
                  dev_d_mat[s, k]    <-  (dev_d$median[counter])
          }
        }
        
        dev_nd_per_study <- rowSums(dev_nd_mat, na.rm = TRUE)
        dev_d_per_study  <- rowSums(dev_d_mat, na.rm = TRUE)
        
        x_hat_nd_mat <- ifelse(x_hat_nd_mat == -1, 0, x_hat_nd_mat)
        x_hat_d_mat  <- ifelse(x_hat_d_mat == -1, 0, x_hat_d_mat)
        
        abs_diff_mtx_nd <-  abs(x_hat_nd_mat - x_nd)
        abs_diff_mtx_d  <-  abs(x_hat_d_mat - x_d)
        
        ## Overall differences in total cell counts:
        message(paste("Overall % difference in total cell counts (D-) = "))
        print(round(100 * sum( abs_diff_mtx_nd ) / sum( x_hat_nd_mat[x_hat_nd_mat != -1] ), 2))
        message(paste("Overall % difference in total cell counts (D+) = "))
        print(round(100 * sum( abs_diff_mtx_d )  / sum( x_hat_d_mat[x_hat_d_mat != -1] ), 2))
        
        rowSums(abs_diff_mtx_nd)
        rowSums(x_hat_d_mat)
        
        ## Study-specific differences in total cell counts:
        ## message(paste("Study-specific % difference in total cell counts (D-) = "))
        ##  print( 100 * abs_diff_mtx_nd / x_hat_nd_mat )
        ## message(paste("Study-specific % difference in total cell counts (D+) = "))
        ## print( 100 * abs_diff_mtx_d  / x_hat_d_mat )
        
        ## Deviance:
        message(paste("Overall Deviance in non-diseased group (D-) = "))
        print( sum(dev_nd_medians, na.rm = TRUE) )
        message(paste("Overall Deviance in non-diseased group (D+) = "))
        print( sum(dev_d_medians,  na.rm = TRUE) )
        
        try({  
        {
            message(paste("Se (diffs) = "))
            ##
            est_Se_OVERALL <- Stan_mod_sample$summary(c("Se")) 
            est_Se_OVERALL_mean <- est_Se_OVERALL$mean*100
            est_Se_OVERALL_median <- est_Se_OVERALL$median*100
            ##
            print( abs(round(est_Se_OVERALL_median, 3) - round(true_Se_OVERALL_weighted, 3)) )
            message(paste("Se (SUM of abs. diffs) = "))
            print( sum(abs(round(est_Se_OVERALL_median, 3) - round(true_Se_OVERALL_weighted, 3))) )
            message(paste("Se (MEAN of abs. diffs) = "))
            print( mean(abs(round(est_Se_OVERALL_median, 3) - round(true_Se_OVERALL_weighted, 3))) )
            message(paste("Se (MEDIAN of abs. diffs) = "))
            print( median(abs(round(est_Se_OVERALL_median, 3) - round(true_Se_OVERALL_weighted, 3))) )
            message(paste("Se (MAX of abs. diffs) = "))
            print( max(abs(round(est_Se_OVERALL_median, 3) - round(true_Se_OVERALL_weighted, 3))) )
            
            message(paste("Sp (diffs) = "))
            ##
            est_Sp_OVERALL <- Stan_mod_sample$summary(c("Sp")) ; round(est_Sp_OVERALL$mean, 3)
            est_Sp_OVERALL_mean <- est_Sp_OVERALL$mean*100
            est_Sp_OVERALL_median <- est_Sp_OVERALL$median*100
            ##
            print( abs(round(est_Sp_OVERALL_median, 3) - round(true_Sp_OVERALL_weighted, 3)) )
            message(paste("Sp (SUM of abs. diffs) = "))
            print( sum(abs(round(est_Sp_OVERALL_median, 3) - round(true_Sp_OVERALL_weighted, 3))) )
            message(paste("Sp (MEAN of abs. diffs) = "))
            print( mean(abs(round(est_Sp_OVERALL_median, 3) - round(true_Sp_OVERALL_weighted, 3))) )
            message(paste("Sp (MEDIAN of abs. diffs) = "))
            print( median(abs(round(est_Sp_OVERALL_median, 3) - round(true_Sp_OVERALL_weighted, 3))) )
            message(paste("Sp (MAX of abs. diffs) = "))
            print( max(abs(round(est_Sp_OVERALL_median, 3) - round(true_Sp_OVERALL_weighted, 3))) )
        }
        })
        
        
        # ## Model_type <- 
        # ## Model_type <- "Cerullo_FIXED_cutpoints"
        # ## Model_type <- "Cerullo_RANDOM_HOMOG_cutpoints"
        # Model_type <- "Cerullo_RANDOM_cutpoints"

        
        if (Model_type == "Jones") { 
              colour <- "red"
              par(mfrow = c(2, 1))
              plot(  log(dev_nd_per_study), ylim = c(-8, 8), col = colour, pch = 19, cex = 2)    ;   abline(h = 0)
              points(log(dev_d_per_study),  ylim = c(-8, 8), col = colour, pch = 17, cex = 2)    ;   abline(h = 0)
        } else if (Model_type == "Cerullo_FIXED_cutpoints") { 
              colour <- "orange"
              points(log(dev_nd_per_study), ylim = c(-8, 8), col = colour, pch = 19, cex = 2)  ;   abline(h = 0)
              points(log(dev_d_per_study),  ylim = c(-8, 8), col = colour, pch = 17, cex = 2)  ;   abline(h = 0)
        } else if (Model_type == "Cerullo_RANDOM_HOMOG_cutpoints") { 
              colour <- "blue"
              points(log(dev_nd_per_study), ylim = c(-8, 8), col = colour, pch = 19, cex = 2)  ;   abline(h = 0)
              points(log(dev_d_per_study),  ylim = c(-8, 8), col = colour, pch = 17, cex = 2)  ;   abline(h = 0)
        } else if (Model_type == "Cerullo_RANDOM_cutpoints") { 
              colour <- "green"
              points(log(dev_nd_per_study), ylim = c(-8, 8), col = colour, pch = 19, cex = 2)  ;   abline(h = 0)
              points(log(dev_d_per_study),  ylim = c(-8, 8), col = colour, pch = 17, cex = 2)  ;   abline(h = 0)
        }
        
        message(paste("Model_type = ", Model_type))

        
       #  # dev.off()
       #  plot( x = 100.0 - est_Sp_OVERALL_mean, y = est_Se_OVERALL_mean, col = "blue")
       #  lines(x = 100.0 - true_Sp_OVERALL_weighted,     y = true_Se_OVERALL_weighted,     col = "green")
       #  # 
       # 
       # 
       #  true_Sp_OVERALL_weighted
       #  ##
       #  true_Se_OVERALL_weighted
       #  try({ 
       #      Stan_mod_sample$summary(c("Se_MU"))  %>% print(n = 100)
       #      Stan_mod_sample$summary(c("Se_MED"))  %>% print(n = 100)
       #  })
       #  try({ 
       #      Stan_mod_sample$summary(c("beta_mu"))  %>% print(n = 100)
       #      Stan_mod_sample$summary(c("beta_SD"))  %>% print(n = 100)
       #  })
       #  
       # ## Stan_mod_sample$summary(c("C"))  %>% print(n = 100)
       #  
       #  try({ 
       #    Stan_mod_sample$summary(c("C_mu"))  %>% print(n = 100)
       #  })
       #  try({ 
       #    Stan_mod_sample$summary(c("C_MU")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("C_MED")) %>% print(n = 100)
       #  })
       #  try({ 
       #    Se_MU <- Stan_mod_sample$summary(c("Se_MU")) ##%>% print(n = 100)
       #    Se_MU <- Se_MU$mean
       #    true_Se_OVERALL_weighted - Se_MU*100
       #    # ##
       #    # Se_MED <- Stan_mod_sample$summary(c("Se_MED")) ##%>% print(n = 100)
       #    # Se_MED <- Se_MED$mean
       #    # true_Se_OVERALL_weighted - Se_MED*100
       #    # ##
       #    Se_EMP <- Stan_mod_sample$summary(c("Se_EMP")) ##%>% print(n = 100)
       #    Se_EMP <- Se_EMP$mean
       #    true_Se_OVERALL_weighted - Se_EMP*100
       #    # ##
       #    # Se_SIM_MED <- Stan_mod_sample$summary(c("Se_SIM_MED")) ##%>% print(n = 100)
       #    # Se_SIM_MED <- Se_SIM_MED$mean
       #    # true_Se_OVERALL_weighted - Se_SIM_MED*100
       #    # ##
       #    # Se_SIM_MU <- Stan_mod_sample$summary(c("Se_SIM_MU")) ##%>% print(n = 100)
       #    # Se_SIM_MU <- Se_SIM_MU$mean
       #    # true_Se_OVERALL_weighted - Se_SIM_MU*100
       #    ##
       #    true_Se_OVERALL_weighted
       #    true_Se_OVERALL_weighted[4:7]
       #  })
       #  ##
       #  try({ 
       #    Stan_mod_sample$summary(c("Sp_MU")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("Sp_MED")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("Sp_EMP")) %>% print(n = 100)
       #    # Stan_mod_sample$summary(c("Sp_SIM")) %>% print(n = 100)
       #    # Stan_mod_sample$summary(c("Sp_SIM_MED")) %>% print(n = 100)
       #    # Stan_mod_sample$summary(c("Sp_SIM_MU")) %>% print(n = 100)
       #    true_Sp_OVERALL_weighted
       #    true_Sp_OVERALL_weighted[4:7]
       #  })
       #  ##
       #  # try({ 
       #  #   Stan_mod_sample$summary(c("C_MU")) %>% print(n = 100)
       #  #   Stan_mod_sample$summary(c("C_SD")) %>% print(n = 100)
       #  #   Stan_mod_sample$summary(c("C_empirical")) %>% print(n = 100)
       #  # })
       #  ##
       #  # try({
       #  #   Stan_mod_sample$summary(c("unc_C_normal_SD_sq"))  %>% print(n = 100)
       #    # Stan_mod_sample$summary(c("raw_MU"))  %>% print(n = 100)
       #    # Stan_mod_sample$summary(c("raw_SD"))  %>% print(n = 100)
       #    # Stan_mod_sample$summary(c("log_increment_SD"))  %>% print(n = 100)
       #    ##
       #  try({
       #    Stan_mod_sample$summary(c("category_means"))  %>% print(n = 100)
       #    Stan_mod_sample$summary(c("category_SDs"))  %>% print(n = 100)
       #    Stan_mod_sample$summary(c("kappa"))  %>% print(n = 100)
       #  })
       #  # })
       #  # try({ 
       #  #   Stan_mod_sample$summary(c("C_mu_empirical"))  %>% print(n = 100)
       #  # })
       #  # try({ 
       #  #   Stan_mod_sample$summary(c("C_mu_medians"))  %>% print(n = 100)
       #  # })
       #  # ##
       #  # # Stan_mod_sample$summary(c("prob_ord_mu"))  %>% print(n = 100)
       #  # # Stan_mod_sample$summary(c("prob_cumul_mu"))  %>% print(n = 100)
       #  # ##
       #  # try({  
       #  # Stan_mod_sample$summary(c("alpha")) %>% print(n = 100)
       #  # alpha <- Stan_mod_sample$summary(c("alpha"),  quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975)))  %>% print(n = 100)
       #  # # alpha_raw <- Stan_mod_sample$summary(c("alpha_raw"))  %>% print(n = 100)
       #  # # 100*c(alpha$mean)
       #  # })
       #  # ##
       #  # try({ 
       #  #   Stan_mod_sample$summary(c("kappa")) %>% print(n = 100)
       #  # })
       #  # try({ 
       #  #   Stan_mod_sample$summary(c("unc_C_normal_SD")) %>% print(n = 100)
       #  #   Stan_mod_sample$summary(c("unc_C_normal_MU")) %>% print(n = 100)
       #  #   Stan_mod_sample$summary(c("unc_C_normal_MED")) %>% print(n = 100)
       #  # })
       #  # 
       #  # print(mean(inf_dir_samples$alpha))
       #  # print(quantile(inf_dir_samples$alpha, c(0.025, 0.50, 0.975)))
       #  # 
       #  # 
       #  # Stan_mod_sample$summary(c("prob_cumul_mu"))  %>% print(n = 100)
       #  try({
       #    Stan_mod_sample$summary(c("location_MU_nd")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("location_MU_d")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("scale_MU_nd")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("scale_MU_d")) %>% print(n = 100)
       #  })
       #    ## 
       #  try({
       #    Stan_mod_sample$summary(c("expected_p_ord")) %>% print(n = 100)
       #    Stan_mod_sample$summary(c("dirichlet_cat_means_phi")) %>% print(n = 100)
       #  })
          ##
        try({
          Stan_mod_sample$summary(c("beta_mu")) %>% print(n = 100)
          Stan_mod_sample$summary(c("beta_SD")) %>% print(n = 100)
          Stan_mod_sample$summary(c("raw_scale_mu")) %>% print(n = 100)
          Stan_mod_sample$summary(c("raw_scale_SD")) %>% print(n = 100)
          ##
          Stan_mod_sample$summary(c("lambda")) %>% print(n = 100)
          
        })
        try({
        Stan_mod_sample$summary(c("Se")) %>% print(n = 100)
        Stan_mod_sample$summary(c("Sp")) %>% print(n = 100)
        })
        try({
          Se <- Stan_mod_sample$summary(c("Se"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          Se_pred <- Stan_mod_sample$summary(c("Se_pred"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          ##
          Sp <- Stan_mod_sample$summary(c("Sp"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          Sp_pred <- Stan_mod_sample$summary(c("Sp_pred"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          ##
          Fp <- Stan_mod_sample$summary(c("Fp"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          Fp_pred <- Stan_mod_sample$summary(c("Fp_pred"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
        })
          
        
        
}


df_true <- tibble(Se_true = true_Se_OVERALL_weighted/100, Sp_true = true_Sp_OVERALL_weighted/100, Fp_true = (100 - true_Sp_OVERALL_weighted)/100)

if (Model_type == "Cerullo_Gat_RANDOM_cutpoints") {
  
        ###
        df_Cerullo_Se <- tibble(Se_median = Se$`50%`, Se_lower = Se$`2.5%`, Se_upper = Se$`97.5%`, Se_pred_lower = Se_pred$`2.5%`, Se_pred_upper = Se_pred$`97.5%`)
        df_Cerullo_Sp <- tibble(Sp_median = Sp$`50%`, Sp_lower = Sp$`2.5%`, Sp_upper = Sp$`97.5%`, Sp_pred_lower = Sp_pred$`2.5%`, Sp_pred_upper = Sp_pred$`97.5%`)
        df_Cerullo_Fp <- tibble(Fp_median = Fp$`50%`, Fp_lower = Fp$`2.5%`, Fp_upper = Fp$`97.5%`, Fp_pred_lower = Fp_pred$`2.5%`, Fp_pred_upper = Fp_pred$`97.5%`)
        df_Cerullo <- tibble(cbind(df_Cerullo_Se, df_Cerullo_Sp, df_Cerullo_Fp)) %>%
          dplyr::mutate(Model = rep("Cerullo", n_thr))
        ##
        df_Cerullo <- tibble(cbind(df_true, df_Cerullo)) ; df_Cerullo
        
} else if (Model_type == "Jones") { 
  
        ###
        df_Jones_Se <- tibble(Se_median = Se$`50%`, Se_lower = Se$`2.5%`, Se_upper = Se$`97.5%`, Se_pred_lower = Se_pred$`2.5%`, Se_pred_upper = Se_pred$`97.5%`)
        df_Jones_Sp <- tibble(Sp_median = Sp$`50%`, Sp_lower = Sp$`2.5%`, Sp_upper = Sp$`97.5%`, Sp_pred_lower = Sp_pred$`2.5%`, Sp_pred_upper = Sp_pred$`97.5%`)
        df_Jones_Fp <- tibble(Fp_median = Fp$`50%`, Fp_lower = Fp$`2.5%`, Fp_upper = Fp$`97.5%`, Fp_pred_lower = Fp_pred$`2.5%`, Fp_pred_upper = Fp_pred$`97.5%`)
        df_Jones <- tibble(cbind(df_Jones_Se, df_Jones_Sp, df_Jones_Fp)) %>%
          dplyr::mutate(Model = rep("Jones", n_thr))
        ##
        df_Jones <- tibble(cbind(df_true, df_Jones)) ; df_Jones
        
}




# Improved function that handles column renaming correctly
create_confidence_polygon <- function(df, 
                                      model_name) {
  
        # Create upper curve (high sensitivity, low false positive)
        upper <- df %>% 
          filter(Model == model_name) %>%
          select(x = Fp_lower, y = Se_upper) %>%
          arrange(x)
        
        # Create lower curve (low sensitivity, high false positive)
        lower <- df %>% 
          filter(Model == model_name) %>%
          select(x = Fp_upper, y = Se_lower) %>%
          arrange(desc(x))
        
        # Combine into a closed polygon
        out <- tibble(bind_rows(upper, lower))
        n_rows <- nrow(out)
        out <- out %>%
          mutate(Model = rep(model_name, n_rows))
        
}




# Create prediction interval polygon function
create_prediction_polygon <- function(df, 
                                      model_name) {
  
        # Create upper curve (high sensitivity, low false positive)
        upper <- df %>% 
          filter(Model == model_name) %>%
          select(x = Fp_pred_lower, 
                 y = Se_pred_upper) %>% 
          arrange(x)
        
        # Create lower curve (low sensitivity, high false positive)
        lower <- df %>% 
          filter(Model == model_name) %>%
          select(x = Fp_pred_upper, 
                 y = Se_pred_lower) %>%
          arrange(desc(x))
        
        # Combine into a closed polygon
        out <- tibble(bind_rows(upper, lower))
        n_rows <- nrow(out)
        out <- out %>%
               mutate(Model = rep(model_name, n_rows))
  
}




Cerullo_polygon_Conf <- create_confidence_polygon(df = df_Cerullo, model_name = "Cerullo") ; Cerullo_polygon_Conf
Jones_polygon_Conf   <- create_confidence_polygon(df = df_Jones,   model_name = "Jones")   ; Jones_polygon_Conf
##

Cerullo_polygon_Pred <- create_prediction_polygon(df = df_Cerullo, model_name = "Cerullo") ; Cerullo_polygon_Pred
Jones_polygon_Pred   <- create_prediction_polygon(df = df_Jones,   model_name = "Jones")   ; Jones_polygon_Pred


df_all <- rbind(df_Jones, df_Cerullo)

df_Jones$Se_median - 
df_Cerullo$Se_median

## --------- Plot 1:
ggplot(data = df_all, 
       mapping = aes(x = Fp_median, y = Se_median, color = Model)) + 
  geom_line(size = 0.5) + 
  geom_point(size = 3) + 
  geom_point(color = "green", size = 3, data = df_all, mapping = aes(x = Fp_true, y = Se_true)) +
  geom_line(color = "green",  size = 0.5, data = df_all, mapping = aes(x = Fp_true, y = Se_true)) + 
  theme_bw(base_size = 16)



## --------- Plot 2:
ggplot(data = df_all, 
       mapping = aes(x = Fp_median, y = Se_median, color = Model)) + 
  # geom_line(size = 0.5) + 
  # geom_point(size = 3) + 
  # geom_point(color = "green", size = 3, data = df_all, mapping = aes(x = Fp_true, y = Se_true)) +
  # geom_line(color = "green",  size = 0.5, data = df_all, mapping = aes(x = Fp_true, y = Se_true)) + 
  theme_bw(base_size = 16) +
  geom_polygon(data = Cerullo_polygon_Conf, aes(x = x, y = y), fill = "red",  alpha = 0.75) +
  geom_polygon(data = Jones_polygon_Conf,   aes(x = x, y = y), fill = "cyan", alpha = 0.75) +
  geom_polygon(data = Cerullo_polygon_Pred, aes(x = x, y = y), fill = "red",  alpha = 0.25) + 
  geom_polygon(data = Jones_polygon_Pred,   aes(x = x, y = y), fill = "cyan", alpha = 0.25)
  

















x_individual_nd <- categorical_to_individual( Stan_data_list$x_cat[[1]], 
                                              binary_disease_indicator = 0, 
                                              missing_indicator = -1)

x_individual_d <- categorical_to_individual(  Stan_data_list$x_cat[[2]], 
                                              binary_disease_indicator = 1, 
                                              missing_indicator = -1)


x_individual_nd ## %>% print(n = 10000)
x_individual_d

x_individual <- tibble(rbind(x_individual_nd, x_individual_d))



require(MASS)
require(gridExtra) 


individual_obs_tibble <- x_individual_nd
group_name <- "non-diseased"
study_index <- 1


plot_distribution_fits(individual_obs_tibble = x_individual_nd,
                       group_name = "non-diseased")
                       



# Function to compare diseased and non-diseased distributions
compare_distributions <- function(x_cat_nd, x_cat_d, study_index = 1) {
  p1 <- plot_distribution_fits(x_cat_nd, study_index, "Non-diseased")
  p2 <- plot_distribution_fits(x_cat_d, study_index, "Diseased")
  
  grid.arrange(p1, p2, ncol = 1)
}

# Example usage for a specific study
# Assuming x_cat_nd is your non-diseased data and x_cat_d is your diseased data
# compare_distributions(x_cat[[1]], x_cat[[2]], 1)

# To examine multiple studies
examine_multiple_studies <- function(x_cat_nd, x_cat_d, study_indices = 1:3) {
  for (i in study_indices) {
    print(paste("Examining Study", i))
    print(compare_distributions(x_cat_nd, x_cat_d, i))
  }
}

# Example usage to examine first 3 studies
# examine_multiple_studies(x_cat[[1]], x_cat[[2]], 1:3)



































# Cerullo_polygon$Fp_lower
# Cerullo_polygon$Fp_upper
# 
# Cerullo_polygon$Se_lower
# Cerullo_polygon$Se_upper
# 
#   # geom_ribbon(data = df_all, mapping  = aes(ymin = Se_lower, ymax = Se_upper, fill = Model)) + 
#   # geom_ribbon(data = df_all, mapping  = aes(xmin = Fp_lower, xmax = Fp_upper, fill = Model))
#   # geom_ribbon(data = df_all, mapping  = aes(ymin = Se_pred_lower, ymax = Se_upper, fill = Model)) + 
#   # geom_ribbon(data = df_all, mapping  = aes(ymin = Fp_pred_lower, ymax = Fp_pred_upper, fill = Model))
# 
# 
#  geom_ribbon()
# 
# 
#   # geom_point(color = "red",   size = 4, data = df_all, mapping = aes(x = Fp_true, y = Se_true)) + 
#   # geom_line(color = "red",    size = 1, data = df_all, mapping = aes(x = Fp_true, y = Se_true)) + 
#   # theme_bw(base_size = 16)
#   
#   
# 
# 
# mu <- 0.50
# sigma <- 1.10
# 
# log_X <- rnorm(n = 100000, mean = mu, sd = sigma)
# 
# X <- exp(log_X)
# 
# plot(density(X), xlim  = c(0, 10))
#      
# 
# mu_X <- exp(mu + 0.5*sigma^2) ; mu_X ; mean(X)
# med_X <- exp(mu) ; med_X ; median(X)
# sd_X <- sqrt((exp(sigma^2) - 1) * exp(2*mu + sigma^2)) ; sd_X ; sd(X)
# 
# quantile_X_lower <- exp(mu + sigma*qnorm(0.16)) ; quantile_X_lower
# quantile_X_upper <- exp(mu + sigma*qnorm(1 - 0.16)) ; quantile_X_upper
# 
# exp(mu + sigma*qnorm(0.50))
# 
# sd_from_q_X <- (log(quantile_X_upper) - log(quantile_X_lower))/2
# 
# sd_from_q_X <- sqrt((quantile_X_upper-med_X)*(med_X-quantile_X_lower)) ; sd_from_q_X
# 
# 
# sd_X_from_quantile <- 
# 
# plot(density(rnorm(n = 10000, mean = med_X, sd = sd_X)))
# 
# true_Se_OVERALL[4:7]
# 
# 
# Stan_mod_sample$summary(c("prob_cumul_mu"))  %>% print(n = 100)
# 
# 
# softmax <- function(x) { 
#    return(log(1 + exp(x)))
# }
# softmax_inv <- function(x) { 
#    return(log(- 1 + exp(x)))
# }
# 
# log1p_exp <- softmax
# 
# 
# 
# median_softmax_X_from_X_normal <- function( unc_C_normal_MU, 
#                                             unc_C_normal_SD) {
#   
#         if (unc_C_normal_MU > 5.0) {
#           
#                     unc_C_normal_MED = unc_C_normal_MU
#               
#         } else if (unc_C_normal_MU < -5.0) { 
#           
#                     unc_C_normal_MED = exp(unc_C_normal_MU);
#               
#         } else { 
#        
#               
#                     unc_C_normal_SD_sq <- unc_C_normal_SD * unc_C_normal_SD;
#                     approx_mean <- log1p_exp(unc_C_normal_MU + 0.5 * unc_C_normal_SD_sq);
#                     
#                     lognormal_adjustment <- log(0.5 * unc_C_normal_SD_sq)
#                     
#                     unc_C_normal_MED <- approx_mean / exp(lognormal_adjustment);
#                     
#                     # base_median <- softmax(unc_C_normal_MU);
#                     # 
#                     # correction <- 0
#                     # if (unc_C_normal_MU >= 0) {
#                     #   
#                     #      ## // For μ ≥ 0, correction is smaller:
#                     #      correction <- 0.2 * exp(-0.5 * unc_C_normal_MU) * unc_C_normal_SD^2;
#                     #   
#                     # } else { 
#                     #   
#                     #      ## // For μ < 0, correction increases
#                     #      correction <- 0.3 * exp(0.2 * abs(unc_C_normal_MU)) * unc_C_normal_SD^2;
#                     #   
#                     # }
#                     # 
#                     # ## // Apply correction (add because softplus underestimates median for moderate μ)
#                     # unc_C_normal_MED <- base_median - 0.5 * correction;
#             
#         }
#   
#         return_list <- list(lognormal_adjustment = lognormal_adjustment,
#                             approx_mean = approx_mean,
#                             unc_C_normal_MED = unc_C_normal_MED)
#         
#         return(return_list)
#   
#   
# }
# 
# 
# 
# unc_C_normal_MU <- 1.5
# unc_C_normal_SD <- 2.5
# unc_C_normal <- rnorm(n = 100000, mean = unc_C_normal_MU, sd = unc_C_normal_SD)
# X <- unc_C_normal
# X_softmax <- softmax(unc_C_normal)
# 
# mean(X_softmax)
# median(X_softmax)
# 
# median_softmax_X_from_X_normal(unc_C_normal_MU = unc_C_normal_MU,
#                                unc_C_normal_SD = unc_C_normal_SD)
# 
# mean(X_softmax)
# median(X_softmax)
# 
# 
# 
# plot(density(X_softmax))
# 
# 
# 
# 
# 
# Y <- exp(unc_C_normal) ## Y is log-normal
# mean(Y)
# median(Y)
# 
# 
# mean_Y <- exp(unc_C_normal_MU + 0.5 * unc_C_normal_SD^2 ) ; mean_Y
# median_Y <- exp(unc_C_normal_MU) ; median_Y
# 
# plot(density(Y))
# 
# plot(density(log(1 + exp(unc_C_normal))))
# 
# Y <- 1 + exp(unc_C_normal)
# mean(Y)
# median(Y)
# Y_m_one <- Y - 1 ## this is log-normal
# mean(Y_m_one)
# median(Y_m_one)
# 
# exp(unc_C_normal_MU) ; median(Y_m_one)
# exp(unc_C_normal_MU + 0.5*unc_C_normal_SD^2) ; mean(Y_m_one)
# 
# 
# 1 + exp(unc_C_normal_MU) ; median(Y)
# 1 + exp(unc_C_normal_MU + 0.5*unc_C_normal_SD^2) ; mean(Y)
# 
# log(1 + exp(unc_C_normal_MU)) ; median(log(1 + exp(unc_C_normal)))
# 
# log(1 + exp(unc_C_normal_MU + 0.5*unc_C_normal_SD^2)) ;   mean(log(1 + exp(unc_C_normal)))
# 
# mean(unc_C_normal)
# median(unc_C_normal)
# 
# 
# 
# mean_Ym1 <- exp(unc_C_normal_MU + 0.5 * unc_C_normal_SD^2 ) ; mean_Ym1
# median_Ym1 <- exp(unc_C_normal_MU) ; median_Ym1
# 
# mean_Y <- 1 + mean_Ym1 ; mean_Y ; mean(Y)
# median_Y <- 1 + median_Ym1 ; median_Y ; median(Y)
# 
# log1p_exp_X <- log1p_exp(X)
# median(log1p_exp_X)
# mean(log1p_exp_X)
# plot(density(log1p_exp_X))
# ##
# plot(density(exp(X)), xlim = c(0, 10000))
# 
# 
# 
# soft_C_normal_MED <- log1p_exp(unc_C_normal_MU) ; soft_C_normal_MED
# median(log1p_exp_X)
# mean(log1p_exp_X)
# median(X)
# 
# 
# 
# log(median_Y)
# log(mean_Y)
# 
# 
# 
# 
#   correction <- -0.5 * unc_C_normal_SD^2;  correction  ##  // Approximate correction
#   
# log(exp(median(X_softmax)) - 1) - correction
# exp(log(exp(median(X_softmax)) - 1) - correction)
#   
# approx_MU <- log(exp(soft_C_normal_MED) - 1) + correction; approx_MU
# 
# 
# 
# 
# 
# # 
# # Stan_mod_sample$summary(c("se"))  %>% print(n = 100)
# # 
# # 
# # 
# # 
# # alpha_array <- array(dim = c(2, n_cat))
# # counter <- 1
# # for (k in 1:n_cat) { 
# # for (c in 1:2) {
# #     alpha_array[c, k] <- alpha$`50%`[counter]
# #     counter <- counter + 1
# #   }
# # }
# #  
# # 
# # p_ord_mu_nd <- alpha_array[1, ] / sum(alpha_array[1,])
# # p_ord_mu_d  <- alpha_array[2, ] / sum(alpha_array[2,])
# # 
# # 
# # 
# # 
# # est_Sp_OVERALL_mean
# # true_Sp_OVERALL
# # 
# # 100 - est_Se_OVERALL_mean
# # true_Se_OVERALL
# # 
# # 
# # 100 - est_Se_OVERALL_mean
# # true_Se_OVERALL
# # 
# # 
# # Stan_mod_sample$summary(c("lambda"))  %>% print(n = 100)
# # 
# # 
# # 
# 
# 
# 
# 
# # 10*c(alpha_raw$mean)
# # 
# # mean(c(alpha_raw$mean))
# # 1/11
# # 
# # 
# # ## Summary estimates:
# # Stan_mod_sample$summary(c("Se"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("Sp"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("se"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("sp"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("beta_mu"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("beta_SD"))  %>% print(n = 100)
# # # ##
# # Stan_mod_sample$summary(c("log_scale_mu"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("log_scale_SD"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("scale"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("alpha"))  %>% print(n = 100)
# # ## Stan_mod_sample$summary(c("log_alpha"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("phi"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("kappa"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("C_nd"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("C_d"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("C_mu"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("C_mu_empirical"))  %>% print(n = 100)
# # ##
# # Stan_mod_sample$summary(c("lambda"))  %>% print(n = 100)
# # ##
# # # Stan_mod_sample$summary(c("phi_d"))  %>% print(n = 100)
# # # Stan_mod_sample$summary(c("phi_nd"))  %>% print(n = 100)
# # 
# # # ## Study-specific estimates:
# # # Stan_mod_sample$summary(c("se"))  %>% print(n = 100)
# # # Stan_mod_sample$summary(c("sp"))  %>% print(n = 100)
# # 
# # Stan_mod_sample$summary(c("cumul_prob"))  %>% print(n = 100)
# # 
# # 
# # 
# # 
# # Stan_mod_sample$summary(c("alpha"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("phi"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("prob_ord_mu"))  %>% print(n = 100)
# # ##
# # df_alpha <- Stan_mod_sample$summary(c("alpha"))  %>% print(n = 100)
# # alpha <- df_alpha$mean
# # alpha_nd <- alpha_d <- c()
# # ##
# # counter <- 0 
# # for (k in 1:n_thr) {
# #   for (c in 1:2) {
# #       counter <- counter + 1
# #       if (c == 1) alpha_nd[k] <- alpha[counter]
# #       if (c == 2) alpha_d[k] <- alpha[counter + 1]
# #   }
# # }
# #   
# #   
# # prob_ord_mu_sim_nd   <- array(NA, dim = c(n_thr, 1000))
# # prob_ord_mu_sim_d    <- array(NA, dim = c(n_thr, 1000))
# # C_mu_sim_nd   <- array(NA, dim = c(n_thr, 1000))
# # C_mu_sim_d    <- array(NA, dim = c(n_thr, 1000))
# # anchor_nd <- 0.0
# # anchor_d  <- 0.0
# # 
# # for (i in 1:1000) {
# #   
# #       prob_ord_mu_sim_nd[,i] <-  c(MCMCpack::rdirichlet(n = 1, alpha_nd))
# #       prob_ord_mu_sim_d[,i]  <-  c(MCMCpack::rdirichlet(n = 1, alpha_d))
# #       
# #       C_mu_sim_nd[1, i] =   anchor_nd - qlogis( plogis(anchor_nd - (-Inf)) -  prob_ord_mu_sim_nd[1 ,i]  )
# #       C_mu_sim_d[1, i]  =   anchor_d  - qlogis( plogis(anchor_d  - (-Inf)) -  prob_ord_mu_sim_d[1, i]  )
# #       for (k in 2:n_thr) {
# #         C_mu_sim_nd[k, i] =  anchor_nd - qlogis( plogis(anchor_nd - C_mu_sim_nd[k - 1, i]) - prob_ord_mu_sim_nd[k, i] );
# #         C_mu_sim_d[k, i]  =  anchor_d  - qlogis( plogis(anchor_d  - C_mu_sim_d[k - 1, i])  - prob_ord_mu_sim_d[k, i] );
# #       }
# #       
# #   
# # }
# # 
# # C_mu_nd <- rowMeans(C_mu_sim_nd)
# # C_mu_d  <- rowMeans(C_mu_sim_d)
# # 
# # 
# # plogis( beta_mu_nd -  C_mu_nd[1:9] )
# # 100.0 - true_Sp_OVERALL
# # 
# # plogis( beta_mu_d  -  C_mu_d[1:9] )
# # true_Se_OVERALL
# # 
# # 
# # ##
# # C_mu <- Stan_mod_sample$summary(c("C_mu"))  %>% print(n = 100)
# # C_mu_means <- C_mu$mean
# # C_mu <- array(dim = c(2, n_thr))
# # counter <- 0 
# # for (k in 1:n_thr) {
# #   for (c in 1:2) {
# #      counter <- counter + 1
# #      C_mu[c, k] <- C_mu_means[counter]
# #   }
# # }
# # C_mu
# # C_mu <- reorder_decreasing(C_mu)
# # ##
# # C_nd_inc <- Stan_mod_sample$summary(c("C_nd_inc"))  %>% print(n = 100)
# # C_d_inc <- Stan_mod_sample$summary(c("C_d_inc"))  %>% print(n = 100)
# # ##
# # beta_mu <- Stan_mod_sample$summary(c("beta_mu"))  %>% print(n = 100)
# # ##
# # plogis( beta_mu_nd - C_nd_inc )
# # 100.0 - true_Sp_OVERALL
# # ##
# #  plogis( - beta_mu_d + C_d_inc$mean)
# # true_Se_OVERALL
# # ##
# # beta_mu <- Stan_mod_sample$summary(c("beta_mu"))  %>% print(n = 100)
# # beta_mu_nd <- beta_mu$mean[1]
# # beta_mu_d  <- beta_mu$mean[2]
# # ##
# # plogis( beta_mu_nd - C_mu[1, ] )
# # 100.0 - true_Sp_OVERALL
# # ##
# # 1.0 - plogis( C_mu[2, ] - beta_mu_d)
# # true_Se_OVERALL
# # ##
# # Stan_mod_sample$summary(c("cutpoints_nd"))  %>% print(n = 100)
# # Stan_mod_sample$summary(c("cutpoints_d"))  %>% print(n = 100)
# # 
# # 
# # 
# # prob_cumul_mu_nd*100
# # true_Sp_OVERALL
# # 
# # 
# # prob_cumul_mu_d*100
# # true_Se_OVERALL
# # 
# # 
# # 
# # Stan_mod_sample$summary(c("beta"))  %>% print(n = 100)
# # 
# # 
# # 
# #  
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# 
# 
# 
# 
# ## prior for softplus raw_scale:
# samps <- rnorm(n = 10000, mu = 0.5, sd = 1.0)
# soft_samps <- log(1 + exp(samps))
# 
# round(quantile(soft_samps, probs = c(0.025, 0.50, 0.975)), 2)
#   
# 
# 
# 
# 
# 
# 
# 
# 
# 
