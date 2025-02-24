


 


setwd("/home/enzocerullo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")

source("R_fn_load_data_ordinal_MA_LC_MVP_sim.R")
source("missing_thr_prep_Stan_data.R")
source("missing_thr_prior_pred_check.R")

#### To overwrite some model options:
{
    ##
    ## Model options for JONES model:
    ##
    Stan_data_list$use_box_cox <- 1
    ##  
    ## Model options for "Cerullo" model:
    ##
    Stan_data_list$estimate_scales <- 0
    Stan_data_list$same_cutpoints_between_groups <- 0
    ##
    Stan_data_list$prior_only <- 0
    ##
    Stan_data_list$prior_kappa_mean <-      rep(0.0, 2)
    Stan_data_list$prior_kappa_SD <-        rep(1.0, 2)
}

 


induced_Dirichlet_ppc_plot(use_alpha_directly = TRUE,
                           prior_sd = Stan_data_list$prior_kappa_SD[1],
                           use_log_kappa = FALSE,
                           n_cat = 12,
                           N = 1000)





##  | ------   Initial values  -------------------------------------------------------------------------
{
  Stan_init_list <- list()
  ##
  Stan_init_list$cutpoints_nd <- t(array(dim = c(n_thr, n_studies), data = seq(from = -2.0, to = 2.0, length = n_thr)))
  Stan_init_list$cutpoints_d <-  t(array(dim = c(n_thr, n_studies), data = seq(from = -2.0, to = 2.0, length = n_thr)))
  ##
  ##
  Stan_init_list$beta_mu <- rep(0.0, 2)
  Stan_init_list$beta_SD <- rep(0.01, 2)
  Stan_init_list$beta_z <-  array(0.0, dim = c(2, n_studies))
  ##
  # Stan_init_list$log_scale_d_mu <- 0.0
  # Stan_init_list$log_scale_d_SD <- 0.01
  # Stan_init_list$log_scale_d_z <- rep(0.0, n_studies)
  ##
  Stan_init_list$kappa <- rep(250, 2)
  Stan_init_list$phi <-  array(dim = c(2, n_cat), data = rep(1/n_cat, n_cat))
  ##
  Stan_init_list$log_scale_mu <- rep(0.0, 2)
  Stan_init_list$log_scale_SD <- rep(0.01, 2)
  Stan_init_list$log_scale_z <-  array(0.0, dim = c(2, n_studies))
  ##
  Stan_init_list$log_scale_d_mu <-  Stan_init_list$log_scale_mu[2]
  Stan_init_list$log_scale_d_SD <-  Stan_init_list$log_scale_SD[2]
  Stan_init_list$log_scale_d_z <-   Stan_init_list$log_scale_z[2, ]
  # ##
  # Stan_init_list$log_alpha <- list(rep(0, n_cat),
  #                                  rep(0, n_cat))
  
  
}





## Model_type <- "Jones"
## Model_type <- "Cerullo_FIXED_cutpoints"
Model_type <- "Cerullo_RANDOM_cutpoints"

## - | ----------  Compile Stan model --------------------------------------------------------------------
{   
  
         if (Model_type == "Jones") {
              file <- file.path(getwd(), "stan_models", "DTA_MA_cts_JONES_BOXCOX_NONsym_sROC.stan")
         } else if (Model_type == "Cerullo_FIXED_cutpoints") { 
              Stan_init_list$cutpoints_nd <-  c(array(dim = c(n_thr, 1), data = seq(from = -2.0, to = 2.0, length = n_thr)))
              Stan_init_list$cutpoints_d <-   c(array(dim = c(n_thr, 1), data = seq(from = -2.0, to = 2.0, length = n_thr)))
              file <- file.path(getwd(), "stan_models", "DTA_MA_ord_FIXEDthr_w_NONsym_sROC.stan")
         } else if (Model_type == "Cerullo_RANDOM_cutpoints") { 
              Stan_init_list$cutpoints_nd <-  t(array(dim = c(n_thr, n_studies), data = seq(from = -2.0, to = 2.0, length = n_thr)))
              Stan_init_list$cutpoints_d <-   t(array(dim = c(n_thr, n_studies), data = seq(from = -2.0, to = 2.0, length = n_thr)))
              file <- file.path(getwd(), "stan_models", "DTA_MA_ord_RANDthr_w_NONsym_sROC.stan")
         }
        
        
          mod <- cmdstan_model(file, 
                             ## force_recompile = TRUE,
                             quiet = FALSE,
                             ## user_header = path_to_cpp_user_header,
                             ## cpp_options = cmdstan_cpp_flags
          )
        
}




## | ------  Run model - using Stan  -------------------------------------------------------------------------
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
                                         refresh = 50,
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




## Summary estimates:
Stan_mod_sample$summary(c("Se"))  %>% print(n = 100)
Stan_mod_sample$summary(c("Sp"))  %>% print(n = 100)
##
Stan_mod_sample$summary(c("se"))  %>% print(n = 100)
Stan_mod_sample$summary(c("sp"))  %>% print(n = 100)
##
Stan_mod_sample$summary(c("beta_mu"))  %>% print(n = 100)
Stan_mod_sample$summary(c("beta_SD"))  %>% print(n = 100)
##
Stan_mod_sample$summary(c("log_scale_mu"))  %>% print(n = 100)
Stan_mod_sample$summary(c("log_scale_SD"))  %>% print(n = 100)
##
Stan_mod_sample$summary(c("cutpoints_nd"))  %>% print(n = 100)
Stan_mod_sample$summary(c("cutpoints_d"))  %>% print(n = 100)
##
Stan_mod_sample$summary(c("phi"))  %>% print(n = 100)
Stan_mod_sample$summary(c("log_kappa"))  %>% print(n = 100)
##
Stan_mod_sample$summary(c("lambda"))  %>% print(n = 100)
##
# Stan_mod_sample$summary(c("phi_d"))  %>% print(n = 100)
# Stan_mod_sample$summary(c("phi_nd"))  %>% print(n = 100)

# ## Study-specific estimates:
# Stan_mod_sample$summary(c("se"))  %>% print(n = 100)
# Stan_mod_sample$summary(c("sp"))  %>% print(n = 100)

Stan_mod_sample$summary(c("cumul_prob"))  %>% print(n = 100)


stan_draws <- Stan_mod_sample$draws()
str(stan_draws)

plot(stan_draws[,,"beta_SD[2]"])






{
        ## Predicted data:
        #### dev    <- Stan_mod_sample$summary(c("dev"))      %>% print(n = 25)
        dev_nd <- Stan_mod_sample$summary(c("dev_nd"))   %>% print(n = 25)
        dev_d  <- Stan_mod_sample$summary(c("dev_d"))    %>% print(n = 25)
        
        
        dev_nd_medians <- ifelse(dev_nd$median == -1, 0, dev_nd$median)
        dev_nd_means   <- ifelse(dev_nd$mean == -1, 0, dev_nd$mean)
        
        dev_d_medians <- ifelse(dev_d$median == -1, 0, dev_d$median)
        dev_d_means   <- ifelse(dev_d$mean == -1, 0, dev_d$mean)
        
        
        sum(dev_nd_medians, na.rm = TRUE)
        sum(dev_nd_means,   na.rm = TRUE)
        
        sum(dev_d_medians, na.rm = TRUE)
        sum(dev_d_means,   na.rm = TRUE)
  
        x_hat_nd <- Stan_mod_sample$summary(c("x_hat_nd")) ##  %>% print(n = 10)
        x_hat_d  <- Stan_mod_sample$summary(c("x_hat_d")) ##   %>% print(n = 10)
    
        x_hat_nd_mat <- round( array(-1, dim = c(n_studies, n_thr)), 0)
        x_hat_d_mat  <- round( array(-1, dim = c(n_studies, n_thr)), 0)
        
        x_nd <-  Stan_data_list$x[[1]]
        x_d  <-  Stan_data_list$x[[2]]
        
        x_nd <- ifelse(x_nd == -1, 0, x_nd)
        x_d  <- ifelse(x_d ==  -1, 0, x_d)
        
        counter <- 0 
        for (k in 1:n_thr) {
          for (s in 1:n_studies) {
                  counter <- counter + 1
                  x_hat_nd_mat[s, k] <- floor(x_hat_nd$median[counter])
                  x_hat_d_mat[s, k]  <- floor(x_hat_d$median[counter])
          }
        }
        
        x_hat_nd_mat <- ifelse(x_hat_nd_mat == -1, 0, x_hat_nd_mat)
        x_hat_d_mat  <- ifelse(x_hat_d_mat == -1, 0, x_hat_d_mat)
        
        abs_diff_mtx_nd <-  abs(x_hat_nd_mat - x_nd)
        abs_diff_mtx_d  <-  abs(x_hat_d_mat - x_d)
        
        ## Overall differences in total cell counts:
        message(paste("Overall % difference in total cell counts (D-) = "))
        print(round(100 * sum( abs_diff_mtx_nd ) / sum( x_hat_nd_mat[x_hat_nd_mat != -1] ), 2))
        message(paste("Overall % difference in total cell counts (D+) = "))
        print(round(100 * sum( abs_diff_mtx_d )  / sum( x_hat_d_mat[x_hat_d_mat != -1] ), 2))
        
        ## Study-specific differences in total cell counts:
        message(paste("Study-specific % difference in total cell counts (D-) = "))
        print( 100 * abs_diff_mtx_nd / x_hat_nd_mat )
        message(paste("Study-specific % difference in total cell counts (D+) = "))
        print( 100 * abs_diff_mtx_d  / x_hat_d_mat )
        
        ## Deviance:
        message(paste("Overall Deviance in non-diseased group (D-) = "))
        print( sum(dev_nd_medians, na.rm = TRUE) )
        message(paste("Overall Deviance in non-diseased group (D+) = "))
        print( sum(dev_d_medians,  na.rm = TRUE) )
        
        ## Model fit / evaluation:
        true_Se_OVERALL <- sim_results$Se_OVERALL_all_tests_all_thresholds[5, ] ; true_Se_OVERALL
        true_Sp_OVERALL <- sim_results$Sp_OVERALL_all_tests_all_thresholds[5, ] ; true_Sp_OVERALL
        
        message(paste("Se (diffs) = "))
        round(true_Se_OVERALL, 3)
        est_Se_OVERALL <- Stan_mod_sample$summary(c("Se")) ; round(est_Se_OVERALL$mean, 3)
        print( abs(round(est_Se_OVERALL$mean, 3) - round(true_Se_OVERALL, 3)) )
        print( sum(abs(round(est_Se_OVERALL$mean, 3) - round(true_Se_OVERALL, 3))) )
        print( max(abs(round(est_Se_OVERALL$mean, 3) - round(true_Se_OVERALL, 3))) )
        
        message(paste("Sp (diffs) = "))
        round(true_Sp_OVERALL, 3)
        est_Sp_OVERALL <- Stan_mod_sample$summary(c("Sp")) ; round(est_Sp_OVERALL$mean, 3)
        print( abs(round(est_Sp_OVERALL$mean, 3) - round(true_Sp_OVERALL, 3)) )
        print( sum(abs(round(est_Sp_OVERALL$mean, 3) - round(true_Sp_OVERALL, 3))) )
        print( max(abs(round(est_Sp_OVERALL$mean, 3) - round(true_Sp_OVERALL, 3))) )
        
}


# abs_diffs <- abs( round(x_hat_nd_mat, 3) - Stan_data_list$x[[1]] )
# abs_diffs_as_pct <- 100*abs_diffs/Stan_data_list$x[[1]]
# 
# 
# abs_diffs_as_pct <- 100 * ( abs_diffs / (Stan_data_list$x[[1]]) )
# rowSums(abs_diffs_as_pct, na.rm = TRUE)

## Model fit / evaluation:
true_Se_OVERALL <- sim_results$Se_OVERALL_all_tests_all_thresholds[5, ]
true_Sp_OVERALL <- sim_results$Sp_OVERALL_all_tests_all_thresholds[5, ]
 

round(true_Se_OVERALL, 3)
est_Se_OVERALL <- Stan_mod_sample$summary(c("Se")) ; round(est_Se_OVERALL$mean, 3)
abs(round(est_Se_OVERALL$mean, 3) - round(true_Se_OVERALL, 3))
sum(abs(round(est_Se_OVERALL$mean, 3) - round(true_Se_OVERALL, 3)))
max(abs(round(est_Se_OVERALL$mean, 3) - round(true_Se_OVERALL, 3)))

round(true_Sp_OVERALL, 3)
est_Sp_OVERALL <- Stan_mod_sample$summary(c("Sp")) ; round(est_Sp_OVERALL$mean, 3)
abs(round(est_Sp_OVERALL$mean, 3) - round(true_Sp_OVERALL, 3))
sum(abs(round(est_Sp_OVERALL$mean, 3) - round(true_Sp_OVERALL, 3)))
max(abs(round(est_Sp_OVERALL$mean, 3) - round(true_Sp_OVERALL, 3)))


Stan_mod_sample$loo()

?loo::loo
loo::loo(x = log_lik)


# x_non_diseased_categorical <- convert_cumulative_to_category(final_data_cumulative$x_non_diseased)
# 
# sum(x_diseased_categorical == 999) / sum(prod(dim(x_diseased_categorical)))
# 
# 
# # 
# 
# 
# 
# 
# y_diseased_df <- categorical_to_individual(  categorical_matrix = x_diseased_categorical, 
#                                              binary_disease_indicator = 1)
# 
# y_non_diseased_df <- categorical_to_individual(  categorical_matrix = x_non_diseased_categorical, 
#                                                  binary_disease_indicator = 0)
# 
# y_df <- tibble(rbind(y_non_diseased_df, y_diseased_df))
# print(y_df)
# 
# 






















Stan_mod_sample$summary(c("beta_SD"))  %>% print(n = 100)