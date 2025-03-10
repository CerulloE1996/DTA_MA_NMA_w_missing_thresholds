


setwd("/home/enzocerullo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")

source("R_fn_load_data_ordinal_MA_LC_MVP_sim.R")

require(dplyr)
require(cmdstanr)

n_studies <- 10
N_per_study_mean <- 5000
N_per_study_SD <- 500
assume_perfect_GS <- 1
seed <- 123

missing_indicator <- -1

options(max.print = 100000)


# Run simulated data - this simulates data from FIVE (5) diagnostic tests (1 BINARY reference test + 4 ORDINAL index tests)
sim_results <- simulate_binary_and_ordinal_MA_LC_MVP_data(n_studies = n_studies,
                                                          N_per_study_mean = N_per_study_mean,
                                                          N_per_study_SD = N_per_study_SD,
                                                          assume_perfect_GS = assume_perfect_GS,
                                                          seed = seed)

y_list <- sim_results$y_list
sim_results$n_total_possible_thresholds_per_test
str(y_list)

true_Se_OVERALL <- sim_results$Se_OVERALL_all_tests_all_thresholds[5, ] ; true_Se_OVERALL
true_Sp_OVERALL <- sim_results$Sp_OVERALL_all_tests_all_thresholds[5, ] ; true_Sp_OVERALL
true_Fp_OVERALL <- 1.0 - true_Sp_OVERALL ; true_Fp_OVERALL

plot(y = true_Se_OVERALL, x = 1 - true_Sp_OVERALL)

sim_results$Se_per_study_all_tests_all_thresholds_list

## Now for the first example we will only take one index test, as initially we are evaluating just a "simple" model
## where there's only a single index test with 12 thresholds (so 13 categories:
n_thr <- 12
n_cat <- n_thr + 1
y_list_example_1 <- list()
for (s in 1:n_studies) { 
    N <- nrow(y_list[[s]])
    y <- array(NA, dim = c(N, 2))
    y[, 1] <- y_list[[s]][, 1]
    y[, 2] <- y_list[[s]][, 5] ## take test #5 as the index test for this example 
    y_list_example_1[[s]] <- y
}

y_list_example_1



# Convert to AGGREGATE data as these more basic NMA/MA models don't use individual-level data:
agg_data_cumulative <- convert_to_aggregate_counts( y_list_example_1, 
                                                    n_studies,
                                                    n_thr)

agg_data_cumulative





## Apply missing thresholds:
agg_data_cumulative_with_missing_thresholds <- agg_data_cumulative



{
      
    ## First 2 studies only report data for the "middle 8" thresholds (starting from threshold 3, so thr = {3, 4, 5, 6, 7, 8, 9, 10})
    studies_subset_vec <- c(1, 2)
    missing_thr_subset_vec <- c(1:2, 11:12)
    agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
                                                                            studies_subset_vec = studies_subset_vec,
                                                                            missing_thr_subset_vec = missing_thr_subset_vec, 
                                                                            missing_indicator = missing_indicator)
    ## The next 2 studies only report at the "middle 4" thresholds (starting at threshold 5): so thr = {5, 6, 7, 8}:
    studies_subset_vec <- c(3, 4)
    missing_thr_subset_vec <- c(1:4, 9:12)
    agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
                                                                            studies_subset_vec = studies_subset_vec,
                                                                            missing_thr_subset_vec = missing_thr_subset_vec,
                                                                            missing_indicator = missing_indicator)
    ## The next two studies only report data at only a single threshold - threshold #7:
    studies_subset_vec <- c(5, 6)
    missing_thr_subset_vec <- c(1:6, 8:12)
    agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
                                                                            studies_subset_vec = studies_subset_vec,
                                                                            missing_thr_subset_vec = missing_thr_subset_vec,
                                                                            missing_indicator = missing_indicator)
    ## The next two studies only report data at 4 non-adjacent thresholds - so thr = {3, 5, 7, 9}:
    studies_subset_vec <- c(7, 8)
    missing_thr_subset_vec <- c(1, 2, 4, 6, 8, 10, 11, 12)
    agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
                                                                            studies_subset_vec = studies_subset_vec,
                                                                            missing_thr_subset_vec = missing_thr_subset_vec,
                                                                            missing_indicator = missing_indicator)
    ## And finally, the last 2 studies report data at ALL thresholds:
    #  -----  Do nothing here - as no missing thresholds for study #10!

}

## Now let's look at the overall % of missing thresholds:
total_missing_thr <- sum(agg_data_cumulative_with_missing_thresholds$x_diseased == 0.999) ; total_missing_thr
prop_missing_thr <- total_missing_thr / (n_studies * n_thr) ; prop_missing_thr
## Just under half the threshold data is missing (46.67%). 

 
 
 
# ## Diseased group:
x_diseased_cumulative <- agg_data_cumulative_with_missing_thresholds$x_diseased
print(paste("x_diseased_cumulative = ")) ; print(x_diseased_cumulative)
## Non-diseased group:
x_non_diseased_cumulative <- agg_data_cumulative_with_missing_thresholds$x_non_diseased
print(paste("x_non_diseased_cumulative = ")) ; print(x_non_diseased_cumulative)





##  | ------   Initial values  -------------------------------------------------------------------------
{
        Stan_init_list <- list()
        ##
        Stan_init_list$cutpoints_nd <- t(array(dim = c(n_thr, 1), data = seq(from = -2.0, to = 2.0, length = n_thr)))
        Stan_init_list$cutpoints_d <-  t(array(dim = c(n_thr, 1), data = seq(from = -2.0, to = 2.0, length = n_thr)))
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
        Stan_init_list$kappa <- rep(10.0, 2)
        Stan_init_list$phi <-  array(dim = c(2, n_cat), data = rep(1/n_cat, n_cat))
        ##
        Stan_init_list$log_scale_mu <- rep(0.0, 2)
        Stan_init_list$log_scale_SD <- rep(0.01, 2)
        Stan_init_list$log_scale_z <-  array(0.0, dim = c(2, n_studies))
}





##  | ------   Stan data   -------------------------------------------------------------------------
{
        Stan_data_list <- list()
        ##
        Stan_data_list$n_studies <-  n_studies
        Stan_data_list$n_thr <-  n_thr
        ##
        Stan_data_list$n_non_diseased <- agg_data_cumulative_with_missing_thresholds$n_non_diseased
        Stan_data_list$n_diseased <-     agg_data_cumulative_with_missing_thresholds$n_diseased
        ##
        Stan_data_list$x_non_diseased <- agg_data_cumulative_with_missing_thresholds$x_non_diseased 
        Stan_data_list$x_diseased <-     agg_data_cumulative_with_missing_thresholds$x_diseased
        ##
        x_with_missings <- list(Stan_data_list$x_non_diseased,
                                Stan_data_list$x_diseased)
        Stan_data_list$x_with_missings <- x_with_missings
        ##
        cutpoint_index <- n <- x <-  list()
        ##
        n_cutpoints <- matrix(nrow = 2, ncol = n_studies)
        ##
        for (c in 1:2) {
            cutpoint_index[[c]] = matrix(-1, nrow = n_studies, ncol = n_thr);
            n[[c]] = matrix(-1, nrow = n_studies, ncol = n_thr);
            x[[c]] = matrix(-1, nrow = n_studies, ncol = n_thr);
        }
        ##
        for (s in 1:n_studies) {
          for (c in 1:2) {
            cutpoint_counter = 0;
            previous_non_missing_index = 1;
            for (k in 1:n_cat) {
              
                  if (k == n_cat) {
                    cond <- TRUE
                  } else { 
                    cond <- (x_with_missings[[c]][s, k] != -1)
                  }
              
                  if (cond == TRUE)   {
                    
                       try({  
                            cutpoint_counter = cutpoint_counter +  1;
                            if (k != n_cat) {
                                cutpoint_index[[c]][s, cutpoint_counter] = k;
                            }
                       })
                        
                        if (cutpoint_counter == 1) {  ## //// First non-missing cutpoint
                              
                              if (c == 1) n[[1]][s, 1] = Stan_data_list$n_non_diseased[s];
                              if (c == 2) n[[2]][s, 1] = Stan_data_list$n_diseased[s];
                          
                        } else if (cutpoint_counter > 1) {
                          
                            try({ 
                              n[[c]][s, cutpoint_counter]         =  x_with_missings[[c]][s, previous_non_missing_index];
                            })
                              x[[c]][s, cutpoint_counter - 1]     =  x_with_missings[[c]][s, previous_non_missing_index];
                          
                        }
                        
                        previous_non_missing_index <- k
                        
                  }
                 ## print(paste("k = ", k))
            }
            n_cutpoints[c, s] <- cutpoint_counter - 1
           ## print(paste("c = ", c))
          }
         ## print(paste("s = ", s))
        }
        ##
        print(paste("n_cutpoints = ")) ; print(n_cutpoints)
        print(paste("cutpoint_index = ")) ; print(cutpoint_index[[1]])
        print(paste("n = ")) ; print(n[[1]])
        print(paste("x = ")) ; print(x[[1]])
        print(paste("x_with_missings = ")) ; print(x_with_missings[[1]])
        ##
        Stan_data_list$n_cutpoints <- n_cutpoints
        Stan_data_list$x <- x
        Stan_data_list$n <- n
        Stan_data_list$x_with_missings <- x_with_missings
        Stan_data_list$cutpoint_index <- cutpoint_index
        ##
        Stan_data_list$cts_thr_values_nd <- seq(from = 1, to = n_thr, by = 1)
        Stan_data_list$cts_thr_values_d  <- seq(from = 1, to = n_thr, by = 1)
        ##
        Stan_data_list$use_box_cox <- 1
        ##
        Stan_data_list$estimate_scales <- 0
        Stan_data_list$same_cutpoints_between_groups <- 0
        Stan_data_list$random_cutpoints <- 0
        ##
        ## Priors:
        ##
        Stan_data_list$prior_alpha    <-        rep(1.0, n_cat)
        Stan_data_list$prior_kappa_mean <-      rep(0.0, 2)
        Stan_data_list$prior_kappa_SD <-        rep(50.0, 2)
        ##
        Stan_data_list$prior_beta_mu_mean <- c(0.0, 0.0)
        Stan_data_list$prior_beta_mu_SD   <- c(1.0, 1.0)
        Stan_data_list$prior_beta_SD_mean <- c(0.0, 0.0)
        Stan_data_list$prior_beta_SD_SD   <- c(0.50, 0.50)
        ##
        Stan_data_list$prior_log_scale_mu_mean <- c(0.0, 0.0)
        Stan_data_list$prior_log_scale_mu_SD   <- c(1.0, 1.0)
        Stan_data_list$prior_log_scale_SD_mean <- c(0.0, 0.0)
        Stan_data_list$prior_log_scale_SD_SD   <- c(0.50, 0.50)
}
 




## - | ----------  Compile Stan model --------------------------------------------------------------------
{
  
        file <- file.path(getwd(), "stan_models", "DTA_MA_ord_RANDthr_w_NONsym_sROC.stan")
        #### file <- file.path(getwd(), "stan_models", "DTA_MA_cts_JONES_BOXCOX_NONsym_sROC.stan")
        
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
          
          n_chains <- 8
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
Stan_mod_sample$summary(c("kappa"))  %>% print(n = 100)
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



## Predicted data:
dev    <- Stan_mod_sample$summary(c("dev"))      %>% print(n = 25)
dev_nd <- Stan_mod_sample$summary(c("dev_nd"))   %>% print(n = 25)
dev_d  <- Stan_mod_sample$summary(c("dev_d"))    %>% print(n = 25)


dev_nd$median
dev_nd$mean

dev_nd$median
dev_nd$mean



{
      x_hat_nd <- Stan_mod_sample$summary(c("x_hat_nd")) ##  %>% print(n = 10)
      x_hat_d  <- Stan_mod_sample$summary(c("x_hat_d")) ##   %>% print(n = 10)
  
      x_hat_nd_array <- array(-1, dim = c(n_studies, n_thr))
      x_hat_d_array <- array(-1, dim = c(n_studies, n_thr))
      
        counter <- 0 
        for (k in 1:n_thr) {
          for (s in 1:n_studies) {
                  counter <- counter + 1
                  x_hat_nd_array[s, k] <- x_hat_nd$median[counter]
                  x_hat_d_array[s, k]  <- x_hat_d$median[counter]
          }
        }
        
        ## Overall differences in total cell counts:
        message(paste("Overall % difference in total cell counts (D-) = "))
        print(round(100 * sum(abs( round(x_hat_nd_array, 3) - Stan_data_list$x[[1]]) ) / sum( Stan_data_list$x[[1]][Stan_data_list$x[[1]] != -1] ), 2))
        message(paste("Overall % difference in total cell counts (D+) = "))
        print(round(100 * sum(abs( round(x_hat_d_array, 3)  - Stan_data_list$x[[2]]) ) / sum( Stan_data_list$x[[2]][Stan_data_list$x[[2]] != -1] ), 2))
        
        ## Study-specific differences in total cell counts:
        message(paste("Overall % difference in total cell counts (D-) = "))
        print(round(100 * (abs( round(x_hat_nd_array, 3) - Stan_data_list$x[[1]] ) ) / ( Stan_data_list$x[[1]][Stan_data_list$x[[1]] != -1] ), 2))
        message(paste("Overall % difference in total cell counts (D+) = "))
        print(round(100 * (abs( round(x_hat_d_array, 3)  - Stan_data_list$x[[2]]) ) / ( Stan_data_list$x[[2]][Stan_data_list$x[[2]] != -1] ), 2))
        
}


abs_diffs <- abs( round(x_hat_nd_array, 3) - Stan_data_list$x[[1]] )
abs_diffs_as_pct <- 100*abs_diffs/Stan_data_list$x[[1]]


abs_diffs_as_pct <- 100 * ( abs_diffs / (Stan_data_list$x[[1]]) )
rowSums(abs_diffs_as_pct, na.rm = TRUE)

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