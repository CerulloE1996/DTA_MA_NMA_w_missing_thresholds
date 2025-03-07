


## setwd("/home/enzocerullo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")

source("R_fn_load_data_ordinal_MA_LC_MVP_sim.R")

require(dplyr)
require(cmdstanr)

n_studies <- 25
N_per_study_mean <- 2500
N_per_study_SD <- 250
assume_perfect_GS <- 1
seed <- 123

set.seed(seed)

missing_indicator <- -1

options(max.print = 100000)


# Run simulated data - this simulates data from FIVE (5) diagnostic tests (1 BINARY reference test + 4 ORDINAL index tests)
sim_results <- simulate_binary_and_ordinal_MA_LC_MVP_data(n_studies = n_studies,
                                                          N_per_study_mean = N_per_study_mean,
                                                          N_per_study_SD = N_per_study_SD,
                                                          assume_perfect_GS = assume_perfect_GS,
                                                          seed = seed)

y_list <- sim_results$y_list
str(y_list)


sim_results$n_total_obs_thr_per_test
sim_results$N_per_study_vec

index_test_chosen_index <- 3
n_thr <- 10


sim_results$Se_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr]

 

# true_Se_OVERALL <- sim_results$Se_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr] ; true_Se_OVERALL
# true_Sp_OVERALL <- sim_results$Sp_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr] ; true_Sp_OVERALL

### plot(y = true_Se_OVERALL, x = 1 - true_Sp_OVERALL)

# sim_results$Se_per_study_all_tests_all_thresholds_list

## Now for the first example we will only take one index test, as initially we are evaluating just a "simple" model
## where there's only a single index test with 12 thresholds (so 13 categories:
##

##
n_cat <- n_thr + 1
y_list_example_1 <- list()
for (s in 1:n_studies) { 
    N <- nrow(y_list[[s]])
    y <- array(NA, dim = c(N, 2))
    y[, 1] <- y_list[[s]][, 1]
    y[, 2] <- y_list[[s]][, index_test_chosen_index] ## take test #5 as the index test for this example 
    y_list_example_1[[s]] <- y
}

y_list_example_1



# Convert to AGGREGATE data as these more basic NMA/MA models don't use individual-level data:
agg_data_cumulative <- convert_to_aggregate_counts( y_list_example_1, 
                                                    n_studies,
                                                    n_thr)

agg_data_cumulative


agg_data_cumulative$Se_per_study
agg_data_cumulative$Sp_per_study
# this is the means of the OBSERVED (i.e. specific to the seed set for the given simulation) study-specific Se's at each threshold:
true_Se_OVERALL <- 100 * colMeans(agg_data_cumulative$Se_per_study) ; true_Se_OVERALL 
true_Sp_OVERALL <- 100 * colMeans(agg_data_cumulative$Sp_per_study) ; true_Sp_OVERALL
# this is equal to: Se @ threshold k  = {number of TOTAL true-positives across all studies} / {TOTAL N across all studies}:
true_Se_OVERALL_weighted <- sim_results$Se_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr]*100 ; true_Se_OVERALL_weighted 
true_Sp_OVERALL_weighted <- sim_results$Sp_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr]*100 ; true_Sp_OVERALL_weighted 
##

true_Se_OVERALL - true_Se_OVERALL_weighted
true_Sp_OVERALL - true_Sp_OVERALL_weighted


## Apply missing thresholds:
agg_data_cumulative_with_missing_thresholds <- agg_data_cumulative



# {
# 
#     ## First 2 studies only report data for the "middle 8" thresholds (starting from threshold 3, so thr = {3, 4, 5, 6, 7, 8, 9, 10})
#     studies_subset_vec <- c(1, 2)
#     missing_thr_subset_vec <- c(1:2, 11:12)
#     agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
#                                                                             studies_subset_vec = studies_subset_vec,
#                                                                             missing_thr_subset_vec = missing_thr_subset_vec,
#                                                                             missing_indicator = missing_indicator)
#     ## The next 2 studies only report at the "middle 4" thresholds (starting at threshold 5): so thr = {5, 6, 7, 8}:
#     studies_subset_vec <- c(3, 4)
#     missing_thr_subset_vec <- c(1:4, 9:12)
#     agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
#                                                                             studies_subset_vec = studies_subset_vec,
#                                                                             missing_thr_subset_vec = missing_thr_subset_vec,
#                                                                             missing_indicator = missing_indicator)
#     ## The next two studies only report data at only a single threshold - threshold #7:
#     studies_subset_vec <- c(5, 6)
#     missing_thr_subset_vec <- c(1:6, 8:12)
#     agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
#                                                                             studies_subset_vec = studies_subset_vec,
#                                                                             missing_thr_subset_vec = missing_thr_subset_vec,
#                                                                             missing_indicator = missing_indicator)
#     ## The next two studies only report data at 4 non-adjacent thresholds - so thr = {3, 5, 7, 9}:
#     studies_subset_vec <- c(7, 8)
#     missing_thr_subset_vec <- c(1, 2, 4, 6, 8, 10, 11, 12)
#     agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
#                                                                             studies_subset_vec = studies_subset_vec,
#                                                                             missing_thr_subset_vec = missing_thr_subset_vec,
#                                                                             missing_indicator = missing_indicator)
#     ## And finally, the last 2 studies report data at ALL thresholds:
#     #  -----  Do nothing here - as no missing thresholds for study #10!
# 
# }
# 
# ## Now let's look at the overall % of missing thresholds:
# total_missing_thr <- sum(agg_data_cumulative_with_missing_thresholds$x_diseased == 0.999) ; total_missing_thr
# prop_missing_thr <- total_missing_thr / (n_studies * n_thr) ; prop_missing_thr
# ## Just under half the threshold data is missing (46.67%).
# 
# #  
#  
 
# ## Diseased group:
x_diseased_cumulative <- agg_data_cumulative_with_missing_thresholds$x_diseased
n_diseased_cumulative <- agg_data_cumulative_with_missing_thresholds$n_diseased
n_total_diseased <- agg_data_cumulative_with_missing_thresholds$n_total_diseased
print(paste("x_diseased_cumulative = ")) ; print(x_diseased_cumulative)
## Non-diseased group:
x_non_diseased_cumulative <- agg_data_cumulative_with_missing_thresholds$x_non_diseased
n_non_diseased_cumulative <- agg_data_cumulative_with_missing_thresholds$n_non_diseased
n_total_non_diseased <- agg_data_cumulative_with_missing_thresholds$n_total_non_diseased
print(paste("x_non_diseased_cumulative = ")) ; print(x_non_diseased_cumulative)


x_nd_cat <- convert_cumulative_to_category( cumulative_matrix = n_non_diseased_cumulative)
x_non_diseased_cumulative
x_nd_cat


x_d_cat <- convert_cumulative_to_category( cumulative_matrix = n_diseased_cumulative)
x_d_cat

x_cat <- list(x_nd_cat, x_d_cat)


# y_nd <- categorical_to_individual(x_cat, missing_indicator = -1, binary_disease_indicator = 1)


##  | ------   Stan data   -------------------------------------------------------------------------
{
        n_cat <- n_thr + 1
        ##
        Stan_data_list <- list()
        ##
        Stan_data_list$n_studies <-  n_studies
        Stan_data_list$n_thr <-  n_thr
        ##
        Stan_data_list$n_non_diseased <- agg_data_cumulative_with_missing_thresholds$n_total_non_diseased
        Stan_data_list$n_diseased <-     agg_data_cumulative_with_missing_thresholds$n_total_diseased
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
        #
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
                  
                      cutpoint_index[[c]][s, 1] <- 1
                      
                      for (k in 2:(n_thr + 1)) {
          
                            if (k == (n_thr + 1)) {
                              cond <- TRUE
                            } else {
                              cond <- (x_with_missings[[c]][s, k] != -1)
                            }
          
                            if (cond == TRUE)   {
                              
                                      cutpoint_counter = cutpoint_counter +  1;
                                      
                                      if (k != (n_thr + 1)) {
                                        cutpoint_index[[c]][s, cutpoint_counter + 1] = k;
                                      }
                                      
                        
                                      # print(paste("cutpoint_counter = ", cutpoint_counter))
                                      # print(paste("previous_non_missing_index = ", previous_non_missing_index))
                                      x[[c]][s, cutpoint_counter] <-   x_with_missings[[c]][s, previous_non_missing_index]
                                      if (cutpoint_counter > 1) {
                                         n[[c]][s, cutpoint_counter - 1] <-   x[[c]][s, cutpoint_counter]
                                      }
                                  
                                      
                                      if (k == n_thr) {
                                        n[[1]][s, n_thr] <- n_total_non_diseased[s]
                                        n[[2]][s, n_thr] <- n_total_diseased[s]
                                      }
                  
                                      previous_non_missing_index <- k
          
                            }
                           ## print(paste("k = ", k))
                      }
                  n_cutpoints[c, s] <- cutpoint_counter
                 ## print(paste("c = ", c))
          }
         ## print(paste("s = ", s))
        }
        x[[1]]
        x_with_missings[[1]]
        n[[1]]
        ##
        # x <- list(x_non_diseased_cumulative, x_diseased_cumulative)
        # n <- list(n_non_diseased_cumulative, n_diseased_cumulative)
        Stan_data_list$n_cutpoints <- n_cutpoints
        Stan_data_list$x <- x
        Stan_data_list$n <- n
        Stan_data_list$x_with_missings <- x_with_missings
        Stan_data_list$cutpoint_index <- cutpoint_index
        ##
        Stan_data_list$cts_thr_values_nd <- seq(from = 1, to = n_thr, by = 1)
        Stan_data_list$cts_thr_values_d  <- seq(from = 1, to = n_thr, by = 1)
        Stan_data_list$cts_thr_values  <- seq(from = 1, to = n_thr, by = 1)
        ##
        Stan_data_list$use_box_cox <- 1
        ##
        Stan_data_list$estimate_scales <- 0
        Stan_data_list$same_cutpoints_between_groups <- 0
        ##
        ## Priors:
        ##
        Stan_data_list$log_alpha <- list()
        Stan_data_list$log_alpha[[1]] <- rep(0.001, n_cat)
        Stan_data_list$log_alpha[[2]] <- rep(0.001, n_cat)
        Stan_data_list$prior_alpha    <-        rep(1.0, n_cat)
        Stan_data_list$prior_kappa_mean <-      rep(0.0, 2)
        Stan_data_list$prior_kappa_SD <-        rep(500.0, 2)
        ##
        Stan_data_list$prior_beta_mu_mean <- rep(0.0, 2)
        Stan_data_list$prior_beta_mu_SD   <- rep(1.0, 2)  ## * 1.702
        Stan_data_list$prior_beta_SD_mean <- rep(0.0, 2)
        Stan_data_list$prior_beta_SD_SD   <- rep(0.50, 2)  ##* 1.702
        ##
        Stan_data_list$prior_raw_scale_mu_mean <- rep(0.50, 2) ## since using: scale = log(1 + exp(raw_scale)) NOT scale = exp(raw_scale). 
        Stan_data_list$prior_raw_scale_mu_SD   <- rep(1.0, 2)  ##* 1.702
        Stan_data_list$prior_raw_scale_SD_mean <- rep(0.0, 2)
        Stan_data_list$prior_raw_scale_SD_SD   <- rep(0.50, 2) ##* 1.702
        ##
        Stan_data_list$prior_only <- 0
        
        ##
        Stan_data_list$x_cat <- x_cat
}
 



# Stan_data_list

# 
# reorder_increasing <- function(matrix_data) {
#   
#       # Apply to each row
#       t(apply(matrix_data, 1, function(row) {
#         # Sort in ascending order
#         sort(row)
#       }))
#   
# }
# 
# 
# reorder_decreasing <- function(matrix_data) {
#   
#   # Apply to each row
#   t(apply(matrix_data, 1, function(row) {
#     # Sort in ascending order
#     sort(row, decreasing = TRUE)
#   }))
#   
# }
# 
# 
# {
#       Stan_data_list$x_diseased     <- reorder_increasing(Stan_data_list$x_diseased)
#       Stan_data_list$x_non_diseased <- reorder_increasing(Stan_data_list$x_non_diseased)
#       
#       for (c in 1:2) { 
#         Stan_data_list$n[[c]] <- reorder_increasing(Stan_data_list$n[[c]])
#         Stan_data_list$x[[c]] <- reorder_increasing(Stan_data_list$x[[c]])
#         Stan_data_list$x_with_missings[[c]] <- reorder_increasing(Stan_data_list$x_with_missings[[c]])
#       }
# }
# 
# 
# 
##
print(paste("n_cutpoints = ")) ; print(Stan_data_list$n_cutpoints)
print(paste("cutpoint_index = ")) ; print(Stan_data_list$cutpoint_index)
print(paste("n = ")) ; print(Stan_data_list$n[[1]])
print(paste("x = ")) ; print(Stan_data_list$x[[1]])
print(paste("x_with_missings = ")) ; print(Stan_data_list$x_with_missings[[1]])
# 
# 
# 
# print(paste("true_Se_OVERALL = ")) ; print(true_Se_OVERALL)
# print(paste("true_Sp_OVERALL = ")) ; print(true_Sp_OVERALL)
# 


