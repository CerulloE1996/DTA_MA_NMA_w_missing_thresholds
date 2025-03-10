


## setwd("/home/enzocerullo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")

source("R_fn_load_data_ordinal_NMA_LC_MVP_sim.R")

require(dplyr)
require(cmdstanr)

n_tests <- 5
n_index_tests <- n_tests - 1
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

sim_results$n_thr_per_test
sim_results$N_per_study_vec

# index_test_chosen_index <- 3
# n_thr <- 10

n_thr <- sim_results$n_thr_per_test ; n_thr
n_cat <- n_thr + 1 ; n_cat
n_thr_max <- max(n_thr) ; n_thr_max
n_cat_max <- max(n_cat) ; n_cat_max

sim_results$Se_OVERALL_all_tests_all_thresholds
# sim_results$Se_OVERALL_all_tests_all_thresholds[1, 1:n_thr[1]]
sim_results$Se_OVERALL_all_tests_all_thresholds[2, 1:n_thr[2]]
 
# true_Se_OVERALL <- sim_results$Se_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr] ; true_Se_OVERALL
# true_Sp_OVERALL <- sim_results$Sp_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr] ; true_Sp_OVERALL

### plot(y = true_Se_OVERALL, x = 1 - true_Sp_OVERALL)

# sim_results$Se_per_study_all_tests_all_thresholds_list

## Now for the first example we will only take one index test, as initially we are evaluating just a "simple" model
## where there's only a single index test with 12 thresholds (so 13 categories:
##

##
y_list_example_1 <- list()
for (s in 1:n_studies) { 
    N <- nrow(y_list[[s]])
    y <- array(NA, dim = c(N, n_tests))
    for (t in 1:n_tests) { 
      y[, t] <- y_list[[s]][, t]
    }
    y_list_example_1[[s]] <- y
}

# y_list_example_1


## Convert to AGGREGATE data as these more basic NMA/MA models don't use individual-level data:
agg_data_cumulative <- convert_to_aggregate_counts( y_list = y_list_example_1, 
                                                    n_studies = n_studies,
                                                    n_tests = n_tests,
                                                    n_thr = n_thr)

agg_data_cumulative$Se_per_study_list


## Get OBSERVED "true" values BEFORE you remove any thresholds and BEFORE you remove any tests for a given study s!!!
Se_per_study_list <- agg_data_cumulative$Se_per_study_list
Sp_per_study_list <- agg_data_cumulative$Sp_per_study_list

true_Se_OVERALL <- true_Sp_OVERALL <- list() 
true_Se_OVERALL_weighted <- true_Sp_OVERALL_weighted <- list() 
true_Se_OVERALL <- list() 
true_Se_OVERALL <- list() 

for (t in 1:n_index_tests) {
  
      ##
      ## this is the means of the OBSERVED (i.e. specific to the seed set for the given simulation) study-specific Se's at each threshold:
      ##
      true_Se_OVERALL[[t]] <- 100 * colMeans(Se_per_study_list[[t]]) ; true_Se_OVERALL[[t]]
      true_Sp_OVERALL[[t]] <- 100 * colMeans(Sp_per_study_list[[t]]) ; true_Sp_OVERALL[[t]]
      ##
      ## this is equal to: Se @ threshold k  = {number of TOTAL true-positives across all studies} / {TOTAL N across all studies}:
      ##
      true_Se_OVERALL_weighted[[t]] <- sim_results$Se_OVERALL_all_tests_all_thresholds[t + 1, 1:n_thr[t + 1]]*100 ; true_Se_OVERALL_weighted[[t]]
      true_Sp_OVERALL_weighted[[t]] <- sim_results$Sp_OVERALL_all_tests_all_thresholds[t + 1, 1:n_thr[t + 1]]*100 ; true_Sp_OVERALL_weighted[[t]]
      ##
      ##
      print(true_Se_OVERALL[[t]] - true_Se_OVERALL_weighted[[t]])
      print(true_Sp_OVERALL[[t]] - true_Sp_OVERALL_weighted[[t]])

}


true_Se_OVERALL[[3]] - true_Se_OVERALL_weighted[[3]]
true_Sp_OVERALL[[3]] - true_Sp_OVERALL_weighted[[3]]


#### ------------------- Apply ** MISSING TESTS ** (hence "NMA" - optional): --------------------------------------------------------------------

outs_NMA <- create_test_in_study_for_NMA( n_studies = n_studies,
                                          n_tests = n_tests,
                                          prob_present = 0.50, 
                                          min_index_tests_per_study = 1)

indicator_index_test_in_study <- outs_NMA$indicator_index_test_in_study ; indicator_index_test_in_study
n_index_tests_per_study      <- outs_NMA$n_index_tests_per_study ; n_index_tests_per_study





#### ------------------- Apply ** MISSING THRESHOLDS ** (optional): -----------------------------------------------------------------------------

agg_data_cumulative_with_missing_thresholds <- agg_data_cumulative

agg_data_cumulative_with_missing_thresholds$x_d_list

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

agg_data_cumulative_with_missing_thresholds$x_nd_list
####
#### Non-diseased group:
####
x_non_diseased_cumulative <- agg_data_cumulative_with_missing_thresholds$x_nd_list
n_non_diseased_cumulative <- agg_data_cumulative_with_missing_thresholds$n_nd_list
n_total_non_diseased      <- agg_data_cumulative_with_missing_thresholds$n_total_nd
print(paste("x_non_diseased_cumulative = ")) ; print(x_non_diseased_cumulative)
####
#### Diseased group:
####
x_diseased_cumulative <- agg_data_cumulative_with_missing_thresholds$x_d_list
n_diseased_cumulative <- agg_data_cumulative_with_missing_thresholds$n_d_list
n_total_diseased      <- agg_data_cumulative_with_missing_thresholds$n_total_d
print(paste("x_diseased_cumulative = ")) ; print(x_diseased_cumulative)



# x_nd_cat <- convert_cumulative_to_category( cumulative_matrix = n_non_diseased_cumulative)
# x_non_diseased_cumulative
# x_nd_cat
# 
# 
# x_d_cat <- convert_cumulative_to_category( cumulative_matrix = n_diseased_cumulative)
# x_d_cat
# 
# 
# x_cat <- list(x_nd_cat, x_d_cat)


# y_nd <- categorical_to_individual(x_cat, missing_indicator = -1, binary_disease_indicator = 1)



##  | ------   Stan data   -------------------------------------------------------------------------
{
        n_cat <- n_thr + 1
        ##
        Stan_data_list <- list()
        ##
        Stan_data_list$n_studies <-  n_studies
        Stan_data_list$n_index_tests   <-  n_index_tests
        ##
        Stan_data_list$n_thr <- n_thr[-c(1)] ## exclude the reference test
        Stan_data_list$n_thr_max <- max(n_thr[-c(1)])
        ##
        Stan_data_list$n_index_tests_per_study <- n_index_tests_per_study
        ##
        ## Index tests only:
        ##
        Stan_data_list$indicator_index_test_in_study <- indicator_index_test_in_study
        ##
        # x_with_missings <- list(x_non_diseased_cumulative, 
        #                         x_diseased_cumulative)
        # Stan_data_list$x_with_missings <- x_with_missings
        ##
        cutpoint_index <- n <- x <- x_with_missings <- list()
        ##
        for (c in 1:2) {
          
              x_with_missings[[c]] <- list()
              x[[c]] <- list()
              n[[c]] <- list()
              cutpoint_index[[c]] <- list()
              ##
              for (t in 1:n_index_tests) {
                  x_with_missings[[c]][[t]] <- array(-1, dim = c(n_studies, n_thr_max)) ## "staggered" array.
                  x[[c]][[t]] <- array(-1, dim = c(n_studies, n_thr_max)) ## "staggered" array.
                  n[[c]][[t]] <- array(-1, dim = c(n_studies, n_thr_max)) ## "staggered" array.
                  cutpoint_index[[c]][[t]] <- array(-1, dim = c(n_studies, n_thr_max)) ## "staggered" array.
              }
        }
        ## Do "x_with_missings" first:
        for (t in 1:n_index_tests) {
            ##
            n_thr_t <- Stan_data_list$n_thr[t]
            for (s in 1:n_studies) {
               ##
               for (k in 1:n_thr_t) {
                  ##
                  for (c in 1:2) {
                     if (c == 1) x_with_missings[[c]][[t]][s, k] <- x_non_diseased_cumulative[[t]][s, k]
                     if (c == 2) x_with_missings[[c]][[t]][s, k] <- x_diseased_cumulative[[t]][s, k]
                  }
                 
               }
              
            }
          
        }
        ##
        Stan_data_list$x_with_missings <- x_with_missings
        ##
        ##
        n_obs_cutpoints <-  matrix(nrow = n_studies, ncol = n_index_tests)
        ##
        for (t in 1:n_index_tests) {
          
                n_thr_t <-  Stan_data_list$n_thr[t]
                ##
                for (s in 1:n_studies) {
                 
                      for (c in 1:2) {
                              
                                  cutpoint_counter = 0;
                                  previous_non_missing_index = 1;
                              
                                  cutpoint_index[[c]][[t]][s, 1] <- 1
                                  
                                  for (k in 2:(n_thr_t + 1)) {
                      
                                        if (k == (n_thr_t + 1)) {
                                          cond <- TRUE
                                        } else {
                                          cond <- (x_with_missings[[c]][[t]][s, k] != -1)
                                        }
                      
                                        if (cond == TRUE)   {
                                          
                                                  cutpoint_counter = cutpoint_counter +  1;
                                                  
                                                  if (k != (n_thr_t + 1)) {
                                                    cutpoint_index[[c]][[t]][s, cutpoint_counter + 1] = k;
                                                  }
                                                  
                                    
                                                  # print(paste("cutpoint_counter = ", cutpoint_counter))
                                                  # print(paste("previous_non_missing_index = ", previous_non_missing_index))
                                                  x[[c]][[t]][s, cutpoint_counter] <-   x_with_missings[[c]][[t]][s, previous_non_missing_index]
                                                  if (cutpoint_counter > 1) {
                                                     n[[c]][[t]][s, cutpoint_counter - 1] <-   x[[c]][[t]][s, cutpoint_counter]
                                                  }
                                              
                                                  
                                                  if (k == n_thr_t) {
                                                    n[[1]][[t]][s, n_thr_t] <- n_total_non_diseased[s]
                                                    n[[2]][[t]][s, n_thr_t] <- n_total_diseased[s]
                                                  }
                              
                                                  previous_non_missing_index <- k
                      
                                        }
                                       ## print(paste("k = ", k))
                                  }
                                  n_obs_cutpoints[s, t] <- cutpoint_counter
                             ## print(paste("c = ", c))
                      }
                     ## print(paste("s = ", s))
                }
        }
        x[[1]]
        x_with_missings[[1]]
        n[[1]]
        cutpoint_index[[1]]
        ##
        # x <- list(x_non_diseased_cumulative, x_diseased_cumulative)
        # n <- list(n_non_diseased_cumulative, n_diseased_cumulative)
        Stan_data_list$n_obs_cutpoints <- n_obs_cutpoints
        ##
        Stan_data_list$x_with_missings <- x_with_missings
        Stan_data_list$x <- x
        Stan_data_list$n <- n
        Stan_data_list$cutpoint_index <- cutpoint_index
        ##
        Stan_data_list$cts_thr_values_nd <- list()
        Stan_data_list$cts_thr_values_d  <- list()
        Stan_data_list$cts_thr_values    <- list()
        for (t in 1:n_index_tests) {
          # n_thr_t <- Stan_data_list$n_thr[t]
          Stan_data_list$cts_thr_values_nd[[t]] <- seq(from = 1, to = n_thr_max, by = 1)
          Stan_data_list$cts_thr_values_d[[t]]  <- seq(from = 1, to = n_thr_max, by = 1)
          Stan_data_list$cts_thr_values[[t]]    <- seq(from = 1, to = n_thr_max, by = 1)
        }
        ##
        Stan_data_list$use_box_cox <- 1
        ##
        Stan_data_list$estimate_scales <- 0
        Stan_data_list$same_cutpoints_between_groups <- 0
        ##
        ## Priors:
        ##
        Stan_data_list$log_alpha <- list()
        Stan_data_list$log_alpha[[1]] <- rep(0.001, n_thr_max + 1)
        Stan_data_list$log_alpha[[2]] <- rep(0.001, n_thr_max + 1)
        Stan_data_list$prior_alpha    <-        rep(1.0, n_thr_max + 1)
        Stan_data_list$prior_kappa_mean <-      rep(0.0, 2)
        Stan_data_list$prior_kappa_SD <-        rep(500.0, 2)
        ##
        Stan_data_list$prior_beta_mu_mean <- rep(0.0,  n_index_tests)
        Stan_data_list$prior_beta_mu_SD   <- rep(1.0,  n_index_tests)  ## * 1.702
        Stan_data_list$prior_beta_SD_mean <- rep(0.0,  n_index_tests)
        Stan_data_list$prior_beta_SD_SD   <- rep(0.50, n_index_tests)  ##* 1.702
        ##
        Stan_data_list$prior_raw_scale_mu_mean <- rep(0.0, n_index_tests)  ## rep(0.50, n_index_tests) ## since using: scale = log(1 + exp(raw_scale)) NOT scale = exp(raw_scale). 
        Stan_data_list$prior_raw_scale_mu_SD   <- rep(1.0,  n_index_tests)  ##* 1.702
        Stan_data_list$prior_raw_scale_SD_mean <- rep(0.0,  n_index_tests)
        Stan_data_list$prior_raw_scale_SD_SD   <- rep(0.50, n_index_tests) ##* 1.702
        ##
        Stan_data_list$prior_only <- 0
        ##
        ## NMA (Nyaga-like) priors:
        ##
        Stan_data_list$prior_beta_mu_mean  <- rep(0.0,  n_index_tests)
        Stan_data_list$prior_beta_mu_SD    <- rep(1.0,  n_index_tests)  ## * 1.702
        Stan_data_list$prior_beta_tau_SD   <- rep(2.5,  n_index_tests)
        Stan_data_list$prior_beta_sigma_SD <- 2.5 # 1.0
        ##
        Stan_data_list$prior_raw_scale_mu_mean  <- rep(0.0,  n_index_tests)
        Stan_data_list$prior_raw_scale_mu_SD    <- rep(1.0,  n_index_tests)  ## * 1.702
        Stan_data_list$prior_raw_scale_tau_SD   <- rep(2.5,  n_index_tests)
        Stan_data_list$prior_raw_scale_sigma_SD <- 2.5 # 1.0
        # ##
        # Stan_data_list$x_cat <- x_cat
        ##
        ## Induced-dirichlet priors:
        ##
        Stan_data_list$prior_alpha <- list()
        for (t in 1:n_index_tests) { 
            Stan_data_list$prior_alpha <- rep(1.0, n_thr_max + 1)
        }
        ##
        ##
        {
                n_total_cat <- 0
                n_total_cutpoints <- 0
                ##
                for (t in 1:n_index_tests) {
                        ##
                        n_thr_t <- Stan_data_list$n_thr[t]
                        n_cat_t <- n_thr_t + 1
                        ##
                        for (s in 1:n_studies) {
                               
                            # if (indicator_index_test_in_study[s, t] == 1) {
                                        ####  for (cut_i in 1:n_obs_cutpoints[s, t]) { //// only loop through the OBSERVED cutpoints for test t in study s
                                        for (k in 1:n_thr_t) {
                                          n_total_cutpoints = n_total_cutpoints + 1;
                                        }
                                        for (k in 1:n_cat_t) {
                                          n_total_cat = n_total_cat + 1;
                                        }
                            # }
                          
                        }
                        
                    
                }
                
                n_total_summary_cat <- 0
                for (t in 1:n_index_tests) {
                  for (k in 1:(Stan_data_list$n_thr[t] + 1)) {
                     n_total_summary_cat = n_total_summary_cat + 1;
                  }
                }
                
 
            ##
            {
              print(paste("n_total_cutpoints = ", n_total_cutpoints))
              print(paste("n_total_cat = ", n_total_cat))
              print(paste("n_total_summary_cat = ", n_total_summary_cat))
            }
            ##
            Stan_data_list$n_total_cutpoints <- n_total_cutpoints
            Stan_data_list$n_total_cat <- n_total_cat
            Stan_data_list$n_total_summary_cat <- n_total_summary_cat
        }
        ####
        #### Convert double-nested lists to 4D arrays:
        ####
        {
          
          if (!(is.array(Stan_data_list$x_with_missings))) {
            
            ## Initialize 4D arrays:
            x_with_missings_array <- array(NA_real_, dim = c(n_index_tests, 2, n_studies, n_thr_max))
            n_array <- array(NA_real_, dim = c(n_index_tests, 2, n_studies, n_thr_max))
            x_array <- array(NA_real_, dim = c(n_index_tests, 2, n_studies, n_thr_max))
            cutpoint_index_array <- array(NA_real_, dim = c(n_index_tests, 2, n_studies, n_thr_max))
            
            
            ## Fill 4D arrays w/ list elements:
            for (t in 1:n_index_tests) {
              for (c in 1:2) {
                # Make sure to check if the data exists
                if (!is.null(Stan_data_list$x_with_missings[[c]][[t]])) {
                  x_with_missings_array[t, c, , ] <- as.matrix(Stan_data_list$x_with_missings[[c]][[t]])
                }
                
                if (!is.null(Stan_data_list$n[[c]][[t]])) {
                  n_array[t, c, , ] <- as.matrix(Stan_data_list$n[[c]][[t]])
                }
                
                if (!is.null(Stan_data_list$x[[c]][[t]])) {
                  x_array[t, c, , ] <- as.matrix(Stan_data_list$x[[c]][[t]])
                }
                
                if (!is.null(Stan_data_list$cutpoint_index[[c]][[t]])) {
                  cutpoint_index_array[t, c, , ] <- as.matrix(Stan_data_list$cutpoint_index[[c]][[t]])
                }
              }
            }
            
            ## Replace list structures with arrays in "Stan_data_list":
            Stan_data_list$x_with_missings <- x_with_missings_array
            Stan_data_list$n <- n_array
            Stan_data_list$x <- x_array
            Stan_data_list$cutpoint_index <- cutpoint_index_array
            
          }
          
        }
        
          
          
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


