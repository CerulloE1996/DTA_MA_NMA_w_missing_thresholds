
 
 
 
 
 
 n_studies <- 10
 N_per_study_mean <- 500
 N_per_study_SD <- 100
 assume_perfect_GS <- 1
 seed <- 123
 
## -| ------------------ R function to simulate a meta-analysis dataset for (binary + ordinal) LC_MVP  --------------------------------------------
 
 
simulate_binary_and_ordinal_MA_LC_MVP_data <- function(   n_studies = 10,
                                                          N_per_study_mean = 500,
                                                          N_per_study_SD = 100,
                                                          assume_perfect_GS = 1,
                                                          seed = 123) {
   
         require(TruncatedNormal)
         
         n_binary_tests <- 1 ## Single binary "reference" test - an imperfect gold standard (GS).
         n_ordinal_tests <- 4 # fixed to 4
         n_tests <- n_binary_tests + n_ordinal_tests
         
         ## set the seed (keep this OUTSIDE the N loop):
         set.seed(seed, kind = "L'Ecuyer-CMRG") 
         
         ## Make storage lists:
         y_ord_and_bin_list <- list()
         Sigma_nd_true_observed_list <- Sigma_d_true_observed_list <- list()
         prev_true_observed_list <- Se_true_observed_list <- Sp_true_observed_list <- list()
         Phi_Se_observed_list <- Phi_Fp_observed_list <- list()
         true_correlations_observed_vec_list <- observed_table_probs_list <- true_estimates_observed_list <- observed_cell_counts_list <- list()
         
         ii_dataset <- 0 ## counter
         
         N_per_study_SD <- 250 ##  N_per_study_mean / 2
         N_per_study_vec <- round(TruncatedNormal::rtnorm(n = n_studies, mu = N_per_study_mean, sd = N_per_study_SD, lb = 100, ub = 2500), 0)
      
         ## True values for the mean Se and Sp across all n_studies studies @threshold of 0:
         Mean_Se_at_threshold_0 <- c(0.60, 0.55, 0.60, 0.65, 0.70)
         Mean_Sp_at_threshold_0 <- c(0.99, 0.95, 0.90, 0.90, 0.85)
         
         if (assume_perfect_GS == 1) { 
           Mean_Se_at_threshold_0[1] <- 0.999999999
           Mean_Sp_at_threshold_0[1] <- 0.999999999
         }
         
         Mean_Fp_at_threshold_0 <- 1 - Mean_Sp_at_threshold_0
         Phi_Mean_Se_at_threshold_0 <- qnorm(Mean_Se_at_threshold_0)
         Phi_Mean_Fp_at_threshold_0 <- qnorm(Mean_Fp_at_threshold_0)
         
         ## Set true values for disease prevelance:
         true_Mean_prev <- 0.20
         true_Mean_probit_prev <- qnorm(true_Mean_prev)
         true_SD_probit_prev <- 0.25
         true_probit_prev_per_study <- rnorm(n = n_studies, mean = true_Mean_probit_prev, sd = true_SD_probit_prev)
         true_prev_per_study <- pnorm(true_probit_prev_per_study)
         
         ## Set true calues for the between-study heterogenity params:
         SD_of_Phi_Se_at_threshold_0 <- 0.50
         SD_of_Phi_Fp_at_threshold_0 <- 0.25
         
         ## Threshold info (test 1 is binary so only 1 threshold @0 on the "latent" scale)::
         n_total_possible_thresholds_per_test <- c(1, 3, 5, 5, 12)
         max_threshold_across_all_tests <- max(n_total_possible_thresholds_per_test)
         threshold_per_study_array <- array(NA, dim = c(n_studies, max_threshold_across_all_tests))
         Mean_of_thresholds_for_all_tests_array <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
         SD_of_thresholds_for_all_tests_array <- array(0.25, dim = c(n_tests, max_threshold_across_all_tests)) ## between-study heterogenity for thresholds 
         
         ## Set the "true values" of the mean threshold between studies for each test (on the * latent * scale):
         ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds! (i.e. the TRUE threshold parameters are "homogenous between classes")
         ## Binary reference test:
         Mean_thresholds_for_test_1 <- rep(NA, n_total_possible_thresholds_per_test[1])
         Mean_thresholds_for_test_1 <- c(0.0)
         Mean_of_thresholds_for_all_tests_array[1, 1:n_total_possible_thresholds_per_test[1]] <- Mean_thresholds_for_test_1
         ## Test 2:
         Mean_thresholds_for_test_2 <- rep(NA, n_total_possible_thresholds_per_test[2])
         Mean_thresholds_for_test_2 <- c(-1.5, -0.50, 0.75)
         Mean_of_thresholds_for_all_tests_array[2, 1:n_total_possible_thresholds_per_test[2]] <- Mean_thresholds_for_test_2
         ## Test 3:
         Mean_thresholds_for_test_3 <- rep(NA, n_total_possible_thresholds_per_test[3])
         Mean_thresholds_for_test_3 <- c(-1.75, -1.00, -0.50, 0.25, 1.25)
         Mean_of_thresholds_for_all_tests_array[3, 1:n_total_possible_thresholds_per_test[3]] <- Mean_thresholds_for_test_3
         ## Test 4:
         Mean_thresholds_for_test_4 <- rep(NA, n_total_possible_thresholds_per_test[4])
         Mean_thresholds_for_test_4 <- c(-1.50, -1.10, -0.75, 0.00, 0.90)
         Mean_of_thresholds_for_all_tests_array[4, 1:n_total_possible_thresholds_per_test[4]] <- Mean_thresholds_for_test_4
         ## Test 5:
         Mean_thresholds_for_test_5 <- rep(NA, n_total_possible_thresholds_per_test[5])
         Mean_thresholds_for_test_5 <- c(-2.50, -2.25, -2.05, -1.75, -1.50, -1.00, -0.50, 0.10, 0.75, 0.95, 1.15, 2.25)
         Mean_of_thresholds_for_all_tests_array[5, 1:n_total_possible_thresholds_per_test[5]] <- Mean_thresholds_for_test_5
         
         s <- 1
         t <- 2
         
         Se_for_current_study_at_threshold_0_list <- Sp_for_current_study_at_threshold_0_list <- Fp_for_current_study_at_threshold_0_list <- list()
         Se_per_study_all_tests_all_thresholds_list <- Sp_per_study_all_tests_all_thresholds_list <- list()
         prev_true_observed_list <- list()
         Se_per_study_ref <- Sp_per_study_ref <- list()
         thresholds_for_all_tests_for_current_study_array_list <- list()
         y_list <- list()
         y_df_list <- list()
         
         for (s in 1:n_studies) {
           
                   # ii_dataset <- ii_dataset + 1 ## increment counter
                   
                   N <- N_per_study_vec[s]
                   true_prev <- true_prev_per_study[s]
                   
                   ## Generate "true" thresholds for study s:
                   ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds! (i.e. the TRUE threshold parameters are "homogenous between classes")
                   thresholds_for_all_tests_for_current_study_array <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
                   for (t in 2:n_tests) {
                     thresholds_for_all_tests_for_current_study_array[t, 1:n_total_possible_thresholds_per_test[t]] <- rnorm(n = n_total_possible_thresholds_per_test[t], 
                                                                                                                             mean = Mean_of_thresholds_for_all_tests_array[t, 1:n_total_possible_thresholds_per_test[t]], 
                                                                                                                             sd = SD_of_thresholds_for_all_tests_array[t, 1:n_total_possible_thresholds_per_test[t]])
                   }
                   
                   thresholds_for_all_tests_for_current_study_array_list[[s]] <- thresholds_for_all_tests_for_current_study_array
                   
                   Se_for_current_study_at_threshold_0 <- pnorm(rnorm(n = n_tests, mean = Phi_Mean_Se_at_threshold_0, sd = SD_of_Phi_Se_at_threshold_0))
                   Fp_for_current_study_at_threshold_0 <- pnorm(rnorm(n = n_tests, mean = Phi_Mean_Fp_at_threshold_0, sd = SD_of_Phi_Fp_at_threshold_0))
                   Sp_for_current_study_at_threshold_0 <- 1 - Fp_for_current_study_at_threshold_0
                   
                   Sigma_highly_varied <- matrix(c(1,  0,     0,        0,        0,
                                                   0,  1,     0.50,     0.25,     0,
                                                   0,  0.50,  1,        0.40,     0.40,
                                                   0,  0.25,  0.40,     1,        0.70,
                                                   0,  0,     0.40,     0.70,     1),
                                                 n_tests, n_tests)
              
                   {
                         true_Fp_at_threshold_0_vec <-        1 - Sp_for_current_study_at_threshold_0
                         true_Se_at_threshold_0_vec <-            Se_for_current_study_at_threshold_0
                         
                         Sigma_d <- Sigma_highly_varied
                         Sigma_nd <-  0.5 * Sigma_highly_varied
                         diag(Sigma_nd) <- rep(1, n_tests)
                   }
                   
                   L_Sigma_d  <-  t(chol(Sigma_d)) # PD check (fails if not PD)
                   L_Sigma_nd  <- t(chol(Sigma_nd)) # PD check (fails if not PD)
                   
                   d_ind <- sort(rbinom(n = N, size = 1, prob = true_prev))
                   n_pos <- sum(d_ind)
                   n_neg <- N - sum(d_ind)
                   
                   ## Simulate the underlying "latent" continuous test responses which are distributed multivariate normal:
                   latent_cts_results_neg <- LaplacesDemon::rmvn(n = n_neg, mu = qnorm(true_Fp_at_threshold_0_vec), Sigma = Sigma_nd)
                   latent_cts_results_pos <- LaplacesDemon::rmvn(n = n_pos, mu = qnorm(true_Se_at_threshold_0_vec), Sigma = Sigma_d)
                   latent_results <- rbind(latent_cts_results_neg, latent_cts_results_pos)
                   
                   ## Now simulate the BINARY data for the reference test (i.e. the "imperfect gold standard" - test # 1):
                   results_pos <- array(NA, dim = c(n_pos, n_tests))
                   results_neg <- array(NA, dim = c(n_neg, n_tests))
                   results <-     array(NA, dim = c(n_neg + n_pos, n_tests))
                   results_neg[, 1] <- ifelse(latent_cts_results_neg[, 1] > 0.0, 1, 0)
                   results_pos[, 1] <- ifelse(latent_cts_results_pos[, 1] > 0.0, 1, 0)
                   
                   ## Now simulate the ORDINAL data for tests 2-5:
                   ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds! (i.e. the TRUE threshold parameters are "homogenous between classes")
                  ##  thresholds_for_all_tests_for_current_study_array[t, 1:n_total_possible_thresholds_per_test[t]]
                 
                   # k <- 1
                   # threshold_lower <- -1000
                   # threshold_upper <- thresholds_for_all_tests_for_current_study_array[2, k]
                   # results_neg[, 2] <- ifelse(latent_cts_results_neg[, 2] %in% c(threshold_lower:threshold_upper), 0, results_neg[, 2])
                   # 
                   # # k <- 2
                   # # t <- 3
                   
                   ## Now simulate the ORDINAL data for tests 2-5:
                   ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds! (i.e. the TRUE threshold parameters are "homogenous between classes")
                   ## Hence we can just use the "latent_results" array directly (no need to split into +'ve and -'ve)
                   for (t in 2:n_tests) {
                     
                           K <- n_total_possible_thresholds_per_test[t] ; K
                           threshold_lower <- -1000
                           threshold_upper <- thresholds_for_all_tests_for_current_study_array[t, 1]
                           results_neg[, t] <- ifelse(latent_cts_results_neg[, t] < threshold_upper, 
                                                      1, 
                                                      results_neg[, t])
                        
                           for (k in 2:(K)) {
                               threshold_lower <- thresholds_for_all_tests_for_current_study_array[t, k - 1] ; threshold_lower
                               threshold_upper <- thresholds_for_all_tests_for_current_study_array[t, k]     ; threshold_upper
                               results_neg[, t] <- ifelse( (latent_cts_results_neg[, t] > threshold_lower) & (latent_cts_results_neg[, t] < threshold_upper), 
                                                             k, 
                                                           results_neg[, t])
                           }
                           
                           threshold_lower <-  thresholds_for_all_tests_for_current_study_array[t, K]
                           threshold_upper <- +1000
                           results_neg[, t] <- ifelse(latent_cts_results_neg[, t] > threshold_lower,  
                                                      K + 1, 
                                                  results_neg[, t])
                       
                   }
                   results_neg
                   for (t in 2:n_tests) {
                     
                     K <- n_total_possible_thresholds_per_test[t] ; K
                     threshold_lower <- -1000
                     threshold_upper <- thresholds_for_all_tests_for_current_study_array[t, 1]
                     results_pos[, t] <- ifelse(latent_cts_results_pos[, t] < threshold_upper, 
                                                1, 
                                                results_pos[, t])
                     
                     for (k in 2:(K)) {
                       threshold_lower <- thresholds_for_all_tests_for_current_study_array[t, k - 1] ; threshold_lower
                       threshold_upper <- thresholds_for_all_tests_for_current_study_array[t, k]     ; threshold_upper
                       results_pos[, t] <- ifelse( (latent_cts_results_pos[, t] > threshold_lower) & (latent_cts_results_pos[, t] < threshold_upper), 
                                                   k, 
                                                   results_pos[, t])
                     }
                     
                     threshold_lower <-  thresholds_for_all_tests_for_current_study_array[t, K]
                     threshold_upper <- +1000
                     results_pos[, t] <- ifelse(latent_cts_results_pos[, t] > threshold_lower,  
                                                K + 1, 
                                                results_pos[, t])
                     
                   }
                   results_pos
                   
 
                   
                   ## Put it all together to get the daaset - y:
                   results <- rbind(results_neg, results_pos)
                   y <- results
                   
                   max(y[,2:n_tests])
                   min(y[,2:n_tests])
                   
                   df <- dplyr::tibble(results, latent_results, d_ind)
                   df_pos <- dplyr::filter(df, d_ind == 1)
                   df_neg <- dplyr::filter(df, d_ind == 0)

                   # prev:
                   prev_true_observed <-  print(round(sum(d_ind)/N, 3))
                   prev_true_observed_list[[s]] <- prev_true_observed
                   
                   # Compute per-study observed TRUE Se:
                   Se_current_study_all_tests_all_thresholds <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
                   for (t in 2:n_tests) {
                       for (k in 1:K) {
                         n_TP_at_current_threshold <- sum(df_pos$results[, t] > k)
                         Se_at_threshold_k <- n_TP_at_current_threshold/n_pos
                         print(paste("Se_at_threshold", k, "for test", t, " = ", round(Se_at_threshold_k, 2)))
                         # Phi_Se_at_threshold_k <- qnorm(Se_at_threshold_k)
                         # print(Se_at_threshold_k)
                         Se_current_study_all_tests_all_thresholds[t, k] <- Se_at_threshold_k
                         Se_per_study_all_tests_all_thresholds_list[[s]] <- Se_current_study_all_tests_all_thresholds
                       }
                   }
             
                   # Compute per-study observed TRUE Fp / Sp:
                   Sp_current_study_all_tests_all_thresholds <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
                   for (t in 2:n_tests) {
                     for (k in 1:K) {
                       n_TN_at_current_threshold <- sum(df_neg$results[, t] > k)
                       Fp_at_threshold_k <- n_TN_at_current_threshold/n_pos
                       # print(paste("Fp_at_threshold_k = ", round(Fp_at_threshold_k, 2)))
                       Sp_at_threshold_k <- 1.0 - Fp_at_threshold_k
                       print(paste("p_at_threshold", k, "for test", t, " = ", round(Sp_at_threshold_k, 2)))
                       # Phi_Sp_at_threshold_k <- qnorm(Sp_at_threshold_k)
                       # print(Sp_at_threshold_k)
                       Sp_current_study_all_tests_all_thresholds[t, k] <- Sp_at_threshold_k
                       Sp_per_study_all_tests_all_thresholds_list[[s]] <- Sp_current_study_all_tests_all_thresholds
                     }
                   }
                   
                     ## For reference test:
                     t <- 1
                     ## Se:
                     Phi_Se_ref <- qnorm(sum(df_pos$results[,t])/nrow(df_pos))
                     Se_ref <- pnorm(Phi_Se_ref)
                     Se_per_study_ref[[s]] <- Se_ref
                     
                     ## Sp:
                     Phi_Fp_ref <-  qnorm( 1.0 - ((nrow(df_neg) - sum(df_neg$results[,t]))/nrow(df_neg))  )
                     Fp_ref <- pnorm(Phi_Fp_ref)
                     Sp_ref <- 1.0 - Fp_ref
                     Sp_per_study_ref[[s]] <- Sp_ref
 
                     
                     y_list[[s]] <- y
                     y_df_list[[s]] <- data.frame(y)
                   
                   # 
                   # Sigma_nd_true_observed <- Sigma_d_true_observed <- array(dim = c(n_tests, n_tests))
                   # observed_correlations <- array(dim = c(n_tests, n_tests))
                   #
                   # for (i in 2:n_tests) {
                   #   for (j in 1:(i-1)) {
                   #     Sigma_nd_true_observed[i, j] <- cor(df_neg$latent_results[,i], df_neg$latent_results[,j])
                   #     Sigma_nd_true_observed[j, i] <-  Sigma_nd_true_observed[i, j]
                   #     Sigma_d_true_observed[i, j] <- cor(df_pos$latent_results[,i], df_pos$latent_results[,j])
                   #     Sigma_d_true_observed[j, i] <-  Sigma_d_true_observed[i, j]
                   #     observed_correlations[i, j] <- cor(y[, i], y[, j])
                   #     observed_correlations[j, i] <-  observed_correlations[i, j]
                   #   }
                   # }
                   # 
                   
                  
                   # true_corrs_observed_vec <- observed_correlations[upper.tri(observed_correlations )]
                   # 
                   # observed_table <- table(y[, 1], y[, 2], y[, 3], y[, 4], y[, 5])
                   # observed_table_probs_vec <- c(unlist(round(prop.table(observed_table), 4)))
                   # 
                   # observed_cell_counts <- observed_table_probs_vec * N
                   # 
                   # true_estimates_observed <-  c(Sigma_nd_true_observed[upper.tri(Sigma_nd_true_observed )],  Sigma_d_true_observed[upper.tri(Sigma_d_true_observed )], Sp_true_observed,  Se_true_observed, prev_true_observed ,
                   #                               true_corrs_observed_vec, observed_table_probs_vec, NA, NA)
                   # true_estimates  <-  c(Sigma_nd[upper.tri(Sigma_nd )],  Sigma_d[upper.tri(Sigma_d )], 1 - true_Fp_vec,  true_Se_vec, true_prev  ,
                   #                       rep(NA, length(true_corrs_observed_vec)), rep(NA, length(observed_table_probs_vec)), NA, NA)
                   #
                   # ## Print some key quantities:
                   # print(paste("N = ", N))
                   # print(paste("prev_true_observed = ", prev_true_observed))
                   # print(paste("Se_true_observed = ", Se_true_observed))
                   # print(paste("Sp_true_observed = ", Sp_true_observed))
                   
                   # ## Populate lists:
                   # y_binary_list[[ii_dataset]] <- y
                   # ##
                   # Sigma_nd_true_observed_list[[ii_dataset]] <- Sigma_nd_true_observed
                   # Sigma_d_true_observed_list[[ii_dataset]] <- Sigma_d_true_observed
                   # ##
                   # Phi_Se_observed_list[[ii_dataset]] <- Phi_Se_observed_vec
                   # Phi_Fp_observed_list[[ii_dataset]] <- Phi_Fp_observed_vec
                   ##
                
                   # ##
                   # true_correlations_observed_vec_list[[ii_dataset]] <- true_corrs_observed_vec
                   # observed_table_probs_list[[ii_dataset]] <- observed_table_probs_vec
                   # true_estimates_observed_list[[ii_dataset]] <- true_estimates_observed
                   # observed_cell_counts_list[[ii_dataset]] <- observed_cell_counts
           
           
         }
         
         y_tibble <- NULL
         y_tibble <- tibble(data.table::rbindlist(y_df_list, idcol = "Study"))
         
         return(list(
           y_list = y_list,
           y_df_list = y_df_list,
           y_tibble = y_tibble,
           ## Between-study true params:
           Mean_Se_at_threshold_0 = Mean_Se_at_threshold_0,
           Mean_Sp_at_threshold_0 = Mean_Sp_at_threshold_0,
           true_Mean_prev = true_Mean_prev,
           ## Between-study true params:
           SD_of_Phi_Se_at_threshold_0 = SD_of_Phi_Se_at_threshold_0,
           SD_of_Phi_Fp_at_threshold_0 = SD_of_Phi_Fp_at_threshold_0,
           ## Between-study true params:
           n_total_possible_thresholds_per_test = n_total_possible_thresholds_per_test,
           ## Within-study true params:
           Se_for_current_study_at_threshold_0_list = Se_for_current_study_at_threshold_0_list,
           Sp_for_current_study_at_threshold_0_list = Sp_for_current_study_at_threshold_0_list,
           thresholds_for_all_tests_for_current_study_array_list = thresholds_for_all_tests_for_current_study_array_list,
           Se_per_study_all_tests_all_thresholds_list = Se_per_study_all_tests_all_thresholds_list,
           Sp_per_study_all_tests_all_thresholds_list = Sp_per_study_all_tests_all_thresholds_list,
           prev_true_observed_list = prev_true_observed_list,
           Se_per_study_ref = Se_per_study_ref,
           Sp_per_study_ref = Sp_per_study_ref
           # ## 
           # Sigma_nd_true_observed_list = Sigma_nd_true_observed_list,
           # Sigma_d_true_observed_list = Sigma_d_true_observed_list,
           # ##
           # Phi_Se_observed_list = Phi_Se_observed_list,
           # Phi_Fp_observed_list = Phi_Fp_observed_list,
           # ##
           # prev_true_observed_list = prev_true_observed_list,
           # Se_true_observed_list = Se_true_observed_list,
           # Sp_true_observed_list = Sp_true_observed_list,
           # ##
           # true_correlations_observed_vec_list = true_correlations_observed_vec_list,
           # observed_table_probs_list = observed_table_probs_list,
           # true_estimates_observed_list = true_estimates_observed_list,
           # observed_cell_counts_list = observed_cell_counts_list
         ))
         
   
 }
 
 
 
 
 
 
  
  
 convert_to_aggregate_counts <- function(y_list, 
                                         n_studies, 
                                         n_thr) {
   
         n_non_diseased <- n_diseased <- numeric(n_studies)
         x_non_diseased <- x_diseased <- matrix(NA, n_studies, n_thr)
         
         for (s in 1:n_studies) {
           
               study_data <- y_list[[s]]
               disease_status <- study_data[,1] # Assuming column 1 is disease status
               test_results <- study_data[,2] # Assuming column 2 is test results
               
               # Get n_diseased and n_non_diseased 
               n_diseased[s] <- sum(disease_status == 1)
               n_non_diseased[s] <- sum(disease_status == 0)
               
               # Get counts for each threshold
               for (k in 1:n_thr) {
                 x_diseased[s, k] <- sum(test_results[disease_status == 1] > k)
                 x_non_diseased[s, k] <- sum(test_results[disease_status == 0] > k)
               }
               
         }
         
         return(list(
           n_diseased = n_diseased,
           n_non_diseased = n_non_diseased, 
           x_diseased = x_diseased,
           x_non_diseased = x_non_diseased
         ))
         
 }
  
  
 
 
 
 




 



apply_thr_missingness <- function(  agg_data_cumulative, 
                                    studies_subset_vec,
                                    missing_thr_subset_vec) { 
  
  agg_data_cumulative$x_diseased[studies_subset_vec, missing_thr_subset_vec] <- 999
  agg_data_cumulative$x_non_diseased[studies_subset_vec, missing_thr_subset_vec] <- 999
  
  return(agg_data_cumulative)
  
}







convert_cumulative_to_category <- function(cumulative_matrix) {
  
      # Get dimensions and create output matrix with one extra column
      n_rows <- nrow(cumulative_matrix)
      n_cols <- ncol(cumulative_matrix)
      category_matrix <- matrix(999, nrow = n_rows, ncol = n_cols + 1)
      
      for (i in 1:n_rows) {
        # Get current row
        row_data <- cumulative_matrix[i, ]
        
        # First category is same as first cumulative count
        category_matrix[i, 1] <- row_data[1]
        
        # Process remaining categories
        for (j in 1:n_cols) {
          current_val <- row_data[j]
          
          # If at last column, this number goes in the last category
          if (j == n_cols) {
            if (current_val != 999) {
              category_matrix[i, n_cols + 1] <- current_val
            }
            next
          }
          
          # Get next value
          next_val <- row_data[j + 1]
          
          # If either current or next value is 999, we can't calculate this category
          if (current_val == 999 || next_val == 999) {
            category_matrix[i, j + 1] <- 999
            next
          }
          
          # Calculate category count as difference between current and next cumulative count
          category_matrix[i, j + 1] <- current_val - next_val
        }
      }
      
      return(category_matrix)
  
}












categorical_to_individual <- function( categorical_matrix, 
                                       binary_disease_indicator) {
  
        # Initialize list to store results
        results <- list()
        n_studies <- nrow(categorical_matrix)
        n_categories <- ncol(categorical_matrix)
        
        for (study in 1:n_studies) {
          # Get counts for current study
          study_counts <- categorical_matrix[study, ]
          
          # Initialize vectors for this study
          study_id <- vector()
          category_values <- vector()
          
          # For each category
          for (cat in 1:n_categories) {
            count <- study_counts[cat]
            
            if (count == 999) {
              # For missing data, create one observation with 999
              study_id <- c(study_id, study)
              category_values <- c(category_values, 999)
            } else if (count > 0) {
              # For non-missing data, replicate the category value count times
              study_id <- c(study_id, rep(study, count))
              category_values <- c(category_values, rep(cat, count))
            }
          }
          
          # Create data frame for this study
          study_data <- data.frame(
            study_id = study_id,
            group = binary_disease_indicator,  # 1 for diseased, 0 for non-diseased
            value = category_values
          )
          
          results[[study]] <- study_data
        }
        
        # Combine all studies into one data frame
        final_data <- do.call(rbind, results)
        rownames(final_data) <- NULL  # Reset row names
        
        return(final_data)
  
}










  
  
  
  