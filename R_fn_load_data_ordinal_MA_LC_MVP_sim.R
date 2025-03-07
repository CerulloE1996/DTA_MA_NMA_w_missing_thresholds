
 
 
 
 
 # 
 # n_studies <- 10
 # N_per_study_mean <- 500
 # N_per_study_SD <- 1000
 # assume_perfect_GS <- 1
 # seed <- 123
 



## -| ------------------ R function to simulate a meta-analysis dataset for (binary + ordinal) LC_MVP  --------------------------------------------
 
 
simulate_binary_and_ordinal_MA_LC_MVP_data <- function(   n_studies = 25,
                                                          N_per_study_mean = 500,
                                                          N_per_study_SD   = 500,
                                                          assume_perfect_GS = 1,
                                                          seed = 123,
                                                          true_Mean_prev = 0.20) {
   
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
     
     # N_per_study_SD <- 250 ##  N_per_study_mean / 2
     N_per_study_vec <- round(TruncatedNormal::rtnorm(n = n_studies, mu = N_per_study_mean, sd = N_per_study_SD, lb = 100), 0)
  
     # ## True values for locations in both groups:
     location_nd <- rep(-1.00, n_tests)
     location_d  <- rep(+1.00, n_tests)
     
     ## For ref test:
     if (assume_perfect_GS == 1) { 
         location_nd[1] <- -100
         location_d[1]  <- +100
     }
     
     ## Set true values for disease prevelance:
     true_Mean_probit_prev <- qnorm(true_Mean_prev)
     true_SD_probit_prev <- 0.25
     true_probit_prev_per_study <- rnorm(n = n_studies, mean = true_Mean_probit_prev, sd = true_SD_probit_prev)
     true_prev_per_study <- pnorm(true_probit_prev_per_study)
     
     ## Set true calues for the between-study heterogenity params:
     scale_nd <- rep(0.25, n_tests)
     scale_d  <- rep(0.50, n_tests)
     
     ## Threshold info (test 1 is binary so only 1 threshold @0 on the "latent" scale)::
     n_total_obs_thr_per_test <- c(1, 5, 10, 25, 50)
     max_threshold_across_all_tests <- max(n_total_obs_thr_per_test)
     threshold_per_study_array <- array(NA, dim = c(n_studies, max_threshold_across_all_tests))
     Mean_of_thr_for_all_tests_array <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
     SD_of_thr_for_all_tests_array <- array(0.125, dim = c(n_tests, max_threshold_across_all_tests)) ## between-study heterogenity for thresholds 
     
     ## Set the "true values" of the mean threshold between studies for each test (on the * latent * scale):
     ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds! 
     ## (i.e. the TRUE threshold parameters are "homogenous between classes")
     ## Binary reference test:
     Mean_thresholds_for_test_1 <- rep(NA, n_total_obs_thr_per_test[1])
     Mean_thresholds_for_test_1 <- c(0.0)
     Mean_of_thr_for_all_tests_array[1, 1:n_total_obs_thr_per_test[1]] <- Mean_thresholds_for_test_1
     ## Test 2:
     Mean_thresholds_for_test_2 <- rep(NA, n_total_obs_thr_per_test[2])
     Mean_thresholds_for_test_2 <- c(-1.75, -1.00, -0.50, 0.25, 1.25)
     length(Mean_thresholds_for_test_2)
     Mean_of_thr_for_all_tests_array[2, 1:n_total_obs_thr_per_test[2]] <- Mean_thresholds_for_test_2
     ## Test 3:
     Mean_thresholds_for_test_3 <- rep(NA, n_total_obs_thr_per_test[3])
     Mean_thresholds_for_test_3 <- c(-1.75, -1.50, -1.00, -0.80, -0.50, 0.25, 1.00, 1.50, 1.80, 2.20)
     length(Mean_thresholds_for_test_3)
     Mean_of_thr_for_all_tests_array[3, 1:n_total_obs_thr_per_test[3]] <- Mean_thresholds_for_test_3
     ## Test 4:
     Mean_thresholds_for_test_4 <- rep(NA, n_total_obs_thr_per_test[4])
     Mean_thresholds_for_test_4 <- c(-2.50, -2.20, -2.00, -1.90, -1.75, -1.50, -1.35, -1.00, -0.80, -0.50, 
                                     -0.25, -0.10, +0.00, +0.10, +0.25, +0.40, +0.80, +1.00, +1.50, +1.80, 
                                     +2.20, +2.40, +2.50, +2.60, +2.70)
     length(Mean_thresholds_for_test_4)
     Mean_of_thr_for_all_tests_array[4, 1:n_total_obs_thr_per_test[4]] <- Mean_thresholds_for_test_4
     ## Test 5:
     Mean_thresholds_for_test_5 <- rep(NA, n_total_obs_thr_per_test[5])
     Mean_thresholds_for_test_5 <- c(-3.00, -2.90, -2.80, -2.70, -2.60, -2.50, -2.35, -2.20, -2.00, -1.90, 
                                     -1.75, -1.50, -1.35, -1.15, -1.00, -0.90, -0.70, -0.60, -0.50, -0.25,
                                     -0.10, +0.00, +0.10, +0.20, +0.35, +0.50, +0.60, +0.70, +0.85, +0.95, 
                                     +1.05, +1.20, +1.30, +1.35, +1.45, +1.65, +1.75, +1.80, +1.90, +2.00, 
                                     +2.10, +2.20, +2.30, +2.40, +2.50, +2.60, +2.70, +2.75, +2.80, +2.85)
     length(Mean_thresholds_for_test_5)
     Mean_of_thr_for_all_tests_array[5, 1:n_total_obs_thr_per_test[5]] <- Mean_thresholds_for_test_5
     
     s <- 1
     t <- 2
     
     100*(1.0 - pnorm((Mean_thresholds_for_test_3 - location_d[3])/exp(scale_nd[3])))
     true_Se_OVERALL
     
     pnorm(location_d[3] - Mean_thresholds_for_test_3 )
     
     Se_for_current_study_at_threshold_0_list <- Sp_for_current_study_at_threshold_0_list <- Fp_for_current_study_at_threshold_0_list <- list()
     Se_per_study_all_tests_all_thresholds_list <- Sp_per_study_all_tests_all_thresholds_list <- list()
     prev_true_observed_list <- list()
     Se_per_study_ref <- Sp_per_study_ref <- list()
     thresholds_for_all_tests_for_current_study_array_list <- list()
     y_list <- list()
     y_df_list <- list()
     ##
     n_TP_at_current_threshold_OVERALL <- array(0.0, dim = c(n_tests, max_threshold_across_all_tests))
     n_FP_at_current_threshold_OVERALL <- array(0.0, dim = c(n_tests, max_threshold_across_all_tests))
     n_pos_OVERALL <- 0 
     n_neg_OVERALL <- 0 
     n_total_OVERALL <- 0
     Se_OVERALL_all_tests_all_thresholds <- array(0.0, dim = c(n_tests, max_threshold_across_all_tests))
     Sp_OVERALL_all_tests_all_thresholds <- array(0.0, dim = c(n_tests, max_threshold_across_all_tests))
     
     
     true_DGP_one_m_Se <- pnorm(location_d - Mean_of_thr_for_all_tests_array)
     true_DGP_Sp       <- pnorm(location_nd - Mean_of_thr_for_all_tests_array)
     
   
   for (s in 1:n_studies) {
     
             # ii_dataset <- ii_dataset + 1 ## increment counter
             
             N <- N_per_study_vec[s]
             true_prev <- true_prev_per_study[s]
             
             ## Generate "true" thresholds for study s:
             ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds!
             ## (i.e. the TRUE threshold parameters are "homogenous between classes")
             thr_for_all_tests_for_current_study_array <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
             for (t in 2:n_tests) {
               thr_for_all_tests_for_current_study_array[t, 1:n_total_obs_thr_per_test[t]] <- rnorm(n = n_total_obs_thr_per_test[t], 
                                                                                                           mean = Mean_of_thr_for_all_tests_array[t, 1:n_total_obs_thr_per_test[t]], 
                                                                                                           sd = SD_of_thr_for_all_tests_array[t, 1:n_total_obs_thr_per_test[t]])
             }
             
             thresholds_for_all_tests_for_current_study_array_list[[s]] <- thr_for_all_tests_for_current_study_array
             
             location_nd_study_s <- rnorm(n = n_tests, mean = location_nd,  sd = scale_nd)
             location_d_study_s  <- rnorm(n = n_tests, mean = location_d,   sd = scale_d)
             
             if (assume_perfect_GS == TRUE) {
                 rho1 <- 0.00
             } else { 
                 rho1 <- 0.20
             }
             
             Sigma_highly_varied <- matrix(c(1,     rho1,  rho1,     rho1,     rho1,
                                             rho1,  1.00,  0.50,     0.20,     0.10,
                                             rho1,  0.50,  1.00,     0.40,     0.40,
                                             rho1,  0.20,  0.40,     1.00,     0.70,
                                             rho1,  0.10,  0.40,     0.70,     1.00),
                                             nrow = n_tests, 
                                             ncol = n_tests)
        
             {
                   Sigma_d <- Sigma_highly_varied
                   Sigma_nd <-  0.5 * Sigma_highly_varied ## Corr(D-) is HALF that of the Corr(D+)
                   diag(Sigma_nd) <- rep(1, n_tests)
             }
             
             L_Sigma_d   <- t(chol(Sigma_d)) # PD check (fails if not PD)
             L_Sigma_nd  <- t(chol(Sigma_nd))   ## BayesMVP:::Rcpp_Chol(Sigma_nd) # PD check (fails if not PD)
             
             d_ind <- sort(rbinom(n = N, size = 1, prob = true_prev))
             ##
             n_pos <- sum(d_ind)
             n_pos_OVERALL <- n_pos_OVERALL + n_pos
             ##
             n_neg <- N - sum(d_ind)
             n_neg_OVERALL <- n_neg_OVERALL + n_neg
             ##
             n_total <- n_pos + n_neg
             n_total_OVERALL <- n_total_OVERALL + n_total
             
             ## Simulate the underlying "latent" continuous test responses which are distributed multivariate normal:
             latent_cts_results_neg <- LaplacesDemon::rmvn(n = n_neg, mu = location_nd_study_s, Sigma = Sigma_nd)
             latent_cts_results_pos <- LaplacesDemon::rmvn(n = n_pos, mu = location_d_study_s,  Sigma = Sigma_d)
             latent_results <- rbind(latent_cts_results_neg, latent_cts_results_pos)
             
             ## Now simulate the BINARY data for the reference test (i.e. the "imperfect gold standard" - test # 1):
             results_pos <- array(NA, dim = c(n_pos, n_tests))
             results_neg <- array(NA, dim = c(n_neg, n_tests))
             y <-     array(NA, dim = c(n_neg + n_pos, n_tests))
             ## BINARY reference test results:
             y[, 1] <- ifelse(latent_results[, 1] > 0.0, 1, 0)
             
             ## Now simulate the ORDINAL data for tests 2-5:
             ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds!
             ## (i.e. the TRUE threshold parameters are "homogenous between classes")
             ## Hence we can just use the "latent_results" array directly (no need to split into +'ve and -'ve)
             ## Negatives / Sp:
             for (t in 2:n_tests) {
               
                       n_thr <- n_total_obs_thr_per_test[t] 
                       
                       threshold_lower <- -1000
                       threshold_upper <- head(thr_for_all_tests_for_current_study_array[t, 1:n_thr], 1)
                       y[, t] <- ifelse(  (latent_results[, t] > threshold_lower) & (latent_results[, t] < threshold_upper),
                                                1, 
                                                y[, t])
                    
                       for (k in 2:(n_thr)) {
                           threshold_lower <- thr_for_all_tests_for_current_study_array[t, k - 1] ; threshold_lower
                           threshold_upper <- thr_for_all_tests_for_current_study_array[t, k]     ; threshold_upper
                           y[, t] <- ifelse(  (latent_results[, t] > threshold_lower) & (latent_results[, t] < threshold_upper), 
                                                     k, 
                                                     y[, t])
                       }
                       
                       threshold_lower <- tail(thr_for_all_tests_for_current_study_array[t, 1:n_thr], 1)  
                       threshold_upper <- +1000
                       y[, t] <- ifelse(     (latent_results[, t] > threshold_lower) & (latent_results[, t] < threshold_upper),
                                                   n_thr + 1, 
                                                   y[, t])
                 
             }
             y
             
             df <- dplyr::tibble(y, latent_results, d_ind)
             df_pos <- dplyr::filter(df, d_ind == 1)
             df_neg <- dplyr::filter(df, d_ind == 0)

             # prev:
             prev_true_observed <-  print(round(sum(d_ind)/N, 3))
             prev_true_observed_list[[s]] <- prev_true_observed

             ## Compute OVERALL observed TRUE Se:
             Se_current_study_all_tests_all_thresholds <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
             for (t in 2:n_tests) {
               n_thr <- n_total_obs_thr_per_test[t] 
               for (k in 1:n_thr) {
                 n_TP_at_current_threshold <- sum(df_pos$y[, t] > k)
                 Se_current_study_all_tests_all_thresholds[t, k] <- n_TP_at_current_threshold/n_pos
                 ## Increment comulative counters:
                 n_TP_at_current_threshold_OVERALL[t, k] <- n_TP_at_current_threshold_OVERALL[t, k] + n_TP_at_current_threshold
               }
             }
             Se_per_study_all_tests_all_thresholds_list[[s]] <- Se_current_study_all_tests_all_thresholds
             
             ## Compute OVERALL observed TRUE Fp / Sp:
             Sp_current_study_all_tests_all_thresholds <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
             for (t in 2:n_tests) {
               n_thr <- n_total_obs_thr_per_test[t] 
               for (k in 1:n_thr) {
                 n_FP_at_current_threshold <- sum(df_neg$y[, t] > k)
                 Fp_at_threshold_k <- n_FP_at_current_threshold/n_neg
                 Sp_at_threshold_k <- 1.0 - Fp_at_threshold_k
                 Sp_current_study_all_tests_all_thresholds[t, k] <- Sp_at_threshold_k
                 ## Increment comulative counters:
                 n_FP_at_current_threshold_OVERALL[t, k] <- n_FP_at_current_threshold_OVERALL[t, k] + n_FP_at_current_threshold
               }
             }
             Sp_per_study_all_tests_all_thresholds_list[[s]] <- Sp_current_study_all_tests_all_thresholds
             
             ## For reference test:
             t <- 1
             ## Se:
             Phi_Se_ref <- qnorm(sum(df_pos$results[ ,t])/nrow(df_pos))
             Se_ref <- pnorm(Phi_Se_ref)
             Se_per_study_ref[[s]] <- Se_ref
             
             ## Sp:
             Phi_Fp_ref <-  qnorm( 1.0 - ((nrow(df_neg) - sum(df_neg$results[ ,t]))/nrow(df_neg))  )
             Fp_ref <- pnorm(Phi_Fp_ref)
             Sp_ref <- 1.0 - Fp_ref
             Sp_per_study_ref[[s]] <- Sp_ref

             y_list[[s]] <- y
             y_df_list[[s]] <- data.frame(y)
             
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
             # true_estimates_observed <-  c(Sigma_nd_true_observed[upper.tri(Sigma_nd_true_observed )], 
             #                             Sigma_d_true_observed[upper.tri(Sigma_d_true_observed )], Sp_true_observed,  Se_true_observed, prev_true_observed ,
             #                               true_corrs_observed_vec, observed_table_probs_vec, NA, NA)
             # true_estimates  <-  c(Sigma_nd[upper.tri(Sigma_nd )],  Sigma_d[upper.tri(Sigma_d )], 1 - true_Fp_vec,  true_Se_vec, true_prev  ,
             #                       rep(NA, length(true_corrs_observed_vec)), rep(NA, length(observed_table_probs_vec)), NA, NA)
             #
             # ## Print some key quantities:
             # print(paste("N = ", N))
             # print(paste("prev_true_observed = ", prev_true_observed))
             # print(paste("Se_true_observed = ", Se_true_observed))
             # print(paste("Sp_true_observed = ", Sp_true_observed))
             #
             # ## Populate lists:
             # y_binary_list[[ii_dataset]] <- y
             # ##
             # Sigma_nd_true_observed_list[[ii_dataset]] <- Sigma_nd_true_observed
             # Sigma_d_true_observed_list[[ii_dataset]] <- Sigma_d_true_observed
             # ##
             # Phi_Se_observed_list[[ii_dataset]] <- Phi_Se_observed_vec
             # Phi_Fp_observed_list[[ii_dataset]] <- Phi_Fp_observed_vec
             #
             # ##
             # true_correlations_observed_vec_list[[ii_dataset]] <- true_corrs_observed_vec
             # observed_table_probs_list[[ii_dataset]] <- observed_table_probs_vec
             # true_estimates_observed_list[[ii_dataset]] <- true_estimates_observed
             # observed_cell_counts_list[[ii_dataset]] <- observed_cell_counts
     
     
   }
   

       ## Compute OVERALL observed TRUE Se:
       for (t in 2:n_tests) {
         n_thr <- n_total_obs_thr_per_test[t]
         for (k in 1:n_thr) {
           Se_OVERALL_all_tests_all_thresholds[t, k] <- n_TP_at_current_threshold_OVERALL[t, k]/n_pos_OVERALL
         }
       }
       # Se_OVERALL_all_tests_all_thresholds <- Reduce(`+`, Se_per_study_all_tests_all_thresholds_list) / length(Se_per_study_all_tests_all_thresholds_list)
      
       
       ## Compute OVERALL observed TRUE Fp / Sp:
       for (t in 2:n_tests) {
         n_thr <- n_total_obs_thr_per_test[t]
         for (k in 1:n_thr) {
           Sp_OVERALL_all_tests_all_thresholds[t, k] <- 1.0 - (n_FP_at_current_threshold_OVERALL[t, k]/n_neg_OVERALL)
         }
       }
       # Sp_OVERALL_all_tests_all_thresholds <- Reduce(`+`, Sp_per_study_all_tests_all_thresholds_list) / length(Sp_per_study_all_tests_all_thresholds_list)
       
       y_tibble <- NULL
       y_tibble <- tibble(data.table::rbindlist(y_df_list, idcol = "Study"))
   
   return(list(
     N_per_study_vec = N_per_study_vec,
     y_list = y_list,
     y_df_list = y_df_list,
     y_tibble = y_tibble,
     ## Between-study true params:
     # Mean_Se_at_threshold_0 = Mean_Se_at_threshold_0,
     # Mean_Sp_at_threshold_0 = Mean_Sp_at_threshold_0,
     true_Mean_prev = true_Mean_prev,
     ## Between-study true params:
     # SD_of_Phi_Se_at_threshold_0 = SD_of_Phi_Se_at_threshold_0,
     # SD_of_Phi_Fp_at_threshold_0 = SD_of_Phi_Fp_at_threshold_0,
     ## Between-study true params:
     n_total_obs_thr_per_test = n_total_obs_thr_per_test,
     ## Within-study true params:
     Se_for_current_study_at_threshold_0_list = Se_for_current_study_at_threshold_0_list,
     Sp_for_current_study_at_threshold_0_list = Sp_for_current_study_at_threshold_0_list,
     thresholds_for_all_tests_for_current_study_array_list = thresholds_for_all_tests_for_current_study_array_list,
     Se_per_study_all_tests_all_thresholds_list = Se_per_study_all_tests_all_thresholds_list,
     Sp_per_study_all_tests_all_thresholds_list = Sp_per_study_all_tests_all_thresholds_list,
     prev_true_observed_list = prev_true_observed_list,
     Se_per_study_ref = Se_per_study_ref,
     Sp_per_study_ref = Sp_per_study_ref,
     Se_OVERALL_all_tests_all_thresholds = Se_OVERALL_all_tests_all_thresholds,
     Sp_OVERALL_all_tests_all_thresholds = Sp_OVERALL_all_tests_all_thresholds,
     true_DGP_one_m_Se = true_DGP_one_m_Se,
     true_DGP_Sp = true_DGP_Sp
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
   
         n_total_non_diseased <- n_total_diseased <- numeric(n_studies)
         x_non_diseased <- x_diseased <- matrix(NA, n_studies, n_thr)
         n_non_diseased <- n_diseased <- matrix(NA, n_studies, n_thr + 1)
         Se_per_study <- Sp_per_study <- matrix(NA, n_studies, n_thr)
         
         for (s in 1:n_studies) {
           
               study_data <- y_list[[s]]
               disease_status <- study_data[,1] # Assuming column 1 is disease status
               test_results   <- study_data[,2] # Assuming column 2 is test results
               
               # Get n_total_diseased and n_total_non_diseased 
               n_total_diseased[s] <- sum(disease_status == 1)
               n_total_non_diseased[s] <- sum(disease_status == 0)
               
               # Get counts for each threshold
               for (k in 1:n_thr) {
                 x_diseased[s, k]     <- sum(test_results[disease_status == 1] <= k) # Pr(testing NEGATIVE at threshold k)
                 n_diseased[s, k] <- x_diseased[s, k]
                 Se_per_study[s, k] <- 1.00 - (x_diseased[s, k] / n_total_diseased[s])
                 ##
                 x_non_diseased[s, k] <- sum(test_results[disease_status == 0] <= k) # Pr(testing NEGATIVE at threshold k)
                 n_non_diseased[s, k] <- x_non_diseased[s, k]
                 Sp_per_study[s, k] <- x_non_diseased[s, k] / n_total_non_diseased[s]
               }
               n_diseased[s, n_thr + 1]     <- n_total_diseased[s]
               n_non_diseased[s, n_thr + 1] <- n_total_non_diseased[s]
               # x_diseased[s, n_thr + 1]     <-  n_total_diseased[s]     -  x_diseased[s, n_thr]
               # x_non_diseased[s, n_thr + 1] <-  n_total_non_diseased[s] -  x_non_diseased[s, n_thr]
         }
         
         return(list(
           n_total_diseased = n_total_diseased,
           n_total_non_diseased = n_total_non_diseased, 
           x_diseased = x_diseased,
           x_non_diseased = x_non_diseased,
           n_diseased = n_diseased,
           n_non_diseased = n_non_diseased,
           Se_per_study = Se_per_study,
           Sp_per_study = Sp_per_study
         ))
         
 }
  
  
 
 
 
 



 
 
 
 
 

 



apply_thr_missingness <- function(  agg_data_cumulative, 
                                    studies_subset_vec,
                                    missing_thr_subset_vec,
                                    missing_indicator = -1) { 
  
  agg_data_cumulative$x_diseased[studies_subset_vec, missing_thr_subset_vec] <- missing_indicator
  agg_data_cumulative$x_non_diseased[studies_subset_vec, missing_thr_subset_vec] <- missing_indicator
  
  return(agg_data_cumulative)
  
}



















# 
# 


convert_cumulative_to_category <- function(cumulative_matrix,
                                           missing_indicator = -1) {

  # Get dimensions and create output matrix with same dimensions plus one column
  n_rows <- nrow(cumulative_matrix)
  n_cols <- ncol(cumulative_matrix)
  category_matrix <- matrix(-1, nrow = n_rows, ncol = n_cols + 1)

  for (i in 1:n_rows) {
    # Get current row
    row_data <- cumulative_matrix[i, ]

    # First category is the first cumulative count
    category_matrix[i, 1] <- row_data[1]

    # Second category is second column minus first column
    if (n_cols >= 2) {
      if (!is.na(row_data[1]) && !is.na(row_data[2])) {
        category_matrix[i, 2] <- row_data[2] - row_data[1]
      }
    }

    # Process remaining categories by taking differences
    for (j in 2:(n_cols-1)) {
      if (!is.na(row_data[j]) && !is.na(row_data[j+1])) {
        category_matrix[i, j+1] <- row_data[j+1] - row_data[j]
      }
    }

    # Last category is the total minus the last cumulative count
    # (Assuming we want all categories to sum to the total)
    if (!is.na(row_data[n_cols])) {
      total <- row_data[n_cols]
      category_matrix[i, n_cols+1] <- 0  # Default to zero

      # If we want a true non-zero last category, we would need the actual total N
      # which isn't directly provided in this data format
    }
  }

  return(category_matrix[, 1:n_cols])

}






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
# categorical_to_individual <- function( categorical_matrix, 
#                                        binary_disease_indicator,
#                                        missing_indicator = -1) {
#   
#         # Initialize list to store results
#         results <- list()
#         n_studies <- nrow(categorical_matrix)
#         n_categories <- ncol(categorical_matrix)
#         
#         for (study in 1:n_studies) {
#           # Get counts for current study
#           study_counts <- categorical_matrix[study, ]
#           
#           # Initialize vectors for this study
#           study_id <- vector()
#           category_values <- vector()
#           
#           # For each category
#           for (cat in 1:n_categories) {
#             count <- study_counts[cat]
#             
#             if (count == missing_indicator) {
#               # For missing data, create one observation with missing_indicator
#               study_id <- c(study_id, study)
#               category_values <- c(category_values, missing_indicator)
#             } else if (count > 0) {
#               # For non-missing data, replicate the category value count times
#               study_id <- c(study_id, rep(study, count))
#               category_values <- c(category_values, rep(cat, count))
#             }
#           }
#           
#           # Create data frame for this study
#           study_data <- data.frame(
#             study_id = study_id,
#             group = binary_disease_indicator,  # 1 for diseased, 0 for non-diseased
#             value = category_values
#           )
#           
#           results[[study]] <- study_data
#         }
#         
#         # Combine all studies into one data frame
#         final_data <- do.call(rbind, results)
#         rownames(final_data) <- NULL  # Reset row names
#         
#         return(final_data)
#   
# }
# 
# 
# 
# 
# 





  
  
  
  