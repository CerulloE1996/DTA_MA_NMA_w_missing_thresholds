
 
 
 
 
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
     n_thr_per_test <- c(1, 5, 10, 25, 50)
     max_threshold_across_all_tests <- max(n_thr_per_test)
     threshold_per_study_array <- array(NA, dim = c(n_studies, max_threshold_across_all_tests))
     Mean_of_thr_for_all_tests_array <- array(NA, dim = c(n_tests, max_threshold_across_all_tests))
     SD_of_thr_for_all_tests_array <- array(0.125, dim = c(n_tests, max_threshold_across_all_tests)) ## between-study heterogenity for thresholds 
     
     ## Set the "true values" of the mean threshold between studies for each test (on the * latent * scale):
     ## NOTE: currently assuming that the diseased and non-diseased latent classes have the SAME set of thresholds! 
     ## (i.e. the TRUE threshold parameters are "homogenous between classes")
     ## Binary reference test:
     Mean_thresholds_for_test_1 <- rep(NA, n_thr_per_test[1])
     Mean_thresholds_for_test_1 <- c(0.0)
     Mean_of_thr_for_all_tests_array[1, 1:n_thr_per_test[1]] <- Mean_thresholds_for_test_1
     ## Test 2:
     Mean_thresholds_for_test_2 <- rep(NA, n_thr_per_test[2])
     Mean_thresholds_for_test_2 <- c(-1.75, -1.00, -0.50, 0.25, 1.25)
     length(Mean_thresholds_for_test_2)
     Mean_of_thr_for_all_tests_array[2, 1:n_thr_per_test[2]] <- Mean_thresholds_for_test_2
     ## Test 3:
     Mean_thresholds_for_test_3 <- rep(NA, n_thr_per_test[3])
     Mean_thresholds_for_test_3 <- c(-1.75, -1.50, -1.00, -0.80, -0.50, 0.25, 1.00, 1.50, 1.80, 2.20)
     length(Mean_thresholds_for_test_3)
     Mean_of_thr_for_all_tests_array[3, 1:n_thr_per_test[3]] <- Mean_thresholds_for_test_3
     ## Test 4:
     Mean_thresholds_for_test_4 <- rep(NA, n_thr_per_test[4])
     Mean_thresholds_for_test_4 <- c(-2.50, -2.20, -2.00, -1.90, -1.75, -1.50, -1.35, -1.00, -0.80, -0.50, 
                                     -0.25, -0.10, +0.00, +0.10, +0.25, +0.40, +0.80, +1.00, +1.50, +1.80, 
                                     +2.20, +2.40, +2.50, +2.60, +2.70)
     length(Mean_thresholds_for_test_4)
     Mean_of_thr_for_all_tests_array[4, 1:n_thr_per_test[4]] <- Mean_thresholds_for_test_4
     ## Test 5:
     Mean_thresholds_for_test_5 <- rep(NA, n_thr_per_test[5])
     Mean_thresholds_for_test_5 <- c(-3.00, -2.90, -2.80, -2.70, -2.60, -2.50, -2.35, -2.20, -2.00, -1.90, 
                                     -1.75, -1.50, -1.35, -1.15, -1.00, -0.90, -0.70, -0.60, -0.50, -0.25,
                                     -0.10, +0.00, +0.10, +0.20, +0.35, +0.50, +0.60, +0.70, +0.85, +0.95, 
                                     +1.05, +1.20, +1.30, +1.35, +1.45, +1.65, +1.75, +1.80, +1.90, +2.00, 
                                     +2.10, +2.20, +2.30, +2.40, +2.50, +2.60, +2.70, +2.75, +2.80, +2.85)
     length(Mean_thresholds_for_test_5)
     Mean_of_thr_for_all_tests_array[5, 1:n_thr_per_test[5]] <- Mean_thresholds_for_test_5
     
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
               thr_for_all_tests_for_current_study_array[t, 1:n_thr_per_test[t]] <- rnorm(n = n_thr_per_test[t], 
                                                                                                           mean = Mean_of_thr_for_all_tests_array[t, 1:n_thr_per_test[t]], 
                                                                                                           sd = SD_of_thr_for_all_tests_array[t, 1:n_thr_per_test[t]])
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
               
                       n_thr <- n_thr_per_test[t] 
                       
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
               n_thr <- n_thr_per_test[t] 
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
               n_thr <- n_thr_per_test[t] 
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
         n_thr <- n_thr_per_test[t]
         for (k in 1:n_thr) {
           Se_OVERALL_all_tests_all_thresholds[t, k] <- n_TP_at_current_threshold_OVERALL[t, k]/n_pos_OVERALL
         }
       }
       # Se_OVERALL_all_tests_all_thresholds <- Reduce(`+`, Se_per_study_all_tests_all_thresholds_list) / length(Se_per_study_all_tests_all_thresholds_list)
      
       
       ## Compute OVERALL observed TRUE Fp / Sp:
       for (t in 2:n_tests) {
         n_thr <- n_thr_per_test[t]
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
     n_thr_per_test = n_thr_per_test,
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
                                         n_tests,
                                         n_thr) {
   
   
           
           n_total_nd <- n_total_d <- numeric(n_studies)
           x_nd_list <- x_d_list <- list()
           n_nd_list <- n_d_list <- list()
           Se_per_study_list <- Sp_per_study_list <- list()
           
           ## Initialise lists:
           n_index_tests <- n_tests - 1
           for (t in 2:n_tests) { ## skip the BINARY reference test (test 1)
             
                   n_thr_t <- n_thr[t]
                   n_cat_t <- n_thr_t + 1
                   ##
                   x_nd_list[[t - 1]] <- matrix(NA, n_studies, n_thr_t)
                   x_d_list[[t - 1]]  <- matrix(NA, n_studies, n_thr_t)
                   ##
                   n_nd_list[[t - 1]] <- matrix(NA, n_studies, n_cat_t)
                   n_d_list[[t - 1]]  <- matrix(NA, n_studies, n_cat_t)
                   ##
                   Se_per_study_list[[t - 1]] <-  matrix(NA, n_studies, n_thr_t)
                   Sp_per_study_list[[t - 1]] <-  matrix(NA, n_studies, n_thr_t)
                   
           }
           
           for (s in 1:n_studies) {
             
                     study_data <- y_list[[s]]
                     disease_status <- study_data[,1] # Col 1 is disease status / reference test
                     
                     ## Get n_total_d and n_total_nd:
                     n_total_d[s]  <- sum(disease_status == 1)
                     n_total_nd[s] <- sum(disease_status == 0)
                   
                     for (t in 2:n_tests) { ## skip the BINARY reference test (test 1)
                       
                             n_thr_t <- n_thr[t]
                             n_cat_t <- n_thr_t + 1
                             test_results   <- study_data[, t]
                             
                             # Get counts for each threshold
                             for (k in 1:n_thr_t) {
                                   ##
                                   ## False-negatives (FN's):
                                   ##
                                   x_d_list[[t - 1]][s, k] <- sum(test_results[disease_status == 1] <= k) # Pr(testing NEGATIVE at threshold k)
                                   n_d_list[[t - 1]][s, k] <- x_d_list[[t - 1]][s, k]
                                   Se_per_study_list[[t - 1]][s, k] <- 1.00 - (x_d_list[[t - 1]][s, k] / n_total_d[s])
                                   ##
                                   ## True negatives (TN's):
                                   ##
                                   x_nd_list[[t - 1]][s, k] <- sum(test_results[disease_status == 0] <= k) # Pr(testing NEGATIVE at threshold k)
                                   n_nd_list[[t - 1]][s, k] <- x_nd_list[[t - 1]][s, k]
                                   Sp_per_study_list[[t - 1]][s, k] <- x_nd_list[[t - 1]][s, k] / n_total_nd[s]
                             }
                             n_d_list[[t - 1]][s,  n_thr_t + 1] <- n_total_d[s]
                             n_nd_list[[t - 1]][s, n_thr_t + 1] <- n_total_nd[s]
                     }
           }
           
         
         return(list(
           n_total_nd = n_total_nd,
           n_total_d  = n_total_d, 
           x_nd_list = x_nd_list,
           x_d_list  = x_d_list,
           n_nd_list = n_nd_list,
           n_d_list  = n_d_list,
           Se_per_study_list = Se_per_study_list,
           Sp_per_study_list = Sp_per_study_list
         ))
         
 }
  
  
 
 
 
 


# Function to create "indicator_test_in_study" indicator matrix
create_test_in_study_for_NMA <- function(   n_studies, 
                                            n_tests, 
                                            prob_present = 0.75, 
                                            min_index_tests_per_study = 1
) {
   
           ## At least 1 index tests per study -> at least 2 total tests per study:
           min_tests_per_study <- min_index_tests_per_study + 1
           
           ## Initialize the matrix:
           indicator_index_test_in_study <- matrix(0, nrow = n_studies, ncol = n_tests)
           
           ## Ref test is ALWAYS included: 
           reference_test_prob <- 1.0
           
           ## For each study, decide which tests are present:
           for (s in 1:n_studies) {
             
                 ## Give reference test higher probability of being included
                 present_probs <- rep(prob_present, n_tests)
                 ## Assuming test 1 is ref test:
                 present_probs[1] <- reference_test_prob  
                 
                 ## Randomly determine which tests are present (1 = present, 0 = absent)
                 test_presence <- rbinom(n_tests, 1, present_probs)
                 
                 ## Ensure at least min_tests_per_study tests are present
                 while (sum(test_presence) < min_tests_per_study) {
                   ## If not enough tests, add one more randomly
                   missing_tests <- which(test_presence == 0)
                   if (length(missing_tests) > 0) {
                     add_test <- sample(missing_tests, 1)
                     test_presence[add_test] <- 1
                   } else {
                     break  # All tests already included
                   }
                 }
                 
                 indicator_index_test_in_study[s, ] <- test_presence
             
           }
           
           ## Optional: Ensure each test appears in at least one study:
           for (t in 1:n_tests) {
               
               if (sum(indicator_index_test_in_study[, t]) == 0) {
                 # If a test doesn't appear in any study, add it to a random study
                 random_study <- sample(1:n_studies, 1)
                 indicator_index_test_in_study[random_study, t] <- 1
               }
               
           }
           
           indicator_index_test_in_study <- indicator_index_test_in_study[1:n_studies, 2:n_tests]
           
           ## Count number of (observed / non-missing): 
           n_index_tests_per_study <- rowSums(indicator_index_test_in_study)
           
           return(list( indicator_index_test_in_study = indicator_index_test_in_study,
                        n_index_tests_per_study = n_index_tests_per_study))
   
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



















categorical_to_individual <- function( categorical_matrix,
                                       binary_disease_indicator,
                                       missing_indicator = -1) {

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

            if (count == missing_indicator) {
              # For missing data, create one observation with missing_indicator
              study_id <- c(study_id, study)
              category_values <- c(category_values, missing_indicator)
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

        return(tibble(final_data))

}









# Create Box-Cox transformed fits for different lambda values
# Function to create Box-Cox grid values
box_cox_grid <- function(x, 
                         lambda) {
  
  if (abs(lambda) < 0.001) {
    return(log(x))
  } else {
    return((x^lambda - 1)/lambda)
  }
  
}




individual_obs_tibble <- x_individual
group_name <- "non-diseased"
study_index <- 1

# 
# 
# # Check what study_id values exist
# unique_studies <- unique(individual_obs_tibble$study_id)
# print("Available study IDs:")
# print(unique_studies)
# 
# # Check what group values exist
# unique_groups <- unique(individual_obs_tibble$group)
# print("Available group values:")
# print(unique_groups)



## Function for log-logistic density:
dloglogistic <- function(x, 
                         alpha, 
                         beta) {
  return((alpha/beta) * (x/beta)^(alpha-1) / (1 + (x/beta)^alpha)^2)
}



## Function for Box-Cox density
box_cox_density <- function(x, 
                            lambda,
                            data) {
  # Only process positive values (Box-Cox requires positive values)
  if(any(x <= 0)) return(rep(0, length(x)))
  
  # Transform the data
  transformed_data <- box_cox_transform(data, lambda)
  
  # Fit normal to transformed data
  mean_transformed <- mean(transformed_data)
  sd_transformed <- sd(transformed_data)
  
  # Calculate density with Jacobian adjustment
  density_values <- dnorm(box_cox_transform(x, lambda), 
                          mean = mean_transformed, 
                          sd = sd_transformed) * x^(lambda-1)
  
  return(density_values)
}




{
      
      ## 
      ## lambda <- -0.5:
      ##
      box_cox_neg05 <- function(x) {
            return((x^(-0.5) - 1.0)/(-0.5))
      }
      bc_neg05_jacobian <- function(x) { 
            return(x^(-1.5))
      }
      bc_neg05_density <- function(x) {
        
            unajusted_dens <- dnorm(box_cox_neg05(x), 
                                    mean = mean(x), 
                                    sd   = sd(x))
            
            return(unajusted_dens * bc_05_jacobian(x))
        
      }
      ## 
      ## lambda <- +0.5:
      ##
      box_cox_05 <- function(x) {
            return((x^0.5 - 1.0)/0.5)
      }
      bc_05_jacobian <- function(x) { 
            return(x^(-0.5))
      }
      bc_05_density <- function(x) {
        
        unajusted_dens <- dnorm(box_cox_05(x), 
                                mean = mean(x), 
                                sd   = sd(x))
        
        return(unajusted_dens * bc_05_jacobian(x))
        
      }
      ## 
      ## lambda > 0:
      ##
      box_cox_pos <- function(x, 
                              lambda) {
            return((x^lambda - 1.0)/lambda)
      }
      bc_pos_jacobian <- function(x,
                                  lambda) { 
            return(x^(1.0 - lambda))
      }
      bc_pos_density <- function(x,
                                 lambda) {
            
            trans_x <- box_cox_pos(x = x, 
                                   lambda = lambda)
            
            unajusted_dens <- dnorm(trans_x, 
                                    mean = mean(x), 
                                    sd   = sd(x))
            
            Jacobian <- bc_pos_jacobian(x = x,
                                        lambda = lambda)
            
            return(unajusted_dens * Jacobian)
        
      }

}




box_cox_transform <- function(x, 
                              lambda) {
  
      if(abs(lambda) < 0.001) return(log(x))
      return((x^lambda - 1)/lambda)
  
}

bc_density <- function(x, 
                       lambda, 
                       observations) {
  
      ## Transform original data:
      transformed_data <- box_cox_transform(observations, 
                                            lambda)
      transformed_mean <- mean(transformed_data)
      transformed_sd   <- sd(transformed_data)
      
      ## Apply Jacobian adjustment:
      jacobian <- x^(lambda-1)
      
      ## Calculate density:
      dnorm(box_cox_transform(x, lambda), 
            mean = transformed_mean, 
            sd = transformed_sd) * jacobian
  
}


study_indexes <- c(1:n_studies)

##
## Function to fit distributions and create plots:
##
plot_distribution_fits <- function(individual_obs_tibble, 
                                   n_studies,
                                   study_indexes = NULL, 
                                   group_name = "non-diseased") {
  
                  if (group_name == "non-diseased") { 
                      group_value <- 0.0
                  } else { 
                      group_value <- 1.0
                  }

                  if (is.null(study_indexes)) { 
                     study_indexes <- seq(from = 1, to = n_studies)
                  }

                  individual_obs_tibble
                  df_filtered <- dplyr::filter(individual_obs_tibble, 
                                               study_id %in% study_indexes, 
                                               group == as.numeric(group_value)) %>% 
                                 print(n = 100)
                    
                  ## Extract test values
                  observations <- df_filtered$value
                  
                  ## Calculate bin width based on data range
                  bin_width <- max(1, (max(observations) - min(observations)) / 30)
                  
                  ## Log-normal (need to handle zeros if present)
                  obs_for_log <- observations
                  if(any(obs_for_log <= 0)) obs_for_log <- observations + 1 # Shift if needed
                  
                  log_normal_fit <- try(MASS::fitdistr(obs_for_log, "lognormal"))
                  if(!inherits(lognormal_fit, "try-error")) {
                    lognormal_meanlog <- lognormal_fit$estimate["meanlog"]
                    lognormal_sdlog <- lognormal_fit$estimate["sdlog"]
                  } else {
                    lognormal_meanlog <- mean(log(obs_for_log))
                    lognormal_sdlog <- sd(log(obs_for_log))
                  }
                  
                  ## Log-logistic (need to handle zeros if present)
                  # Fit log-logistic parameters
                  # For log-logistic, we'll estimate shape and scale parameters
                  # using method of moments from log-transformed data
                  log_data <- log(obs_for_log)
                  mu_log   <- mean(log_data)
                  sd_log   <- sd(log_data)
                  
                  ## Convert to log-logistic parameters (approximation)
                  ## Shape parameter (alpha) related to the variability
                  alpha <- pi/(sd_log * sqrt(3))
                  ## Scale parameter (beta) related to the median
                  beta <- exp(mu_log)
                  # ##
                  # ## Create box-cox densities:
                  # ##
                  # bc_neg1_density   <- bc_density(lambda = -1.0, x = observations)
                  # bc_neg05_density  <- bc_density(lambda = -0.5, x = observations)
                  # ## For lambda > 0:
                  # bc_05_density  <- bc_density(lambda = 0.5, x = observations)
                  # bc_1_density   <- bc_density(lambda = 1.0, x = observations)
                  # bc_15_density  <- bc_density(lambda = 1.5, x = observations)
                  # bc_2_density   <- bc_density(lambda = 2.0, x = observations)
                  # # bc_25_density  <- bc_density(lambda = 2.5, x = observations)
                  # bc_3_density   <- bc_density(lambda = 3.0, x = observations)
                  ##
                  ## Plot title:
                  ##
                  if (length(study_indexes == n_studies)) { 
                        plot_title <- paste("Distribution of observed", group_name, "data for ALL", n_studies, "studies")
                  } else { 
                        plot_title <- paste("Distribution of observed", group_name, "data for all studies", study_index)
                  }
                  ##
                  ## Create plot:
                  ##
                  plot <- ggplot(data.frame(value = observations), aes(x = value)) +
                          ##
                          ## Histogram of OBSERVED data:
                          ##
                          geom_histogram(aes(y = after_stat(density)), binwidth = bin_width, 
                                         fill = "lightblue", color = "darkblue", alpha = 0.8) +
                          ##
                          ## Add density estimate:
                          ##
                          geom_density(color = "black", size = 2.0) +
                          ##
                          ## Add log-normal density (if data is positive):
                          ##
                          stat_function(fun = function(x) dlnorm(x,
                                                                 meanlog = lognormal_meanlog, 
                                                                 sdlog   = lognormal_sdlog),
                                        color = "green", size = 1.0, linetype = "dashed") + 
                          ##
                          ## Add log-logistic density:
                          ##
                          stat_function(fun = function(x) dloglogistic(x, 
                                                                       alpha, 
                                                                       beta),
                                        color = "blue", size = 1.0, linetype = "dashed") + 
                          ##
                          ## Add Box-Cox densities:
                          ##
                          stat_function(fun = function(x) bc_density(x, lambda = -1.0, observations = observations), color = "darkblue", size = 1.5) +
                          stat_function(fun = function(x) bc_density(x, lambda = -0.5, observations = observations), color = "darkred",  size = 1.5) +
                          ##
                          stat_function(fun = function(x) bc_density(x, lambda = +0.5,  observations = observations), color = "orange", size = 1.5) +
                          stat_function(fun = function(x) bc_density(x, lambda = +1.0,  observations = observations), color = "green",  size = 1.5) +
                          stat_function(fun = function(x) bc_density(x, lambda = +1.5,  observations = observations), color = "purple", size = 1.5) +
                          stat_function(fun = function(x) bc_density(x, lambda = +2.0,  observations = observations), color = "red",    size = 1.5) +
                          ## 
                          ## Formatting:
                          ##
                          labs(title = plot_title,
                               subtitle = paste("n =", length(observations)),
                               x = "Test Value", 
                               y = "Density") +
                          theme_minimal() +
                          ##
                          ## Add legends:
                          ##
                          annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.90, label = "Log-normal (dashed line)",   color = "green") +
                          annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.85, label = "Log-logistic (dashed line)", color = "blue") + 
                          ## Box-cox:
                          annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.75, label = "Box-Cox (λ=-1.0)", color = "darkblue") +
                          annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.70, label = "Box-Cox (λ=-0.5)", color = "darkred") + 
                          ## Box-cox:
                          annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.65, label = "Box-Cox (λ=+0.5)", color = "orange") + 
                          annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.60, label = "Box-Cox (λ=+1.0)", color = "green") + 
                          annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.55, label = "Box-Cox (λ=+1.5)", color = "purple") + 
                          annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.50, label = "Box-Cox (λ=+2.0)", color = "red")
                
                    print(plot)
                        
                    return(plot)
  
}


























  
  
  
  
