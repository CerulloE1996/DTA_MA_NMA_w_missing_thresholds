





## read in data
data_Baron <- read.csv(file.path(getwd(), "Baron_et_al_data", "Baron et al validation data_2.csv"))

data <- data_Baron %>% 
  mutate(id = 1:nrow(data)) %>%
  dplyr::select(id,
                age,
                currentminid, 
                language, # covariate (accuracy may differ between groups as did in original study)
                children, # use this as "population"/"herd" variable 
                cesd1,
                cesd2,
                cesd3,
                cesd4,
                cesd5,
                cesd6,
                cesd7,
                cesd8,
                cesd9,
                cesd10,
                totalcesd,
                phq1,
                phq2,
                phq3,
                phq4,
                phq5,
                phq6,
                phq7,
                phq8,
                phq9,
                totalphq)


data <- tibble(data)


data

data$id

#View(data)
# check variables 
table(data$language) # no missings for language covariate 
# (author said language, ethnicity or religion would affect accuracy - but only looking at language for this analysis)
table(data$children) # no missings
table(data$currentminid) # binary reference test, no missings 
table(data$cesd1)
table(data$cesd2)
table(data$cesd3)
table(data$cesd4) # 1 missing
table(data$cesd5)
table(data$cesd6)
table(data$cesd7) # 1 missing
table(data$cesd8)
table(data$cesd9) 
table(data$cesd10) # 2 missing
filter(data, cesd4=="")
filter(data, cesd7=="")
filter(data, cesd10=="")
# 3 total with missing CES-D items
# all CES-D item's 
table(data$totalcesd) # no empty categories for CES-D (some with only 1-5 patients, more severe scores)
table(data$phq1) # 1 missing
table(data$phq2)
table(data$phq3)
table(data$phq4)
table(data$phq5)
table(data$phq6) # 3 missing
table(data$phq7)
table(data$phq8) # 1 missing
table(data$phq9) # 1 missing
filter(data, phq1=="")
filter(data, phq6=="")
filter(data, phq8=="")
filter(data, phq9=="")
# 6 total with missing PHQ-9 items 
table(data$totalphq) # no empty categories for PHQ-9 (some with only 1-5 patients, more severe scores)

# see how many patients have missing items for PHQ-9 or CES-D
filter(data, cesd1 == "" | 
         cesd2 == "" | 
         cesd3 == "" | 
         cesd4 == "" | 
         cesd5 == "" |
         cesd6 == "" | 
         cesd7 == "" | 
         cesd8 == "" | 
         cesd9 == "" | 
         cesd10 == "" |
         phq1 == "" |
         phq2 == "" |
         phq3 == "" |
         phq4 == "" |
         phq5 == "" |
         phq6 == "" |
         phq7 == "" |
         phq8 == "" |
         phq9 == "" )
# 9 in total with missing items 

# code variables 
data_2 <- data %>% 
  mutate(phq1_num = case_when(phq1 == "Not at all" ~ 1,
                              phq1 == "Several days" ~ 2,
                              phq1 == "More than half the days" ~ 3,
                              phq1 == "Nearly every day" ~ 4,
                              phq1 == "" ~ 999), 
         phq2_num = case_when(phq2 == "Not at all" ~ 1,
                              phq2 == "Several days" ~ 2,
                              phq2 == "More than half the days" ~ 3,
                              phq2 == "Nearly every day" ~ 4,
                              phq2 == "" ~ 999), 
         phq3_num = case_when(phq3 == "Not at all" ~ 1,
                              phq3 == "Several days" ~ 2,
                              phq3 == "More than half the days" ~ 3,
                              phq3 == "Nearly every day" ~ 4,
                              phq3 == "" ~ 999),
         phq4_num = case_when(phq4 == "Not at all" ~ 1,
                              phq4 == "Several days" ~ 2,
                              phq4 == "More than half the days" ~ 3,
                              phq4 == "Nearly every day" ~ 4,
                              phq4 == "" ~ 999),
         phq5_num = case_when(phq5 == "Not at all" ~ 1,
                              phq5 == "Several days" ~ 2,
                              phq5 == "More than half the days" ~ 3,
                              phq5 == "Nearly every day" ~ 4,
                              phq5 == "" ~ 999),
         phq6_num = case_when(phq6 == "Not at all" ~ 1,
                              phq6 == "Several days" ~ 2,
                              phq6 == "More than half the days" ~ 3,
                              phq6 == "Nearly every day" ~ 4,
                              phq6 == "" ~ 999),
         phq7_num = case_when(phq7 == "Not at all" ~ 1,
                              phq7 == "Several days" ~ 2,
                              phq7 == "More than half the days" ~ 3,
                              phq7 == "Nearly every day" ~ 4,
                              phq7 == "" ~ 999),
         phq8_num = case_when(phq8 == "Not at all" ~ 1,
                              phq8 == "Several days" ~ 2,
                              phq8 == "More than half the days" ~ 3,
                              phq8 == "Nearly every day" ~ 4,
                              phq8 == "" ~ 999),
         phq9_num = case_when(phq9 == "Not at all" ~ 1,
                              phq9 == "Several days" ~ 2,
                              phq9 == "More than half the days" ~ 3,
                              phq9 == "Nearly every day" ~ 4,
                              phq9 == "" ~ 999)) %>% 
  mutate(cesd1_num = case_when(cesd1 == "Rarely or none of the time (less than 1 day)" ~ 1,
                               cesd1 == "Some or little of the time (1-2 days)" ~ 2,
                               cesd1 == "Occasionally or a moderate amount of time (3-4 days)" ~ 3,
                               cesd1 == "All of the time (5-7 days)" ~ 4, 
                               cesd1 == "" ~ 999),
         cesd2_num = case_when(cesd2 == "Rarely or none of the time (less than 1 day)" ~ 1,
                               cesd2 == "Some or little of the time (1-2 days)" ~ 2,
                               cesd2 == "Occasionally or a moderate amount of time (3-4 days)" ~ 3,
                               cesd2 == "All of the time (5-7 days)" ~ 4,
                               cesd2 == "" ~ 999),
         cesd3_num = case_when(cesd3 == "Rarely or none of the time (less than 1 day)" ~ 1,
                               cesd3 == "Some or little of the time (1-2 days)" ~ 2,
                               cesd3 == "Occasionally or a moderate amount of time (3-4 days)" ~ 3,
                               cesd3 == "All of the time (5-7 days)" ~ 4,
                               cesd3 == "" ~ 999),
         cesd4_num = case_when(cesd4 == "Rarely or none of the time (less than 1 day)" ~ 1,
                               cesd4 == "Some or little of the time (1-2 days)" ~ 2,
                               cesd4 == "Occasionally or a moderate amount of time (3-4 days)" ~ 3,
                               cesd4 == "All of the time (5-7 days)" ~ 4,
                               cesd4 == "" ~ 999),
         cesd5_num = case_when(cesd5 == "Rarely or none of the time (less than 1 day)" ~ 1,
                               cesd5 == "Some or little of the time (1-2 days)" ~ 2,
                               cesd5 == "Occasionally or a moderate amount of time (3-4 days)" ~ 3,
                               cesd5 == "Occacionally or a moderate amount of time (3-4 days)" ~ 3,
                               cesd5 == "All of the time (5-7 days)" ~ 4,
                               cesd5 == "" ~ 999),
         cesd6_num = case_when(cesd6 == "Rarely or none of the time (less than 1 day)" ~ 1,
                               cesd6 == "Some or little of the time (1-2 days)" ~ 2,
                               cesd6 == "Occasionally or a moderate amount of time (3-4 days)" ~ 3,
                               cesd6 == "All of the time (5-7 days)" ~ 4,
                               cesd6 == "" ~ 999),
         cesd7_num = case_when(cesd7 == "Rarely or none of the time (less than 1 day)" ~ 1,
                               cesd7 == "Some or little of the time (1-2 days)" ~ 2,
                               cesd7 == "Occasionally or a moderate amount of time (3-4 days)" ~ 3,
                               cesd7 == "All of the time (5-7 days)" ~ 4,
                               cesd7 == "" ~ 999),
         cesd8_num = case_when(cesd8 == "Rarely or none of the time (less than 1 day)" ~ 1,
                               cesd8 == "Some or little of the time (1-2 days)" ~ 2,
                               cesd8 == "Occasionally or a moderate amount of time (3-4 days)" ~ 3,
                               cesd8 == "Occacionally or a moderate amount of time (3-4 days)" ~ 3,
                               cesd8 == "All of the time (5-7 days)" ~ 4,
                               cesd8 == "" ~ 999),
         cesd9_num = case_when(cesd9 == "Rarely or none of the time (less than 1 day)" ~ 1,
                               cesd9 == "Some or little of the time (1-2 days)" ~ 2,
                               cesd9 == "Occasionally or a moderate amount of time (3-4 days)" ~ 3,
                               cesd9 == "All of the time (5-7 days)" ~ 4,
                               cesd9 == "" ~ 999),
         cesd10_num = case_when(cesd10 == "Rarely or none of the time (less than 1 day)" ~ 1,
                                cesd10 == "Some or little of the time (1-2 days)" ~ 2,
                                cesd10 == "Occasionally or a moderate amount of time (3-4 days)" ~ 3,
                                cesd10 == "All of the time (5-7 days)" ~ 4,
                                cesd10 == "" ~ 999)) %>%
  mutate(totalcesd_plus_1 = totalcesd + 1, 
         totalphq_plus_1 = totalphq + 1,
         children_num = case_when(children == "No" ~ 0, 
                                  children == "Yes" ~ 1), 
         language_num = case_when(language == "Afrikaans" ~ 1,
                                  language == "Xhosa" ~ 2,
                                  language == "Zulu" ~ 3), 
         MINI_num = case_when(currentminid == "No depression" ~ 0, 
                              currentminid == "Current depressive episode" ~ 1))



# CES-D
max(data_2$totalcesd_plus_1)
length(unique(data_2$totalcesd_plus_1)) # no empty categories (except for scores > 30 / 31 for plus_1)
## table(data_2$totalcesd_plus_1)

data_2$totalcesd_plus_1

# PHQ-9
max(data_2$totalphq_plus_1)
length(unique(data_2$totalphq_plus_1)) # we have 2 empty categories here 
data_2$totalphq_plus_1
## table(data_2$totalphq_plus_1) # categories 26 and 27 are empty but 28 is not -
# we will merge 25 and 28 together and consider this category as "24 (25 for plus_1) or more"
# no need to create an empty category for scores > 27 (or 28 for plus_1)
# therefore we will have 25 - 1 = 24 thresholds for the PHQ-9 to estimate



# response / test data 
y     <- array(c(data_2$MINI_num,
                 data_2$totalphq_plus_1,
                 data_2$totalcesd_plus_1), 
               dim = c(nrow(data_2),3))


y_tibble <- tibble(REF_MINI = y[,1],
                   PHQ_9    = y[,2],
                   CES_D    = y[,3])

y_tibble


y_tibble_nd <- dplyr::filter(y_tibble, REF_MINI == 0)
y_tibble_d  <- dplyr::filter(y_tibble, REF_MINI == 1)
##
## SAVE:
##
saveRDS(object = y_tibble, file = "tibble_Baron_et_al.R")
saveRDS(object = y_tibble_nd, file = "tibble_Baron_et_al_nd.R")
saveRDS(object = y_tibble_d, file = "tibble_Baron_et_al_d.R")


##
## Cumulative Baron PHQ-9 data:
##
{
  
        # phq_9_nd <- y_tibble_nd$PHQ_9
        # phq_9_d  <- y_tibble_d$PHQ_9
        
        
        disease_status <- y_tibble$REF_MINI 
        
        ## Get n_total_d and n_total_nd:
        n_total_d   <- sum(disease_status == 1)
        n_total_nd  <- sum(disease_status == 0)
        
        x_nd <- x_d <- n_nd <- n_d <- c()
        Se_true <- Sp_true <- c()
        
        {
          
              n_thr_t <- 26
              n_cat_t <- n_thr_t + 1
              test_results   <- y_tibble$PHQ_9
              
              # Get counts for each threshold
              for (k in 1:n_thr_t) {
                ##
                ## False-negatives (FN's):
                ##
                x_d[k] <- sum(test_results[disease_status == 1] <= k) # Pr(testing NEGATIVE at threshold k)
                n_d[k] <- x_d[k]
                Se_true[k] <- 1.00 - (x_d[k] / n_total_d)
                ##
                ## True negatives (TN's):
                ##
                x_nd[k] <- sum(test_results[disease_status == 0] <= k) # Pr(testing NEGATIVE at threshold k)
                n_nd[k] <- x_nd[k]
                Sp_true[k] <- x_nd[k] / n_total_nd
              }
              
              n_d[n_thr_t + 1]   <- n_total_d
              n_nd[ n_thr_t + 1] <- n_total_nd
          
        }
        
        {
            n_nd_new <- n_nd[2:(n_thr_t + 1)]
            n_d_new  <- n_d[2:(n_thr_t + 1)]
            n_nd <- n_nd_new
            n_d <- n_d_new
        }
        
        Baron_data_list <- list()
        ##
        Baron_data_list$y_tibble <- y_tibble
        ##
        Baron_data_list$x_nd <- x_nd
        Baron_data_list$x_d  <- x_d
        ##
        Baron_data_list$n_nd <- n_nd
        Baron_data_list$n_d  <- n_d
        ##
        Baron_data_list$n_total_nd <- n_total_nd
        Baron_data_list$n_total_d  <- n_total_d
        ##
        Baron_data_list$n_total_nd <- n_total_nd
        Baron_data_list$n_total_nd <- n_total_nd
        ##
        Baron_data_list$Se_true_PHQ_9 <- Se_true
        Baron_data_list$Sp_true_PHQ_9 <- Sp_true
  
}
##
## SAVE:
##
saveRDS(object = Baron_data_list, file = "Baron_data_list.RDS")

c <- 1
s <- 1

{
        Baron_dummy_data_list <- list()
        ##
        Baron_dummy_data_list$n_studies <- 10
        Baron_dummy_data_list$n_thr <- 26
        Baron_dummy_data_list$n_obs_cutpoints <- rep(26, Baron_dummy_data_list$n_studies)
        ##
        ## PRIORS:
        ##
        SD <- 1
        exp(-0.5*(0 - 2*SD)) ; exp(-0.5*(0 + 2*SD))
        exp(+0.5*(0 - 2*SD)) ; exp(+0.5*(0 + 2*SD))
        ##
        Baron_dummy_data_list$prior_beta_mu_mean <- 0.0
        Baron_dummy_data_list$prior_beta_mu_SD   <- 2.0
        Baron_dummy_data_list$prior_beta_SD_mean <- 0.0
        Baron_dummy_data_list$prior_beta_SD_SD   <- 1.0
        ##
        Baron_dummy_data_list$prior_raw_scale_mu_mean <- 0.0
        Baron_dummy_data_list$prior_raw_scale_mu_SD   <- 1.0
        Baron_dummy_data_list$prior_raw_scale_SD_mean <- 0.0
        Baron_dummy_data_list$prior_raw_scale_SD_SD   <- 1.0
        ##
        Baron_dummy_data_list$prior_alpha <- rep(1.0, 27)
        ##
        Baron_dummy_data_list$prior_dirichlet_cat_means_alpha <- rep(1.0, 27)
        Baron_dummy_data_list$prior_dirichlet_cat_SDs_mean    <- rep(0.0, 27)
        Baron_dummy_data_list$prior_dirichlet_cat_SDs_SD      <- rep(0.10, 27)
        ##
        Baron_dummy_data_list$prior_only <- 0
        Baron_dummy_data_list$kappa_lb <- 1.0
        
        {
          
          x_with_missings <-  x <- x <- cutpoint_index <- list()
          for (c in 1:2) {
            x_with_missings[[c]] <- array(-1, dim = c(Baron_dummy_data_list$n_studies, Baron_dummy_data_list$n_thr))
            n[[c]] <- array(-1, dim = c(Baron_dummy_data_list$n_studies, Baron_dummy_data_list$n_thr))
            x[[c]] <- array(-1, dim = c(Baron_dummy_data_list$n_studies, Baron_dummy_data_list$n_thr))
            cutpoint_index[[c]] <- array(-1, dim = c(Baron_dummy_data_list$n_studies, Baron_dummy_data_list$n_thr))
          }
          ##
          for (c in 1:2) {
            for (s in 1:Baron_dummy_data_list$n_studies) {
               cutpoint_index[[c]][s, 1:Baron_dummy_data_list$n_thr] <- seq(from = 1, to = Baron_dummy_data_list$n_thr, by = 1)
            }
          }
          ##
          for (s in 1:Baron_dummy_data_list$n_studies) {
            x_with_missings[[1]][s, 1:Baron_dummy_data_list$n_thr] <- Baron_data_list$x_nd
            x_with_missings[[2]][s, 1:Baron_dummy_data_list$n_thr] <- Baron_data_list$x_d
            ##
            x[[1]][s, 1:Baron_dummy_data_list$n_thr] <- Baron_data_list$x_nd
            x[[2]][s, 1:Baron_dummy_data_list$n_thr] <- Baron_data_list$x_d
            ##
            n[[1]][s, ] <- Baron_data_list$n_nd
            n[[2]][s, 1:Baron_dummy_data_list$n_thr] <- Baron_data_list$n_d
          }
          ##
          Baron_dummy_data_list$x_with_missings <- x_with_missings
          Baron_dummy_data_list$x <- x
          Baron_dummy_data_list$n <- n
          Baron_dummy_data_list$cutpoint_index <- cutpoint_index
        }
}

## Check that n >= x:
x_nd <= n_nd
x_d <= n_d

 


 

 

##
## --------------------  Run Stan model to get "true" parameter values for the PHQ-9 test based on real-world data: -------------------------------------
##


##
## COMPILE Stan model:
##
{
  
      file <- file.path(getwd(), "stan_models", "DTA_MA_Gat_FIXEDthr.stan")
      
      outs <- R_fn_compile_Stan_model(Stan_model_file_path = file, 
                                      CCACHE_PATH = "/usr/bin/ccache")
      
      mod <- outs$mod

}


##
## Initial values:
##
{
      Baron_dummy_init_list <- list()
      ##
      ##
      Baron_dummy_init_list$beta_mu <- 0.0
      Baron_dummy_init_list$beta_SD <- 0.001
      Baron_dummy_init_list$beta_z <-  rep(0.001, Baron_dummy_data_list$n_studies)
      ##
      Baron_dummy_init_list$raw_scale_mu <-  0.0
      Baron_dummy_init_list$raw_scale_SD <-  0.001
      Baron_dummy_init_list$raw_scale_z  <-  rep(0.001, Baron_dummy_data_list$n_studies)
      ##
      Baron_dummy_init_list$C <-   seq(from = -2.0, to = 2.0, length = Baron_dummy_data_list$n_thr)
}


##
## Sample the model:
##
{
        
        seed <- 123
        
        n_chains <- 4
        init_lists_per_chain <- rep(list(Baron_dummy_init_list), n_chains) 
        
        n_burnin <- 500
        n_iter   <- 1000
        
        tictoc::tic()
        
        Stan_mod_sample <- mod$sample( seed = seed,
                                       data = Baron_dummy_data_list,
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



##
## Outputs to obtain "true" raw parameter values for simulation study:
##
{
        try({
          Stan_mod_sample$summary(c("beta_mu")) %>% print(n = 100)
          Stan_mod_sample$summary(c("beta_SD")) %>% print(n = 100)
          ##
          Stan_mod_sample$summary(c("raw_scale_mu")) %>% print(n = 100)
          Stan_mod_sample$summary(c("raw_scale_SD")) %>% print(n = 100)
          ##
          Stan_mod_sample$summary(c("lambda")) %>% print(n = 100)
          
        })
        try({
          Stan_mod_sample$summary(c("C")) %>% print(n = 100)
          Stan_mod_sample$summary(c("C_MU")) %>% print(n = 100)
        })
        # try({
        #   Stan_mod_sample$summary(c("Se")) %>% print(n = 100)
        #   Stan_mod_sample$summary(c("Sp")) %>% print(n = 100)
        # })
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




Se_true_PHQ_9 <- Baron_data_list$Se_true_PHQ_9
Sp_true_PHQ_9 <- Baron_data_list$Sp_true_PHQ_9

Se_medians <- Se$`50%`
Sp_medians <- Sp$`50%`


Se_true_PHQ_9
Se_medians
Se_true_PHQ_9 - Se_medians
avg_abs_diff_Se <- sum(abs(Se_true_PHQ_9 - Se_medians)) / length(Se_medians)
avg_abs_diff_Se

Sp_true_PHQ_9
Sp_medians
Sp_true_PHQ_9 - Sp_medians
avg_abs_diff_Sp <- sum(abs(Sp_true_PHQ_9 - Sp_medians)) / length(Sp_medians)
avg_abs_diff_Sp


avg_abs_diff_Se*100
avg_abs_diff_Sp*100

## Plot true vs. estimates ROC curve to make sure model params "match" the true data:
plot( x = 1.0 - Sp_medians,   y = Se_medians, col = "blue")
lines(x = 1.0 - Sp_true_PHQ_9,   y = Se_true_PHQ_9,     col = "green")









##
## Store + save "true" parameter values to use for simulation:
##
{
    true_mean_PHQ_9_Cerullo_Gat_params_list <- list()
    ##
    true_mean_PHQ_9_Cerullo_Gat_params_list$beta_mu <-  -1.67 
    true_mean_PHQ_9_Cerullo_Gat_params_list$raw_scale_mu  <- -0.0610 
    true_mean_PHQ_9_Cerullo_Gat_params_list$C <- c(-2.00,  -1.53, -1.16,  -0.859, -0.578, 
                                              -0.330, -0.102, 0.0908, 0.264,  0.459, 
                                              0.660,   0.824, 0.944,  1.02,   1.25,
                                              1.33,    1.52,  1.67,   1.98,   2.11,
                                              2.20,    2.45,  2.54,   2.80,   3.05, 3.08)
    ##
    ## SAVE file:
    ##
    saveRDS(object = true_mean_PHQ_9_Cerullo_Gat_params_list, file = "true_mean_PHQ_9_Cerullo_Gat_params_list.RDS")
    
}

# true_mean_PHQ_9_Cerullo_Gat_params_list <- readRDS(file = "true_mean_PHQ_9_Cerullo_Gat_params_list.RDS")


exp(-0.5*(-0.0610))
exp(+0.5*(-0.0610))



ariable  mean median     sd    mad    q5   q95  rhat ess_bulk ess_tail
<chr>    <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
  1 beta_mu  -1.67  -1.67 0.0389 0.0381 -1.73 -1.60  1.00    2559.    3112.
# A tibble: 1 × 10
variable   mean median     sd    mad      q5    q95  rhat ess_bulk ess_tail
<chr>     <dbl>  <dbl>  <dbl>  <dbl>   <dbl>  <dbl> <dbl>    <dbl>    <dbl>
  1 beta_SD  0.0213 0.0170 0.0181 0.0153 0.00147 0.0549  1.00    3135.    1825.
# A tibble: 1 × 10
variable        mean  median     sd    mad     q5     q95  rhat ess_bulk ess_tail
<chr>          <dbl>   <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl>    <dbl>    <dbl>
  1 raw_scale_mu -0.0610 -0.0611 0.0312 0.0314 -0.112 -0.0102  1.00    1621.    2391.
# A tibble: 1 × 10
variable       mean median     sd    mad      q5    q95  rhat ess_bulk ess_tail
<chr>         <dbl>  <dbl>  <dbl>  <dbl>   <dbl>  <dbl> <dbl>    <dbl>    <dbl>
  1 raw_scale_SD 0.0158 0.0126 0.0135 0.0115 0.00108 0.0424  1.00    2796.    1915.
Error in h(simpleError(msg, call)) : 
  error in evaluating the argument 'x' in selecting a method for function 'print': Can't find the following variable(s) in the output: lambda
In addition: There were 50 or more warnings (use warnings() to see the first 50)
# A tibble: 26 × 10
   variable    mean  median     sd    mad      q5     q95  rhat ess_bulk ess_tail
   <chr>      <dbl>   <dbl>  <dbl>  <dbl>   <dbl>   <dbl> <dbl>    <dbl>    <dbl>
 1 C[1]     -2.00   -2.00   0.0335 0.0331 -2.05   -1.94    1.00    1483.    1904.
 2 C[2]     -1.53   -1.53   0.0270 0.0271 -1.58   -1.49    1.00    1743.    2528.
 3 C[3]     -1.16   -1.16   0.0232 0.0232 -1.19   -1.12    1.00    2263.    3011.
 4 C[4]     -0.859  -0.859  0.0214 0.0215 -0.894  -0.824   1.00    3004.    3165.
 5 C[5]     -0.578  -0.578  0.0202 0.0202 -0.612  -0.545   1.00    3896.    3465.
 6 C[6]     -0.330  -0.330  0.0197 0.0198 -0.362  -0.298   1.00    5394.    3686.
 7 C[7]     -0.102  -0.102  0.0202 0.0201 -0.136  -0.0692  1.00    5537.    3390.
 8 C[8]      0.0908  0.0912 0.0209 0.0209  0.0566  0.125   1.00    5237.    3589.
 9 C[9]      0.264   0.264  0.0219 0.0216  0.227   0.299   1.00    4386.    3245.
10 C[10]     0.459   0.459  0.0228 0.0231  0.421   0.496   1.00    3957.    3268.
11 C[11]     0.660   0.660  0.0240 0.0241  0.620   0.699   1.00    3565.    3514.
12 C[12]     0.824   0.824  0.0251 0.0255  0.783   0.865   1.00    3409.    3175.
13 C[13]     0.944   0.944  0.0258 0.0261  0.902   0.986   1.00    3471.    3170.
14 C[14]     1.02    1.02   0.0265 0.0269  0.974   1.06    1.00    3446.    3029.
15 C[15]     1.25    1.25   0.0284 0.0282  1.20    1.30    1.00    4074.    3356.
16 C[16]     1.33    1.33   0.0292 0.0293  1.28    1.38    1.00    4096.    3426.
17 C[17]     1.52    1.52   0.0316 0.0321  1.47    1.58    1.00    4769.    3346.
18 C[18]     1.67    1.67   0.0333 0.0341  1.62    1.72    1.00    5437.    3486.
19 C[19]     1.98    1.98   0.0392 0.0391  1.92    2.05    1.00    7174.    3629.
20 C[20]     2.11    2.11   0.0432 0.0443  2.04    2.18    1.00    7159.    3601.
21 C[21]     2.20    2.20   0.0461 0.0469  2.13    2.28    1.00    7518.    3455.
22 C[22]     2.45    2.45   0.0557 0.0550  2.36    2.55    1.00    7031.    3402.
23 C[23]     2.54    2.54   0.0598 0.0601  2.45    2.64    1.00    7091.    3260.
24 C[24]     2.80    2.80   0.0769 0.0787  2.68    2.93    1.00    7813.    3583.
25 C[25]     3.05    3.04   0.104  0.108   2.88    3.22    1.00    7090.    3492.
26 C[26]     3.08    3.08   0.107  0.110   2.91    3.26    1.00    6777.    3571.
Error in h(simpleError(msg, call)) : 



# 
# 
# 
# 
# 
# 
# ## Dummy identital MA dataset:
# 
# 
# 
# # merge 25 and 28 categories together to make "25 or greater" category
# data_3_ungrouped <- data_2
# data_3_grouped <- mutate(data_2, 
#                          totalphq_plus_1 = case_when(totalphq_plus_1 == 28 ~ 25,
#                                                      totalphq_plus_1 == 25 ~ 25,
#                                                      totalphq_plus_1 %in% c(1:24) ~ totalphq_plus_1)) %>%
#   mutate(totalphq_plus_1 = case_when(totalphq_plus_1 %in% c(1:5) ~ 1,
#                                      totalphq_plus_1 %in% c(6:20) ~  totalphq_plus_1 - 4,
#                                      totalphq_plus_1 %in% c(21:25) ~ 17)) %>%
#   mutate(totalcesd_plus_1 = case_when(totalcesd_plus_1 %in% c(1:5) ~ 1,
#                                       totalcesd_plus_1 %in% c(6:21) ~  totalcesd_plus_1 - 4,
#                                       totalcesd_plus_1 %in% c(22:31) ~ 18))
# 
# 
# # response / test data 
# nrow(data_3_grouped)
# y     <- array(c(data_3_grouped$MINI_num,
#                  data_3_grouped$totalphq_plus_1,
#                  data_3_grouped$totalcesd_plus_1), 
#                dim = c(nrow(data_3_grouped),3))
# 
# # saveRDS(y, file = "depression_data_ordinal_no_covariates.R")
# # 
# # # categorical covariates 
# # X_cat <- array(c(data_3_grouped$children_num, 
# #                  data_3_grouped$language_num), 
# #                dim = c(941, 2))
# # 
# # saveRDS(X_cat, file = "depression_data_ordinal_cat_covariates.R")
# # 
# # # continuous covariate 
# # X_cts <- array(c(data_3_grouped$age),
# #                dim = c(941, 1))
# # 
# # saveRDS(X_cts, file = "depression_data_ordinal_cts_covariates.R")
# 
# data_3_grouped$totalphq_plus_1
# table(data_3_grouped$totalphq_plus_1) #25 categories now 
# length(unique(data_3_grouped$totalphq_plus_1))
# 
# length(unique(data_3_grouped$totalcesd_plus_1))
# 
# 
# 
# data <- tibble(data.frame(cbind(y, X_cat[,2]))) %>%
#   dplyr::rename(MINI = X1, PHQ_9 = X2, CES_D = X3, Grp = X4)
# 
# data_grp3 <- filter(data, Grp == 3)
# 
# # PHQ-9 
# # Se's
# for (i in 1:max(data_grp3$PHQ_9)) {
#   print(round(sum((data_grp3$PHQ_9 > i) & (data_grp3$MINI == 1))/ sum((data_grp3$MINI == 1)),2))
# }
# # Sp's
# for (i in 1:max(data_grp3$PHQ_9)) {
#   print( 1 -  round(sum((data_grp3$PHQ_9 > i) & (data_grp3$MINI == 0)) / sum((data_grp3$MINI == 0)) ,2))
# }
# 
# data_grp2 <- filter(data, Grp == 2)
# 
# for (i in 1:max(data_grp2$PHQ_9)) {
#   print(round(sum((data_grp2$PHQ_9 > i) & (data_grp2$MINI == 1))/ sum((data_grp2$MINI == 1)),2))
# }
# # Sp's
# for (i in 1:max(data_grp2$PHQ_9)) {
#   print( 1 -  round(sum((data_grp2$PHQ_9 > i) & (data_grp2$MINI == 0)) / sum((data_grp2$MINI == 0)) ,2))
# }
# 
# data_grp1 <- filter(data, Grp == 1)
# 
# for (i in 1:max(data_grp1$PHQ_9)) {
#   print(round(sum((data_grp1$PHQ_9 > i) & (data_grp1$MINI == 1))/ sum((data_grp1$MINI == 1)),2))
# }
# # Sp's
# for (i in 1:max(data_grp1$PHQ_9)) {
#   print( 1 -  round(sum((data_grp1$PHQ_9 > i) & (data_grp1$MINI == 0)) / sum((data_grp1$MINI == 0)) ,2))
# }
# 
# 
# 
# 
# 
