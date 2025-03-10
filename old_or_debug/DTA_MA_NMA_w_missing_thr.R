


setwd("/home/enzocerullo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")

source("R_fn_load_data_ordinal_MA_LC_MVP_sim.R")

require(dplyr)


n_studies <- 10
N_per_study_mean <- 1000
N_per_study_SD <- 500
assume_perfect_GS <- 1
seed <- 123


# Run simulated data - this simulates data from FIVE (5) diagnostic tests (1 BINARY reference test + 4 ORDINAL index tests)
sim_results <- simulate_binary_and_ordinal_MA_LC_MVP_data(n_studies = n_studies,
                                                          N_per_study_mean = N_per_study_mean,
                                                          N_per_study_SD = N_per_study_SD,
                                                          assume_perfect_GS = assume_perfect_GS,
                                                          seed = seed)

y_list <- sim_results$y_list
sim_results$n_total_possible_thresholds_per_test
str(y_list)

## Now for the first example we will only take one index test, as initially we are evaluating just a "simple" model
## where there's only a single index test with 12 thresholds (so 13 categories:
n_thr <- 12
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



  
  
## First 2 studies only report data for the "middle 8" thresholds (starting from threshold 3, so thr = {3, 4, 5, 6, 7, 8, 9, 10})
studies_subset_vec <- c(1, 2)
missing_thr_subset_vec <- c(1:2, 11:12)
agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
                                                                        studies_subset_vec = studies_subset_vec,
                                                                        missing_thr_subset_vec = missing_thr_subset_vec)
## The next 2 studies only report at the "middle 4" thresholds (starting at threshold 5): so thr = {5, 6, 7, 8}:
studies_subset_vec <- c(3, 4)
missing_thr_subset_vec <- c(1:4, 9:12)
agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
                                                                        studies_subset_vec = studies_subset_vec,
                                                                        missing_thr_subset_vec = missing_thr_subset_vec)
## The next two studies only report data at only a single threshold - threshold #7:
studies_subset_vec <- c(5, 6)
missing_thr_subset_vec <- c(1:6, 8:12)
agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
                                                                        studies_subset_vec = studies_subset_vec,
                                                                        missing_thr_subset_vec = missing_thr_subset_vec)
## The next two studies only report data at 4 non-adjacent thresholds - so thr = {3, 5, 7, 9}:
studies_subset_vec <- c(7, 8)
missing_thr_subset_vec <- c(1, 2, 4, 6, 8, 10, 12)
agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
                                                                        studies_subset_vec = studies_subset_vec,
                                                                        missing_thr_subset_vec = missing_thr_subset_vec)
## And finally, the last 2 studies report data at ALL thresholds:
#  -----  Do nothing here - as no missing thresholds for study #10!

## Now let's look at the overall % of missing thresholds:
total_missing_thr <- sum(agg_data_cumulative_with_missing_thresholds$x_diseased == 999) ; total_missing_thr
prop_missing_thr <- total_missing_thr / (n_studies * n_thr) ; prop_missing_thr
## Just under half the threshold data is missing (46.67%). 

 
# # Now final_data contains the format needed for the Stan model w/ some thresholds marked as 999 (missing):
# final_data_cumulative$n_diseased
# final_data_cumulative$n_non_diseased
# ## Data:
# final_data_cumulative$x_diseased
# final_data_cumulative$x_non_diseased
# 
# 
# final_data_cumulative <- 


# ## Diseased group:
# x_diseased_cumulative <- agg_data_cumulative_with_missing_thresholds$x_diseased
# print(paste("x_diseased_cumulative = ")) ; print(x_diseased_cumulative)
# x_diseased_categorical <- convert_cumulative_to_category(x_diseased_cumulative)
# print(paste("x_diseased_categorical = ")) ; print(x_diseased_categorical)
# ## Non-diseased group:
# x_non_diseased_cumulative <- agg_data_cumulative_with_missing_thresholds$x_non_diseased
# print(paste("x_non_diseased_cumulative = ")) ; print(x_non_diseased_cumulative)
# x_non_diseased_categorical <- convert_cumulative_to_category(x_non_diseased_cumulative)
# print(paste("x_non_diseased_categorical = ")) ; print(x_non_diseased_categorical)


x_non_diseased_categorical <- convert_cumulative_to_category(final_data_cumulative$x_non_diseased)

sum(x_diseased_categorical == 999) / sum(prod(dim(x_diseased_categorical)))







y_diseased_df <- categorical_to_individual(  categorical_matrix = x_diseased_categorical, 
                                             binary_disease_indicator = 1)

y_non_diseased_df <- categorical_to_individual(  categorical_matrix = x_non_diseased_categorical, 
                                                 binary_disease_indicator = 0)

y_df <- tibble(rbind(y_non_diseased_df, y_diseased_df))
print(y_df)
























