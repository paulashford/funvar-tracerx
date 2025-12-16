# script/benchmarking/calc_class_metrics.R
# Calculate classification metrics inc MCC, precision, recall, F1, accuracy and others
# uses tidymodels/yardstick library [https://yardstick.tidymodels.org/reference]
# 24 05 2024
# P. Ashford
get_binary_classification_metrics <- function(df_test, truth_col = "truth", predicted_class_col = "predicted") {
    require(dplyr)
    require(yardstick)
    return( 
        bind_rows(
            # accuracy
            accuracy(df_test, 
                truth = {{ truth_col }}, 
                estimate = {{ predicted_class_col }}),
            # precision and recall
            precision(df_test, 
                truth = {{ truth_col }}, 
                estimate = {{ predicted_class_col }}),
            recall(df_test, 
                truth = {{ truth_col }}, 
                estimate = {{ predicted_class_col }}),
            # f1 and mcc
            f_meas(df_test, 
                truth = {{ truth_col }}, 
                estimate = {{ predicted_class_col }}),
            mcc(df_test, 
                truth = {{ truth_col }}, 
                estimate = {{ predicted_class_col }}),
            # PPV and NPV
            ppv(df_test, 
                truth = {{ truth_col }}, 
                estimate = {{ predicted_class_col }}),
            npv(df_test, 
                truth = {{ truth_col }}, 
                estimate = {{ predicted_class_col }})
        )
    )
}

# Sample function for slicing equal samples 
fn_equal_slice <- function(df, group_var, n_size) {
    require(dplyr)
    df_out  <- df %>%
        group_by( {{group_var}} ) %>%
        slice_sample(n = n_size) %>% 
        ungroup()
  return(df_out)
}

# create sampling dataframe
fn_create_sampling_df <- function(df, true_col_var, num_samples = 100, slice_prop = 0.25) {
    require(dplyr)
    n_slice <- round( slice_prop * nrow( filter( df, {{ true_col_var }} == 1 ) ), 0 )
  
    # Initialise sampling df
    df_samples <-   slice_sample( 
                        df %>% mutate(slice_num=0, .before = "gene"),
						n = 0
    	            )
    # Take num_samples balanced samples of df of size n_slice
	for(this_slice in 1:num_samples) {
        # Sample slice
	    this_sample = fn_equal_slice(df, {{ true_col_var }}, n_slice) %>%
	        mutate(slice_num = this_slice)
	    # Add to sampling df
        df_samples <- bind_rows(df_samples, this_sample)
	}
    # Return grouped by slice_num for calculation of sampling statistics
    return( df_samples %>% group_by(slice_num) )
}

# summary stats
get_summary_stats <- function(df, group_col, est_col) {
    df_calc_stat <- df %>%
        group_by( {{ group_col }} ) %>%
            summarise(	
                est_mean = mean( {{ est_col }} ), 
                est_sd = sd( {{ est_col }} ), 
                est_meansd_l = mean( {{ est_col }} ) - sd( {{ est_col }} ),
                est_meansd_h = mean( {{ est_col }} ) + sd( {{ est_col }} ),
                est_median = median( {{ est_col }} ),
                n = n()
        )
    return(df_calc_stat)
}

