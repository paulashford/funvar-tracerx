# script/diversity_analysis/create_gene_fie_count_table.R 
# Nov 23

# This is for plotting long-tail distribution.
# Script simplifies/groups/counts the diversity table for use in:
# 	plot/plot_fie_count_gene_dist_long_tail.ipynb

library( tidyverse )

# diversity functions module
source('diversity_functions.R')

# data for analysis
data_dir <- '../../data'

# analysis path
an_dir  <- 'working'

# read gene table
data_file_name <- 'diversity_yj07f_full_table'
data_file_type <- 'gene'  # info used elsewhere!
# data_file_type <- 'funfam'  # info used elsewhere!
data_file_name <- paste0( data_file_name, '_', data_file_type, '.csv' )
data_imported <- read_delim( file.path( data_dir, data_file_name ), delim = ',' ) 

# filter
txps = c( 'N' )  # Tx only
cancers <- c( 'LUSC', 'LUAD' )
atyp <- 'timing'
timings <- c( 'early', 'late', 'unknown'  )   # aka pre- and post-duplication and unk
typedats <- c( 'FIE' )
var_classes <- c( 'Missense_Mutation') 
clonal_filter <- 'none'

# delete any existing analysis groups 
if ( exists('dfa_groups') ){ rm(dfa_groups) }
# create analysis groups for FIE and/or background muts
for ( typedat in typedats ){
  	# Apply filter function
  	dfa_filt <- filter_analysis( data_imported, 
                            	type_dat = typedat,
                            	variant_class = var_classes,
                            	txp = txps, 
								timing_in = timings,
                            	cancer_in = cancers,
                            	remove_cols = c( 'VARIANT_CLASS', 'SCORE_THRESHOLD', 'DRIVERMUT', 'MUTANT_EXPRESSED', 'NUM_ROWS' ),
                            	simpler_labels = TRUE,
                            	pyclone_filt = clonal_filter
                        )
  	# group according to analysis type (timing etc)
  	dfa_group <- group_analysis( dfa_filt, analysis_type = atyp )
  
  	# is this 1st time in loop?
  	if ( !exists('dfa_groups') ) { 
      dfa_groups <- dfa_group 
    } else {
      dfa_groups <- bind_rows( dfa_groups, dfa_group )
    }
}

# write table for plots
write_delim( dfa_groups, file.path( an_dir, 'fies_by_gene_tx_luad_lusc_pre_post_unk.tsv' ), delim = '\t' ) 
