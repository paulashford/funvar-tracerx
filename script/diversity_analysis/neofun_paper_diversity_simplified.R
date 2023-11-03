# Diversity scoring with iNEXT and iNEXT.3D
# 23 02 2023
# Simpler files based on script/v3.1/neofun_v3.1_paper_counts_yj_diversity.sql and diversity VIEW
# filtering is done here, not in various Oracle exports
# See: script/v3.1/neofun_v3.1_paper_counts_yj_diversity.sql
# Also: neofun/script/sql/vw/vw_diversity.sql
# 22 09 23 collated dir see: script/v3.2/diversity_analysis/README_diversity.txt
# + changed var name: root_dir to data_dir
# + changed path locations here
# 01 11 23 collated for submission to git/funvar-tracerx/script/diversity_analysis

# NOTE to summarise groups:
# test <- df_analysis %>%
# 	# group_by(SPECIES_TYPE, DATA_SOURCE, TXP, CANCER_TYPE, OVERALL_TIMING)  %>%
#     group_by(SPECIES_TYPE, DATA_SOURCE, TXP, CANCER_TYPE, OVERALL_TIMING)  %>%
# group_data(test)
# # ----------------------

library( readr )
library( tidyverse )
library( dplyr )
library( tidyr )
library( iNEXT.3D )
library( ggplot2 )

# diversity functions module
source('diversity_functions.R')

# data for analysis
data_dir <- '../../data'

# analysis path
an_dir  <- 'working'

# MAIN DIVERSITY ANALSYSIS SOURCE FILES - fies and missense/syn mutations by gene or funfam (require filtering for FIE, background etc)
data_file_name <- 'diversity_yj07f_full_table'
data_file_type <- 'gene'  # info used elsewhere!
# data_file_type <- 'funfam'  # info used elsewhere!
data_file_name <- paste0( data_file_name, '_', data_file_type, '.csv' )
data_imported <- read_delim( file.path( data_dir, data_file_name ), delim = ',' ) 

# Working copy of analysis df
df_analysis   <- data_imported

# filter variables also used elsewhere
cancers <- c( 'LUAD' )
# cancers <- c( 'LUSC' )
# cancers <- c( 'NSCLC' )  # c( 'LUSC', 'LUAD' )
txps = c( 'N' )  # Tx only
# txps <- c( 'N', 'Y' )

# Simplifies the SPECIES_TYPE field to make plot keys easier to read
simplify_lables <- TRUE

# timing only or duplicated-clonality
# timing
# atyp <- 'timing'
# timings <- c( 'early', 'late'  )  

# duplicated-clonality
atyp <- 'duplicated-clonality'; 
timings <- c( 'early', 'late', 'unknown' )   # NOTE: need unknown in tests with ANY_REGION_DUPLICATRED...

# timing-clonality
# atyp <- 'timing-clonality'
# timings <- c( 'early', 'late'  )  

# FIEs or bckground muts or both
# nfd = not a FIE or a known driver
typedats <- c( 'FIE' )
# typedats <- c( 'nfd' )
# typedats <- c( 'FIE', 'nfd' )

# variant types 
# var_classes <- c( 'Silent', 'Missense_Mutation')
var_classes <- c( 'Missense_Mutation') 


# C/S analysis
# clonal_filter <- 'C'
# clonal_filter <- 'S'
clonal_filter <- 'none'


# delete any existing analysis groups 
if ( exists('dfa_groups') ){ rm(dfa_groups) }
# create analysis groups for FIE and/or background muts
for ( typedat in typedats ){
  # Apply filter function
  dfa_filt <- filter_analysis( df_analysis, 
                            	type_dat = typedat,
                            	variant_class = var_classes,
                            	txp = txps, 
								              timing_in = timings,
                            	cancer_in = cancers,
                            	remove_cols = c( 'VARIANT_CLASS', 'SCORE_THRESHOLD', 'DRIVERMUT', 'MUTANT_EXPRESSED', 'NUM_ROWS' ),
                            	simpler_labels = simplify_lables,
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

# dataset label for types
label_info <- paste0( unlist(typedats), collapse="_and_" )

# pivot for iNEXT
df_pivs <- pivot_analysis( dfa = dfa_groups, analysis_type = atyp, names_separator = '-' )

# ---------------------
# iNEXT runs
# ---------------------
# bootstrap parameters and diversity orders to include
knots <- 200; nboot <- 5000; 
conf <- 0.95
q <- c( 1 );  # Shannon diversity

# convert to (simpler) iNEXT-suitable format
df <- convert_for_iNEXT( df_pivs )

# file naming
label_core <- paste0( 'diversity_', 
                        data_file_name, '_', 
                        atyp, 
                        '_q_', paste( unlist(q), collapse='' ), 
                        '_knots_', knots, 
                        '_nboot_', nboot ) 

# Add cancer type and type of run (just FIEs, FIEs and background etc) and simpler labels
dat_sources <- paste0(txps,collapse=""); if (dat_sources == 'NY'){dat_sources<-'Tx_TCGA'}

label_full <- paste0( label_core, 
                        '_simpler_labels_', paste0( simplify_lables ), '_',
                        label_info, '_', 
                        paste0( cancers, collapse = '-' ), '_',
                        'txp_', dat_sources,
                        '_conf_', conf
                    )

# iNEXT diversity tables & save RDA 
idst <- iNEXT3D( df, 
                    diversity = 'TD', 
                    q = q, 
                    datatype = 'abundance', 
                    knots = knots, 
                    nboot = nboot, 
                    conf = conf
					          , endpoint = 200 )
save( idst, file = file.path( an_dir, paste0( label_full, '.rda' ) ) )

# run/save plots
facetvar = "Order.q" 
# facetvar = "None" 
# log2 plot
logplot <- TRUE; log_type <- 'log2'
# normal plot
logplot <- FALSE; log_type <- ''

inext_plot( idst, 
            atyp = atyp, 
            data_file_type=data_file_type, 
            an_dir=an_dir, 
            label_full=label_full, 
            logplot=logplot, 
            log_type=log_type, 
            facetvar=facetvar
          )



