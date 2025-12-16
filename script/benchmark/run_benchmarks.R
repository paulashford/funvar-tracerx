# run_benchmarks.R
# Sept 2025
# P.Ashford
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(yardstick))
# options(tibble.print_max = 15, tibble.print_min = 15)
# options(tibble.width = Inf)

# bench_dir	<- 'data/resources/benchmark/original'
bench_dir	<- 'data/resources/benchmark/minimal'
# bench_ver			<- 'v02'

# Classification calculation functions
source('script/benchmark/bench_functions.R')
source('script/benchmark/calc_class_metrics.R')
# Load required datasets
source('script/benchmark/bench_data_paths.R')

cat( paste0("FunVar benchmarking - please wait loading datasets...\n") )
source('script/benchmark/load_bench_datasets.R')

# ---------------------------------------------------
# MAKE CORE DF FOR ALL GENES AND MUTS IN SCOPE
# ---------------------------------------------------
# Using clust res (eg EGFR 858)
cat( paste0("Creating benchmark table...\n") )
df_core_resno  	<- bind_rows(
						unique( select( df_cosmic123_basic_hv, gene, mutation, mutation_res_num) ),
						unique( select( df_clinvarb_basic_hv, gene, mutation, mutation_res_num) ),
						unique( select( df_dbsnp_basic_hv, gene, mutation, mutation_res_num) ),
					)  %>% unique() %>% arrange(gene, mutation_res_num, mutation)

# Check HUGO for coreresno (should - and is - fine given indidividual datasets checked)
df_core_resno_hv	<-	add_hugo_valid(df_core_resno, .check_uprot = FALSE) %>%
					select( gene, mutation, mutation_res_num, external_uniprot_id)

# FunFam annot for all df_core_resno_hv
df_ff_core			<- left_join( df_funfam_ffr,
									df_core_resno_hv,
									by = c('uniprot_acc' = 'external_uniprot_id'),
									keep = FALSE, 
									na_matches = 'never', 
									relationship = 'many-to-many' 
								) %>%
						filter( mutation_res_num >= AA_RANGE_LOW,
								mutation_res_num <= AA_RANGE_HIGH ) 
# filter to just 3D PDB rep available
df_ff_core_cath		<- filter(df_ff_core, REP_SOURCE_ID == 'cath') %>%
							select( gene,
									mutation_res_num,
									uniprot_acc,
									SUPERFAMILY_ID, 
									FUNFAM_NUMBER,
									REP_ID,
									funfam_name = NAME
							) 									
									
# ---------------------------------------------------
# MAKE BENCH DF FOR ALL REFS AND TEST DATASETS
# ---------------------------------------------------
join_columns = c('gene', 'mutation_res_num')
include_cols = NULL
# Create ref set 
df_bench_rs <- 	create_bench_refset( df_ff_core_cath, 
										join_columns = join_columns, 
										include_cols = include_cols, 
										join_relationship = 'many-to-many' ) 

group_cols 		<- c('gene', 'mutation_res_num','uniprot_acc', 'SUPERFAMILY_ID', 'FUNFAM_NUMBER', 'REP_ID', 'funfam_name')
flag_cols 		<- setdiff(colnames(df_bench_rs), group_cols)
# Ensure rows unique within group cols
df_bench_rsu	<- enforce_unique_bench_dataset_with_group_max(df_bench_rs, group_cols, flag_cols)

# ---------------------------------------------------
# ADD 3D PREDICTION FLAGS USING RESIDUE NUMBERS
# Specific mutations are not available for all 3D tools
# Tests use cluster residue predictions matching residues 
# with a colocated CMC or benign variant.
# ---------------------------------------------------
cat( paste0("Adding binary FIE & 3D method flags...\n") )
df_bm <- df_bench_rsu
join_columns <- c('gene', 'mutation_res_num')

# FIEs (TCGA pancancer)
df_bm 	<- add_flag( df_bm,
					unique( select(df_tc_mc3_basic, all_of( c(join_columns, include_cols) ) ) ), 
					join_cols = join_columns,
					flag_name = 'tcga_pan_fie3',
					relationship_join = 'many-to-many'
					)

# GDC: other 3D prediction tools
# GDC: HOTSPOT3D
df_bm	<- 	add_flag( df_bm,
					unique( select(df_gdc_hotspot3d, all_of( c(join_columns, include_cols) ) ) ), 
					join_cols = join_columns,
					flag_name = 'hotspot3d',
					relationship_join = 'many-to-many'
					)
# 3D HOTSPOTS
df_bm	<- 	add_flag( df_bm,
					unique( select(df_gdc_3dhotspots, all_of( c(join_columns, include_cols) ) ) ), 
					join_cols = join_columns,
					flag_name = 'threedhotspots',
					relationship_join = 'many-to-many'
					)

# HOTMAPS
df_bm	<- 	add_flag( df_bm,
					unique( select(df_gdc_hotmaps, all_of( c(join_columns, include_cols) ) ) ), 
					join_cols = join_columns,
					flag_name = 'hotmaps',
					relationship_join = 'many-to-many'
					)

# --------------------------------------------------------
# ADD ACTUAL POS AND ACTUAL NEG COLS; RUN BENCHMARKS
# --------------------------------------------------------
# TRUE FLAG / FALSE FLAG
# 1. Using POS: cosmic CMC; NEG: clinvar-benign or dbSNP-benign
cat( paste0("Adding ACTUAL POS (COSMIC CMC) and NEG (Clinvar benign & dbSNP)...\n") )
df_bm 	<- df_bm %>% 
					mutate( ACTUAL_NEG = ifelse( ( clinvar_benign == 1 | dbSNP_benign == 1), 1, 0) ) %>%
					mutate( ACTUAL_POS = ifelse( cosmic_cmc == 1, 1, 0 ) ) %>% 
					mutate( ACTUAL_NEG = ifelse( ACTUAL_POS == 1 , 0, ACTUAL_NEG ) ) %>% 
					filter( ACTUAL_NEG == 1 | ACTUAL_POS == 1 ) 

# Get all flag columns now on df_bm
flag_cols_bm <- setdiff( colnames(df_bm), group_cols )

# Create sampling df comprised num_samples balanced samples of ACTUAL_POS/ACTUAL_NEG
# each slice has size 25% of total ACTUAL_POS rows
n_samples 			<- 1000
slice_prop_act_pos 	<- 0.25
cat( paste0("Running balanced sampling for mean/sd...(~2 mins)\n") )
df_samples 	<- fn_create_sampling_df(
									mutate( df_bm, across( all_of(flag_cols_bm), ~ as_factor(.x) ) ),
									true_col_var = ACTUAL_POS, 
									num_samples = n_samples, 
									slice_prop = slice_prop_act_pos 
								)

# predicted_class_col = "tcga_pan_fie3"
cat( paste0("\n\n\nBinary classification metrics for 3D methods.\n") )
cat( paste0("Mean/SD estimates from ", 
				as.character(n_samples), " balanced random samples of ", 
				as.character(100 * slice_prop_act_pos), "% of benchmarking table entries marked ACTUAL_POS from COSMIC CMC.\n\n") )

# FIEs TCGA PANCANCER
cat( paste0("TCGA pancancer FIEs (FIE score >=3) [tcga_pan_fie3]\n") )
df_calc_tcga_pan_fie3 		<- get_binary_classification_metrics(df_samples, truth_col = "ACTUAL_POS", predicted_class_col = "tcga_pan_fie3")
get_summary_stats(df_calc_tcga_pan_fie3, .metric, .estimate)
			
# HOTSPOT3D
cat( paste0("\nHotspot3D [hotspot3d]\n") )
df_calc_hotspot3d 		<- get_binary_classification_metrics(df_samples, truth_col = "ACTUAL_POS", predicted_class_col = "hotspot3d")
get_summary_stats(df_calc_hotspot3d, .metric, .estimate)

# HOTMAPS
cat( paste0("\nHotMAPS [hotmaps]\n") )
df_calc_hotmaps 		<- get_binary_classification_metrics(df_samples, truth_col = "ACTUAL_POS", predicted_class_col = "hotmaps")
get_summary_stats(df_calc_hotmaps, .metric, .estimate)

# 3DHOTSPOTS
cat( paste0("\n3D Hotspots [threedhotspots]\n") )
df_calc_threedhotspots 		<- get_binary_classification_metrics(df_samples, truth_col = "ACTUAL_POS", predicted_class_col = "threedhotspots")
get_summary_stats(df_calc_threedhotspots, .metric, .estimate)

cat("\nFinished.\n")



