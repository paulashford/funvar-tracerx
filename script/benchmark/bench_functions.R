# bench_functions.R
# P. Ashford
# Sept 2025
# benchmarking utility functions - see also: script/benchmark/calc_class_metrics.R

# Read as tibble from flat file with specific defaults
read_wrap <- function(dataset_name, data_path, delim = '\t', colnames = TRUE, coltypes = NULL, trimws = TRUE, lazy = FALSE, guess_max_rows = 100000, show_col_types = TRUE){
	require(readr)
	file_name <- file.path(data_path, dataset_name)
	
	if (delim == ','){
		return( read_csv( file_name, 
							col_names = colnames, 
							col_types = coltypes, 
							trim_ws = trimws, 
							lazy = lazy, 
							guess_max = guess_max_rows,
							show_col_types = FALSE 
						) 
				) 
	} else {
		return( read_delim( file_name, 
							col_names = colnames, 
							col_types = coltypes, 
							delim = delim, 
							trim_ws = trimws, 
							lazy = lazy, 
							guess_max = guess_max_rows,
							show_col_types = FALSE 
						) 
				) 
	}
}

# Adds mutation_res_num column (eg 224) given a mutation_res (eg A224)
fn_add_res_num	<- function(df){
	require(dplyr)
	require(tidyr)
	df	<- df %>%
		separate_wider_regex( 	cols = 'mutation_res', 
								c( mutation_res_letter = "^[A-Z]{1}", mutation_res_num = "[0-9]+" ), 
								too_few = "align_start" , cols_remove = FALSE ) %>% 
				select(-mutation_res_letter) 
	
	df$mutation_res_num <- as.numeric(df$mutation_res_num)		
	return(df)	
}

# Create new benchmark reference set by adding binary flags for each reference data type.
# many-many joins are to deal with either mutation or mutation residue number joins (with non-redundancy then enforced elsewhere)
create_bench_refset <- function( df_core_grn, 
										join_columns = c('gene', 'mutation_res_num'), 
										include_cols = c('mutation'), 
										join_relationship = 'many-to-many' ) {
		require(dplyr)
		require(tidyr)
		
		# Initialise tibble
		df_bench_refset <- {{ df_core_grn }}
	
		# Flag all input entries 1 (for convenience/validation checks)
		df_b0 	<- add_flag( df_bench_refset, 
									unique( select( df_bench_refset, all_of(c(join_columns, include_cols)) ) ), 
									join_cols = join_columns,
									flag_name = 'bench_refset',
									relationship_join = join_relationship
							)

		# clinvar-benign
		df_b1 	<- add_flag( df_bench_refset, 
									unique( select( df_clinvarb_basic, all_of(c(join_columns, include_cols)) ) ), 
									join_cols = join_columns,
									flag_name = 'clinvar_benign',
									relationship_join = join_relationship
							) %>%
					filter(clinvar_benign == 1)

		# COSMIC CMC (somatic mutations in cancer)
		df_b3 	<- add_flag( df_bench_refset,
									unique( select(df_cosmic123_basic, all_of(c(join_columns, include_cols)) ) ), 
									join_cols = join_columns,
									flag_name = 'cosmic_cmc',
									relationship_join = join_relationship
							) %>%
					filter(cosmic_cmc == 1)

		# dbSNP benign
		df_b5 	<- add_flag( df_bench_refset,
									unique( select(df_dbsnp_basic, all_of(c(join_columns, include_cols)) ) ), 
									join_cols = join_columns,
									flag_name = 'dbSNP_benign',
									relationship_join = join_relationship
							) %>%
					filter(dbSNP_benign == 1)

		# Now bind all rows together and clean NAs -> 0 
		df_bench_flags	<- 	bind_rows( 	df_b0,	# bench_refset
										df_b1,	# clinvarb
										df_b3,	# COSMIC CMC
										df_b5	# dbSNP benign
									) %>%
							arrange(gene, mutation_res_num) %>%
							replace_na( 
								list(bench_refset=0, clinvar_benign = 0, cosmic_cmc = 0, dbSNP_benign = 0) 
								)
	return(df_bench_flags)
}

# Add binary flag col on tbl_main (default 0). If matching geme+mutation exists in tbl_flag set flag 1.
add_flag <- function(tbl_main, tbl_flag, join_cols = c('gene', 'mutation'), flag_name = 'flag', relationship_join = 'one-to-one') {
    # Check required join cols
    if ( !all( join_cols %in% names(tbl_main) ) || !all( join_cols %in% names(tbl_flag) ) ) {
        stop( "Both tables must contain all join_cols: ", paste(join_cols, collapse = ", ") )
    }
    
	# Add binary flag column to tbl_main based on matching rows in tbl_flag
    df_result <- tbl_main %>%
        left_join( 
			tbl_flag %>% mutate(match = 1), 
			by = join_cols, 
			relationship = {{ relationship_join }},
			keep = FALSE, 
			na_matches = 'never' 
		) %>%
        mutate( {{flag_name}} := if_else( is.na(match), 0, 1) ) %>%
        select(-match)
  
    return(df_result)
}

# LEFT JOIN convenience fn which automatically converts NAs
fn_left_join <- function(dfx, dfy, join_columns, join_relationship) {
	require(dplyr)
	require(tidyr)
	return(
			left_join( {{ dfx }}, {{ dfy }},
						by = join_columns, 
						relationship = join_relationship,
						copy = FALSE, 
						keep = FALSE, 
						na_matches = 'na' 
					) 		
					%>% mutate( across( !where(is.double), ~ replace_na(.x, '-') ) )
					%>% mutate( across( where(is.double),  ~ replace_na(.x, 0) ) ) 
	)		
}

# HUGO gene check helper function
add_hugo_valid	<- function(df, df_hugo_exp_up_map = df_hugo_core_exp, .check_uprot = TRUE){
	require(dplyr)
	return( 
			gene_to_uniprot( {{ df }}, {{ df_hugo_exp_up_map }} ) 
				%>% gene_add_hugo_valid( .check_uprot = {{.check_uprot}} ) 
				%>% select( -c(ensembl_gene_id, hugo_gene_approved) ) 
				%>% rename( gene_orig = gene, gene = gene_hugo_valid )
	) 
}


# add column: hugo_gene_valid which is the correct HUGO gene name 
# expects a df_map of form returned from gene_to_uniprot()
gene_add_hugo_valid <- function(df_map, drop_alias_names = TRUE, drop_alt_reading_frames = TRUE, .check_uprot = TRUE) {
	require(dplyr)

	# In df_map these are valid HUGO gene names (hugo_gene_map_type from HGNC mapping dataset)
	df_map_hugo 	<- filter( df_map, hugo_gene_map_type == 'hugo_gene' ) %>% 
						mutate( gene_hugo_valid = hugo_gene_approved ) 
							
	# Rows with these gene names are not valid HUGO
	df_map_nothugo 	<- filter( df_map, hugo_gene_map_type != 'hugo_gene' )
	df_map_nothugo 	<- df_map_nothugo %>% 
						filter( !gene %in% pull( df_map_hugo, gene ) ) %>%
						mutate( gene_hugo_valid = hugo_gene_approved ) 

	df_out <- bind_rows(df_map_hugo, df_map_nothugo) %>%
					relocate( hugo_gene_map_type, .after = gene ) %>%
					relocate( gene_hugo_valid, .after = hugo_gene_map_type ) %>%
					arrange( gene, gene_hugo_valid ) 

	if (drop_alias_names){
		df_out <- filter( df_out, hugo_gene_map_type != 'alias_symbols' )
	}

	# A workaround...
	if (drop_alt_reading_frames & .check_uprot){
		df_out <- filter(df_out, !(gene == 'CDKN2A' & uniprot_acc == 'Q8N726'))
	}

	return(df_out)
}

# HUGO gene to uniprot
# NOTE: uses HUGO approved and all alias and previous symbols to match
# empty uniprots will be returned by default if just using a gene alias/prev -> approved map!
gene_to_uniprot <- function(df, df_hugo_exp_map, return_all = TRUE) {
	require(dplyr)
	df_map 	<- 	left_join(	df,
							df_hugo_exp_map,
							# filter(df_hugo_exp_map, hugo_gene_map_type == 'hugo_gene'),
							by = c('gene' = 'hugo_gene_map_name'),
							keep = TRUE,
							na_matches = 'na',
							relationship = 'many-to-many' 
					) 
					# %>% rename(uniprot_acc = uniprot_acc.x)
	
	if (!return_all){
		df_map 	<- 	filter( df_map, !is.na(uniprot_acc) )
	}
	return(df_map)
}

# Ensure only single row with flags per grouping cols (gene, mutation_res_no, mutation,...)
# This will take a df with a set of flags and ensure only one row per group of the grouping cols
enforce_unique_bench_dataset_with_group_max <- function(df, 
														group_cols = c('gene', 'mutation_res_num', 'mutation', 'uniprot_acc'), 
														flag_cols = c('bench_refset', 'cosmic_cmc', 'gdc_ctat', 'clinvar_benign', 'clinvar_pathogenic',  'dbSNP_benign')
														)
{
	require(dplyr)
	require(tidyr)
	
	df_out <- df %>%
		group_by(across(all_of(group_cols))) %>%
		mutate(across(all_of(flag_cols), ~ if_else(sum(.x) > 0, 1, 0))) %>%
		ungroup()
	
	return(unique(df_out))

}

# Function to filter rows based on a column and a list of values; returns columns in col_list
filter_rows <- function(df, filter_col, filter_vals, col_list){
	df_filt <- df %>% filter( {{filter_col}} %in% filter_vals ) %>% 
		select(all_of(col_list)) %>% unique()
	return(df_filt)
}

# Function to count number of unique residues in a column
# filter_col is the column to filter on with filter_vals list of values
count_unique_mut_residue <- function(df, filter_col, filter_vals){
	df_filt <- df %>% filter( {{filter_col}} %in% filter_vals ) %>% select(gene, uniprot_acc, mutation_res_num) %>% unique()
	return(nrow(df_filt))	
}
