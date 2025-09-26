# script/ces/cesa_funfam_tcga_pancancer.R
# CES Plot cancer effect size analysis for TCGA pan-cancer
# NOTE: run cesa_analysis_create_working_datasets.R first!
# Aug 2025

# Info for cancereffectsizeR package
# https://townsend-lab-yale.github.io/cancereffectsizeR/articles/cancereffectsizeR.html

# NOTE: in case of cryptic error messages, suggestion is to upgrade BioConductor e.g.
# BiocManager::install(version = "3.21") # or newer, when available

# Some dependencies are large, so we increase the download time limit to be safe
options(timeout = 6000)
library(cancereffectsizeR)
library(data.table)
library(fst)
library(tidyverse)

resources_dir 	<- file.path('data/resources')
benchmark_dir	<- file.path(resources_dir, 'benchmarking')
dat_fies		<- 'data/mutations_inc_fies_tcga_pancancer_funfam_pdb_res.tsv'

# max display length for FunFam names
max_name_len	<- 40

run_name  <- 'FunVar'
# Set histology filter and load a working set
histo_filter	<- 'PANCAN'
run_name <- str_c('TCGA', histo_filter, run_name, sep = "-")

# Load working set created by: cesa_analysis_create_working_datasets.R
working_set_rda <- paste0('data/rda_ces/working_set_', run_name, '.rda')
# ***
# working_set_rda <- paste0('private/data/rda_ces/ces_tcga_pancancer_funvar_v3p_dnds', '.rda')

if (file.exists(working_set_rda)){
	load(working_set_rda)	
} else{
	stop( paste0('Working set ', working_set_rda, ' not found - run cesa_analysis_create_working_datasets.R for TRACERx ', histo_filter, ' first.') )
}

# source analysis functions script
source('script/ces/cesa_functions.R')

# FunFam general info from CATH
dt_funfams 	<- read_delim( file.path(resources_dir, 'cath/cath_funfam_all_v4_2.tsv') ) %>%
				filter(REP_SOURCE_ID == 'cath')

# TCGA mutation and FIE data with FunFam_PDB_RES col
dt_tcga_mut_fie   <- as.data.table( read_delim(dat_fies, delim = '\t', guess_max = Inf) )
# simplify table format
dt_tcga_mut_fie[ , CHR := as.character(CHR)] 
dts_ff_muts_pancan 	<- simplify_funfam_table(dt_tcga_mut_fie)
# get FIEs and drivers and other muts
dt_ff_fies 		<- dt_tcga_mut_fie %>% filter(SCORE_THRESHOLD == 'MC3')

# load dbSNP 'high quality' (most certain) benign SNPs
dt_dbsnp 		<- read_delim(file.path(benchmark_dir, 'dbsnp_pre_processed.tsv'), delim = '\t', col_names = TRUE) %>%
					mutate(aa_change = paste0(mut_from, '/', mut_to))
dt_ff_dbsnp		<-  dt_tcga_mut_fie %>%
						inner_join(dt_dbsnp,
							by = join_by('SOURCE_HUGO_SYMBOL' == 'gene',
											'VM_SEQ_NO' == 'mut_res',
											'VM_AA_CHANGE' == 'aa_change'
								)
							)

# --------------------------------------------------
# calc effects for FIEs pancancer
# --------------------------------------------------
# We select ces variants using gene names. However, the hg19 reference gene names have differences from our later funvar data
# so we need to join on genomic positions and use the hg19 name rto select variants using the ces select_variants() function
tib_cesa_ff <- as_tibble(dts_ff_muts_pancan %>% 
							inner_join(cesa$variants,
								by = join_by('chr' == 'chr', 
									'mut_pos' >= 'start', 
									'mut_pos' <= 'end',
									'mut_aa_change' == 'aachange'
								)
							)
)		 
tib_ffsu <- create_ff_groups(tib_cesa_ff)

# Filter analysis to more sensible numbers for figures...
# FunFam rep domains for those uniquely identified by FunVar in Supplemental Table 11
st11_doms <- c('3n5nX02', '2vn9A02', '1b55A00', '4p7aA02', '1fbvA02', '1m14A01', '3e04A02', '3v00A02', '1nf1A02', '1zkcA00', '4yonA02', '2y3aA02', '1a07A00', '2jwaA00', '1m14A02', '3my0A01', '1dbhA01', '1vr2A02', '1w9cA00')
dt_ff_fies_st11 <- dt_ff_fies[REP_ID %in% st11_doms]
# Calc cancer effect size for analysis FIEs
df_out_pancan_st11 <- run_ces_analysis(cesa, tib_ffsu = tib_ffsu, df_to_analyse = dt_ff_fies_st11)

# Diverse funfams/FIEs ST5
st5_doms <- c('1cf7A00','4nftA00','2cv5A00','2zv2A02','3omvA02','1b6cB02','4pypA00','1a5eA00','1n8zC02','1q46A01','1jhjA00','2di7A01','1f05A00','3rgfA01','1uwhA01','3eqcA01','5b8aA00','1u0eA02','3bfkA00','3mwyW04','5jcpA02','1kaoA00','2wtbA02','5ax8A02','2hr7A03','1l8kA00','1d5rA01','2ffuA01')
dt_ff_fies_st5 <- dt_ff_fies[REP_ID %in% st5_doms]
# Calc cancer effect size for analysis FIEs
df_out_pancan_st5 <- run_ces_analysis(cesa, tib_ffsu = tib_ffsu, df_to_analyse = dt_ff_fies_st5)

# dbSNP benign
# Calc cancer effect size for dbSNPs is pancancer muts
df_out_dbsnp <- run_ces_analysis(cesa, tib_ffsu = tib_ffsu, df_to_analyse = dt_ff_dbsnp)

#########################################################
# select run for compound variant based ces analysis
# ** Note: choose ST11, ST5, or dnSNP below **
#########################################################
# Plot FIE sites
# dt_out <- df_out_pancan_st11_sites_fies
# run_name_compound = paste0(run_name, '-compound-FunFam-PDB-res-pancan-ST11-sites')

# dt_out <- df_out_pancan_st5_sites_fies
# run_name_compound = paste0(run_name, '-compound-FunFam-PDB-res-pancan-ST5-sites')

# Plot dbSNP sites
dt_out <- df_out_dbsnp
run_name_compound = paste0(run_name, '-compound-FunFam-PDB-res-dbSNPb')

#########################################################
# Run compound variant based ces analysis
#########################################################
# unique variant_ids
unique_variants <- unique(dt_out, by = c('variant_id'))

compound_vars <- define_compound_variants(cesa, variant_table = unique_variants, by = 'FUNFAM_PDB_RES', merge_distance = Inf)
cesa <- ces_variant(cesa, variants = compound_vars, run_name = run_name_compound) 

# Join COMPOUND VARS to FunVar table
dt_cesa_sel_ff <- arrange(cesa$selection[[run_name_compound]], -selection_intensity) %>%
  mutate(
    	split_variant_name = str_split(variant_name, "_"),
		sf_id 	= map_chr(split_variant_name, ~ paste(.x[1])),
		ff_id 	= as.numeric( map_chr(split_variant_name, ~ paste(.x[2])) ),
		rep_id	= map_chr(split_variant_name, ~ paste(.x[3])),
    	funfam 	= map_chr(split_variant_name, ~ paste( .x[ 1:min( 3, length(.x) ) ], collapse = "_" ) ),
		funfam_variant = str_sub(variant_name, 1, nchar(variant_name) -2 )
		
		) %>% select(-split_variant_name)
# View(dt_cesa_sel_ff)

# CES analysis selection by FunFam joined with mutation/FIE table and FunFam info (for the full FunFam names)
dt_cesa_sel_ff_muts <- dt_cesa_sel_ff %>%
							inner_join( dt_tcga_mut_fie,
								by = join_by('funfam_variant' == 'FUNFAM_PDB_RES'
										),
								relationship = "one-to-many"
							) %>%
							inner_join(dt_funfams,
								by = join_by('sf_id' == 'SUPERFAMILY_ID',
											'ff_id' == 'FUNFAM_NUMBER',
											'rep_id' == 'REP_ID'),
								relationship = "many-to-one"
							) %>%
							# Create better names for funfam and variant labels
							# 3.40.190.10_203444_1b0lA01 + cDNA FLJ53691, highly similar to Serotransferrin
							mutate(
								funfam_rep_domain_id = str_c(str_trunc(NAME, width =  max_name_len),
									' (', rep_id, ')')
							)
# TEST check unique mutations
if ( length(unique(dt_cesa_sel_ff_muts$MUTATION_ID)) != length(dt_cesa_sel_ff_muts$MUTATION_ID) ){
	stop('Mutation IDs are not unique on dt_cesa_sel_ff_muts')
}else {
   print('PASS: All mutation IDs on dt_cesa_sel_ff_muts are unique.')
}
# View(dt_cesa_sel_ff_muts)

# Simplify 
dt_cesa_sel_ff_muts_sml <- dt_cesa_sel_ff_muts %>% 
	select( colnames(cesa$selection[[run_name_compound]] ),
			sf_id, ff_id, rep_id, funfam, funfam_rep_domain_id, funfam_variant,
			CANCER_TYPE, FIE_threshold = SCORE_THRESHOLD, TRACERx_driver = DRIVERMUT, 
			source_gene = SOURCE_HUGO_SYMBOL, UniProt_ACC = VM_UNIPROT_ACCESSION, seq_num = VM_SEQ_NO, aa_change = VM_AA_CHANGE 
			) %>% 
			dplyr::rename(variant_name_orig = variant_name) %>%
			dplyr::rename(variant_name = funfam_variant) %>%
			relocate(all_of( c('funfam_rep_domain_id', 'source_gene','UniProt_ACC', 'seq_num', 'aa_change') ), .before = 'variant_type' ) %>%
			relocate('variant_name_orig', .after = 'TRACERx_driver') %>%
			select(-shared_cov) %>%
			unique() %>%
			arrange(funfam_rep_domain_id, -selection_intensity) 

# TEST check unique mutations
if ( length(unique(dt_cesa_sel_ff_muts_sml$MUTATION_ID)) != length(dt_cesa_sel_ff_muts_sml$MUTATION_ID) ){
	stop('Mutation IDs are not unique on dt_cesa_sel_ff_muts_sml')
}else {
   print('PASS: All mutation IDs on dt_cesa_sel_ff_muts_sml are unique.')
}
# View(dt_cesa_sel_ff_muts_sml)

# Plot cancer effects
plot( format_plots(
		plot_effects(dt_cesa_sel_ff_muts_sml, 
						group_by = "funfam_rep_domain_id", 
						topn = 200, 
						label_individual_variants = FALSE
						
				) )
)

