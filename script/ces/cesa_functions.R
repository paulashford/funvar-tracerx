# ces cancer effect size analysis functions
# P. Ashford
# Aug/Sept 2025

# Info for cancereffectsizeR:
# https://townsend-lab-yale.github.io/cancereffectsizeR/articles/cancereffectsizeR.html

# Get CES variants for gene list (aac: amino acid coding only)
get_ces_variants <- function(cesa, gene_list, dt_var_pos = NULL, aac = TRUE){
	require(cancereffectsizeR)
	vars <- select_variants(
		cesa,
		genes = gene_list,
		min_freq = 0, 
		variant_ids = NULL, 
		gr = NULL, 
		variant_position_table = dt_var_pos, 
		include_subvariants = FALSE, 
		padding = 0, 
		collapse_lists = FALSE, 
		remove_secondary_aac = TRUE
	)
	if (!is.null(vars) & aac == TRUE) {vars <- vars[variant_type == 'aac']}
	return(vars)
}

# Ces get_variants() wrapper for FunFam bases analysis
ces_var_wrapper <- function(funfam_resnum, cesa, dt_vars){
	require(cancereffectsizeR)

	cesa_varsf <- NULL
	# get cesa variants per gene(s)
	cesa_vars <- get_ces_variants(cesa, dt_vars[FUNFAM_PDB_RES == funfam_resnum, genes][[1]])
	if (!is.null(cesa_vars)){
			cesa_varsf <- inner_join(dt_vars[FUNFAM_PDB_RES == funfam_resnum], 
									cesa_vars, 
									by = join_by(
										'chr' == 'chr', 
										'mut_pos' >= 'start', 
										'mut_pos' <= 'end',
										'mut_aa_change' == 'aachange'))
									

			# filter variants based on funfam relevant
			# cesa_varsf <- cesa_vars[center_nt_pos %in% dt_vars[FUNFAM_PDB_RES == funfam_resnum, mut_pos]]
			# cesa_varsf <- cesa_vars[ (start <= dt_vars[FUNFAM_PDB_RES == funfam_resnum, mut_pos]) & (end_pos >= dt_vars[FUNFAM_PDB_RES == funfam_resnum, mut_pos]) ]
			# cesa_varsf[, funfam := funfam_resnum]			
			# Ensure integer maf Processed 4373 groups out of 22226. 20% done. Time elapsed: 4583s. ETA: 18712s.Error in `[.data.table`(dt_ff_analyse, , ces_var_wrapper(FUNFAM_PDB_RES,  : 	#   Column 25 of result for group 4374 is type 'double' but expecting type 'integer'. Column types must be consistent for each group.
			cesa_varsf$maf_prevalence <- as.integer(cesa_varsf$maf_prevalence)
	}
return(cesa_varsf)
}


# Run Cancer effect calcs for funfams
# cesa: pre calc cancer effects analuysis set
# tib_ffsu: Simplified and unique on FUNFAM_PDB_RES, center_nt_pos and genes (this is just used to query the MAF via get_ces_variants() while preserving FIE/FunFam info)
# 		: eg tib_ffsu <- select(tib_ffs, FUNFAM_PDB_RES, chr, mut_pos, mut_aa_change, genes) %>% unique() %>% arrange(FUNFAM_PDB_RES, mut_pos)
# df_to_analyse: FIEs? non-drivers? (funfam/cesa table)
run_ces_analysis <- function(cesa, tib_ffsu, df_to_analyse){
	require(cancereffectsizeR)
	require(data.table)
	require(dplyr)
	dt_ff_analyse <- as.data.table(
					tib_ffsu %>% 
						filter(FUNFAM_PDB_RES %in% df_to_analyse$FUNFAM_PDB_RES)
					)
	return(
		dt_ff_analyse[, ces_var_wrapper(FUNFAM_PDB_RES, cesa, dt_ff_analyse), by = .I]
	)
}

# 
simplify_funfam_table <- function(dt_funfam){
	require(dplyr)
	dt_ff <- as_tibble(dt_funfam)
	dt_ff <- dt_ff %>%
			select(CHR, START_POSITION, REF_ALLELE, VAR_ALLELE, 
			SOURCE_HUGO_SYMBOL, VM_UNIPROT_ACCESSION, VM_SEQ_NO, VM_AA_CHANGE,
			CANCER_TYPE, DATA_SOURCE, TXP,  TUMOUR_ID, MUTATION_ID, FUNFAM_PDB_RES) %>%
			dplyr::rename('chr' = CHR, 
							mut_pos = START_POSITION,
							ref = REF_ALLELE,
							alt = VAR_ALLELE,
							source_gene = SOURCE_HUGO_SYMBOL) %>%
			mutate( split_aa_change = str_split(VM_AA_CHANGE, '/'),
					aa_from = map_chr(split_aa_change, ~ .x[1]),
					aa_to = map_chr(split_aa_change, ~ .x[2]),
					mut_aa_change = paste0(aa_from, VM_SEQ_NO, aa_to)
			)
	 
	return(dt_ff)
}

# create a simplified variants and gene-list table for using with cancereffectsize function get_ces_variants()
# Simplified and unique on FUNFAM_PDB_RES, center_nt_pos and genes (this is just used to query the MAF via get_ces_variants() while preserving FIE/FunFam info)
create_ff_groups <- function(tib_cesa){
	require(dplyr)
	# FunFam-based 'compound variants' analysis
	tib_ffg <- tib_cesa %>% 
		group_by(FUNFAM_PDB_RES) %>% 
		summarise(gene_list = list(gene)) %>%
		# ensure unique gene list
		mutate(genes = map(gene_list, ~ unique(.))) %>%
		select(FUNFAM_PDB_RES, genes)

	tib_ffsu <- tib_cesa %>%
		select(FUNFAM_PDB_RES, chr, mut_pos, ref.x, alt.x, mut_aa_change, TUMOUR_ID, MUTATION_ID) %>%
		inner_join(tib_ffg, by = 'FUNFAM_PDB_RES') %>%
		select(FUNFAM_PDB_RES, chr, mut_pos, ref.x, alt.x, mut_aa_change, genes) %>% 
			unique() %>% 
			arrange(FUNFAM_PDB_RES, mut_pos)
	
	return(tib_ffsu)

}



# Standard formatting of Supplemental Figure plots for ces figs
format_plots <- function(ggp, fs_title = 21, fs_axis = 12 ){
	return(
		ggp +	theme(panel.grid.major = element_line(colour = "#d9d6d6")) +
				theme(axis.line = element_line(linewidth = 1, colour = "#9b9898" )) +
				theme(axis.ticks = element_line(linewidth = 1, colour = "#000000" )) +
				theme(axis.title = element_text(size = fs_title)) +
				theme(axis.text.x = element_text(angle = 45, hjust = 1, size = fs_axis)) + 
				theme(axis.text.y = element_text(size = fs_axis)) +
				theme(legend.position = 'inside',
        				legend.byrow = TRUE,
 						legend.position.inside = c(.97, 0.02),
 						legend.justification = c("right","bottom"),
        				legend.background = element_rect(colour = "black"),
						legend.title = element_text(face = "bold", size = fs_axis),
						legend.title.position = "bottom",
				        legend.key.spacing.x = unit(17.5, 'cm'),
						legend.margin = margin(2, 10, 2, 10)) # upper, right, lower and left padding
        			# legend.title = element_blank())
	)
}

