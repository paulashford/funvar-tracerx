# script/v3.1/script/neofun_v3.1_paper_stat_functions.R
# 29 09 2022
# Nov '23 submission - tidy/remove any hard-coded paths etc

##################################################################################
# out_summary(df) 
# Summarise counts as per groups by fies, tumours etc
# df should be _yfa (ie FIE table) format
##################################################################################
    # GENERAL OUTPUT as per groups		
	out_summary <- function( dtf ){
		return( dtf %>%
			summarise(
                    num_fies = n_distinct( m_id ), 
                    # num_fies = n_distinct( MUTATION_ID ),
					num_genes = n_distinct( SOURCE_HUGO_SYMBOL ),
					num_tumours = n_distinct( t_id ),
                    # num_tumours = n_distinct( TUMOUR_ID ),
					num_ff = n_distinct( SF_ID, REP_ID ),
					fies_per_tumour = num_fies/num_tumours,
					fies_per_gene = num_fies/num_genes,
					num_rows = n(),
					.groups = "keep"
			)
		)
	}
    # Add group percentage col to summary df
    # eg would give total FIEs per TXP and cancer_type so can use for %ages and apply also
    # for % in timing classes etc
    out_summary_perc <- function( dtf ){
        return(
            dtf %>%
                out_summary() %>%
                    mutate( fie_perc_group = round(100 * num_fies/tot_fies_group, 1) ) 
        )
    }

##################################################################################
# st_fisher run Fisher's exact test 
#  test  <- c(208, 88, 5351, 7450)
#           c(fie_early, fie_late, mut_early, mut_late)
##################################################################################
fie_fischers <- function(test_list, test_type ="fisher", conf = 0.95) {
    # Construct contingency table
    mat_test <- matrix(test_list,
                    nrow = 2,
                    dimnames = list(row_cat = c("Early", "Late"),
                                    col_cat = c("FIEs", "Other mutations") 
                                    )
                    )
    print(mat_test)
    # choose & do test
    if (test_type == "fisher") {
        test <- fisher.test(mat_test,
                            alternative = "two.sided",
                            conf.level = conf)
    } else if (test_type == "chisq") {
        test <- chisq.test(mat_test)
    } else {
        print("Invalid test type choose fisher or chisq ")
        return(-1)
    }
    return(test)
}

##################################################################################
# st01 Fischer's test early v late FIEs and mutations
# Note yj04 dt_tot_fies doesn't include muts, so best to group yj07 counts instead
##################################################################################
stat_st01 <- function(dt_gene_fie, dt_gene_SYN, txp='N', cancer_type='LUAD'){
    fies <- dt_gene_fie %>%
        group_by( TXP, CANCER_TYPE, OVERALL_TIMING ) %>%
        filter( CANCER_TYPE == cancer_type & TXP == txp )
        # summarise(tot_fies = sum(NUM_FIES))
    syns <- dt_gene_SYN %>%
        group_by( TXP, CANCER_TYPE, OVERALL_TIMING ) %>%
        filter( CANCER_TYPE == cancer_type & TXP == txp )
        # summarise(tot_syns = sum(NUM_SYNONYMOUS_MUTS))
        
    # Fisher's input:
    # filter for correct data source & cancer type
    # fies <- filter(fie_group, CANCER_TYPE == cancer_type & TXP == txp)
    # syns <- filter(syn_group, CANCER_TYPE == cancer_type & TXP == txp)
    # Make input matrix
    test_mat <- c(  fies[fies$OVERALL_TIMING == 'early', ]$tot_fies,
                    fies[fies$OVERALL_TIMING == 'late', ]$tot_fies,
                    syns[syns$OVERALL_TIMING == 'early', ]$tot_syns,
                    syns[syns$OVERALL_TIMING == 'late', ]$tot_syns
                )

    return( fie_fischers(test_mat) )
}
