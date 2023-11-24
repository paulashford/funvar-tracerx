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

##################################################################################
# st02 Basic percentages of FIEs/muts by tumours by timing etc
# updaed 10 08 2023
##################################################################################
# timed df for muts or fiew eg: dt_tot_timed_fies  (from _yj05 in stats_load.R) 
stat_st02 <- function( dt_timed,  txps=c('N','Y'), cancer_types=c('LUAD', 'LUSC'), score='MC3', timings=c('early', 'late', 'unknown' ) ){
    df_filt <- filter(    dt_timed,
                            TXP %in% txps &
                            CANCER_TYPE %in% cancer_types &
                            SCORE_THRESHOLD == score &
                            OVERALL_TIMING_MOD %in% timings
                        )

    # #remove "NA" literal in Tx
    # if (txp == 'N'){
    #     df_filt <- filter(df_filt, OVERALL_TIMING_MOD !='NA')
    # }
    
    # if (remove_unknown=='Y'){
    #     df_filt <- filter(df_filt, OVERALL_TIMING_MOD !='unknown')
    # }

    # Total timed FIEs
    # tot_timed <-  sum(df_filt$NUM_FIES)
    tot_timed <-  sum(df_filt$NUM_MUTS)

    # Calc % FIE per time
    # df_filt <- df_filt %>%   
    #                 rowwise() %>%
    #                 # select( TXP, CANCER_TYPE, OVERALL_TIMING_MOD, NUM_MUTS ) %>% 
    #                 mutate( fie_perc_timed = 100 * ( NUM_MUTS/tot_timed ) )
    
    return(df_filt)

    # print(paste(c('woof', tot_timed)))   

}

##################################################################################
# st03 Percentage of FIEs/genes that are CGC
# use script/v3.1/scrivener_mgmt/scriv014*
##################################################################################
stat_st03 <- function(dt_cgc_count,  txp='N', cancer_type='LUAD'){
    # If combo then no need to filter cancer type
    # if (cancer_type == 'COMBO'){
    #     cgc_counts <- filter(   dt_cgc_count,
    #                             TXP == txp
    #                       )
    # } else {
    #    cgc_counts <- filter(    dt_cgc_count,
    #                             TXP == txp,
    #                             CANCER_TYPE == cancer_type
    #                       )
    # }
    # tot_cgc_genes <-  sum(cgc_counts$NUM_GENES)

    # scriv014_COMBO.tsv
    # DATA_SOURCE	TXP	CGC_GENE	NUM_GENES	NUM_FIES
    # Tx	N	N	85	89
    # Tx	N	Y	57	261
    # Tx	Y	N	217	226
    # Tx	Y	Y	69	282
    
    # Tx FIEs 
    89 + 261 = 350
    # FIE CGC 
    261 / 350 = 0.7457143
    # FIE Novel 
    89 / 350 =  0.2542857
    
    # Tx genes
    85 + 57 = 142
    # CGC
    57 / 142 = 0.4014085
    # Novel
    85 / 142 = 0.5985915
    
    # FIEs per gene CGC
    261 / 57
    4.578947
    # FIEs per gene novel
    89 / 85
    1.047059

    # scriv014.tsv
    # DATA_SOURCE	TXP	CANCER_TYPE	CGC_GENE	NUM_GENES	NUM_FIES
    # Tx	N	LUAD	N	58	62
    # Tx	N	LUAD	Y	36	197
    # Tx	N	LUSC	N	27	27
    # Tx	N	LUSC	Y	21	64
    # Tx	Y	LUAD	N	116	121
    # Tx	Y	LUAD	Y	39	136
    # Tx	Y	LUSC	N	101	105
    # Tx	Y	LUSC	Y	30	146

}

##################################################################################
# st04 Percentage FIEs that are driverMuts
# use yj11 unfiltered
# ONLY WORKS FOR TX
##################################################################################
stat_st04 <- function(dt_yj11){
    dt_yj11 %>%
        group_by(TXP, CANCER_TYPE, TX_DRIVERMUT) %>%
        summarise(n = n())

    # Groups:   TXP, CANCER_TYPE [4]
    #   TXP   CANCER_TYPE TX_DRIVERMUT      n
    #   <chr> <chr>       <chr>         <int>
    # 1 N     LUAD        FALSE            70
    # 2 N     LUAD        TRUE            189
    # 3 N     LUSC        FALSE            29
    # 4 N     LUSC        TRUE             62
    # 5 Y     LUAD        NotApplicable   257
    # 6 Y     LUSC        NotApplicable   251    
    
# [ins] r$> 251/350
# [1] 0.7171429

}

##################################################################################
# st05 Percentage EARLT FIEs that are driverMuts or CGC
# use yj11 unfiltered
##################################################################################
stat_st05 <- function(dt_yj11){
    dt_yj11 %>%
        # group_by(TXP, CANCER_TYPE, OVERALL_TIMING_MOD,SOURCE_HUGO_SYMBOL, CGC_GENE) %>%
        # group_by(TXP, OVERALL_TIMING_MOD, CGC_GENE) %>%
        group_by(TXP, OVERALL_TIMING_MOD, TX_DRIVERMUT) %>%
        summarise(n = n())

    # Groups:   TXP, CANCER_TYPE [4]
    #   TXP   CANCER_TYPE TX_DRIVERMUT      n
    #   <chr> <chr>       <chr>         <int>
    # 1 N     LUAD        FALSE            70
    # 2 N     LUAD        TRUE            189
    # 3 N     LUSC        FALSE            29
    # 4 N     LUSC        TRUE             62
    # 5 Y     LUAD        NotApplicable   257
    # 6 Y     LUSC        NotApplicable   251    
    
# [ins] r$> 251/350
# [1] 0.7171429

}

##################################################################################
# st5 to st12 got here...
##################################################################################

##################################################################################
# st13 FIE scores and CGC distribution
# 17 12 2022
# Supp Fig _yg
# data: _yj13b
# plot:Dropbox/bioinf/neofun/paper/script/v3.1/plot/nfe_plot8_CGC_and_Tx_driver_distributions.ipynb
##################################################################################
# Chi-squared for the discrete distribution of NFE scores depending
# https://stats.stackexchange.com/questions/95559/testing-if-two-non-normal-distributions-are-significantly-different-k-s-or-wilc
# TCGA pancancer (not TxP)
# /Users/ash/Dropbox/bioinf/neofun/paper/calc/_yj13/_yj13_tcga_only_count_fies_by_score_and_cgc.tsv
# DATA_SOURCE	TXP	SCORE_THRESHOLD	CGC	NFE_SCORE_D_MF_H	NUM_FIES	NUM_ROWS
# TCGA	N	ANY	N	1	1698	1698
# TCGA	N	ANY	N	2	2757	2757
# TCGA	N	ANY	N	3	2308	2308
# TCGA	N	ANY	N	4	1480	1480
# TCGA	N	ANY	N	5	608	608
# TCGA	N	ANY	N	6	154	154
# TCGA	N	ANY	N	7	39	39

# TCGA	N	ANY	Y	1	172	172
# TCGA	N	ANY	Y	2	470	470
# TCGA	N	ANY	Y	3	1113	1113
# TCGA	N	ANY	Y	4	1043	1043
# TCGA	N	ANY	Y	5	1618	1618
# TCGA	N	ANY	Y	6	948	948
# TCGA	N	ANY	Y	7	329	329
# TCGA	N	ANY	Y	8	35	35
# horrid - this should be table read!!
stat_st13 <- function(){
    tcga=matrix(
        c(  1698,
            2757,
        	2308,
        	1480,
        	608,
        	154,
        	39,
            0,
    	172,
    	470,
    	1113,
    	1043,
    	1618,
    	948,
    	329,
    	35
        ), 8, 2,
    dimnames=list( fie_score = c(1, 2, 3, 4, 5, 6, 7, 8),
                    cgc_gene = c("N", "Y")
                )
    )
# tcga
#          cgc_gene
# fie_score    N    Y
#         1 1698  172
#         2 2757  470
#         3 2308 1113
#         4 1480 1043
#         5  608 1618
#         6  154  948
#         7   39  329
#         8    0   35

chisq.test(tcga)
# v2    Pearson's Chi-squared test: X-squared = 6011.4, df = 8, p-value < 2.2e-16
# v3.1  Pearson's Chi-squared test: X-squared = 4116.1, df = 7, p-value < 2.2e-16  
}