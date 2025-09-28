# Dropbox/bioinf/neofun/paper/script/v3.2/submission/neofun_paper_stats_calc.R
# Compute various stats used in paper
# Nov '23 submission - tidy/remove any hard-coded paths etc

# Init
options("width"=300)   # screen width - was 80  chars
options(tibble.width = Inf)   # tibble width in terms of num cols
options(tibble.print_max = 14, tibble.print_min = 7)

# load STs
source('script/neofun_paper_stats_load.R')

# functions
source('script/neofun_paper_stat_functions.R')

# §1 FIE counts by Txp/cancer_type/CGC/driver
# -------------------------------------------
dt_fiet <- dt_hfiet

	# (a) FIEs by Tx/Txp
	dtf <- dt_fiet %>%
		group_by( TXP ) 
	out_summary(dtf)	
	
	# (b) FIEs by Tx/Txp v LUAD/LUSC
	dtf <- dt_fiet %>%
		group_by( TXP, CANCER_TYPE ) 
	out_summary(dtf)	

	# (c) FIEs by Tx/Txp - CGC
	dtf <- dt_fiet %>%
		group_by( TXP, CGC ) 
	out_summary(dtf)	

	# (d) FIEs by Tx/Txp - CGC and Tx drivermut
	dtf <- dt_fiet %>%
		group_by( TXP, CGC, TX_DRIVERMUT ) 
	out_summary(dtf)	

	# (e) FIEs by Tx/Txp - Tx drivermut
	dtf <- dt_fiet %>%
		group_by( TXP, TX_DRIVERMUT ) 
	out_summary( dtf )

	# (f) FIEs by Tx/Txp - Tx drivermut and LUAD/LUSC
	dtf <- dt_fiet %>%
		group_by( TXP, TX_DRIVERMUT, CANCER_TYPE )
	out_summary( dtf )

	# (g) FIE-genes ident in >5 tumours [if inc cancer_type note filter will miss SMARCA4 ie is not >5 in both LUAD and LUSC]
	dt_fiet %>%
		group_by( TXP, CGC,  SOURCE_HUGO_SYMBOL ) %>%
		# group_by( TXP, CGC, CANCER_TYPE, SOURCE_HUGO_SYMBOL ) %>%
		summarise( num_fies = n_distinct( m_id ), .groups = "keep" ) %>%
		filter(num_fies>5) 

# §1S: FIE SUMMARY TABLES [full outputs in script/neofun_paper_stats_outputs.txt]
# -------------------------------------------------------------------------------
	# ** (a) TXP **
		#   TXP   num_fies num_genes num_tumours num_ff num_rows
		# 1 N          355       109         224     95      355
		# 2 Y          446       216         337    198      446
	
	# ** (b) TXP/CANCER_TYPE **
		#   TXP   CANCER_TYPE num_fies num_genes num_tumours num_ff num_rows
		# 1 N     LUAD             264        80         164     69      264
		# 2 N     LUSC              91        37          60     36       91
		# 3 Y     LUAD             212       121         161    119      212
		# 4 Y     LUSC             234       112         176    103      234
	
	# ** (c) TXP/CGC ** 
		#   TXP   CGC   num_fies num_genes num_tumours num_ff num_rows
		# 1 N     N           91        84          69     75       91
		# 2 N     Y          264        25         197     27      264
		# 3 Y     N          188       171         152    158      188
		# 4 Y     Y          258        45         238     49      258

	# ** (d) TXP/CGC/TX DRIVERMUT **
		#   TXP   CGC   TX_DRIVERMUT  num_fies num_genes num_tumours num_ff num_rows
		# 1 N     N     FALSE               87        82          65     74       87
		# 2 N     N     TRUE                 4         2           4      2        4
		# 3 N     Y     FALSE               14         9          14      8       14
		# 4 N     Y     TRUE               250        18         191     21      250
		# 5 Y     N     NotApplicable      188       171         152    158      188
		# 6 Y     Y     NotApplicable      258        45         238     49      258

	# ** (e) TXP/TX DRIVERMUT (TXP 'N' only **
		#   TXP   TX_DRIVERMUT  num_fies num_genes num_tumours num_ff num_rows
		# 1 N     FALSE              101        91          72     80      101
		# 2 N     TRUE               254        20         192     23      254

	# ** (f) TXP/TX DRIVERMUT + LUAD/LUSC (TXP 'N' only) **
		#   TXP   TX_DRIVERMUT  CANCER_TYPE num_fies num_genes num_tumours num_ff num_rows
		# 1 N     FALSE         LUAD              72        64          50     57       72
		# 2 N     FALSE         LUSC              29        28          22     27       29
		# 3 N     TRUE          LUAD             192        17         144     20      192
		# 4 N     TRUE          LUSC              62        10          48     12       62

# §1C: FIE SPECIFIC CALCS
# -------------------------
	# ** (ai) Tx NSCLC driver muts **
		#  [Abstract]
		# From TRACERx lung FIEs, 71.5% (n=254 FIEs in 20 genes) were known NSCLC drivers and 28.5% (n=101 FIEs in 91 genes) were other mutations that could be novel drivers acting via their impacts on protein interface regions or binding sites, including those for ligands and nucleic acids.
		# NSCLC drivers: TX_DRIVERMUT =='Y' [§1e] 254 FIEs in 20 genes  [254/]
		round( 100 * (254 / 355) , 1 ) 
		# 71.5%
		# 'other muts': TX_DRIVERMUT =='N' [§1e]  101 FIEs in 91 genes
		round( 100 * (101 / 355) , 1 ) 
		# 28.5%
		
	# ** (aii) FIEs also known Tx driver muts LUAD LUSC **
		# "Using driver mutations classified with the TRACERx pipeline26, 254 FIEs (71.5%) were classed as TRACERx drivers and with similar frequency in LUAD and LUSC (~1.3 FIEs per tumour)"
		# Use §1f [added fies_per_tumour calc at this point and altered display/col width options for script too]
		#   TXP   TX_DRIVERMUT  CANCER_TYPE num_fies num_genes num_tumours num_ff fies_per_tumour num_rows
		# 1 N     FALSE         LUAD              72        64          50     57            1.44       72
		# 2 N     FALSE         LUSC              29        28          22     27            1.32       29
		# 3 N     TRUE          LUAD             192        17         144     20            1.33      192
		# 4 N     TRUE          LUSC              62        10          48     12            1.29       62
	
	# ** (aiii) FIEs not known Tx driver muts or CGC **
		# In total 91 FIEs in 69 tumours occurred in 84 genes that were not known cancer genes (4 of these FIEs in 2 genes were TRACERx driver mutations), corresponding to just over 1 FIE per gene. 
		# Use §1c [added fies per gene calc at this point too]
		# 	  TXP   CGC   num_fies num_genes num_tumours num_ff fies_per_tumour fies_per_gene num_rows
		# 1 N     N           91        84          69     75            1.32          1.08       91
		# 2 N     Y          264        25         197     27            1.34         10.6       264
		# 3 Y     N          188       171         152    158            1.24          1.10      188
		# 4 Y     Y          258        45         238     49            1.08          5.73      258
		# Use §1d for CGC/Tx driver correspndence
		# 	TXP   CGC   TX_DRIVERMUT  num_fies num_genes num_tumours num_ff fies_per_tumour fies_per_gene num_rows
		# 1 N     N     FALSE               87        82          65     74            1.34          1.06       87
		# 2 N     N     TRUE                 4         2           4      2            1             2           4
		# 3 N     Y     FALSE               14         9          14      8            1             1.56       14
		# 4 N     Y     TRUE               250        18         191     21            1.31         13.9       250

	# ** (b) Enzyme functions of TRACERx FIE genes
		# EC from UniProt [ST_yu [ST3]] searched using all FIE gene UniProt IDs for Tx [ST_yfa [ST2]]
		# Tx FIE genes having an EC number [ST_yu [ST3]]: 62
		# Number of distinct FIE-genes in Tx: TXP=='N' [§1a]: 109
		round( 100 * (62 / 109) , 0 ) 
		# 57%

	# ** (c) Tx NSCLC CGC **
		# "Overall, 264 TRACERx FIEs (74.4% of FIEs identified) were found in 25 known cancer genes, an average of over 10 FIEs per known cancer gene"
		# TRACERx FIEs/genes: TXP=='N' & CGC=='Y' [§1c] / 264 FIEs in 25 CGC genes
		# total Tx FIEs: TXP=='N' [§1c] 264 + 91 or TXP=='N' [§1a] 355
		round( 100 * (264 / 355) , 1 ) 
		# 74.4%
		# "In total 91, FIEs in 69 tumours occurred in genes that were not known cancer genes, with on average just over 1 FIE per gene."

	# ** (d) Genes with >5 FIEs **
		# Only 7 genes had FIEs identified in more than 5 tumours and all of these were known cancer genes: TP53 (n=74 FIEs), PIK3CA (28), CDKN2A (13), KRAS (90), SMARCA4 (7) plus LUAD-specific EGFR (10) and BRAF (6). 
		# Use §1g (it creates output directly - see note there re LUAD specific EGFR BRAF)
		#    TXP   CGC   SOURCE_HUGO_SYMBOL num_fies
		#  1 N     Y     BRAF                      8
		#  2 N     Y     CDKN2A                   13
		#  3 N     Y     EGFR                     10
		#  4 N     Y     KRAS                     90
		#  5 N     Y     PIK3CA                   28
		#  6 N     Y     SMARCA4                   7
		#  7 N     Y     TP53                     74

# §2 COHORT counts _yj03
# ----------------------
	# note: already a count table - just filter if nec
	dt_yj03 %>%	
		filter( DATA_SOURCE == 'Tx' & CANCER_TYPE %in% c( 'LUAD', 'LUSC' ) )
	# 	DATA_SOURCE TXP   CANCER_TYPE NUM_PATIENTS NUM_TUMOURS
	# 1 Tx          N     LUAD                 235         239	[Tx421]
	# 2 Tx          N     LUSC                 134         134	[Tx421]
	# 3 Tx          Y     LUAD                 387         387	[TCGA-lung]
	# 4 Tx          Y     LUSC                 342         342	[TCGA-lung]


# §2C COHORT count calcs
# ----------------------
	# ora_exported_SI_yj03.tsv NOTE: v03_22
	# "Using these tuneable sites, we then applied FunVar to missense SNVs identified from the TRACERx NSCLC study using tumours classified as either 
	# lung adenocarcinoma (LUAD, n=239 tumours) or lung squamous cell carcinoma (LUSC, n=134 tumours). "
	# DATA_SOURCE	TXP	CANCER_TYPE	NUM_PATIENTS	NUM_TUMOURS
	# Tx	N	LUAD	235	239
	# Tx	N	LUSC	134	134
	# Tx	N	default	46	46
	# Tx	Y	LUAD	387	387
	# Tx	Y	LUSC	342	342


# §3 FIE counts by timing PRE-DUPLICATION('early') and POST-DUPLICATION('late')
# -----------------------------------------------------------------------------
	# FIE count by timing groups / percs
		# (a) Timing FIE count by TXP 
			dat <- dt_fiet %>%
				group_by( TXP ) %>%
				mutate( tot_fies_group = n_distinct(m_id) ) %>%
				ungroup() 
			# summarise by timing with % in each (by txp/cancer type)
			dats <- dat %>%
				group_by( TXP, OVERALL_TIMING, tot_fies_group )
			out_summary_perc( dats )
		
		# (b) Timing FIE count by TXP and cancer type
			dbtc <- dt_fiet %>%
				group_by( TXP, CANCER_TYPE ) %>%
				mutate( tot_fies_group = n_distinct(m_id) ) %>%
				ungroup() 
			# summarise by timing with % in each (by txp/cancer type)
			dbtcs <- dbtc %>%
				group_by( TXP, CANCER_TYPE, OVERALL_TIMING, tot_fies_group )
			out_summary_perc( dbtcs )
	
		# (c) Timing FIE count by TXP and CGC
			dctc <- dt_fiet %>%
				group_by( TXP, CGC ) %>%
				mutate( tot_fies_group = n_distinct(m_id) ) %>%
				ungroup() 
			# summarise by timing with % in each (by txp/cancer type)
			dctcs <- dctc %>%
				group_by( TXP, CGC, OVERALL_TIMING, tot_fies_group )
			out_summary_perc( dctcs )
	
# §3S: FIE counts by timing PRE- POST-DUPLICATION SUMMARY TABLES
# ---------------------------------------------------------------
	# ** (a) Timing FIE count by TXP ** 
	   	#   TXP   OVERALL_TIMING tot_fies_group num_fies num_genes num_tumours num_ff fies_per_tumour fies_per_gene num_rows fie_perc_group
		# 1 N     early                     355      191        41         148     38            1.29          4.66      191           53.8
		# 2 N     late                      355       62        43          55     39            1.13          1.44       62           17.5
		# 3 N     unknown                   355      102        50          73     47            1.40          2.04      102           28.7
		# 4 Y     early                     446      131        75         114     73            1.15          1.75      131           29.4
		# 5 Y     late                      446       19        19          19     19            1             1          19            4.3
		# 6 Y     unknown                   446      296       140         251    133            1.18          2.11      296           66.4

	# ** (b) Timing FIE count by TXP and cancer type **
		#    TXP   CANCER_TYPE OVERALL_TIMING tot_fies_group num_fies num_genes num_tumours num_ff fies_per_tumour fies_per_gene num_rows fie_perc_group
		#  1 N     LUAD        early                     264      140        29         108     27            1.30          4.83      140           53  
		#  2 N     LUAD        late                      264       38        27          33     26            1.15          1.41       38           14.4
		#  3 N     LUAD        unknown                   264       86        44          61     40            1.41          1.95       86           32.6
		#  4 N     LUSC        early                      91       51        18          40     17            1.27          2.83       51           56  
		#  5 N     LUSC        late                       91       24        18          22     18            1.09          1.33       24           26.4
		#  6 N     LUSC        unknown                    91       16        12          12     14            1.33          1.33       16           17.6
		#  7 Y     LUAD        early                     212       74        43          64     43            1.16          1.72       74           34.9
		#  8 Y     LUAD        late                      212       10        10          10     10            1             1          10            4.7
		#  9 Y     LUAD        unknown                   212      128        75         114     76            1.12          1.71      128           60.4
		# 10 Y     LUSC        early                     234       57        35          50     35            1.14          1.63       57           24.4
		# 11 Y     LUSC        late                      234        9         9           9      9            1             1           9            3.8
		# 12 Y     LUSC        unknown                   234      168        74         137     69            1.23          2.27      168           71.8

	# ** (c) Timing FIE count by TXP and CGC **
		#    TXP   CGC   OVERALL_TIMING tot_fies_group num_fies num_genes num_tumours num_ff fies_per_tumour fies_per_gene num_rows fie_perc_group
		#  1 N     N     early                      91       23        23          20     20            1.15          1          23           25.3
		#  2 N     N     late                       91       33        33          29     30            1.14          1          33           36.3
		#  3 N     N     unknown                    91       35        34          30     32            1.17          1.03       35           38.5
		#  4 N     Y     early                     264      168        18         140     20            1.2           9.33      168           63.6
		#  5 N     Y     late                      264       29        10          27     11            1.07          2.9        29           11  
		#  6 N     Y     unknown                   264       67        16          53     17            1.26          4.19       67           25.4
		#  7 Y     N     early                     188       57        56          51     52            1.12          1.02       57           30.3
		#  8 Y     N     late                      188       12        12          12     12            1             1          12            6.4
		#  9 Y     N     unknown                   188      119       110         101    106            1.18          1.08      119           63.3
		# 10 Y     Y     early                     258       74        19          71     21            1.04          3.89       74           28.7
		# 11 Y     Y     late                      258        7         7           7      7            1             1           7            2.7
		# 12 Y     Y     unknown                   258      177        30         174     32            1.02          5.9       177           68.6

# §3S:FIE counts by timing PRE- POST-DUPLICATION SUMMARY TABLES SPECIFIC CALCS
# ----------------------------------------------------------------------------
	# (a) % early FIEs that are CGC
		d3sa <- dctc %>%
				filter( OVERALL_TIMING=='early' ) %>%
				filter( TXP == 'N' ) %>%
				group_by( TXP, CGC, tot_fies_group )
		
		# early FIE count Tx 191
		tmp <- out_summary_perc(d3sa)
		tot_early <- sum(tmp$num_fies)
		d3sa <- mutate( tmp, perc_fie_early = round(100 * (num_fies / tot_early), 1 ) )
		d3sa
 
	# (b) % early / late overall by gene (TX421)
		d3sb 	<- dt_fiet %>%
			filter( TXP == 'N' ) %>%
			group_by( OVERALL_TIMING ) %>%
			mutate( tot_fies_group = n_distinct(m_id) )
		
		d3sbs <- d3sb %>%
			group_by( TXP, OVERALL_TIMING, SOURCE_HUGO_SYMBOL, tot_fies_group ) %>%
			out_summary_perc()
		
	# (c) Tx genes early / late / (+ unknown timing) with >4 FIEs per gene
		d3sc <- d3sb %>%
			group_by( OVERALL_TIMING, SOURCE_HUGO_SYMBOL, tot_fies_group ) %>%
			out_summary_perc( ) %>%
			filter( fies_per_gene > 4 )
			

############################################################################################
# st01a FIEs highly enriched early
# uses diversity table [loaded as dt_div_gene]
# which has mutations (missense and synonymous) and FIEs for Tx and TCGA-lung
############################################################################################
# FIEs grouped by data source, cancer type and timing
dt_gene_fie <- dt_div_gene %>%
	filter( SCORE_THRESHOLD == 'MC3' ) %>%
	group_by( SPECIES_TYPE, DATA_SOURCE, TXP, CANCER_TYPE, OVERALL_TIMING ) %>%
	mutate( tot_fies = sum( NUM_MUTS ) ) %>%
	select( -c( 'VARIANT_CLASS', 'PYCLONECLONAL_SC', 'ANY_REGION_DUPLICATED', 'DRIVERMUT', 'MUTANT_EXPRESSED', 'SPECIES', 'NUM_MUTS', 'NUM_ROWS' )) %>%
	distinct()

# synonymous variants in mapped genes grouped by data source, cancer type and timing
# these are genes in which FIEs could in theory be found (3D available, has func site etc)
dt_gene_SYN <- dt_div_gene %>%
	filter( is.na( SCORE_THRESHOLD ) ) %>%
	filter( VARIANT_CLASS == 'Silent' ) %>%
	group_by( SPECIES_TYPE, DATA_SOURCE, TXP, CANCER_TYPE, OVERALL_TIMING ) %>%
	mutate( tot_syns = sum( NUM_MUTS ) ) %>%
	select( -c( 'VARIANT_CLASS', 'PYCLONECLONAL_SC', 'ANY_REGION_DUPLICATED', 'DRIVERMUT', 'MUTANT_EXPRESSED', 'SPECIES', 'NUM_MUTS', 'NUM_ROWS' ) ) %>%
	distinct( )

# Tx LUAD
stat_st01(dt_gene_fie, dt_gene_SYN, txp='N', cancer_type = 'LUAD' )
# row_cat FIEs Other mutations
#   Early  136             678
#   Late    48            1214

#         Fisher's Exact Test for Count Data
# data:  mat_test
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  3.569965 7.298476
# sample estimates:
# odds ratio 
#   5.069004 

# Tx LUSC
stat_st01(dt_gene_fie, dt_gene_SYN, txp='N', cancer_type = 'LUSC' )
#        col_cat
# row_cat FIEs Other mutations
#   Early   51             535
#   Late    24             533

#         Fisher's Exact Test for Count Data

# data:  mat_test
# p-value = 0.002744
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.257426 3.649651
# sample estimates:
# odds ratio 
#    2.11573 


