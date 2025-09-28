# script/v3.2/script/neofun_v3.2_paper_stats_load.R
# Aug 2023: note for v3.2_0526 just use the supplementary tables folder as root (not calc/_yj01/v01 etc)
# Nov '23 submission - tidy/remove any hard-coded paths etc
# loads datasets
library(tidyverse)
		 
# source( 'neofun_paper_stat_functions.R' )

st_dir <- 'data'

# FIEs TRACERx and TCGA-lung (score >=3)
# st_yfa / _yj07d_0919.tsv 
# dt_fiet <- read_delim( file.path( st_dir, 'st_yfa_yj07d_0919.tsv' ), delim = '\t' )
dt_hfiet <- read_delim( file.path( st_dir, 'hst_yfa_yj07d_0919.tsv' ), delim = '\t' )

# st_yj03
dt_yj03 <- read_delim( file.path( st_dir, 'st_yj03_0919.tsv' ), delim = '\t' )

# diversity table (has mutations (missense and synonymous) and FIEs for Tx and TCGA-lung)
dt_div_gene <- read_csv( file.path( st_dir, 'diversity_yj07f_full_table_gene.csv' )  )

