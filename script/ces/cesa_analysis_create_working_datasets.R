# cesa_analysis_create_working_datasets.R
# CES cancer effect size analysis 
# Create working set for CES via MAF for TRACERx LUAD, LUSC, or TCGA pan-cancer mutations
# P.Ashford Aug 2025
# Info for cancereffectsizeR (The "Prepare CES data" section below is based on this guide)
# https://townsend-lab-yale.github.io/cancereffectsizeR/articles/cancereffectsizeR.html

# NOTE: in case of cryptic error messages, suggestion is to upgrade BioConductor e.g.
# BiocManager::install(version = "3.21") # or newer, when available

options(timeout = 6000)
library(cancereffectsizeR)
library(data.table)
library(fst)
library(tidyverse)

# source analysis functions script
source('script/ces/cesa_functions.R')


# Parameters for type of working set
run_name  <- 'FunVar'
# Compute all variants? (we are using compound variant sets for functional families)
compute_all_vars <- FALSE

# ** TRACERx LUAD or LUSC **
# data_source <- 'Tx'
# histo_filter	<- 'LUAD'
# histo_filter	<- 'LUSC'

# ** TCGA pan-cancer **
data_source <- 'TCGA' 
histo_filter	<- 'PANCAN'

# CES tissue-specific covariate sets
if ( data_source == 'Tx' & histo_filter %in% c('LUAD','LUSC') ){
  
  # TRACERx mutation and FIE tables [via Oracle export script/sql/mut_tx_pass_with_cancer_type.sql]
  dt_tx_ora     <- as.data.table( read.delim('data/mut_tx_with_cancer_type_idh.tsv') )

  # CES tissue-specific covariate set
  covariate_set <- 'lung'

} else if (data_source == 'TCGA' & histo_filter == 'PANCAN') {
  # TCGA pancancer mutation and FIE tables [via Oracle export script/sql/mut_tcga_pass_with_cancer_type.sql]
  dt_tcga_ora   <- as.data.table( read.delim('data/mut_tcga_with_cancer_type.tsv') )

  # not tissues specific [https://townsend-lab-yale.github.io/cancereffectsizeR/reference/gene_mutation_rates.html]
  covariate_set <- 'hg19' 

} else{
  stop(paste0('Data source ', data_source, ' with histology_filter ', histo_filter, ' is not supported.'))
}

run_name <- str_c(data_source, histo_filter, run_name, sep = "-")

# ----------------------------------
# Create MAF format for CES
# ----------------------------------
# TRACERx mutations (oracle)
if (data_source == 'Tx'){
  dt_mttxo  	<- dt_tx_ora[VARIANT_TYPE == 'SNP']
  dt_mttxo[, chrom := sapply( CHR_TEXT, function(x){ sub('chr', '', x, fixed = TRUE) } )]
  dt_mttxo[PASS == 'TRUE', FILTER := 'PASS']
  dt_mttxo[, source := 'funvar_export_mut_tx']
  if (histo_filter != 'PANCAN') {dt_mttxo <- dt_mttxo[CANCER_TYPE == histo_filter] } 
  dt_txvcf <- dt_mttxo[, c('chrom', 'START_POSITION', 'TCGABARCODE_TUMOUR_ID', 'REF_ALLELE', 'VAR_ALLELE', 'FILTER', 'source')]

  # Create a VCF suitable for CES analysis
  setnames(dt_txvcf, c('chrom', 'START_POSITION', 'REF_ALLELE', 'VAR_ALLELE', 'TCGABARCODE_TUMOUR_ID'), c('#CHROM', 'POS', 'REF', 'ALT', 'tumour_id')) 
  dt_txvcf <- dt_txvcf[, c('#CHROM', 'POS', 'tumour_id', 'REF', 'ALT', 'FILTER', 'source')]    
  list_of_tables <- split(dt_txvcf, by = "tumour_id")
  dt_tx_maf <- vcfs_to_maf_table( list_of_tables, sample_ids = names(list_of_tables) )
  dt_maf <- dt_tx_maf
}

if (data_source == 'TCGA'){
  # TCGA mutations (oracle)
  # dt_mttcga  	<- dt_tcga_ora[VARIANT_CLASS %in% c('Missense_Mutation', 'Silent')]
  dt_mttcga  	<- dt_tcga_ora[VARIANT_TYPE == 'SNP']
  dt_mttcga[, chrom := as.character(CHROMOSOME)]
  dt_mttcga[STRAND == '+', strandn := 1]
  dt_mttcga[STRAND == '-', strandn := -1]
  dt_mttcga[PASS == 'PASS', FILTER := 'PASS']
  dt_mttcga[, source := 'mut_tcga_with_cancer_type.tsv']
  if (histo_filter != 'PANCAN') {dt_mttcga <- dt_mttcga[CANCER_TYPE == histo_filter] } 
  dt_tcga_vcf <- dt_mttcga[, c('chrom', 'START_POSITION', 'TCGABARCODE_TUMOUR_ID', 'REF_ALLELE', 'VAR_ALLELE', 'FILTER', 'source')]

  # Create a VCF suitable for CES analysis
  setnames(dt_tcga_vcf, c('chrom', 'START_POSITION', 'REF_ALLELE', 'VAR_ALLELE', 'TCGABARCODE_TUMOUR_ID'), c('#CHROM', 'POS', 'REF', 'ALT', 'tumour_id')) 
  dt_tcga_vcf <- dt_tcga_vcf[, c('#CHROM', 'POS', 'tumour_id', 'REF', 'ALT', 'FILTER', 'source')]    
  list_of_tables <- split(dt_tcga_vcf, by = "tumour_id")
  dt_tcga_maf <- vcfs_to_maf_table( list_of_tables, sample_ids = names(list_of_tables) )
  dt_maf <- dt_tcga_maf
}

# ----------------------------------
# Prepare CES data
# ----------------------------------
# Prepare data
maf <- preload_maf(maf = dt_maf, refset = "ces.refset.hg19")
# Create cancereffectsizeR analysis and load data
cesa <- CESAnalysis(refset = "ces.refset.hg19")
cesa <- load_maf(cesa = cesa, maf = dt_maf, maf_name = run_name)

# Run cesa trinucleotide, gene mutation rates and ces_variant for TRACERx or TCGA-pancancer datasets
#   note:list_ces_signature_sets() # ces.refset.hg19: COSMIC_v3, COSMIC_v3.1, COSMIC_v3.2 # ces.refset.hg38: COSMIC_v3.4, COSMIC_v3.2
# TRACERx
if ( data_source == 'Tx' & histo_filter %in% c('LUAD', 'LUSC') ) {  
  print( paste0('Setting up CESA for ', data_source, ' histology type ', histo_filter) )
  
  # Infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis 
  signature_exclusions	<- suggest_cosmic_signature_exclusions(cancer_type = histo_filter, treatment_naive = TRUE)
  cesa <- trinuc_mutation_rates(
            cesa = cesa, 
            signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
            signature_exclusions = signature_exclusions
          )
  
  # dNdSCv with lung-specififc covariates
  if (covariate_set=='lung') {
      print( paste0('COVARATE SET: LUNG') )
      cesa <- gene_mutation_rates(cesa, covariates = ces.refset.hg19$covariates$lung)
    } else {
      print( paste0('COVARATE SET hg19') )
       cesa <- gene_mutation_rates(cesa, covariates = ces.refset.hg19$covariates$hg19)
    }
  
  # Only necessary if want to calc all variants (we use compund variant sets instead - see next section)
  # Infer scaled selection coefficients under the default model of clonal selection. (All variants, not just MAF >1)
  if (compute_all_vars){ 
     print( paste0('CALCULATING ALL VARIANTS...') )
    cesa <- ces_variant(cesa, variants = cesa$variants, run_name = run_name) 
  }
    
# Pancancer (TCGA)                         
} else if (data_source == 'TCGA' & histo_filter == "PANCAN"){
  print( paste0('Setting up CESA for ', data_source, ' histology type ', histo_filter) )
  
  # signature_exclusions	<- NA
  cesa <- trinuc_mutation_rates(
            cesa = cesa, 
            signature_set = ces.refset.hg19$signatures$COSMIC_v3.2 
          ) 

  # dNdSCv with hg19 
  if (covariate_set=='hg19') {
    print( paste0('COVARATE SET hg19') )
    cesa <- gene_mutation_rates(cesa, covariates = ces.refset.hg19$covariates$hg19)
  } else {
     stop( paste0('Covarate set ', covariate_set, ' is not currently supported for pan-cancer, consult CES documentation.') )
  }

  # ----------------------------------------------------------------------------------------------------------------
  # Infer scaled selection coefficients under the default model of clonal selection. 
  # Only necessary if want to calc all variants (we use compund variant sets instead)
  # Note if just want MAF>1 use:  cesa <- ces_variant(cesa, run_name = run_name)
  if (compute_all_vars){ 
    print( paste0('CALCULATING ALL VARIANTS...') )
    cesa <- ces_variant(cesa, variants = cesa$variants, run_name = run_name) 
    }
  # ----------------------------------------------------------------------------------------------------------------

} else {
  stop("FunVar / cesa: histology type currently unsupported")
}

# SAVE WORKING SET FOR USE BY EITHER:
#   script/ces/cesa_funfam_tracerx_luad_lusc.R [For TRACERx]
#   script/ces/cesa_funfam_tcga_pancancer.R [For TCGA pan-cancer]
save.image( paste0('data/rda_ces/working_set_', run_name, '.rda') )

