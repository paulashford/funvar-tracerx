# script/benchmark/bench_data_paths.R
# Sept 2025
# P. Ashford
# Defines paths for current resource datasets used in benchmarking
# and required by: script/benchmark/load_bench_datasets.R
# Caller: script/benchmark/run_benchmarks.R

# TCGA pancancer FIEs
simp_mt_tcga_missense 	<- 	'simple_muttab_tcga_missense.tsv'
simp_mt_tcga_mc3		<-	'simple_muttab_tcga_mc3.tsv'
simp_mt_tcga_fie_scores	<-	'simple_muttab_tcga_fie_score.tsv'

# Simplified formats of mutation/SNP datasets
# clinvar
simp_mt_clinvarb		<- 'simple_muttab_clinvar_benign.tsv'
# dbSNP
simp_mt_dbsnp			<- 'simple_muttab_dbsnp_benign.tsv'
# COSMIC cancer mutation census (CMC) with clinical significance level 1, 2, or 3.
simp_mt_cosmic123		<- 'simple_muttab_cosmic123.tsv'

# Gene/protein cross-refs
hg_dl					<- '20240731'
hg_gene_map_plus		<- 	paste0('hugo_id_map_plus_', hg_dl, '.tsv')
hg_core					<- 	paste0('hugo_core_', hg_dl, '.tsv')
# expanded version to allow mapping to / from any alias or previous name
hg_core_exp				<- 	paste0('hugo_core_exp_', hg_dl, '.tsv')

# Simplified GDC datasets
# hotspot3d
simp_mt_gdc_hotspot3d 	<- 'simple_muttab_gdc_hotspot3d.tsv'
# hotmaps
simp_mt_gdc_hotmaps 	<- 'simple_muttab_gdc_hotmaps.tsv'
# 3Dhotspots
simp_mt_gdc_3dhotspots 	<- 'simple_muttab_gdc_3dhotspots.tsv'

# FunFams (v4.2), CATH (v4.2) representative domains and amino acid range boundaries
funfam_ranges		<- 'human_funfam_member_ranges_v4_2.tsv'
funfam_info			<- 'cath_funfam_all_v4_2.tsv'

