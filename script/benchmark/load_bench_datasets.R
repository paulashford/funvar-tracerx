#script/benchmark/load_bench_datasets.R
# Sept 2025
# P. Ashford
# Load benchmarking datasets
# called from script/benchmark/run_benchmarks.R

# HUGO 
df_hugo_core		<- arrange( read_wrap(hg_core, bench_dir), hugo_gene_approved )
df_hugo_core_exp	<- arrange( read_wrap(hg_core_exp, bench_dir), hugo_gene_approved )

# TCGA pancancer FIEs
# TCGA pancancer score >=3
df_tc_mc3_basic		<- arrange( read_wrap(simp_mt_tcga_mc3, bench_dir), gene, mutation )
df_tc_mc3_basic		<- fn_add_res_num(df_tc_mc3_basic)
df_tc_mc3_basic_hv	<- add_hugo_valid(df_tc_mc3_basic, .check_uprot = FALSE) 

# TCGA pancanc missense muts (FF mapped 3D)
df_tc_mis_basic		<- arrange( read_wrap(simp_mt_tcga_missense, bench_dir), gene, mutation )
df_tc_mis_basic		<- fn_add_res_num(df_tc_mis_basic)
df_tc_mis_basic_hv	<- add_hugo_valid(df_tc_mis_basic, .check_uprot = FALSE) 

# Clinvar
# benign
df_clinvarb_basic 		<- arrange( read_wrap(simp_mt_clinvarb, bench_dir), gene, mutation )
df_clinvarb_basic_hv	<-	add_hugo_valid(df_clinvarb_basic, .check_uprot = FALSE)

# dbSNP benign
df_dbsnp_basic 		<- arrange( read_wrap(simp_mt_dbsnp, bench_dir), gene )	
df_dbsnp_basic		<- fn_add_res_num(df_dbsnp_basic) %>% filter(gene != 'N/A')
df_dbsnp_basic_hv	<- 	add_hugo_valid(df_dbsnp_basic, .check_uprot = FALSE)

# COSMIC Cancer Mutation Census (CMC) - clinically sig only
df_cosmic123_basic		<- arrange( read_wrap(simp_mt_cosmic123, bench_dir), gene, mutation_res_num )
df_cosmic123_basic_hv	<- 	add_hugo_valid(df_cosmic123_basic, .check_uprot = FALSE)

# GDAC pancancer 3D tools for comparison
cancer_type		 	<- 'PANCAN'

# HotSpot3D
df_gdc_hotspot3d		<- arrange( read_wrap(simp_mt_gdc_hotspot3d, bench_dir), gene, mutation_res_num )
df_gdc_hotspot3d_hv 	<- add_hugo_valid(df_gdc_hotspot3d, .check_uprot = FALSE)

#HOTMAPS
df_gdc_hotmaps			<- arrange( read_wrap(simp_mt_gdc_hotmaps, bench_dir), gene, mutation_res_num )
df_gdc_hotmaps_hv 		<- add_hugo_valid(df_gdc_hotmaps, .check_uprot = FALSE)

# 3d hotspots
df_gdc_3dhotspots		<- arrange( read_wrap(simp_mt_gdc_3dhotspots, bench_dir), gene, mutation_res_num )
df_gdc_3dhotspots_hv 	<- add_hugo_valid(df_gdc_3dhotspots, .check_uprot = FALSE)

# FunFam domain ranges on uniprots
df_ffr		<-  read_wrap(funfam_ranges, bench_dir) %>%
					select( -c(MEMBER_ID, GENE_ID, GENE_NAME, SWISSPROT) ) %>%
					rename( uniprot_acc = UNIPROT_ACC ) %>%
					arrange(SUPERFAMILY_ID, FUNFAM_NUMBER, AA_RANGE) %>%
					separate_wider_delim( 	cols = 'AA_RANGE', 
											delim = '-', 
											names = c( 'AA_RANGE_LOW', 'AA_RANGE_HIGH' ), 
											cols_remove = FALSE 
										) 
df_ffr$AA_RANGE_LOW = as.double(df_ffr$AA_RANGE_LOW)
df_ffr$AA_RANGE_HIGH = as.double(df_ffr$AA_RANGE_HIGH)
									
# FunFam info including CATH repid
df_funfam	<-  read_wrap(funfam_info, bench_dir) %>%
					select( -c(INCLUSION_E_VALUE,INCLUSION_BITSCORE,NUM_REP_ANNOTATIONS,NUM_MEMBERS_IN_SEED_ALN) ) 
						# %>% filter(REP_SOURCE_ID == 'cath') 

# FunFam annotation on ffr
df_funfam_ffr		<- left_join( df_ffr,
							df_funfam,
							by = c('SUPERFAMILY_ID' = 'SUPERFAMILY_ID', 
									'FUNFAM_NUMBER' = 'FUNFAM_NUMBER'),
							keep = FALSE, na_matches = 'never', relationship = 'many-to-one'
						) 

