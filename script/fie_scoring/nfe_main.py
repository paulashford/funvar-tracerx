# nfe_main.py 
# 24/06/2019 and Dec '19, Jan '20...
# Main call-point for FunVar NFE (FIE) score calculations

# Nov '23 version for publication 
#   - remove any hard coded paths)
#   - create data export from database.py 'nfe_pfh_v090' SQL ->data/nfe_score_database_py_funvar_archive_nfe_pfh_mutfam_090_exported.tsv
#   - uses HID mutation id - see data/nfe_score_database_py_export.sql and class files

#-----------------------------------------------------------------------
# (1) Imports / initialise
#-----------------------------------------------------------------------
import sys
import os
import pandas as pd
import math
# import altair as alt
from IPython.display import display

# Useful!
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 100)

# McLachlan + Grantham index data
AA_INDEX_PATH = '../../data/aaindex'

# data out dir for scored FIEs
nfe_dir = 'working'

# append packages for classes
sys.path.append('./packages')

# Database N/A now a flat file
# from database import OracleConnection as DbOracleConnection
# from database import SelectStatements as DbSelectStatements

# Mutations
from mutations import Mutation as Mutation
# Scoring FIEs
from neofunevent import score_nfe_by_mutation

# -----------------------------------------------------------------
# # (2) Oracle connection - connect to orengodev1 (Oracle 12c)
# #-----------------------------------------------------------------------
# N/A for submission version

#-----------------------------------------------------------------------
# (3) SQL select for NFEs
#-----------------------------------------------------------------------
# SQL - N/A for submission version

# Use exported file
sql_data_exported='../../data/nfe_score_database_py_funvar_archive_nfe_pfh_mutfam_090_exported_no_mut_id.tsv'

pd_nfe = pd.read_table(sql_data_exported)
# pd_nfe.head()
# pd_nfe.info()

#-----------------------------------------------------------------------
# (4) FUNCTION: Iterate mutations  
#-----------------------------------------------------------------------
# This will calculate McLachlan, delta size and polymorphic properties.
# Threshold for deciding polymorphism
def iterate_mutations( pd_nfe, poly_thresh = 0.00001 ):
#POLY_CUTOFF = 0.00001
    # Get mutation list (includes PFH and MutClust from NFE_v6)
    mutations = [ Mutation( 
                    nfe_type = row.NFE_TYPE,
                    data_source = row.DATA_SOURCE,
                    # mutation_id = row.MUTATION_ID,
                    mutation_id = row.HID,
                    hid = row.HID,
                    cancer_type = row.CANCER_TYPE,
                    variant_type = row.VARIANT_TYPE,
                    variant_class = row.VARIANT_CLASS,
                    hugo_symbol = row.SOURCE_HUGO_SYMBOL,
                    uniprot_acc = row.VM_UNIPROT_ACCESSION,
                    uniprot_seq_no = row.VM_SEQ_NO,
                    uniprot_aa_change = row.VM_AA_CHANGE,
                    num_patients = row.NUM_PATIENTS,
                    mut_count_sum_p_corr = row.MUT_COUNT_SUM_P_CORR,
                    weighted_mut_sum_p_corr = row.WEIGHTED_MUT_SUM_P_CORR,
                    gnomad_af = row.GNOMAD_AF_CONVERTED,
                    diseases = row.DISEASES,
                    aa_index_path = AA_INDEX_PATH,
                    poly_threshold = poly_thresh
                    ) 
                    for row in pd_nfe.itertuples()
                ]     
    return mutations
    
#-----------------------------------------------------------------------
# (5) FUNCTION: Collate Mutations with calculated McLachlan etc properties
#-----------------------------------------------------------------------
def mutation_properties( mutations ):
    mutation_props = [ {
                        'nfe_type':     mutation.nfe_type,
                        'data_source':  mutation.data_source,
                        # 'mutation_id':  mutation.mutation_id,
                        'mutation_id':  mutation.hid,
                        'hid':          mutation.hid,
                        'cancer_type':  mutation.cancer_type,
                        'variant_type': mutation.variant_type,
                        'variant_class': mutation.variant_class,
                        'hugo_symbol':  mutation.hugo_symbol,
                        'uniprot_acc':  mutation.uniprot_acc,
                        'uniprot_seq_no': mutation.uniprot_seq_no,
                        'uniprot_aa_change': mutation.uniprot_aa_change,
                        'aa_from':      mutation.uniprot_aa_change[0],
                        'aa_to':        mutation.uniprot_aa_change[2],
                        'num_patients': mutation.num_patients,
                        'mut_count_sum_p_corr': mutation.mut_count_sum_p_corr,
                        'weighted_mut_sum_p_corr': mutation.weighted_mut_sum_p_corr,
                        'gnomad_af':    mutation.gnomad_af,
                        'diseases':     mutation.diseases,
                        'mclachlan':    mutation.mcLachlan, 
                        'grantham':     mutation.grantham,
                        'deltasize':    mutation.aminoSizeChg,
                        'ispoly':       mutation.isPolymorphic
                        } for mutation in mutations
                    ]
    return mutation_props

#-----------------------------------------------------------------------
# (6) Calculate mutation properties for all mutations & 
#       merge with original data frame
#-----------------------------------------------------------------------
mutations = iterate_mutations( pd_nfe )
mutation_props = mutation_properties( mutations )

pd_nfe_props = pd.merge(    
                            pd_nfe,
                            pd.DataFrame( mutation_props ),
                            how = 'inner', 
                            left_on = 'HID',
                            right_on='mutation_id' 
                        )

#-----------------------------------------------------------------------
# (7) Calculate NFE scores using the full pandas df of mutations and their properties
#-----------------------------------------------------------------------
# What is Q3 size change for all mutations in this set?  (70.9 at time of running - NFE_v6)
delta = pd_nfe_props.deltasize.abs().quantile(0.75)

# 01 09 2022 Grantham scores *include size change calculation*
#x-scrivener-item:///Users/ash/Dropbox/bioinf/neofun/scrivener/neofun_v3.1_sectioning.scriv?id=9ACA4C6E-B6FA-4387-8F2A-F7A313A8641B
pd_nfe_props[['mclachlan','grantham']].describe()
# mclachlan	grantham
# count	16756.000000	16756.000000
# mean	0.508236	81.687873
# std	0.951459	50.787531
# min	0.000000	5.000000
# 25%	0.000000	43.000000
# 50%	0.000000	64.000000
# 75%	1.000000	109.000000
# max	3.000000	215.000000
# Grantham based scores will ignore amino size change and calc using:
# gs > median           -> score+1
# gs > 75 percentile    -> score+2
g_vals = [pd_nfe_props.grantham.abs().quantile(0.50), pd_nfe_props.grantham.abs().quantile(0.75)]

pd_nfe_props['nfe_score'] = 0
pd_nfe_props['nfe_score_d'] = 0
pd_nfe_props['nfe_score_d_mf'] = 0
pd_nfe_props['nfe_score_d_mf_h'] = 0

# GRANTHAM-based (pass grantham_vals = [median, 75%] )
for nfe_row in pd_nfe_props.itertuples():
    # Unbiased score (simple score)
    nfe_row_score_simple = score_nfe_by_mutation( nfe_row, deltas=delta, score_type='simple', mc_p_cutoff=0.001, grantham_vals=g_vals )
    pd_nfe_props.loc[pd_nfe_props.mutation_id == nfe_row.mutation_id ,'nfe_score'] = nfe_row_score_simple
    
    # + diseases
    nfe_row_score_d = score_nfe_by_mutation( nfe_row, deltas=delta, score_type='s_d', mc_p_cutoff=0.001, grantham_vals=g_vals )
    pd_nfe_props.loc[pd_nfe_props.mutation_id == nfe_row.mutation_id ,'nfe_score_d'] = nfe_row_score_d
    
    # + diseases + MutFam
    nfe_row_score_d_mf = score_nfe_by_mutation( nfe_row, deltas=delta, score_type='s_d_mf', mc_p_cutoff=0.001, grantham_vals=g_vals )
    pd_nfe_props.loc[pd_nfe_props.mutation_id == nfe_row.mutation_id ,'nfe_score_d_mf'] = nfe_row_score_d_mf
    
    # + diseases + MutFam + hotspot
    nfe_row_score_d_mf_h = score_nfe_by_mutation( nfe_row, deltas=delta, score_type='s_d_mf_h', mc_p_cutoff=0.001, grantham_vals=g_vals )
    pd_nfe_props.loc[pd_nfe_props.mutation_id == nfe_row.mutation_id ,'nfe_score_d_mf_h'] = nfe_row_score_d_mf_h

#-----------------------------------------------------------------------
# (8) Write scored NFE's to pickle & Excel
#-----------------------------------------------------------------------
# Float->int format changes
pd_nfe_props[['NEAR_SCONS_90', 'NEAR_MCSA', 'NEAR_NUC', 'NEAR_LIG']] = pd_nfe_props[['NEAR_SCONS_90', 'NEAR_MCSA', 'NEAR_NUC', 'NEAR_LIG']].fillna(value=0)
pd_nfe_props = pd_nfe_props.astype( {       "mclachlan": int, 
                                            "grantham":int,
                                            "nfe_score":int,
                                            "nfe_score_d":int,
                                            "nfe_score_d_mf":int,
                                            "nfe_score_d_mf_h":int,  
                                            "NEAR_ANGSTROMS":int,
                                            "NEAR_SCONS_90":int,
                                            "NEAR_MCSA":int,
                                            "NEAR_NUC":int,
                                            "NEAR_LIG":int 
                                            })
# Float format 10f (fine for gnomAD and MutClust P vals!)
pd.options.display.float_format = '{:.10f}'.format
pd_nfe_props.dtypes

# Save FIE scores file
fie_scored_file = 'funvar_fie_scores_v090_nmid.tsv'
# EXPORT tsv
pd_nfe_props.to_csv( 
                    os.path.join( nfe_dir, fie_scored_file ), 
                    sep='\t', 
                    index=False, 
                    float_format='%.10f' 
                    )
# Save Pickle too
fie_scored_pickle_file = 'funvar_fie_scores_v090_nmid.pkl'
pd_nfe_props.to_pickle( os.path.join( nfe_dir, fie_scored_pickle_file ) )
