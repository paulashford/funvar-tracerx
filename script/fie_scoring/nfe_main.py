# nfe_main.py 
'''
Paul Ashford  
24/06/2019

 24/06/2019 and Dec '19, Jan '20...
 Main call-point for FunVar NFE (FIE) score calculations

 Nov '23 version for publication 
   - remove any hard coded paths)
   - create data export from database.py 'nfe_pfh_v090' SQL ->data/nfe_score_database_py_funvar_archive_nfe_pfh_mutfam_090_exported.tsv
   - uses HID mutation id - see data/nfe_score_database_py_export.sql and class files
 Nov '25 tidy imports and refactor 
'''

# Imports / initialise
import sys
import os
import pandas as pd
# root package path
sys.path.append('/Users/ash/git/funvar-tracerx/script')
# Mutation class
from fie_scoring.packages.mutations import Mutation as Mutation
# FIE scoring
from fie_scoring.packages.neofunevent import score_nfe_by_mutation
# McLachlan + Grantham index data
AA_INDEX_PATH = 'data/aaindex'
# data out dir for scored FIEs
nfe_dir = 'script/fie_scoring/working'

def main():
    # SQL - N/A for submission version
    print('FunVar-FIE scoring')
    print('Loading annotated mutations (includes clusters and functional site proximity)... ')
    # Use exported file
    sql_data_exported='data/nfe_score_database_py_funvar_archive_nfe_pfh_mutfam_090_exported_no_mut_id.tsv'
    # sql_data_exported='../../data/nfe_pfh_mutfam_090_af_for_fie_scoring_INC_MUTID.tsv'
    pd_nfe = pd.read_table(sql_data_exported)

    # Calculate mutation properties for all mutations & merge with original data frame
    print('Calculating mutation properties...')
    mutations = iterate_mutations( pd_nfe )
    mutation_props = mutation_properties( mutations )
    pd_nfe_props = pd.merge(    
                            pd_nfe,
                            pd.DataFrame( mutation_props ),
                            how = 'inner', 
                            left_on = 'HID',
                            right_on='mutation_id' 
                        )
    # Calculate FIE scores using the full df of mutations and properties
    print('Calculating FIE-scores...')

    # What is Q3 size change for all mutations in this set?  (70.9 at time of running - NFE_v6)
    delta = pd_nfe_props.deltasize.abs().quantile(0.75)

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

    # Write scored FIEs to tab-sep file
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

    # Output tsv
    fie_scored_file = 'funvar_fie_scores_output.tsv'
    print('Writing tab-separated file: ' + os.path.join( nfe_dir, fie_scored_file ))
    pd_nfe_props.to_csv( 
                        os.path.join( nfe_dir, fie_scored_file ), 
                        sep='\t', 
                        index=False, 
                        float_format='%.10f' 
                        )
    

# Function: Iterate mutations  
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
    
# Function: Collate Mutations with calculated McLachlan etc properties
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


if __name__ == '__main__':
    main()
