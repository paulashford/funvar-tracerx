# neofunevent.py
'''
Functional Impact Events (FIEs) - score based on Mutation properties
such as gnomAD polymorhisms, McLachlan score and size change
10/12/2019

01 09 2022
Grantham scores (see nfe_main §7)
For gratham-based scoring amino size change is not used, instead if 
gs > median -> score+1
gs > 75 percentile -> score+2
'''
#import operator

'''
FIE score function should return identical scores for mutations in different patients/tumours that share:
- same gene (uniprot)
- same sequence position 
- same amino acid change
These will also therefore share:
- same FunFam
- same functional site proximities

They may not share the same gnomAD score (triplet codon redundancy), but in practice (by inspection of scores in hotspots) they seem to.

As a guiding principle this means that aggregate calculations, such as hotspot mutation frequency, 
should not be included, as then the scores would vary across datasets.

Highly significant clusters....
As MutClust residues are filtered with this (from VW_NFE_MUTCLUSTS_ONLY in JOIN .mutclust_analysis  ):
    mut_count_sum_p_corr < 0.05 OR weighted_mut_sum_p_corr < 0.05

A sensible approach is to use the weighted corrected p-val (has literature support)
with a stringent cut-off (p<0.001) - most of the unweighted p-values are then <0.05 too (5 not in LUAD+LUSC). E.g.:
    
'''
import math

# Load pre-merged nfe data (without any scoring)
def load_nfe_merge_data(nfe_pickle_file):
    import pandas as pd
    try:
        pd_nfe_calc = pd.read_pickle(nfe_pickle_file)
    except Exception:
        print( "Unable to open pickle file: " + nfe_pickle_file )
    else:
        return pd_nfe_calc

# deltas: deltaSize threshold 
#   defaults to a pre-calculated Q3 default, for a previous version of NFE data, assuming it doesn't change much with more data.
def score_nfe_by_mutation(nfe_row, deltas = 70.9, score_type = 'simple', mc_p_cutoff = 0.005, grantham_vals = [999,999]):
    # NOTE grantham = [median grantham of all FIEs, 75% grantham of all FIEs]
    # Simple score.
    # Is Grantham median and 75% passed in?
    # ... no
    if grantham_vals == [999,999]:
        tot_nfe_by_mutation = (
                            __score_mclachlan(nfe_row) + 
                            __score_deltasize(nfe_row, deltas) + 
                            __score_fs_on_pred(nfe_row) + 
                            __score_fs_on_known(nfe_row) + 
                            __score_fs_near_pred(nfe_row) +
                            __score_fs_near_known(nfe_row) +
                            __highly_sig_mutclust(nfe_row, mc_p_cutoff)
                            #__score_poly(nfe_row) +  # This is accounted for in McLach and deltaS
                            )
    else:
    # yes
        tot_nfe_by_mutation = (
                            __score_grantham(nfe_row, grantham_vals) +
                            # __score_mclachlan(nfe_row) + 
                            # __score_deltasize(nfe_row, deltas) + 
                            __score_fs_on_pred(nfe_row) + 
                            __score_fs_on_known(nfe_row) + 
                            __score_fs_near_pred(nfe_row) +
                            __score_fs_near_known(nfe_row) +
                            __highly_sig_mutclust(nfe_row, mc_p_cutoff)
                            #__score_poly(nfe_row) +  # This is accounted for in McLach and deltaS
                            )
    # Add in other terms if unbiased_score == False (see script notes)
    if score_type == 'simple':
        pass
    # simple plus diseases
    elif score_type == 's_d':
        tot_nfe_by_mutation += __score_diseases(nfe_row)
                              
    # simple plus diseases + mutfam
    elif score_type == 's_d_mf':
        tot_nfe_by_mutation += (
                            __score_diseases(nfe_row) +
                            __score_mutfam(nfe_row)
                            )    
    
    # simple plus diseases + mutfam + hotspot
    # These are needed if taking MAX score per amino/amino_change on a rep_id
    # Essentially, MutClust and PFH *may* be double counting, in cases where MutClust has already found 
    # a 'hotspot'...
    elif score_type == 's_d_mf_h':
        tot_nfe_by_mutation += (
                            __score_diseases(nfe_row) +
                            __score_mutfam(nfe_row) +
                            __score_pfh(nfe_row)
                            )            
    else:
        pass
                
    return tot_nfe_by_mutation

'''
POSITIVE TERMS
'''
# NFEs (of either type "MC" or "PFH")

# ON known func sites
def __score_fs_on_known(nfe_row):
    s_fs_on_known = 0
    if nfe_row.ON_MCSA or nfe_row.ON_NUC or nfe_row.ON_LIG or nfe_row.ON_PPI:
        s_fs_on_known = 1
    return s_fs_on_known
# NEAR known func sites (doesn't include PPI)
def __score_fs_near_known(nfe_row):
    s_fs_near_known = 0
    if nfe_row.NEAR_MCSA or nfe_row.NEAR_NUC or nfe_row.NEAR_LIG:
        s_fs_near_known = 1
    return s_fs_near_known

# ON predicted func site
def __score_fs_on_pred(nfe_row):
    s_fs_on_pred = 0
    if nfe_row.ON_SCONS_90:
        s_fs_on_pred = 1
    return s_fs_on_pred
# NEAR predicted func site (doesn't include PPI)
def __score_fs_near_pred(nfe_row):
    s_fs_near_pred = 0
    if nfe_row.NEAR_SCONS_90:
        s_fs_near_pred = 1
    return s_fs_near_pred

# PFHs (these are reported as NFEs of type 'PFH')
def __score_pfh(nfe_row):   
    s_pfh = 0
    if nfe_row.NFE_TYPE == 'PFH':
        s_pfh = 1
    return s_pfh

# Physico-chemical properties
# If polymorphic, this score is not set
def __score_mclachlan(nfe_row):
    s_mclachlan = 0
    if (nfe_row.mclachlan == 0) and (not nfe_row.ispoly):
        s_mclachlan = 1
    return s_mclachlan

# Grantham (see notes at top - note overrides mclachlan and amino size chg if used)
# If polymorphic, this score is not set
def __score_grantham(nfe_row, grantham_vals):
    s_grantham = 0
    # is grantham > 75% percentile of all putative-FIE granthams?
    if (nfe_row.grantham > grantham_vals[1]) and (not nfe_row.ispoly):
        s_grantham = 2
    # is grantham > median of all putative-FIE granthams?
    if ( (nfe_row.grantham > grantham_vals[0]) and (nfe_row.grantham <= grantham_vals[1]) ) and (not nfe_row.ispoly):
        s_grantham = 1
    return s_grantham

# If polymorphic, this score is not set
# Size change can be negative, so check absolute value >= delta (based on Q3 of abs values for all NFEs)
def __score_deltasize(nfe_row, deltas):
    s_deltasize = 0
    if (abs(nfe_row.deltasize) >= deltas) and (not nfe_row.ispoly):
        s_deltasize = 1
    return s_deltasize

# Is this a *highly significant* MutClust?
# See script comments
def __highly_sig_mutclust(nfe_row, p_cutoff):
    s_high_sig_mc = 0
    if (nfe_row.NFE_TYPE == 'MC') and (nfe_row.weighted_mut_sum_p_corr < p_cutoff):
        s_high_sig_mc = 1
    return s_high_sig_mc

# Other diseases noted?
def __score_diseases(nfe_row):
    s_diseases = 0
    if nfe_row.DISEASES != '-':
        s_diseases = 1
    return s_diseases

# Is MutFam?
def __score_mutfam(nfe_row):
    s_mutfam = 0
    if nfe_row.MUTFAM == 'Y':
        s_mutfam = 1
    return s_mutfam

'''
NEGATIVE TERMS
'''
# def __score_poly(nfe_row):
#     s_poly = 0
#     if nfe_row.ispoly:
#         s_poly = -1
#     return s_poly


