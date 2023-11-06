# mutations.py
'''
Mutation classes including basic Mutation, those mapped to CATH in MutationMap
14/06/2019 and 10/12/2019

23/01/2020: Simplified Mutation class to remove GD status and timing, but add nfe_type and num_patients.
PFHs contain multiple mutations with same AA change at same seq position
(i.e. from different patients). However, certain properties will be the same for each mutation,
so to keep simple Mutation class has nfe_type: {'MC', 'PFH'} and num_patients. 
For 'MC' (MutClust), num_patients = 1.  GD status and mutation timing will vary depending on patient
tumour, so have been removed, rather than the more 'proper' solution of creating collections of mutations
for PFHs, where 100s of instances would essentially contain the same info...

01 09 2022
Added Grantham similarity in addition to McLachlan - copied to covid_funvar code
git/funvar/covid_funvar_pipeline/py/funvar_scoring/packages/funvar_covid/mutations.py
Note: amino acid mutation scoring matrices are here:
https://www.genome.jp/aaindex/AAindex/list_of_matrices
'''
import operator
from resources import aaindex2b as aaindex

class Mutation():
    '''
    Mutation core properties such as IDs, variant type, AA change...
    '''
    NFE_TYPES = ('Mutation','MC', 'PFH')
    DATA_SOURCES = ('Tx', 'TCGA')
    VARIANT_TYPES = ('SNP')
    VARIANT_CLASSES = ('Missense_Mutation', 'Silent')
    
    # McLachlan similarity / size change indexes and Grantham
    # https://www.genome.jp/aaindex/AAindex/list_of_matrices
    # 16 03 2022 added Grantham (GRAR740104 Chemical distance (Grantham, 1974))
    MCLACH_INDEX_TYPES = ('MCLA710101','MCLA720101', 'GRAR740104') 
    
    # Threshold for polymorphisms
    POLY_THRESHOLD = 0.0000001 
    
    def __init__(self, 
                 nfe_type, 
                 data_source, 
                 mutation_id, 
                 hid, 
                 cancer_type, 
                 variant_type, 
                 variant_class,
                 hugo_symbol, 
                 uniprot_acc, 
                 uniprot_seq_no, 
                 uniprot_aa_change, 
                 num_patients,
                 mut_count_sum_p_corr,
                 weighted_mut_sum_p_corr,
                 gnomad_af, 
                 diseases, 
                 aa_index_path, 
                 poly_threshold = POLY_THRESHOLD, 
                 mclach_index_types = MCLACH_INDEX_TYPES):
           
        ''' Create new Mutation '''
        self.nfe_type = nfe_type
        self.data_source = data_source
        self.mutation_id = mutation_id
        self.hid = hid
        self.cancer_type = cancer_type
        self.variant_type = variant_type        # SNP
        self.variant_class = variant_class      # Missense_Mutation / Silent
        self.hugo_symbol = hugo_symbol
        self.uniprot_acc = uniprot_acc
        self.uniprot_seq_no = uniprot_seq_no
        self.uniprot_aa_change = uniprot_aa_change
        self.num_patients = num_patients
        self.mut_count_sum_p_corr = mut_count_sum_p_corr
        self.weighted_mut_sum_p_corr = weighted_mut_sum_p_corr
        self.gnomad_af = gnomad_af
        self.diseases = diseases
        self.aa_index_path = aa_index_path

        # Initialise McLach/Grantham scores
        self.mcLachlan = self.__calculateMcLach()
        self.grantham = self.__calculateGrantham()
        # Initialise polymorphic variant flag
        self.poly_threshold = poly_threshold
        self.isPolymorphic = self.__calculate_polymorphic()
        # Initialise Amino size change
        self.aminoSizeChg = self.__calculateAminoSizeChg()
    
    @property
    def nfe_type(self):
        return self.__nfe_type
    @nfe_type.setter
    def nfe_type(self,x):
        if x not in self.NFE_TYPES:
            raise ValueError("%s is not a valied NFE TYPE." % x)
        else:
            self.__nfe_type = x

    @property
    def data_source(self):
        return self.__data_source
    @data_source.setter 
    def data_source(self, x):
        if x not in self.DATA_SOURCES:
            raise ValueError("%s is not a valid mutation data source." % x)
        else:
            self.__data_source = x
        
    @property
    def variant_type(self):
        return self.__variant_type
    @variant_type.setter 
    def variant_type(self, x):
        if x not in self.VARIANT_TYPES:
            raise ValueError("%s is not a valid mutation variant type." % x)
        else:
            self.__variant_type = x
    
    @property
    def variant_class(self):
        return self.__variant_class
    @variant_class.setter 
    def variant_class(self, x):
        if x not in self.VARIANT_CLASSES:
            raise ValueError("%s is not a valid mutation variant class." % x)
        else:
            self.__variant_class = x

    def __calculate_polymorphic(self):
        if self.gnomad_af >= self.poly_threshold:
            return True
        else:
            return False

    def __calculateMcLach(self, indexType = 'MCLA720101'):
        if indexType not in self.MCLACH_INDEX_TYPES:
            raise ValueError( "%s is not a valid McLachlan index type." % indexType)

        from resources import aaindex2b as aaindex
        aaindex.init(path = self.aa_index_path, index = '2')
        mcIndex = aaindex.get('MCLA720101')  # later version of McLach scores
        
        # Calculate McLach score for mutation
        #self.__mcLachlan  = mcIndex.get( self.uniprot_aa_change[0], self.uniprot_aa_change[2] )
        return mcIndex.get( self.uniprot_aa_change[0], self.uniprot_aa_change[2] )
    
    def __calculateGrantham(self, indexType = 'GRAR740104'):
        if indexType not in self.MCLACH_INDEX_TYPES:
            raise ValueError( "%s is not a valid amino acid index type." % indexType)

        # from resources import aaindex2b as aaindex
        aaindex.init(path = self.aa_index_path, index = '2')
        mcIndex = aaindex.get(indexType)  
        
        # Calculate Grantham scores for mutation
        return mcIndex.get( self.uniprot_aa_change[0], self.uniprot_aa_change[2] )

    def __calculateAminoSizeChg(self):
        from resources import aaindex2b as aaindex
        # Volume
        aaindex.init(path = self.aa_index_path, index = '1')
        volIndex = aaindex.get('PONJ960101')
        #print ( 'AminoVolChg params: ', volIndex, self.aaChangeVarMap[0], self.aaChangeVarMap[2] )
        aminoFromSize = volIndex.get(self.uniprot_aa_change[0])
        aminoToSize = volIndex.get(self.uniprot_aa_change[2])
        return aminoToSize - aminoFromSize
