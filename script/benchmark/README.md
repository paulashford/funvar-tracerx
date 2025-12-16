# FunVar-FIE protocol benchmarking

README.md
P. Ashford
November 2025

FunVar protocol benchmarking of pancancer FIEs against predicted driver mutations obtained from other 3D cancer driver prediction tools.
FIEs and other driver predictions are all referenced against independent annotations:
*Actual Positives*: Clincally significant COSMIC Cancer Mutation Census (CMC)[^1]
*Actual Negatives*: ClinVar[^2] (benign variants) and dbSNP[^3] (benign variants).

Minimal versions of all datasets needed to run benchmarking are included. Please see associated references for full source data and descriptions.

---

## Run benchmark

1. **Clone repository to local directory**
```bash
cd path/to/repos
git clone git@github.com:paulashford/funvar-tracerx.git
cd funvar-tracerx
```

2. **Install pre-requisite R packages (if necessary)**
[tidyverse](https://tidyverse.tidyverse.org/)
[yardstick](https://yardstick.tidymodels.org/)

3. **Run benchmarks**

```bash
cd path/to/repos/funvar-tracerx
Rscript script/benchmark/run_benchmarks.R
```

## Reference datasets

**COSMIC Cancer Mutation Census (CMC)**
Source download v100 on 30/05/2024
https://cancer.sanger.ac.uk/cosmic/download/cancer-mutation-census/v100/alldata-cmc
Filtered for `MUTATION_SIGNIFICANCE_TIER` 1, 2, or 3 and `Mutation_Description_AA` "Substitution - Missense".

**ClinVar**
Source download on 30/05/2024
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
Filtered for benign SNVs where `ReviewStatus` "reviewed by expert panel".

**dbSNP**
Source download on 25/04/2024
https://www.ncbi.nlm.nih.gov/snp/
Filtered for `Function Class` "missense variant" and `Clinical Significance` "benign" and `Organism` "humans".

**HUGO Gene Nomenclature Committee**
Downloaded 31/07/2024 using custom column query as per [link](https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_pub_ensembl_id&col=md_prot_id&col=md_ensembl_id&col=md_eg_id&col=gd_name_aliases&col=md_refseq_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit)
Note an additional "expanded" form of the table (hugo_core_exp_20240731) was also created from the original by pivoting all alias gene names and previous gene names to long-table format to allow for simple cross-checking against all of these for any given gene.

**CATH FunFams**
Functional Families v4.2[^4] from the CATH database ([CATH-db.info](https://www.cathdb.info))
https://www.cathdb.info/wiki/doku/?id=release_notes#cath-plus_version_42



**HotSpot3D[^5]**, **HotMAPS[^6]**, and **3dHotSpots[^7]** obtained via PanCan Atlas at the [Genomic Data Commons (GDC)]((https://gdc.cancer.gov/about-data/publications/pancan-driver))
The PanCan Atlas provides a compendium driver prediction tools and pre-computed predictions for TCGA cancer types and pancancer as described in Bailey *et al* 2018[^8]. We used significance thresholds for each og these 3D methods as described in the GDC compendium.

[^1]: Tate JG, Bamford S, Jubb HC, Sondka Z, Beare DM, Bindal N, Boutselakis H, Cole CG, Creatore C, Dawson E, et al. 2019. COSMIC: the Catalogue Of Somatic Mutations In Cancer. Nucleic Acids Res 47: D941–D947
[^2]: Landrum MJ, Lee JM, Riley GR, Jang W, Rubinstein WS, Church DM, Maglott DR. 2014. ClinVar: public archive of relationships among sequence variation and human phenotype. Nucleic Acids Res 42: D980-985.
[^3]: Sherry ST, Ward M-H, Kholodov M, Baker J, Phan L, Smigielski EM, Sirotkin K. 2001. dbSNP: the NCBI database of genetic variation. Nucleic Acids Res 29: 308–311.
[^4]: Dawson NL, Orengo CA. 2015a. Diversity in protein domain superfamilies. Curr Opin Genet Dev 35: 40–49. Lee D, Sillitoe I, Dawson NL, Lees JG, Orengo CA. 2015b. Functional classification of CATH superfamilies: a domain-based approach for protein function annotation. Bioinformatics 31: 3460–3467
[^5]: Niu B, Scott AD, Sengupta S, Bailey MH, Batra P, Ning J, Wyczalkowski MA, Liang W-W,Zhang Q, McLellan MD, et al. 2016. Protein-structure-guided discovery of functional mutations across 19 cancer types. Nat Genet 48: 827–837.
[^6]: Tokheim C, Bhattacharya R, Niknafs N, Gygax DM, Kim R, Ryan M, Masica DL, Karchin R. 2016. Exome-scale discovery of hotspot mutation regions in human cancer using 3D protein structure. Cancer Res 76: 3719–3731.
[^7]: Gao J, Chang MT, Johnsen HC, Gao SP, Sylvester BE, Sumer SO, Zhang H, Solit DB, Taylor BS, Schultz N, et al. 2017. 3D clusters of somatic mutations in cancer reveal numerous rare mutations as functional targets. Genome Med 9: 4.
[^8]: Bailey MH, Tokheim C, Porta-Pardo E, Sengupta S, Bertrand D, Weerasinghe A, Colaprico A, Wendl MC, Kim J, Reardon B, et al. 2018. Comprehensive Characterization of Cancer Driver Genes and Mutations. Cell 173: 371-385.e18.