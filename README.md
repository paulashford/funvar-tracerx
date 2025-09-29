## funvar-tracerx 
28/09/2025

**This repository contains scripts relevant to submission of FunVar-tracerx manuscript and includes code-generated figures, FunVar scoring algorithm and stats calculations.**

**FunVar makes use of mutation clusters from MutClust https://github.com/paulashford/mutclust**

**FIE scores include a scoring term from Mutationally ernriched functional Families (MutFams) https://github.com/paulashford/mutfam**

**Note some large files (>50Mb) are zipped in the data directory and need extracting for certain scripts to work correctly (e.g. for ./script/ces/*.R)**


| path  				| | Description |
| ----------- 			| ----------- | ----------- |
| `./script` 			| 						| Collated scripts in submission form (ie remove hard-coded paths and any non-relevant code) |
|  						| `/fie_scoring` 		| Run the FunVar scoring using a flat file export[1]  |
|  						| `/fie_scoring/packages` | Classes for mutations, FIEs and Grantham scoring (AAindex)[2]  |
|  						| `/fie_scoring/working`| Datasets generated from FIE scoring |
|						| `/diversity_analysis`	| Hill-Shannon diversity calculations & plots |
|						| `/ces`	| Cancer Effect Size calculations for FIEs grouped by FunFam and alignment residue using cancereffectsizeR |
|	`./data`			|						| Datasets for figures and calculations incl db exports so can run locally |
|	`./plot`			|						| Plots (except diversity) and plot data & working directory |


[1] FunVar-FIE scoring can be run interactively from: 
`script/fie_scoring/nfe_main.py`

- **using input data (mutations, clusters, sites etc):** `./data/nfe_score_database_py_funvar_archive_nfe_pfh_mutfam_090_exported_no_mut_id.tsv`

- **input data was generated from Oracle db by:**
`./data/nfe_score_database_py_export.sql`

[![CC BY 4.0][cc-by-shield]][cc-by]

This work (excluding AAindex[2]) is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]
---
[2] AAindex is provided by Kawashima et al 
Kawashima, S. and Kanehisa, M.; AAindex: amino acid index database. Nucleic Acids Res. 28, 374 (2000). [PMID:10592278]
https://www.genome.jp/aaindex/aaindex_help.html


[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg




