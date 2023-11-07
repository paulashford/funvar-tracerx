## funvar-tracerx 
31/10/2023

**This repository contains scripts relevant to submission of FunVar-tracerx manuscript and includes code-generated figures, FunVar scoring algorithm and stats calculations.**

| Syntax 				| Description | Description |
| ----------- 			| ----------- | ----------- |
| `./script` 			| 						| Collated scripts in submission form (ie remove hard-coded paths and any non-relevant code) |
|  						| `/fie_scoring` 		| Run the FunVar scoring using a flat file export[1]  |
|  						| `/fie_scoring/packages` | Classes for mutations, FIEs and Grantham scoring (aaindex)  |
|  						| `/fie_scoring/working`| Datasets generated from FIE scoring |
|						| `/diversity_analysis`	| Hill-Shannon diversity calculations & plots |
|	`./data`			|						| Datasets for figures and calculations incl db exports so can run locally |
|	`./plot`			|						| Plots (except diversity) and plot data & working directory |


[1] Scoring can be run interactively from: script/fie_scoring/nfe_main.py
**Input data (mutations, clusters, sites etc):**
`./data/nfe_score_database_py_funvar_archive_nfe_pfh_mutfam_090_exported_no_mut_id.tsv`
**Generated from Oracle db by:**
`./data/nfe_score_database_py_export.sql`


Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg




