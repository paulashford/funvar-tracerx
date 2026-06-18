# FunVar Functional Impact Events (FIEs) in cancer

## FunVar-TRACERx 

**"Gene duplication is associated with gene diversification and potential neofunctionalisation in lung cancer evolution" (Ashford *et al* 2025 )[^1]**

Repository includes FunVar-FIE scoring algorithm, FIE-gene/FIE-FunFam diversity scoring,
Cancer Effect Size (CES) calculations, benchmarking, and associated datasets.

In addition, FunVar-FIE scoring relies on previously published methods:

- [MutClust](https://github.com/paulashford/mutclust) - Calculation of significant cancer mutation clusters on protein structures by permutation testing

- [MutFam](https://github.com/paulashford/mutfam) - Mutationally enriched CATH Functional Families (used as a component of FIE-scores)

Benchmarking and resources are described in [/benchmark/README.md](./script/benchmark/README.md).

*Note: some large files (>50Mb) are zipped in the data directory and require extracting prior to running CES.*

## Paths and descriptions

| path  				| 							| Description |
| ----------- 			| ----------- 				| ----------- |
| `./script` 			| 							| Collated scripts including summary calculations and statistics|
|  						| `/fie_scoring` 			| Run FunVar-FIE scoring on a pre-processed mutation dataset|
|  						| `/fie_scoring/packages` 	| Classes for mutations, FIEs and FunVar scoring|
|  						| `/fie_scoring/resources` 	| aaindex package|
|  						| `/fie_scoring/working`	| FIE score outputs|
|						| `/diversity_analysis`		| Hill-Shannon diversity calculations & plots |
|						| `/ces`					| Cancer Effect Size calculations with cancereffectsizeR for FIEs grouped by FunFam and alignment residue number|
|						| `/benchmark`				| Binary classification benchmarking (see: [/benchmark/README.md](./script/benchmark/README.md)) |
| `./data`				|							| General datasets for figures, calculations, and FIE scoring|
|						| `/cath`					| CATH v4.2 FunFams and UniProt mapping |
|						| `/benchmark/minimal`		| Datasets for benchmarking and comparison with other structure-based driver prediction algorithms|
| `./plot`				|							| Plots for FunVar scores by CGC gene category, long-tail of FIEs per gene, and CES outputs |

## Running FunVar-FIE scoring

```shell
# Create a virtual environment and install libraries
# cd <git-cloned-dir>
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

# Run FIE-scoring with dataset of mutations annotated 
# with clusters and functional sites 
python script/fie_scoring/nfe_main.py
# Output written to script/fie_scoring/working/
```

[![CC BY 4.0][cc-by-image]][cc-by]
This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-shield]][cc-by]

AAindex data is provided by Kawashima *et al* [^2].
aaindex.py licensed under a [BSD-2-Clause license](BSD).
Further info: https://pymolwiki.org/AAindex

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

[^1]: Ashford P, Frankell AM, Piszka Z, Pang CSM, Abbasian M, Bakir MA, Jamal-Hanjani M, McGranahan N, Swanton C, Orengo CA. Gene duplication is associated with gene diversification and potential neofunctionalization in lung cancer evolution. Genome Res. 2026 Mar 2;36(3):561-577. doi: 10.1101/gr.278663.123. PMID: 41714147; PMCID: PMC12951968. 
https://doi.org/10.1101/gr.278663.123
[^2]: Kawashima, S. and Kanehisa, M.; AAindex: amino acid index database. Nucleic Acids Res. 28, 374 (2000). [PMID:10592278]
https://www.genome.jp/aaindex/aaindex_help.html
