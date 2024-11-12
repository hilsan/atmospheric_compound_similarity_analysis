# Similarity analysis for atmospheric compounds 

Published version at  https://doi.org/10.5281/zenodo.14007731 .  


A compound similarity analysis of atmospheric compounds conducted on atmospheric datasets from Wang et al. (https://doi.org/10.5194/acp-17-7529-2017), , Isaacman-VanVertz and Aumont (https://doi.org/10.5194/acp-21-6541-2021), Tabor et al. (https://doi.org/10.1039/C9TA03219C ) datasets, and non atmospheric datasets from https://massbank.eu/MassBank/, https://mona.fiehnlab.ucdavis.edu/, Ramakrishnan et al. (10.1038/sdata.2014.22) and Krabrov et al. (DOI: 10.1039/d2cp03966d). The analysis is conducted in three parts. 1) analysis of presence of atmospherically relevant fragments, molecular size and atomic ratios. 2) Comparison of molecular fingerprints (MACCS and RDKit topological fingerprint) through TSNE clustering. 3) Comparison of molecular fingerprint by calculating pairwise similarities between all molecules in two datasets. This project requires a version of aprl_ssp (https://github.com/stakahama/aprl-ssp/tree/master/SMARTSpatterns) and my other project MolAnalysis (https://gitlab.com/hilsan/MolAnalysis.git, https://doi.org/10.5281/zenodo.14007835). Note that there are hard coded paths in some scripts that will need updating depending on your setup and file system.

# Analysis and data processing:
# src:
 1. data:
        get_massbankeu_smiles.sh -bash script to preprocess massbankeu msp file (used in make_dataset)
        get_nabla_dft_smiles.sh - bash script to preprocess nabla_dft msp file (used in make_dataset)
        get_smiles_from_mb_eu_msp.sh - bash script to preprocess mb_eu msp file (used in make_dataset)
        get_smiles_from_mona_msp.sh - bash script to preprocess mona msp file (used in make_dataset)
        make_dataset.py - read raw compound lists and filters them. Remove error producing SMILES and duplicates or sulphur 
 2. features:
        build_descriptors.py - create  features, fingerprints, and fungroups scripts to 
        build_features.py - functions to calculate atomic rations, elemental composition and molecular size.
        build_fingerprints.py - create molecular fingerprints (MACCS and topological fingerprints)
        build_FunGroups.py - calculate functional group statistics. 
        load_datasets.py
 3. fingerprint_comparison:
        calc_tsne.py  - perform a t-SNE analysis of two datasets 
        compare_fps.py - perform a tanimoto similarity analysis of two dataset files. 
        sum_bins.py - post analysis of the tanimoto similarity statistics in a directory.


# Project run and visualization 
# notebook:
 'notebook/run_project.ipynb' - Outlines the sequence of commands to project analysis. This assumes that there exists a data/raw folder with subfolders: gecko(txt with a list of SMILES), mb_eu (msp file), mona/all/ (msp file), nablaDFT (csv file), qm9/xyz/ (xyz file for each molecule), quinones (csv file), wang (with Supporting information CSV.csv). For the manuscript, all datasets were collected 8th of January, 2024.

 'notebook/for_manuscript.ipynb' - Includes code to create figures for manuscript from.
 'notebook/psat.ipynb' - Code to create plot comparing GeckoQ vapor pressures to those from Handbook of Chemistry and Physics

 data and reports empty folders. First for putting the raw data (see above for instructions) and second where to write figures and analysis.
 



Developed in Python ver. 3.10.8. See requirements.txt for relevant library versions.
