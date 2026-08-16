# High-Throughput Virtual Screening (HTVS) Pipeline

Created by: Justin Keele and created under Dr. Jagdish Patel, University of Idaho  
Supported by the College of Engineering Professional Skills Academy, the Department of Chemical and Biological Engineering, and funding from the DeVlieg Foundation.

## Overview
This repository contains the foundational architecture and scripts for a modular, automated High-Throughput Virtual Screening (HTVS) pipeline. The pipeline is meant to handle large-scale molecular docking. The scripts have been successfully run on libraries of 870,000 molecules downloaded raw from the ZINC database for our lab's high throughput virtual screening campaign against YeF3, using the 7B7D and 2IW3 targets. This process started with the raw SMILES strings downloads, then filtered and processed with OpenBabel/RDKit and finally run with the AutoDock-Vina software. We are currently working on developing tools for GNINA, CmDock/RxDock, and potentially others. Dock6 is included but is very fragile since the program is no longer maintained. 

**Status: Alpha / Active Development**  
*Note: The scripts in this repository are currently set up with hardcoded paths specific to a local environment and the U of I HPC cluster. Users cloning or forking this repository would need to manually configure the directory paths in the shell scripts to match their local environment and the directory structure. This pipeline is currently intended for my laptop and HPC account to be used on internal Patel lab development and benchmarking; it will likely have many breaking changes*

### Current Directory Structure
```
HTVS_pipeline/
├── docs/                  # Documentation and generated benchmarking graphs
├── environments/          # Conda env required files
├── libraries/             # Aggregated SMILES subsets (e.g., ZINC, Enamine)
├── modules/               # Independent execution modules for docking/processing
│   ├── analysis/          # Log parsing and metric extraction scripts
│   ├── dock6/             # DOCK6 execution scripts
│   ├── gnina/             # GNINA execution scripts
│   ├── rdkit_prep/        # 3D coordinate generation and filtering scripts
│   ├── rxdock/            # CmDock/RxDock execution scripts
│   ├── vina/              # AutoDock Vina execution scripts
│   └── ZINC_download/     # ZINC database retrieval and combination scripts
├── results/               
│   └── yef3/
│       └── gnina/         # Contains parsed results, SDFs, and Meeko temp files for specific runs
└── targets/               # Target protein data
    └── yef3/
        ├── configs/       # Vina/GNINA grid box parameters
        ├── ligand_batches/# ligand batches ready for processing/docking
        └── receptors/     # Cleaned receptor files
```

# Future Development
This pipeline is currently being expanded to include:

- Automated ZINC database aggregation and filtering directly within the pipeline. Additional options to generate tautomers and isomers, pH range. 
- Additional docking modules (CmDock, ArtiDock, ProLIF, REINVENT4).
- A possible consensus scoring model to cross-reference hits across multiple scoring platforms.
- Dynamic configuration files to eliminate hardcoded paths and allow for user-defined hyperparameter tuning (e.g., variable exhaustiveness) via command-line arguments.
- Partially automated receptor preparation
- An API for users to easily select and define runs

# License & Credit
This code is open for academic and research use. If you fork, modify, or integrate this architecture into your own workflows or machine learning models, please cite this repository and its developer.

# Acknowledgements​ and Citations

Computational resources were provided in part by Research Computing and Data Services (RCDS). We also thank the Drug Discovery Working Group of the
Institute for Modeling Collaboration and Innovation (IMCI) for providing additional support. 

Tools and servers included are: ZINC_database, RDKit, OpenBabel, dock6 (legacy), gnina, Rxdock, CmDock, AutoDock-Vina, Conda, and various Python libraries. Credit goes to each of the creators of these tools. 

Notable papers used in this work include:

Condon DE, Schroeder BK, Rowley PA, Ytreberg FM (2025) Discovery of novel
targets for important human and plant fungal pathogens via an automated
computational pipeline HitList. PLoS One 20(6): e0323991.​

Andersen, Christian B F et al. “Structure of eEF3 and the mechanism of transfer RNA
release from the E-site.” Nature vol. 443,7112 (2006): 663-8. ​

Ranjan, N., Pochopien, A.A., Chih‐Chien Wu, C. et al. Yeast translation elongation
factor eEF3 promotes late stages of tRNA translocation. EMBO J 40,
EMBJ2020106449 (2021). ​

Dever, Thomas E, and Rachel Green. “The elongation, termination, and recycling
phases of translation in eukaryotes.” Cold Spring Harbor perspectives in biology vol.
4,7 a013706. 1 Jul. 2012.​
