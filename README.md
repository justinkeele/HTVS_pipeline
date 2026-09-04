# High-Throughput Virtual Screening (HTVS) Pipeline
Created by: Justin Keele (Biology, back-end modules scripts) and London Beglich-Bal (API/CLI Development, C++ optimization), University of Idaho  
Supported by the College of Engineering Professional Skills Academy, the Department of Chemical and Biological Engineering, and funding from the DeVlieg Foundation.

## Overview
This repository contains the foundational architecture and scripts for a modular, automated High-Throughput Virtual Screening (HTVS) pipeline. The pipeline is meant to handle large-scale molecular docking. The scripts have been successfully run on libraries of 870,000 molecules downloaded raw from the ZINC database for a inital high throughput virtual screening campaign against YeF3, using the 7B7D and 2IW3 targets. This process started with the raw SMILES strings downloads, then filtered and processed with OpenBabel/RDKit and finally run with the AutoDock-Vina software. We are currently working on developing tools for GNINA, CmDock/RxDock, and potentially others. Dock6 is included but is very fragile since the program is no longer maintained. 

The ultimate goal of this project is a user-friendly API and GUI for a biologist and/or student who is unfamiliar with bash and python coding languages. While tools like CHARMM-GUI, DockM8 and Schrodinger Suite attempt to do this, they all have significant downsides. CHARMM-GUI is mainly focused on molecular dynamics and requires a lot of command-line usage, making it both non-applicable to HTVS and complex for the average user. Webservers like DockM8 are very simplified, preventing rigorous parameter tweaking. Schrodinger is an incredibly expensive Enterprise-grade software suite. It does not explain the biology and parameters behind it, making it a black-box for anyone learning HTVS. It is also, as noted prohibitively expensive. The goal of this project is to guide a user through HTVS, allowing default parameters that may be tweaked as a user gains experience, teaching the parameters and programs along the way. 

**Status: Alpha / Active Development**  
*Note: The scripts in this repository are currently set up with hardcoded paths specific to a local environment and the U of I HPC cluster. Users cloning or forking this repository would need to manually configure the directory paths in the shell scripts to match their local environment and the directory structure. This pipeline is currently being transitioned from use with my private laptop and HPC account to be usable for other people; it has many hardcoded paths and likely won't work on another device or HPC account. There is a secondary branch I am working on for analysis and benchmarking.

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

- Automated ZINC database aggregation and filtering directly within the pipeline.
- Additional options in RDKit to generate tautomers and isomers, pH range. 
- Additional docking/fingerprinting modules (CmDock, ArtiDock, ProLIF). REINVENT4 for generating new hits from a parent molecule, with a target molecule design. 
- Consensus scoring models to cross-reference hits across multiple scoring platforms.
- Eliminate hardcoded paths and allow user-defined hyperparameter tuning (e.g., variable exhaustiveness) via command-line arguments.
- An API for users to easily select and define runs

 In the future:
- Partially automated receptor preparation
- Test across multiple devices
- A complete install.sh script for easy setup
- A web-based GUI
- Dynamic HPC-scheduler setup

# License & Credit
This code is open for academic and research use and licensed under the GNU General Public License v3.0 (GPLv3).

This guarantees the code remains open-source. Under the terms of the GPLv3, any individual or laboratory that forks, modifies, or integrates any part of this architecture into their own workflows or publications must:
1) Make their resulting project entirely open-source.
2) Explicitly cite this repository and its original authors (Justin Keele & London Beglich-Bal).

Failure to adhere to these terms constitutes a violation of the open-source license.

# Acknowledgements​ and Citations

Computational resources were provided in part by Research Computing and Data Services (RCDS). Justin Keele also thanks the Drug Discovery Working Group of the
Institute for Modeling Collaboration and Innovation (IMCI) for providing initial additional support in learning about HTVS architecture and designing a pipeline to target YeF3. 

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
