
# ont-rna-modifications

Master's Thesis Project — Friedrich Schiller University Jena

This project investigates how RNA modifications, specifically **5-methylcytosine (m5C)** and **pseudouridine (Ψ or psU)**, influence the raw signal characteristics of Oxford Nanopore Technologies (ONT) sequencing, particularly as modified 5-mers pass through the R9 pore model.


## Overview

- The raw RNA sequencing data was generated using the ONT R9 pore chemistry.
- Separate datasets were used for each modification type (m5C and Ψ) to ensure accurate comparative analysis.
- Various statistical and signal-processing techniques were applied to characterize how modifications alter current signal profiles.


## Project Structure

- **`baches/`**  
  Contains shell scripts for running various tools, including:
  - Converting POD5 files to SLOW5 format  
  - Running the ONT basecaller  
  - Aligning ONT signal events to the reference genome using `f5c eventalign`

- **`distribution_analysis/`**  
  Includes Jupyter notebooks for each RNA modification (e.g., m5C, psU). These notebooks call functions from utility folders to perform statistical analysis and generate plots.

- **`distribution_analysis_utils/`** and **`mapping_utils/`**  
  Contain Python modules with utility functions for:
  - Statistical analysis  
  - Machine learning workflows  
  - Visualization of signal data

- **`mapping_utils/`**  
  Also includes scripts for:
  - Parsing basecalling outputs  
  - Mapping signal-level events to the reference genome  
  - Extracting Gaussian features from signal segments



## Data Availability

Due to privacy constraints and the large volume of raw data, **sequencing reads and processed datasets are not included** in this repository.

If you are interested in data access or collaboration, please contact me directly.


## Contact

Feel free to reach out with any questions or collaboration ideas:

- 📧 mohammad.noori.vareno@uni-jena.de  
- 📧 hadivareno@gmail.com
