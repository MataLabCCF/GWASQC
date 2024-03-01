# GWASQC

Genotyping Quality Control focused on admixed populations. This pipeline is part of a LARGE-PD variant calling and QC 
procedure, specifically designed for LARGE-PD Phase 2 data.

## Options

```
Genotyping QC to GWAS

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --inputFile INPUTFILE
                        Input file without suffix
  -o OUTPUTNAME, --outputName OUTPUTNAME
                        Output name
  -O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
                        Output folder
  -I INFO, --info INFO  File with covar information. The pipeline requires at least three columns (ID, Sex, Phenotype)

Optional arguments:
  -e, --erase           Set to exclude all temporary files
  -p PHENONAME, --phenoName PHENONAME
                        Name of the phenotype column in the info file (default: PHENO)
  -s STEPS, --steps STEPS
                        File with the order of the steps (one step per line)
  -P POPFILE, --popFile POPFILE
                        File with two columns: Individual ID and Population ID.This approach is important
                        forrelationship control that should be done per population
  -S SAVEPERPOP, --savePerPop SAVEPERPOP
                        Save the QCed data with all pops and per population

Programs arguments:
  --plink2 PLINK2       Path to Plink 2 (default: plink2)
  --plink1 PLINK1       Path to Plink 1.9 (default: plink)
  --NAToRA NATORA       Path to NAToRA (default: NAToRA.py)
  --python PYTHON       Python > 3 with networkX instaled (default: python)
  --king KING           Path to king program (default: king)

```

## Input files

- **Input**: Data with genotyped information. It can be VCF(vcf or vcf.gz), Plink Bfile (default to plink 1.90) or 
Plink Pfile (default to plink 2)
- **Steps**: List of QC steps to be performed and the parameters. If you did provide any file, the program will run the 
default QC with default parameters. See QC Steps Section.
- **Pop file**: File with two columns. The first column is the individual ID and the second is the population ID. This 
file is used by relationship control and to the flag -S (to be implemented)
- **Info**: Provide covariate information in a file with at least 3 columns and a header:
  - **ID:** Sample ID as in the input data. Possible header names: ID, ind, sample_external_id.
  - **SEX:** Sex information for the sample. Possible header names: Sex, gender.
    - We accept four different values for SEX:
      - *Males:*
        - `0` (0/1 codification)
        - `1` (1/2 codification, default to Plink)
        - `M` or `Male`
      - *Females:*
        - `1` (0/1 codification)
        - `2` (1/2 codification, default to Plink)
        - `F` or `Female`

  - **STATUS:** Status information for the sample. Possible header names: Pheno, status, affection_status.
    - *Case:*
      - `0` (0/1 codification)
      - `1` (1/2 codification, default to Plink)
      - `case`, `pd case`, `affected`

    - *Control:*
      - `0` (0/1 codification)
      - `1` (1/2 codification, default to Plink)
      - `control`, `pd control`, `unaffected`, `non-affected`, `non affected`, `not affected`

#### Example for info file:

```csv
ID,SEX,STATUS
Sample1,M,Case
Sample2,F,Control
Sample3,Male,Control
Sample4,Female,Affected
```

Our code is not case-sensitive, and it accommodates different conventions for header names. 
If your file has headers with names like "id" instead of "ID" our pipeline will handle it seamlessly.

## QC Steps

In this section, we provide a brief overview of the quality control (QC) steps performed by the pipeline. 
The steps are presented in the default order.

### Removal of Samples without Sex and Status (Mandatory)

- **Description:**
  - We prioritize data integrity by requiring sex and status information for all samples. Samples lacking this essential information are removed from the analysis.
  
- **Implementation:**
  - We remove all samples with missing data in the Sex and Status columns of the covariate info file.

### Sex-Check (Optional)

- **Description:**
  - This step involves checking the consistency between the declared sex (provided in the info file) and the inferred 
  sex based on the sexual chromosome information using Plink1.90.
  
- **Implementation:**
  - We use Plink1.90 with specific thresholds to determine sex discrepancies.
    - Females: F < 0.5
    - Males: F > 0.8

  - **Note:**
    - By default, Plink uses F < 0.2 to female. In admixed populations, females might commonly have F statistics higher 
    than 0.2, what can cause an unnecessary sample loss.

### Redundant variants (Optional)
- **Description:**
  - The genotyping platform used in LARGE-PD Phase 2 sometimes reports genotypes for variants on the "+" strand
  and other times on "-" strand. This can lead to issues, especially for variants involving A/T and C/G, 
  potentially introducing bias into the analysis.
- **Implementation:**
  - We remove all variants that are A/T or C/G.

### Sample missing data (Optional)
- **Description:**
  - High missing data in samples may indicate issues such as problems with DNA extraction, degraded DNA, or 
  low DNA concentration, which can impact the reliability of genotyping data.

- **Implementation:**
  - To ensure data quality, we exclude samples with missing data exceeding a specified cutoff.
    - Default cutoff: 0.05 (5%)

### Variant missing data (Optional)
- **Implementation:**
  - To ensure data quality, we exclude variants with missing data exceeding a specified cutoff.
    - Default cutoff: 0.05 (5%)

### Duplicate variants (Optional)
- **Description:**
  - Certain analyses may require the absence of duplicate variants to avoid redundancy and potential biases in 
  downstream results.

- **Implementation:**
  - We identify and compare duplicated variants, calculating the missing data for each duplicate. 
  To maintain data quality, we keep the variant with the least missing data.

### Heterozigosity (Optional)
- **Description:**
  - Deviations in heterozygosity can serve as indicators of sample contamination or inbreeding within the dataset.

- **Implementation:**
  - We calculate the sample heterozygosity and provide an optional step to remove samples that deviate from a specified cutoff.
    - Default cutoff: 3 (sd).

- **Note:**
  - While some studies remove samples based on Plink F statistics, it's important to note that this approach in admixed 
  populations can result in unnecessary sample loss.

### HWE in Controls (Optional)
- **Description:**
  - Deviations from Hardy Weinberg Equilibrium (HWE) can indicate issues in the genotyping process. 
  It is a common practice to set different HWE cutoffs for cases and controls.

- **Implementation:**
  - Variants with p-values below the specified cutoff are removed.
    - Default cutoff: 1e-6.

### HWE in Cases (Optional)
- **Description:**
  - Deviations from Hardy Weinberg Equilibrium (HWE) can indicate issues in the genotyping process. 
  It is a common practice to set different HWE cutoffs for cases and controls.

- **Implementation:**
  - Variants with p-values below the specified cutoff are removed.
    - Default cutoff: 1e-10.

### Relationship Control (Optional)
- **Description:**
  - The presence of related samples can cause bias in genetic analysis. 

- **Implementation:**
  - We use NAToRA software to remove related samples. We consider related samples which kinship coefficient > cutoff
    - Default cutoff: 0.0884 (second degree).
- **Note:**
  - If your data have samples from different population, the best standard is to calculate the relationship by 
  population


### Example file to --steps following the default order and cutoffs:

```tsv
sex-check
ATCG
mind	0.05
geno	0.05
duplicate
heterozygosity	3
HWE control	1e-6
HWE case	1e-10
relationship	0.0884
```

## Programs used

In this pipeline we use:

- Plink1
- Plink2
- Python3
  - Libraries required:
    - NetworkX *
    - Numpy
- NAToRA ([Leal et al. 2022](https://doi.org/10.1016/j.csbj.2022.04.009))*

\* Used on relationship control

## Common Problems

-  **PLINK2 errors**:
   - Error: Unrecognized flag ('--het').
     - PLINK2 is a tool in development, so we strongly recommend download the latest version and use it.

- **NAToRA errors**:
  - ImportError: No module named networkx
    - Install networkX in your python.
      - *pip install networkx* or *python -m pip install networkx*
    - If the problem persists
      - Some computer environments have different python installed. To fix this
        - type *which python* and copy the path and add on the --python in your command line
      


## Acknowledgements

This work is supported by NIH Grant R01 1R01NS112499-01A1, MJFF Grant ID: 18298, ASAP-GP2 and Parkinson's Foundation

We also want to thank all users that report problems/bugs.

 ### Contact 

 Developer: Thiago Peixoto Leal. PhD (PEIXOTT@ccf.org or thpeixotol@hotmail.com)


