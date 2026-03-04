# IR multiGWAS Pipeline

Welcome! This repository contains the full analytical pipeline used in our publication:  
**_Genome-Wide Discovery Reveals Adipose-Specific and Systemic Regulators of Insulin Resistance_**

This project implements a comprehensive **multi-trait GWAS framework** to identify genetic regulators of insulin resistance (IR), leveraging complementary statistical approaches for robust discovery.

Full GWAS summary statistics produced with G-SEM can be found in the following zenodo repository:
Independent IR loci and their G-SEM and CPASSOC association can be found in Supplementary Table 2 of the manuscript

**NOTE** this repository is still under construction. Our publication will soon be sent for reviewing 

---

## What This Pipeline Does

We compute multi-trait IR associations using two complementary methods:

- **CPASSOC** – cross-phenotype association testing  
- **G-SEM (Genomic Structural Equation Modeling)** – multivariate genetic modeling framework  

Together, these approaches provide a powerful strategy to detect both adipose-specific and systemic regulators of IR.

---

## Repository Structure & Workflow

The code is organized as a **sequential analysis pipeline**, designed for full reproducibility.

### 🔹 Step 1 — Data Preparation  
Download the data provided in **Supplementary Table 1** and place it in the `raw_data/` folder.  
This is required to replicate the analyses.

### 🔹 Step 2 — GWAS Curation  
The first analytical step harmonizes and curates all required GWAS summary statistics.

### 🔹 Step 3 — run the analyses!!

### 🔹 Outputs  
All results are written to the `output/` directory.

We provide the full output folder structure in this repository so users can:
- Reproduce all intermediate steps  
- Regenerate intermediary files  
- Fully replicate the published results  

---

## Software Requirements

Most analyses were performed using: R/4.3.1
Other softwares that have been utilized and should be installed include:

```
CPASSOC - version
LDSC - version
DEPICT - version
GoShifter - version
CHEERS - version
HOMER - version
```


