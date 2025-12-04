# Long-Read Haplotype Typer for Forensic Sequencing

## Overview
This tool performs haplotype genotyping from long-read sequencing data (e.g., ONT, PacBio) for a **single microhaplotype locus**. It is designed for forensic genetics research where long-fragment microhaplotypes (like MiniHaps) are used.

![Workflow Diagram](workflow_diagram.png)

*Figure: Schematic workflow of the genotyping logic.*

## Input Files

### 1. SAM File
- **Description**: Alignment file where reads are mapped to **a single reference sequence**. The reference should be a fragment encompassing the entire target haplotype locus (e.g., MH-001).
- **Format**: Standard SAM format.
- **Note**: If you have a BAM file, convert it using `samtools view -h input.bam > output.sam`.

### 2. SNP Positions Excel File (`input_snp_positions.xlsx`)
- **Description**: Defines the **relative physical positions** of SNPs within the haplotypes for one or more loci.
- **Required Format**:

| Locus  | SNP_positions                    |
|--------|----------------------------------|
| MH-001 | 107,158,199,378,487,829,843,871  |
| ...    | ...                              |

- **Note**: The SNP positions are 1-based coordinates relative to the **provided reference sequence** in the SAM file.

## Output Files

### 1. Haplotype Ranking (`output_01_haplotype_ranking.xlsx`)
- **Description**: Provides detailed statistics and a ranked list of all complete haplotypes observed in the data.
- **Format**:

| Site   | Information          | Value                               |
|--------|----------------------|-------------------------------------|
| MH-001 | SNP positions        | [100, 155, 427, 501, 707, 728, 794] |
| MH-001 | Total reads         | 34729                               |
| MH-001 | Covering reads      | 23132                               |
| MH-001 | Coverage rate       | 66.61%                              |
| ...    | ...                  | ...                                 |
| MH-001 | Haplotype pattern   | Count                               |
| MH-001 | GCCTGGA             | 8482 (36.67%)                       |
| MH-001 | TCCGCAA             | 7560 (32.68%)                       |
| ...    | ...                  | ...                                 |

*This file is useful for manual inspection and quality control.*

### 2. Genotype Call (`output_02_genotype_call.xlsx`)
- **Description**: Contains the final genotype call (heterozygous or homozygous) for the individual at the target locus, based on the Allelic Depth Ratio (ADR) calculation.
- **Format**:

| Site   | Genotype     | ADR    | Haplotype1 | Count1 | Haplotype2 | Count2 |
|--------|--------------|--------|------------|--------|------------|--------|
| MH-001 | heterozygous | 0.8913 | GCCTGGA    | 8482   | TCCGCAA    | 7560   |

## Genotyping Logic
1.  **Parse SAM**: Extract read info and CIGAR strings.
2.  **Map SNPs**: For each read, use the CIGAR string to find bases at the SNP positions defined in the Excel file.
3.  **Haplotype Construction**: Construct a haplotype string for each read that covers all target SNPs.
4.  **Rank Haplotypes**: Count and rank all complete haplotypes → Output File 1.
5.  **Genotype Call**:
    - Calculate **ADR** = (Count of 2nd-ranked haplotype) / (Count of top-ranked haplotype).
    - **Heterozygous Call (0.1 ≤ ADR ≤ 1)**: Output top two haplotypes.
    - **Homozygous Call (0 ≤ ADR < 0.1)**: Output top haplotype. Each allele count = (Top count) / 2.
    - → Output File 2. *The ADR threshold is adjustable in the code.*

## Usage
1.  Prepare your SAM file and SNP positions Excel file.
2.  Run the script:
    ```bash
    python haplotype_typer.py your_input.sam your_snp_positions.xlsx
    ```
3.  Find the results in the generated Excel files.

## Notes & Limitations
- **Proof-of-Concept**: This script is designed for a **single sample and single locus**. It serves as a core module for building batch-processing pipelines.
- **SAM Format Required**: Convert BAM to SAM if necessary.
- **Parameter Adjustment**: Users may need to adjust paths and the ADR threshold based on their specific data characteristics.