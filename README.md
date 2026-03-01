# DIVERS
Genome-wide detection of <ins>D</ins>eep <ins>I</ins>ntronic <ins>V</ins>ariants with <ins>E</ins>ffect on <ins>R</ins>ecursive <ins>S</ins>plicing

## Introduction
- Human transcripts with long introns can undergo recursive splicing (RS) - a stepwise splicing mechanism that seamlessly removes long intronic sequences by segmenting them into smaller fragments. 

- Deep intronic variants that disrupt RS could alter splicing, leading to intronic retentions with/without frameshifts and premature stops, thus resulting in deleterious and potentially pathogenic gene products, contributing to human diseases and traits.

- DIVERS is a genome-wide computational method to systematically detect and interpret deep intronic variants that may impact RS. DIVERS identifies five types of RS-variants: "RS_AGGT_loss", "RS_BP_loss", "RS_3SS_gain", "RS_5SS_gain", and "RS_CRYP".

- This standalone version can be easily implemented by a one-line command, for large/small-scale analyses. We also provide the [webserver version](https://hgidsoft.rockefeller.edu/DIVERS) with a user-friendly interface, for small-scale analyses.

## News
- Apr 2025: DIVERS webserver & GitHub were launched.  
- Sep 2024: DIVERS was applied to multiple databases with deep analyses.
- Apr 2024: DIVERS prototype was completed.
- Jul 2023: DIVERS idea was conceived.

## Usage 
Current version: version-1

### Dependency
The code is written in [python3](https://www.python.org/downloads/), and requires [bedtools](https://bedtools.readthedocs.io/en/latest/) installed.

### Reference datasets
Due to the file size limit, please download them from [DIVERS_Detection.bed](http://hgidsoft.rockefeller.edu/DIVERS/standalone.html) and put them into your DIVERS folder.

### File Format
**Input:** Variants in VCF format (GRCh38/hg38), with 5 mandatory tab-delimited fields (CHROM, POS, ID, REF, ALT).
  - Check the example: sample_variants_DIVERS.vcf
  - DIVERS will append all other annotation fields from your input to the end of its output

**Output:** DIVERS-detected variants will be output in CSV format, with the following annotations.
  - **SAMPLE**: sample name (only for DIVERS-batch.py)
  - **CHROM, POS, ID, REF, ALT**: (exactly the same as input)
  - **STRAND**: +/-
  - **GENE**: gene symbol
  - **TRANSCRIPT**: transcript ID (e.g. ENST123456789)
  - **IVS#**: rank of the intron in this gene (e.g., IVS1, IVS2, IVS3)
  - **RS#**: rank of the RS in this intron (e.g., RS1, RS2, RS3)
  - **RS_CONSEQ**: predicted consequences (RS_AGGT_loss, RS_BP_loss, RS_3SS_gain, RS_5SS_gain, RS_CRYP), where Xnt suggests the potential size of intronic insertion
  - **RS_SCORE**: weighted confidence score of the RS-site (1-9, the higher the better)
  - **IVS_SIZE**: intron size
  - **BP_POS**: BP position of the RS-site
  - **PPT**: pyrimidine content (%) in the [BP+2, RS-5] region
  - **RS_START**: start position of the essential RS-site AGGT
  - **RS_END**: end position of the essential RS-site AGGT
  - **RNASEQ**: support from RNA-seq splice junction reads? (Y/N)
  - **CLIP**: support from U2AF eCLIP peaks? (Y/N)
  - **RNALM**: support from RNA language model prediction? (Y/N)
  - **PHYLOP**: all four essential RS-site nucleotides show positive phyloP conservation scores? (Y/N)
  - **PHASTCONS**: RS-site overlaps a conserved element with PhastCons score > 0.5? (Y/N)
  - **RARE**: no human variants with MAF > 0.01% at all four essential RS-site nucleotides? (Y/N)
  - All other annotation fields from the input data (Note: for DIVERS-batch.py, all input VCFs must contain the same annotation fields)

### Command & Parameters (DIVERS_VCF.py)
```
python DIVERS-single.py -i variants.vcf
```

Parameter | Type | Description | Default
----------|------|-------------|--------------
*-i*|file|variants in VCF format, with 5 fields (CHROM, POS, ID, REF, ALT)|N.A.

### Command & Parameters (DIVERS_VCF_batch.py)
```
python DIVERS-batch.py -d foldername/ -s samplelist.txt -o output.csv
```

Parameter | Type | Description | Default
----------|------|-------------|--------------
*-d*|str|directory of VCF files|N.A.
*-s*|file|sample list (without .vcf extension) in the above directory|N.A.
*-o*|str|output CSV filename|N.A.


## Reference
- ***Zhang P. et al.*** Genome-wide detection of human deep intronic variants that affect recursive splicing. 2026.

## Contact
> **Developer:** Peng Zhang, Ph.D.

> **Email:** pzhang@rockefeller.edu

> **Laboratory:** Laboratory of Human Genetics of Infectious Diseases

> **Institution:** The Rockefeller University, New York, NY, USA
