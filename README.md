# DIVERS
Genome-wide detection of <ins>D</ins>eep <ins>I</ins>ntronic <ins>V</ins>ariants with <ins>E</ins>ffect on <ins>R</ins>ecursive <ins>S</ins>plicing

## Introduction
- Human transcripts of long introns can be removed from pre-mRNA by recursive splicing (RS), a multi-step splicing mechanism that seamlessly converts a long intronic sequence into smaller segments for removal. 

- Deep intronic variants that disrupt RS could alter splicing, leading to intronic retentions with/without frameshifts and premature stops, thus resulting in deleterious and pathogenic gene products, underlying human diseases and traits. 

- DIVERS is a genome-wide computational method for detecting deep intronic variants that may impact RS, systematically and efficiently. DIVERS currently identifies five types of variants: "RS-AGGT", "RS-BP/BP2", "RS-AGAIN", "RS-DW5SS", and "CRYPRS".

- This standalone version can be easily implemented by a one-line command, for large-scale analyses. We also provide the [webserver version](https://hgidsoft.rockefeller.edu/DIVERS) with a user-friendly interface, for small-scale analyses.

## News
- Apr 2025: DIVERS webserver & GitHub were launched.  
- Sep 2024: DIVERS was applied to multiple databases with deep analyses.
- Apr 2024: DIVERS prototype was completed.
- Jul 2023: DIVERS idea was conceived.

## Usage 
Current version: version-1

### Dependency
The code is written in [python3](https://www.python.org/downloads/), and requires [bedtools](https://bedtools.readthedocs.io/en/latest/) installed.

### Reference dataset
Due to the file size limit, please download the [DIVERS_Detection.bed](http://hgidsoft.rockefeller.edu/DIVERS/standalone.html) and put it into your DIVERS folder.

### File Format
**Input:** Variants in VCF format (GRCh38/hg38), with 5 mandatory tab-delimited fields (CHROM, POS, ID, REF, ALT).
  - Check the example: sample_variants_DIVERS.vcf
  - DIVERS will append all other annotation fields from your input to the end of its output

**Output:** DIVERS-detected variants will be output in CSV format, with the following annotations.
  - SAMPLE: sample name (only for DIVERS_VCF_batch.py)
  - CHROM, POS, ID, REF, ALT: (exactly the same as input)
  - STRAND: the strand +/- where the variant found affecting RS
  - GENE: gene symbol
  - TRANSCRIPT: transcript ID (e.g. ENST123456789)
  - IVS#: the ranking number of the intron in the gene (e.g. IVS1, IVS2, IVS3)
  - IVS_SIZE: the size of the intron
  - RS#: the ranking number of the RS in this intron (e.g. RS1, RS2, RS3)
  - RS_CONSEQ: the predicted consequences (RS-AGGT, RS-BP/BP2, RS-AGAIN, RS-DW5SS-xnt, CRYPRS-DW5SS/UP3SS-xnt), where xnt suggesting the size between the paired cryptic splice sites
  - RS_SCORE: the weighted confidence score of RS-site (1-5, the higher the better)
  - RS_POS: the first position of the essential RS-site AGGT
  - BP_POS: the BP position of the RS-site
  - PPT: the pyrimidine content in the PPT region
  - CLIP: if the RS-site is supported by eCLIP-U2AF data (Y/N)
  - RARE: if the RS-site is absent of common human variants (Y/N)
  - CONSERV: if the RS-site is positively scored by conservation score phyloP (numerical)
  - RNALM: if the RS-site is predicted by RNA language model (0-1)
  - All other annotation fields from the input data: Note: if you are using DIVERS_VCF_batch.py, all VCF files should have the same annotation fields.

### Command & Parameters (DIVERS_VCF.py)
```
python DIVERS_VCFP_one.py -i variants.vcf
```

Parameter | Type | Description | Default
----------|------|-------------|--------------
*-i*|file|variants in VCF format, with 5 fields (CHROM, POS, ID, REF, ALT)|N.A.

### Command & Parameters (DIVERS_VCF_batch.py)
```
python DIVERS_VCF_batch.py -d ./foldername/ -s samplelist.txt -o output.csv
```

Parameter | Type | Description | Default
----------|------|-------------|--------------
*-d*|str|directory of VCF files|N.A.
*-s*|file|sample list (without .vcf extension) in the above directory|N.A.
*-o*|str|output CSV filename|N.A.


## Reference
- ***Zhang P. et al.*** Genome-wide detection of human deep intronic variants that affect recursive splicing. 2025.

## Contact
> **Developer:** Peng Zhang, Ph.D.

> **Email:** pzhang@rockefeller.edu

> **Laboratory:** Laboratory of Human Genetics of Infectious Diseases

> **Institution:** The Rockefeller University, New York, NY, USA
