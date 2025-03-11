# DIVERS
Genome-wide detection of <ins>D</ins>eep <ins>I</ins>ntronic <ins>V</ins>ariants with <ins>E</ins>ffect on <ins>R</ins>ecursive <ins>S</ins>plicing

## Introduction
- Human transcripts of long introns can be removed from pre-mRNA by recursive splicing (RS), a multi-step splicing mechanism that seamlessly converts a long intronic sequence into smaller segments for precise removal. 

- Deep intronic variants that disrupt RS could in principle alter splicing, leading to intronic retentions with/without frameshifts and premature stops, thus resulting in deleterious and pathogenic gene products, underlying human diseases and traits. 

- DIVERS is a genome-wide computational approach for detecting deep intronic variants that may impact RS, systematically and efficiently. DIVERS currently identifies five types of variants: "RS_AGGT", "RS_BP", "RS_AGAIN", "RS_DW5SS", and "CRYP_RS".

- This standalone version can be easily implemented by a one-line command, for large-scale analysis. We also provide the [webserver version](https://hgidsoft.rockefeller.edu/DIVERS) with user-friendly interface, for a small-scale analysis.

## News
- Mar 2025: DIVERS webserver & GitHub was launched.  
- Apr 2024: DIVERS prototype was completed.
- Feb 2023: DIVERS idea was conceived.

## Usage 
Current version: version-1

### Dependency
The code is written in [python3](https://www.python.org/downloads/), and requires [bedtools](https://bedtools.readthedocs.io/en/latest/) installed.

### Reference dataset
Due to the file size limit, please download the [DIVERS reference dataset](http://hgidsoft.rockefeller.edu/DIVERS/standalone.html) and put it into your DIVERS folder.

### File Format
**Input:** Variants in VCF format (GRCh38/hg38), with 5 mandatory fields (CHROM, POS, ID, REF, ALT) tab-delimited.
  - sample_variants_DIVERS.vcf

**Output:** DIVERS-detected variants will be output in CSV format, with the following annotations.
  - SAMPLE: sample name (only for DIVERS_VCF_batch.py)
  - CHROM, POS, ID, REF, ALT: (exactly the same as input)
  - STRAND: the strand +/- where the variant found affecting RS
  - GENE: gene symbol
  - TRANSCRIPT: transcript ID (e.g. ENST123456789)
  - IVS#: the ranking number of the intron in the gene (e.g. IVS1, IVS2, IVS3)
  - IVS_SIZE: the size of the intron
  - RS_COUNT: the total count of RS in this intron
  - RS#: the ranking number of the RS in this intron (e.g. RS1, RS2, RS3)
  - RS_POS: the first position of the essential RS-site AGGT
  - RS_UP_DW: the distance to its upstream and downstream RS-sites (e.g. -1111_+2222)
  - CLIP: if the RS-site is supported by eCLIP-U2AF data (Y/N)
  - RS_CONSEQ: the predicted consequences (RS_AGGT, RS_BP, RS_BP2, RS_AGAIN, RS_DW5SS_x-nt, CRYP_RS_DW5SS_x-nt, CRYP_RS_UP3SS_x-nt), where x-nt suggesting the size between the paired cryptic splice sites


### Command & Parameters (BPHunter_VCF.py)
```
python DIVERS_VCF.py -i variants.vcf
```

Parameter | Type | Description | Default
----------|------|-------------|--------------
*-i*|file|variants in VCF format, with 5 fields (CHROM, POS, ID, REF, ALT)|N.A.

### Command & Parameters (BPHunter_VCF_batch.py)
```
python DIVERS_VCF_batch.py -d directory -s samplelist.txt -o output.txt
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
