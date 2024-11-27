# HCsig
High-grade serous ovarian carcinoma-derived copy number signatures in diagnostic and archival cervical samples
## Introduction
Ovarian cancer remains a leading cause of mortality, with nearly two-thirds of patients succumbing to the disease due to late-stage diagnosis and a lack of effective screening methods. Early detection is crucial, as women diagnosed with localized disease have a survival rate exceeding 90%. Most high-grade serous tubo-ovarian carcinomas (HGSC) originate from epithelial cells in the fallopian tubes, with precursor lesions such as serous tubal intraepithelial carcinomas (STICs) and earlier dysplastic cells harbouring TP53 mutations. These findings suggest that early tumor-driving mutations may be detectable in asymptomatic women. Using shallow whole-genome sequencing (sWGS), we investigated CNAs in cervical samples from HGSC patients and controls. We applied HGSC- and pan-cancer-derived copy number signatures and developed novel signatures (HCsig) using unsupervised clustering. A centroid-based prediction model based on HCsig achieved promising sensitivity and specificity in an external validation cohort. These findings highlight the potential of sWGS on standard cervical samples for non-invasive early detection of HGSC.

### Bioinformatic workflow
[BINP52_CNA_Framework](https://github.com/IngridHLab/BINP52_CNA_Framework), a pipeline to generate copy number profiles and detect copy number signatures from shallow whole genome sequencing (sWGS) samples.

The version of tools and packages to be used will be specified in each step (see Chapter 3). The scripts within the pipeline are based on Snakemake (6.15.1), Python (v3.11.6), R (v4.3.2), and Java (JDK 11).
- (1) Preprocessing. This step includes quality assessment and quality trimming on the raw reads. (`Fastp` will be used for QC and trimming, together with `fastqc` and `multiQC` to generate the QC reports.)
- (2) Alignment. The human reference genome will be indexed. The reads will be mapped to the reference genome. (`BWA` will be used for both indexing and alignment.)
- (3) Clean-up. After alignment, the SAM files will be sorted and the PCR duplicates will be marked and removed. Also, the .sorted.deduplicated.sam will be converted to BAM files. The BAM files will be indexed for later analysis. (`Picard` will be used for sorting SAM, marking duplicates, removing duplicates and converting SAM to BAM. `samtools` will be used for generating the clean_up stats and for indexing the BAM files.)
- (4) Relative copy number profile. The BAM files will be analyzed through fixed-size binning, filtering, correction, and normalization to generate the read counts per bin. This data will then be used for the segmentation of bins and for generating the relative copy number profile. (`QDNAseq` will be used for this step.)
- (5) Ploidy and cellularity solutions. The output file from `QDNAseq` contains a relative copy number, and we need to estimate ploidy and cellularity in our samples to generate our final absolute copy number profile for comparison. (`Rascal` will be used for this step to find the solutions that best fit our study samples.)
- (6) Absolute copy number profile. We will further use other information (such as TP53 allele frequency) to infer the tumour fraction to select the best ploidy and cellularity solution. We apply this best solution to our relative copy number profile and generate the final absolute copy number profile for each sample. (`Rascal` will be used for this step.)
- (7) Comparison with the recent HGSC signatures (n=7). The functions should be loaded from the github repository: https://bitbucket.org/britroc/cnsignatures.git .
- (8) Comparison with the Pan-Cancer signatures (n=17). The package `CINSignatureQuantification` will be used to generate the samply-by-component matrix for the Pan-Cancer chromosomal instability signatures.
- (9) Comparison with the panConusig signatures (n=25). Tools including `Battenberg` (`alleleCounter`, `impute2` and `beagle5` were included in this package), `ASCAT.sc` and `panConusig` will be used in this step.  
- (10) Generate and validate HCsig signatures (n=5) and the HGSC prediction model. Tools including [SRIQ clustering](https://github.com/sunnyveerla/SRIQ), and [CentriodClassification](https://github.com/sunnyveerla/CentriodClassification) will be used in this step.

### Analysis pipeline execution steps:
#### (1) Preprocessing, (2) Alignment, (3) Clean-up.
```bash
```
#### (4) Relative copy number profile, (5) Ploidy and cellularity solutions, (6) Absolute copy number profile.
```bash
```
#### (7) Generate sample-by-component matrix using the GitHub package [CNsignatures](https://bitbucket.org/britroc/cnsignatures.git)
```bash
```
#### (10) Predict HCsig signatures (n=5) using [CentriodClassification](https://github.com/sunnyveerla/CentriodClassification).
```bash
```
