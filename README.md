# Assessing Pathogenicity in Genetic Variants
![images](https://github.com/Gho-Ost/pathogenicity-assessment/assets/100212265/f5215263-d2a4-4920-9f4c-de4b7a59d365)

## Introduction and Problem Description

The problem of assessing pathogenicity in genetic variants is a problem of classification. The classification system was created by The American College of Medical Genetics and Genomics and the Association for Molecular Pathology (ACMG-AMP). 5 classes can be distinguished: benign, likely benign, variant of unknown significance (VUS), likely pathogenic, and pathogenic.

---

## Data

Raw data format (vcf) specification: https://samtools.github.io/hts-specs/VCFv4.1.pdf

VEP (CSQ) outputs: http://www.ensembl.org/info/docs/tools/vep/vep_formats.html#vcfout

To read .csv converted data: 

```py
from utils.utils import get_data

EE_015 = get_data("data/EE_015/", ["default", "csq"])
EE_050 = get_data("data/EE_050/", ["default", "csq"])
EE_069 = get_data("data/EE_069/", ["default", "csq"])

df = pd.concat([EE_015, EE_050, EE_069], ignore_index=True, axis=0)
```

---

## File structure

```

└───archive
└───data
    ├───EE_sample*
    ├───EE_015
    ├───EE_050
    └───EE_069
          ├───EE_069.vcf.gz
          ├───EE_069_default.csv.gz
          ├───EE_069_genotype.csv.gz
          └───EE_069_csq.csv.gz
```

*EE_sample contains uncompressed files

---

[Working Documentation](https://docs.google.com/document/d/1QrPL4XlauwmgChU2wR5oaxm3lQT9XguRHlArkw-dGnk/edit?fbclid=IwAR0bLvaZl5aDMawowjTp23NeM8kCLT2UjOY_lNQLWdC-6atJqYklR94vMzc)
