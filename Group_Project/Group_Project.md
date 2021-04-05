Group Project: Investigating the validity of *Mnemiopsis leidyi* as a
model organism for the study of multicellularity
================
Charlotte Barclay and Gabriel d’Alba
31/03/2021

  - [Introduction](#introduction)
  - [Method](#method)
  - [Results](#results)
  - [Conclusion](#conclusion)

``` r
pscp -P 22 C:\Users\cbarc\OneDrive\Documents\BMEG591E_Genome_Informatics\Group_Project\Mnemiopsis_leidyi.MneLei_Aug2011.dna.nonchromosomal.fa.gz cbarcl01@gi-edu-sv4.bme.ubc.ca:/home/cbarcl01/Group_Project

pscp -P 22 C:\Users\cbarc\OneDrive\Documents\BMEG591E_Genome_Informatics\Group_Project\AGCP01.1.fsa_nt.gz cbarcl01@gi-edu-sv4.bme.ubc.ca:/home/cbarcl01/Group_Project
```

## Introduction

overview of the original study, the scope of your re-analysis, and why
you chose it

![Schema of phylogenetic position of Mnemiopsis
leidyi](C:\\Users\\cbarc\\OneDrive\\Desktop\\git_temp\\CB_BMEG591E-repository\\Group_Project\\Mle.jpg)

## Method

integrated data processing, QC, analysis, results, graphs, and other
data, as well as written explanations for what is being done and why,
and interpretations of results. This should flow in chronological order
(e.g. starting with fastqs and ending with the last graph).

**Create an Index**

``` bash
bowtie2-build /home/cbarcl01/Group_Project/Mnemiopsis_leidyi.MneLei_Aug2011.dna.nonchromosomal.fa.gz MleIdx
```

**Install BLAT**

``` bash
conda install -c bioconda ucsc-blat
```

**Alignment**

``` bash
#Unzip fas_nt.ga format
gunzip /home/cbarcl01/Group_Project/AGCP01.1.fsa_nt.gz > /home/cbarcl01/Group_Project/AGCP01_extract.fasta #note I tried this without the output defined first

bowtie2 -x /home/cbarcl01/Group_Project/MleIdx \ -U /home/cbarcl01/Group_Project/AGCP01.1.fsa_nt  \ -S /home/cbarcl01/Group_Project/TEST.sam
```

## Results

integrated data processing, QC, analysis, results, graphs, and other
data, as well as written explanations for what is being done and why,
and interpretations of results. This should flow in chronological order
(e.g. starting with fastqs and ending with the last graph).

## Conclusion

Summarize your findings and contrast this with what the original study
found. Remark on anything surprising or anything you would do
differently next time.
