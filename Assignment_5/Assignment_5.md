Assignment 5
================

  - [Assignment Overview](#assignment-overview)
      - [0. Getting Started](#0-getting-started)
      - [1. ChIP signal tracks](#1-chip-signal-tracks)
          - [1. a) convert to bigwig](#1-a-convert-to-bigwig)
          - [1. b) Copy to computer and load to
            IGV](#1-b-copy-to-computer-and-load-to-igv)
      - [2. Narrow vs Broad peaks](#2-narrow-vs-broad-peaks)
          - [2. a) computeMatrix](#2-a-computematrix)
          - [2. b) Heatmap using
            Deeptools](#2-b-heatmap-using-deeptools)
      - [3. Peak calling](#3-peak-calling)
          - [a. Peak calling with macs2](#a-peak-calling-with-macs2)
          - [b. Understanding the peaks](#b-understanding-the-peaks)
      - [c. Peak enrichments](#c-peak-enrichments)

# Assignment Overview

By now you must have become vaguely familiar with ChIP-seq data but
might be a little confused about what to do after the alignment. Well,
this assignment’s aim is to walk you through a ChIP-seq pipeline
post-alignment. We will be analyzing 3 different histone modification
marks (H3K27me3, H3K4me3 and H3K27ac). In order to identify the
enrichments for each epigenetic mark, we also need to use the *input*
which represents the DNA content of the sheared chromatin sample prior
to immunoprecipitation. All the files can be found under the following
path: **/usr/local/share/data/assignment\_5/** .

  - H3K27me3 (H3K27me3\_chr3\_subset.bam)

  - H3K4me3 (H3K4me3\_chr3\_subset.bam)

  - H3K27ac (H3K27ac\_chr3\_subset.bam)

  - input (input\_chr3\_subset.bam)

A couple of things to remember:

  - When uploading your completed assignment to your GitHub directory,
    remember to specify a **github\_document** instead of the default
    *html\_document* on your .Rmd file.

  - Double check that all the files have been uploaded to your
    repository and that you are able to visualize the html-like version
    on GitHub.

  - Be mindful of your space on the server\! Delete ALL unnecessary
    files and keep a clean environment.

## 0\. Getting Started

We will be using a couple of new tools this time. Before we move on to
the practical part, make sure you have them all installed.

  - Integrative Genomics Viewer (IGV): Interactive tool to visualize
    different data types of genetic information (e.g. bam, bed files).
    You will install this tool to your **local computer**. To visualize
    where the reads of our ChIP analysis mapped in the genome. To
    install it, follow the instructions on this website:
    *<https://software.broadinstitute.org/software/igv/home>*

  - Deeptools
    (<https://deeptools.readthedocs.io/en/develop/index.html>): Software
    to analyze high-throughput data that allows to create easy to
    visualize figures. This will be installed on the server as part of
    your conda environment.

  - macs2
    (<https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md>):
    Tool to capture enrichments of ChIP-seq analysis. This will be
    installed on the server as part of your conda environment.

<!-- end list -->

``` bash

#?# Add macs2 and deeptools to your environment created on A1 - 1 pt

conda activate Gnme_Assignment_1 
conda install -c bioconda macs2
conda install -c bioconda deeptools


## Install IGV to your local computer after downloading it from: https://software.broadinstitute.org/software/igv/home
```

## 1\. ChIP signal tracks

ChIP-seq experiments require a control against which to compare the ChIP
signnal. Often this is the “input DNA”, sheared chromatin that has not
been immmunoprecipitated. For this assignment we are using an input and
three different epigenetic marks. These histone modifications mark
states of active (H3K27ac, H3K4me3) or inactive (H3K27me3) gene
transcription and have different coverage of the genomic region where
they are located. To better visualize the differences, we will create
bigWig files from previously aligned, filtered and indexed bam files.
BigWig files are indexed, compressed files that can be used to visualize
signals across the genome. Here, we will be using it to graph the
coverage across the genome.

### 1\. a) convert to bigwig

``` bash

## Use the bamCoverage command (included in deepTools) to convert the bam files outlined below (located in the /usr/local/share/data/assignment_5/ directory) into bigwig files that represent the read coverage (e.g. density of reads across the genome). Each file has already been indexed using sambamba and hence, has it's bam index (.bai) file associated that you will need to run bamCoverage. 
## Tip: Remember to add the .bw extension to your specified output! 

#?# Type the commands you use for converting each of the four bam files to bigwig files - 2 pts

mkdir ./Assignment_5
cd ./Assignment_5 ##making and moving to new subdirectory for ease of file organisation


bamCoverage -b /usr/local/share/data/assignment_5/H3K27ac_chr3_subset.bam -o ./H3K27ac_chr3_subset.bw
bamCoverage -b /usr/local/share/data/assignment_5/H3K27me3_chr3_subset.bam -o ./H3K27me3_chr3_subset.bw
bamCoverage -b /usr/local/share/data/assignment_5/H3K4me3_chr3_subset.bam -o ./H3K4me3_chr3_subset.bw
bamCoverage -b /usr/local/share/data/assignment_5/input_chr3_subset.bam -o ./input_chr3_subset.bw


### resulting files:

ls

H3K27ac_chr3_subset.bw   H3K4me3_chr3_subset.bw
H3K27me3_chr3_subset.bw  input_chr3_subset.bw
```

### 1\. b) Copy to computer and load to IGV

``` bash

### Follow this steps: 

## 1. Copy the bigwig files from the previous steps to your computer 


pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/home/cbarcl01/Assignment_5/H3K27ac_chr3_subset.bw C:\Users\cbarc\OneDrive\Desktop
pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/home/cbarcl01/Assignment_5/H3K4me3_chr3_subset.bw C:\Users\cbarc\OneDrive\Desktop
pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/home/cbarcl01/Assignment_5/H3K27me3_chr3_subset.bw C:\Users\cbarc\OneDrive\Desktop
pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/home/cbarcl01/Assignment_5/input_chr3_subset.bw C:\Users\cbarc\OneDrive\Desktop


## 2.Load all the files bigwig signal track files onto IGV on your local computer, select the "autoscale" option for each file on their individual tracks
## Tip, use the "File" tab to "load from file" option to choose the files from your computer directories

## 3. Change the visualization of the files to see the following position: ** chr3:93,470,124-93,471,058 **

#?# Include a screenshot of your IGV session right after this code chunk (using the Rmarkdown syntax)  - 2 pt

See below

#?# Explore this region by zooming in and out and navigating around. What do you see? Is there something similar across all the different files that stands out on this region? Is there anything peculiar about the DNA sequence at this locus?- 3 pt

This region exibits very similar expression across both the active and inactive markers



## Tip: sometimes track signal will be truncated to a pre-set maximum. If you right-click the track label (left) and select "autoscale", it will automatically adjust the scales of the track to the range of your data. Try playing around with the track settings to see their effect.

## 4. This file (/usr/local/share/data/assignment_5/hg38_blacklisted_regions.bed) contains the hg38 blacklisted regions. Load it into IGV along with all the other files. 

### 4.1 first I copied to my local computer then I loaded it into IGV 

pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/usr/local/share/data/assignment_5/hg38.blacklist.bed C:\Users\cbarc\OneDrive\Desktop

## 5. Look at the following region again (chr3:93,470,124-93,471,058). 
#?# What can you say now about the similarities between the files at this region? 

The similarities between these files align with the blacklisted region, potentially therefore what we are seeing is due to the noise at the blacklisted region rather than true signal reads.

#?#In your own words, explain what a blacklisted region is and if you think this region should be excluded a ChIP-seq analysis. - 1.5 pt


Blacklisted regions are regions within the genome that consistently have high (or "erroneous") signal/read counts, regardless of the ChIP-seq experiment/cell type being examined. These regions can impact any data normalisation or lead to false positive peaks so should be excluded. 
```

**Add screenshot of your IGV session here**

![IGV\_screenshot1](https://github.com/cbarcl01/CB_BMEG591E-repository/blob/master/Assignment_5/IGV_Screenshot.png)

## 2\. Narrow vs Broad peaks

While exploring the bigwig files of the epigenetic marks on IGV, you
probably noticed that they can look very different from each other and
some of them ressemble the input more closely than others. Different
epigenetic marks can have very different signatures based on their
distribution patterns across the genome. When calling peaks, we often
lump these into two different categories of reads: broad and narrow.
Active transcription marks (H3K4me3 and H3K27ac) tend to form a sharper
coverage peaks at transcription start sites (H3K27ac is also at
enhancers), while repression marks cover a broader area (H3K27me3).
Here, we’re going to inspect their distributions relative to genes.

### 2\. a) computeMatrix

``` bash


## Here, we have created three bigWig track files, one for each epigenetic mark, which show the read coverage normalized using the input. They are found here: /usr/local/share/data/assignment_5/
# H3K4me3_norm.bw
# H3K27me3_norm.bw
# H3K27ac_norm.bw


## The deepTools ** computeMatrix reference-point ** command calculates scores to represent the reads mapping to specified regions of the genome across different files. 
## Use computeMatrix to compute a matrix for the signal tracks for each histone modification outlined above (which we will use to create a plot in the following step), with the following criteria: 

## - We will use the regions in reference_genes.bed located under the /usr/local/share/data/assignment_5/ directory as the basis for the plot.
## - Include the surrounding 1kb
## - Use all 3 input-normalized bigWig files (H3K4me3, H3K27ac, H3K27me3) as signal tracks
#?# Write the command you used to run it below: - 1.5 pt

computeMatrix scale-regions -S /usr/local/share/data/assignment_5/H3K4me3_norm.bw /usr/local/share/data/assignment_5/H3K27ac_norm.bw /usr/local/share/data/assignment_5/H3K27me3_norm.bw -R /usr/local/share/data/assignment_5/reference_genes.bed -a 1000 -b 1000 -o ./computematrixRP.gz
```

### 2\. b) Heatmap using Deeptools

``` bash

## Now that the scores matrix has been computed, we can use it to create a heatmap to provide a better visual representation of the reads distrubution across our reference genes (provided in the reference_gened.bed file)
## Use the deepTools ** plotHeatmap ** function to create a heatmap following this criteria: 
## - Use the matrix from the previous point
## - Use the Blues colormap
## - Create 3 clusters within the heatmap according to the patterns of the reads distrubution across the files using heirarchical clustering
#?# Type the command you used to run it below: - 1.5

plotHeatmap -m ./computematrixRP.gz \
     -out Heatmap1.png \
     --colorMap 'Blues' --kmeans 3
     
##then used the below code to copy to local comp for screenshot     
pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/home/cbarcl01/Assignment_5/Heatmap1.png C:\Users\cbarc\OneDrive\Desktop 

#?# Add a screenshot of the plot right after this code chunk using Rmarkdown syntaxis - 1 pt 


#?# Explain what you are looking at (Axes, colours, curves). Where are the marks located? What are the differences between the clusters? - 3 pts

The heat maps display data in a grid where each row represents a gene and each column represents a sample. The colour and intensity of the boxes is used to represent changes (not absolute values) of gene expression. The bottom axes displays gene distance in base pairs with the Transcription Start Site and Transcription End Site marked. The 3 clusters show high/medium/low expression regions of the genome.  

For the line graphs above, the height of the line reflects the presence of epigenetic marks.

For the 2 active epigenetics marks (H3K27ac, H3K4me3), H3K4me3 seems to be more prevalent. The epigenetic marks appear to be found near the Transcription Start Site

Cluster 1 seems to have the hughest itesnity out of the the three. Cluster 1 also seems to show a peak in intensity immediately after the TSS, with decrease in intensity suggesting silencing at TES for the active H3K4me3 and H3K27ac. For both active marks the highest expression (the darkest colour) in Cluster 2 is located before the TSS with values decreasing immediately after the TSS. Cluster 3 seems to be more homogenous/constant across the markers, with a narrower eak in expression at the TSS in H3K4me3.
```

**Add screenshot here**

![Heatmap for normalised bigwig
files](https://github.com/cbarcl01/CB_BMEG591E-repository/blob/master/Assignment_5/Heatmap1.png)

``` bash
## Now the above heatmap was made with the ratio of ChIP to input. Repeat the process above, but this time using the raw bigwig files (not input-normalized). 
#?# Type the commands you used for this below - 1 pt (H3K4me3, H3K27ac, H3K27me3)

computeMatrix scale-regions -S  /home/cbarcl01/Assignment_5/H3K4me3_chr3_subset.bw  /home/cbarcl01/Assignment_5/H3K27ac_chr3_subset.bw /home/cbarcl01/Assignment_5/H3K27me3_chr3_subset.bw -R /usr/local/share/data/assignment_5/reference_genes.bed -a 1000 -b 1000 -o ./computematrix2.gz

plotHeatmap -m ./computematrix2.gz \
     -out Heatmap2.png \
     --colorMap 'Blues' --kmeans 3
     
pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/home/cbarcl01/Assignment_5/Heatmap2.png C:\Users\cbarc\OneDrive\Desktop 


#?# Include a screenshot of this analysis, below this code block. - 1 pt

#?# How does this compare to the input-normalized data? Why do you think this is? - 1 pt

Because the data is not input normalised, the signals appear stronger so the intensities can be more clearly defined. However this is not necessarily reflective of the relative difference in expression.
```

**Add screenshot here**

![Heatmap with raw bigwig
files](https://github.com/cbarcl01/CB_BMEG591E-repository/blob/master/Assignment_5/Heatmap2.png)

***How does this compare to the input-normalized data? Why? 1 pt***

## 3\. Peak calling

Now we want to identify enriched regions of the genome for each of our
three histone marks. In order to get the enrichments, we will run the
**macs2** program to call the peaks for each epigenetic mark.

### a. Peak calling with macs2

``` bash

## Tip: Make sure to read the documentation (using the -h flag) for the *masc2 callpeak* command

## Run the callpeak command of the macs2 tool, once for each of H3K27ac, H3K27Me3, H3K4Me3
## Each time, use the input as the control 
## For H3K27Me3, you should call "broad" peaks because these signals tend to form longer domains and not acute peaks.

#?# Type the commands you used below: - 1.5 pt

macs2 callpeak -t /usr/local/share/data/assignment_5/H3K27ac_chr3_subset.bam -c /usr/local/share/data/assignment_5/input_chr3_subset.bam -n H3K27ac_callpeaks

macs2 callpeak -t /usr/local/share/data/assignment_5/H3K4me3_chr3_subset.bam -c /usr/local/share/data/assignment_5/input_chr3_subset.bam -n H3K4me3_callpeaks

macs2 callpeak -t /usr/local/share/data/assignment_5/H3K27me3_chr3_subset.bam -c /usr/local/share/data/assignment_5/input_chr3_subset.bam --broad -n H3K27me3_callpeaks
```

### b. Understanding the peaks

Macs2 calls the peaks by analyzing the enrichments of sequences in the
genome compared to the input. The algorithm uses different measures to
ensure the enrichments are statistically significant and make biological
sense based on the sequence length, the expected peaks (narrow or broad)
and other parameters. We will focus on the H3K4me3 mark to visualize how
are the peaks are identified.

``` bash

## 1. Copy the H3K4me3 .narrowPeak file to your local computer

pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/home/cbarcl01/Assignment_5/H3K4me3_callpeaks_peaks.narrowPeak C:\Users\cbarc\OneDrive\Desktop 


## 2. Open IGV with and load the hg38 reference 

Done

## 3. Load the following files:

### a. H3K4me3 bigwig file that you created in section 1 
### b. Input  bigwig file that you created in section 1 
## Note: remember to autoscale each file track!

## 4. Go to this position: * chr3:44,608,952-44,670,849 *

Done

#?# Compare the input and H3K4me3 signal tracks (a and b), are there regions that you consider to have an enriched signal in the H3K4me3 track compared to the input? If yes, how many?- 0.5 pt 

Yes, two regions of 'peaks' can be seen in the H3K4me3 signal track.

#?# Why do you consider those regions to have an enrichment? Would you need to confirm that the enrichment is statistically significant? If so, what do you think would be a way to test for the significance (you can mention the tools used on this assignment or explain conceptually what would you do to confirm an enrichment)? - 1 pt

These regions could have enrichment because there appears to be an increase in coverage in these two regions. However this is just an observation, so in order to confirm if the enrichment is statistically significant, we would need to calculate the p-value which is the probability that the result observed could occur under the null hypothesis and test for this.We can use Macs2 call peaks feature to analyzing the enrichments of sequences compared to the input and test for statistical significance.

## 5. Load into the same session the following files:

### Copy across to local

pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/usr/local/share/data/assignment_5/H3K4me3_norm.bw C:\Users\cbarc\OneDrive\Desktop 

pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/usr/local/share/data/assignment_5/H3K4me3_callpeaks_peaks.narrowPeak C:\Users\cbarc\OneDrive\Desktop 

### c. H3K4me3 normalized using the input bigWig file that was pre-computed and you used on section 2
### d. H3K4me3 .narrowPeak file that you created in section 3a
## Note: remember to autoscale each file track!

## 6. Go to the same position as in step 4

#?# Does the region/s you thought had an enrichment show differences between the input and the H3K4me3 mark in the bamCompare (file c)? - 0.5 pt

Yes the two peaks identified in the previous question differ from the input file and the H3K4me3_callpeaks_peaks.narrowpeak confirms this.


#?# Are these visually enriched regions you selected, found in file d (narrowPeak)? What does that mean? - 1 pt

Yes this is supported by the presence of the bar at the same regions in the H3K4me3_callpeaks_peaks.narrowpeak track. This means that this region is statistically signifiant in H3K4me3, compared to the input file (based on the algorithim run in the call peaks function).

#?# Add a screenshot of your IGV session right after this chunk, using Rmarkdown syntax - 1pt
```

***SCREENSHOT OF IGV***

![IGV\_Screenshot2](https://github.com/cbarcl01/CB_BMEG591E-repository/blob/master/Assignment_5/IGV_Screenshot2.png)

## c. Peak enrichments

For this assignment, we are working with 3 different epigenetic marks:
H3K4me3, H3K27me3 and H3K27ac. Each of these marks a different
transcription state (i.e., activation or repression). Thus, to better
visualize the differences in the called peaks for each of these
epigenetic marks you will create a heatmap plot (like the one you
created on part 2) using their called peaks files (.narrowPeak or
.broadPeak).

``` bash


### Create 3 heatmaps following the specifications you used on part 2 (one for the peaks called for each epigenetic mark, bot containing data from all three tracks). Use the peak files of each of the epigenetic marks as reference files. Use ONLY the non-input normalized files: H3K27ac_norm.bw H3K27me3_norm.bw H3K4me3_norm.bw

#?# Write the commands you used to compute the matrices and create the heatmaps below: - 3 pt

computeMatrix scale-regions -S /usr/local/share/data/assignment_5/H3K4me3_norm.bw /usr/local/share/data/assignment_5/H3K27ac_norm.bw /usr/local/share/data/assignment_5/H3K27me3_norm.bw -R /home/cbarcl01/Assignment_5/H3K27ac_callpeaks_peaks.narrowPeak -a 1000 -b 1000 -o ./computematrixH3K27ac.gz

plotHeatmap -m ./computematrixH3K27ac.gz \
     -out HeatmapH3K27ac.png \
     --colorMap 'Blues' --kmeans 3

pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/home/cbarcl01/Assignment_5/HeatmapH3K27ac.png C:\Users\cbarc\OneDrive\Desktop 


computeMatrix scale-regions -S /usr/local/share/data/assignment_5/H3K4me3_norm.bw /usr/local/share/data/assignment_5/H3K27ac_norm.bw /usr/local/share/data/assignment_5/H3K27me3_norm.bw -R /home/cbarcl01/Assignment_5/H3K4me3_callpeaks_peaks.narrowPeak -a 1000 -b 1000 -o ./computematrixH3K4me3.gz

plotHeatmap -m ./computematrixH3K4me3.gz \
     -out HeatmapH3K4me3.png \
     --colorMap 'Blues' --kmeans 3
     
pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/home/cbarcl01/Assignment_5/HeatmapH3K4me3.png C:\Users\cbarc\OneDrive\Desktop 


computeMatrix scale-regions -S /usr/local/share/data/assignment_5/H3K4me3_norm.bw /usr/local/share/data/assignment_5/H3K27ac_norm.bw /usr/local/share/data/assignment_5/H3K27me3_norm.bw -R /home/cbarcl01/Assignment_5/H3K27me3_callpeaks_peaks.broadPeak  -a 1000 -b 1000 -o ./computematrixH3K27me3.gz

plotHeatmap -m ./computematrixH3K27me3.gz \
     -out HeatmapH3K27me3.png \
     --colorMap 'Blues' --kmeans 3
     
pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/home/cbarcl01/Assignment_5/HeatmapH3K27me3.png C:\Users\cbarc\OneDrive\Desktop 


#?# Add screenshots of the 3 heatmaps you got using the epigenetic marks' peak files as reference files. Add them after this code chunk in the following order: H3K4me3, H3K27ac, H3K27me3 - 1.5 pt
```

***H3K4me3 screenshot***

![Heatmap for
H3K4me3](https://github.com/cbarcl01/CB_BMEG591E-repository/blob/master/Assignment_5/HeatmapH3K4me3.png)

***H3K27ac screenshot***

![Heatmap for
H3K27ac](https://github.com/cbarcl01/CB_BMEG591E-repository/blob/master/Assignment_5/HeatmapH3K27ac.png)

***H3K27me3 screenshot***

![Heatmap for
H3K27me3](https://github.com/cbarcl01/CB_BMEG591E-repository/blob/master/Assignment_5/HeatmapH3K27me3.png)

``` bash

#?# Do you see an overlap between the peaks of different epigenetic marks? Which epigenetic marks? - 1 pt

Yes, there is an overlap between H3K27ac and H3K4me3

#?# Why do you think these epigenetic marks overlap? - 1 pt

These are both describing active transcription state, so would expect there to be overlap

#?# Explain the pattern you see for the peaks of H3K27me3, do they look similar or different to the other two marks? Why do you think this happens? - 1pt



#?# Why do you think the borders of the elements have such clearly-defined borders with respect to the ChIP signal? -1 pt

As it is comparing the to the call to peak statistical anaylsis, it is comparing to a binary yes/no option which makes the borders much more clearly defined, with respect to the TSS and TSE.
```
