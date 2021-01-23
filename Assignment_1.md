---
title: "BMEG 400E: Assignment 1"
output:
  html_document:
    toc: true
    toc_depth: 4
    toc_float: true
    code_folding: hide
---

This is the homework submission for assignment 1 BMEG591E for Charlotte Barclay.

# Actions
## 0. Starting Github

 a. Create a private repository on Github - 1pt -- DONE
   
      https://github.com/cbarcl01/CB_BMEG591E-repository.git

  b. Add **BMEGGenoInfo** as collaborator - 1pt -- DONE
    

  c. Clone the repository on your local computer - 1pt --DONE
  
## 2. Getting Familiar with the Server

**Windows system:** 

  a. Install a terminal emulator like Putty (*https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html*). 
  b. This will allow you to do a SSH connection using the server's address and your credentials. 
  
**Student Answer**

I have downloaded Putty to use the SSH connection and access the server

### a. Creating a directory 

**Student Answer**

```{bash, eval=False}
mkdir Assignment_1
```

### b. Check a directory's path 

When using terminal, paths work as the addresses of the files and directories you want to access and easily move between them.


**Student Answer**

```{bash, eval=False}
pwd
/home/cbarcl01
```

### c. Moving within directories

Access the newly created directory


How would you move back to your home directory (*/home/your_user_name*)?


**Student Answer**

***Move to newly created directory***

```{bash, eval=False}
cd /home/cbarcl01/Assignment_1
```

***To get to the home directory using the complete directory path***

```{bash, eval=False}
cd /home/cbarcl01
```

***To get to the home directory using the shortcut***

```{bash, eval=False}
cd ..
```

### d. Explore the dataset 

The sequencing data that we will be using is paired-end. This means that each DNA fragment in the sequencer has been sequenced twice, one on each end (5' and 3'). Choose one of the reads files (1 or 2) for the following exercises.


**Student Answer**

***To view the first 5 rows of the dataset***

```{bash, eval=False}
head -5 ./SRR12506919_1_subset.fastq
```

*Output*

@SRR12506919.1 1 length=151
GNTTCTCAGAGGCCCGGAGCGAGGCGAAGGCGTCGCCTTCTCGGAACCTGCCTCTGCCGGGCGCTCGCCGCTTCCCAGCGCGGGCTCGGGCCGCCTGCGTCCCCTAGCTAGGCGGGGTCGGTGCCCCGAGGGAGGGCCGGGGAGCGGGGCA
+SRR12506919.1 1 length=151
F#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF::FFFFFFFFFFFFFFFFFF:FFF:FFFF,FFFFFFFFFFFF,
@SRR12506919.2 2 length=151


***To view the last 6 rows of the dataset***
```{bash, eval=False}
tail -6 ./SRR12506919_1_subset.fastq
```

*Output*

+SRR12506919.667556 667556 length=151
FFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@SRR12506919.667557 667557 length=151
CGCAACAGCCGCGGATCCACCCCGGGGCGCGGGGGCCCAGCCAGGACGGCGGGAGACGGAGGCCCAGAGACCCCCCCTCACCCGCAGCCCGGACACTGCAGCCGGACCGGGTCCCTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
+SRR12506919.667557 667557 length=151
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF::FFFFFFF:FFFFFFFFFFFF:FFFF:FF:FFFF:FFF,:FFFFF


***Explore the file with the less command***
```{bash, eval=False}
less ./SRR12506919_1_subset.fastq
```
used q to exit


### e. Piping 

Because this is a very large dataset we will proceed to subset it in order to make it more manageable for the assignment. Using the commands that you learned above:


**Student Answer**

***Total number of lines of the file***

```{bash, eval=FALSE}
wc -l SRR12506919_1_subset.fastq
```

*Output*: 2670228


***Select only the id lines***

```{bash, eval=FALSE}
grep @SRR12506919. SRR12506919_1_subset.fastq
```

*Output*: 667557 lines of id reads, last line reads<br> @SRR12506919.667557 667557 length=151

***How many reads are in the file***

```{bash, eval=FALSE}
grep @SRR12506919. SRR12506919_1_subset.fastq | wc -l
```

*Output*: 667557


***Select only the reads id***

```{bash, eval=FALSE}
grep @SRR12506919. SRR12506919_1_subset.fastq | cut -d ' ' -f 1
```


### f. Saving an output 


**Student answer**

***Save a file that contains only the read ids***

```{bash, eval=FALSE}
grep @SRR12506919. SRR12506919_1_subset.fastq | cut -d ' ' -f 1 > /home/cbarcl01/Assignment_1/read_ids.txt
```

***Now list all the files in a directory***

```{bash, eval=FALSE}
cd /home/cbarcl01/Assignment_1
ls
```

*Output*: 

read_ids.txt


***What do you see? Was the subset file created correctly?***

Yes, see result in example above


### g. Creating a backup 

There will be times where you will want to save a copy of a dataset or file before modifying in case something goes wrong. For those cases, you can create a copy of any file or directory using the "copy" command


**Student answer**

***Create a copy of the reads ids file***

```{bash, eval=FALSE}
cp read_ids.txt backup_read_ids.txt
```

***Change the name of the backup reads id file***

```{bash, eval=FALSE}
mv backup_read_ids.txt copy_read_ids.txt
```

***Delete the copy***

```{bash, eval=FALSE}
rm copy_read_ids.txt
```

## 3. Managing Conda Environments

### a. Create a conda environment

In order to run the reads alignments against the human genome, there are a few tools that we will need:

  - fastQC (*https://www.bioinformatics.babraham.ac.uk/projects/fastqc/*): comprehensive quality control measures of sequencing data.
  
  
  - bowtie2 (*http://bowtie-bio.sourceforge.net/bowtie2/index.shtml*): alignments of sequencing data. 
  
  - samtools (*http://www.htslib.org/doc/samtools.html*): set of tools to manage sam/bam alignment files 
  
  
To install them, we will be making use of the conda environments. Conda allows you to create and manage environments for any programming language. Managing this environments mean that you can work with specific versions of different programs at any given moment, simply by loading the desired environment. You can find more information about this resource here: *https://docs.conda.io/en/latest/* . 


**Student answer**

```{bash, eval=FALSE}
conda create -n Gnme_Assignment_1
```

*Output*:

Collecting package metadata (current_repodata.json): done <br>
Solving environment: done


### b. Add programs to your conda environment

Now that the environment has been created, its time to add the packages that we will need. Conda has an active community and a great documentation (*https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html*) which you can use throughout the course to help answer any questions you may have. 


**Student answer**

***Add fastqc and bowtie2***

```{bash, eval=FALSE}
conda install fastqc
conda install bowtie2
```

***Install samtools***

```{bash, eval=FALSE}
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -n Gnme_Assignment_1 samtools==1.11
```

*Output*:

Collecting package metadata (current_repodata.json): done
Solving environment: done


## 4. Performing Alignments 

### a. Data quality check

We will use the widely used fastQC software to do a quick inspection of the data quality. Once it has run, it will give you an html report of the quality of the data, that you can open using a web browser. More information on how to read the output can be found here: *https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf* and on the tool's website.

 
**Student answer**

***Activate your conda environment***

```{bash, eval=FALSE}
conda activate Gnme_Assignment_1
```

***Run fastqc on fastq files***

```{bash, eval=FALSE}
fastqc *fastq -o /home/cbarcl01
```

***Copy files to your computer***

```{bash, eval=FALSE}
pscp -P 22 cbarcl01@gi-edu-sv4.bme.ubc.ca:/home/cbarcl01/SRR12506919_1_subset_fastqc.zip C:\Users\cbarc\OneDrive\Desktop
```

***Comments on data quality***

***1. SRR12506919_1_subset***

*Summary Output* <br>
PASS	Basic Statistics	SRR12506919_1_subset.fastq <br>
PASS	Per base sequence quality	SRR12506919_1_subset.fastq <br>
PASS	Per sequence quality scores	SRR12506919_1_subset.fastq <br>
PASS	Per base sequence content	SRR12506919_1_subset.fastq <br>
WARN	Per sequence GC content	SRR12506919_1_subset.fastq <br>
PASS	Per base N content	SRR12506919_1_subset.fastq <br>
PASS	Sequence Length Distribution	SRR12506919_1_subset.fastq <br>
PASS	Sequence Duplication Levels	SRR12506919_1_subset.fastq <br>
FAIL	Overrepresented sequences	SRR12506919_1_subset.fastq <br>
PASS	Adapter Content	SRR12506919_1_subset.fastq <br>

*FAIL for overrepresented sequences if any sequence is found to represent more than 1% of the total*


***2. SRR12506919_2_subset***

*Summary Output* <br>
PASS	Basic Statistics	SRR12506919_2_subset.fastq <br>
PASS	Per base sequence quality	SRR12506919_2_subset.fastq <br>
PASS	Per sequence quality scores	SRR12506919_2_subset.fastq <br>
PASS	Per base sequence content	SRR12506919_2_subset.fastq <br>
FAIL	Per sequence GC content	SRR12506919_2_subset.fastq <br>
PASS	Per base N content	SRR12506919_2_subset.fastq <br>
PASS	Sequence Length Distribution	SRR12506919_2_subset.fastq <br>
PASS	Sequence Duplication Levels	SRR12506919_2_subset.fastq <br>
PASS	Overrepresented sequences	SRR12506919_2_subset.fastq <br>
PASS	Adapter Content	SRR12506919_2_subset.fastq <br>

*FAIL per sequence GC if the sum of the deviations from the normal distribution represents more than 30% of the reads*
  

### b. Running a process on the background: screen 

The processes that we are about to run, can take a long time to finish. Thus, we will make use of the *screen* tool, which allows you to run a process in the background while continue using your server session without the risk of ending the process due to a bad internet connection, accidentally closing the tab or other random circumstances. 


**Student Answer**

***Start a background process***

```{bash, eval=FALSE}
screen -S background_CB 
```

***Run a command***

```{bash, eval=FALSE}
wc -l usr/local/share/data/assignment_1/SRR12506919_1_subset.fastq
```

***Get out of background screen***

```{bash, eval=FALSE}
ctrl + A + D
```

***return to background screen***

```{bash, eval=FALSE}
screen -r background_CB
```


### c. Checking the server capacity and your current runs: htop

Another way to check on the progress of your runs and the server status is with the *htop* command. This will open a screen showing all the processes that are being currently being run in the server. Our server only has 2 CPUs/cores, the green bar next to each code, represents how much of that node it is currently in use. Please be considerate when running jobs and do not monopolize the server that we all share.



**student answer**

***check server capacity and current runs***

```{bash, eval=FALSE}
htop
```
*Output*: 

[![htop output](C:\Users\cbarc\OneDrive\Pictures\htop_output.jpg)]


### d. Genome indexing - bowtie2

Now, we will need to create an index of the human genome that bowtie2 will use to map our sequences in the genome. In order to do this, you will need to use the previously downloaded files of the human genome with the desired build (e.g. hg19,hg38), you can find those files within the server here: */usr/local/share/human/*

Indexes take a lot of computational resources and time to run. The one we need for the human genome would take around 3 hours to be done.  For this question, you can just enter the command here. You can run the command to see that Bowtie2 starts to build an index, and kill it with `ctrl+c` once you have confirmed it is working. 


### e. Alignment

We are working with paired-end data. Thus, you will need to make sure to use both fastq files to align the sequences to the genome.
**IMPORTANT:** Run with *default parameters* DO NOT specify any non-essential paramaters.

**Time flag**: This step will take up to 30 mins 


**Student answer**

```{bash, eval=FALSE}
bowtie2 -x /usr/local/share/indexes/hg38_bowtie2_index \ -1 /usr/local/share/data/assignment_1/SRR12506919_1_subset.fastq \ -2 /usr/local/share/data/assignment_1/SRR12506919_2_subset.fastq \ -S /home/cbarcl01/eg3.sam
```

*Output*

667557 reads; of these:
  667557 (100.00%) were paired; of these:
    97993 (14.68%) aligned concordantly 0 times
    466591 (69.90%) aligned concordantly exactly 1 time
    102973 (15.43%) aligned concordantly >1 times
    ----
    97993 pairs aligned concordantly 0 times; of these:
      26971 (27.52%) aligned discordantly 1 time
    ----
    71022 pairs aligned 0 times concordantly or discordantly; of these:
      142044 mates make up the pairs; of these:
        111490 (78.49%) aligned 0 times
        11802 (8.31%) aligned exactly 1 time
        18752 (13.20%) aligned >1 times
91.65% overall alignment rate


### f. Viewing the alignments

Now, we will make use of **samtools** to review some basic features of the alignment. For more information about this tool: *http://www.htslib.org/doc/*


**Student answer**

***Mapped reads***

611812

***Unmapped reads***

55745

***Command used***

I ran the following commands (-F for mapped, -f for unmapped), then divided the output by 2 to account for forward and reverse strands

```{bash, EVAL=false}
samtools view -c -F 4 eg3.sam
samtools view -c -f 4 eg3.sam
```




