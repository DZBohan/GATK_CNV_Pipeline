# GATK CNV Pipeline

#### Bohan Zhang September 13th

## Index

- [1. Introduction of the Pipeline](#1)

	- [1.1 Preparation step](#1.1)

	- [1.2 CNV analysis step](#1.2)

	- [1.3 Gene level CNV analysis](#1.3)

- [2. Required Files List](#2)

	- [2.1 Project level files](#2.1)

	- [2.2 Sample level files](#2.2)

- [3. Files Editing](#3)

	- [3.1 Bam files viewing](#3.1)

	- [3.2 Intervals.bed](#3.2)

	- [3.3 Reference.dict](#3.3)

	- [3.4 Geneinfo.txt](#3.4)

	- [3.5 Prepare the bamIdsUniq.txt](#3.5)

- [4. Environment Setting](#4)

	- [4.1 Installation of GATK](#4.1)

	- [4.2 R environment setting](#4.2)

	- [4.3 HPC modules loading](#4.3)

- [5. Introduction of Scripts](#5)

	- [5.1 Preparation step](#5.1)

	- [5.2 Main step](#5.2)

- [6. Introduction of Results](#6)

	- [6.1 CNV plots](#6.1)

	- [6.2 Segment level CNV table](#6.2)

	- [6.3 Gene level CNV table](#6.3)


## <h2 id="1">1. Introduction of the Pipeline</h2>

This pipeline uses eight GATK tools and an R package called CNTools. Here are the GATK tools used.

* PreprocessIntervals
* AnnotateIntervals
* CollectReadCounts
* CreatReadCountPanelOfNormals
* DenoiseReadCounts
* CollectAllelicCounts
* ModelSegments
* PlotModeledSegments

Here are the official web pages introducing the [GATK CNV analysis steps](https://gatk.broadinstitute.org/hc/en-us/articles/360035531092) and the [CNTools package](https://bioconductor.org/packages/release/bioc/html/CNTools.html).

The pipeline is theoretically divided into three steps, the `preparation step`, the `CNV analysis step`, and the `gene level CNV analysis step`, but practically, the pipeline includes two scripts, a preparation script and a main script containing both the `CNV analysis step`, and the `gene level CNV analysis step`.

Now, I am going to have a brief introduction to the three steps first.

### <h2 id="1.1">1.1 Preparation step</h2>

Here is a flowchart of the preparation step. The wave shape stands for files, and the rectangle stands for tools.

![Prestep_flow](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/preflow.png?raw=true)

This step aims to generate the interval files, including a `.interval_list` file and a `intervals.tsv` file, which are needed in the following step.

### <h2 id="1.2">1.2 CNV analysis step</h2>

Here is also a flowchart of the CNV analysis step.

![CNA_flow](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/mainflow.png?raw=true)

In this step, you are first collecting the read counts of both the tumor bam file and the normal bam file and collecting allelic counts of both the bam files. Then, you are supposed to generate a panel using the normal bam file.

The denoised step is required before doing the CNV analysis. By doing the denoised step, you will get a more straightforward plot, which can help you know the CNV of a sample easier. The denoised step is also one of the advantages of the GATK CNV analysis process.

![Denoised](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/denoised.png?raw=true)

The next thing is the CNV analysis. You will use the ModelSegments tool to get the CNV results tables and the PlotModeledSegments tool to get the CNV results plots. The results tables and plots will be introduced in chapter 6.

### <h2 id="1.3">1.3 Gene level CNV analysis</h2>

GATK tools provide the CNV analysis on the segment level. However, if we want to know the CNV of a specific gene, we need to introduce another tool called `CNTools`. With the help of genomic information, CNTools can convert the output table of GATK CNV pipeline from the segment level to the gene level. 


## <h2 id="2">2. Required Files List</h2>

### <h2 id="2.1">2.1 Project level files</h2>

These six files are the same among different samples in one project.

* intervals.bed (exist but need editing)
* reference.fa (exist)
* reference.fai (along with reference.fa)
* reference.dict (exist but need editing)
* bamIdsUniq.txt (need to be generated yourself)
* geneinfo.txt (can be downloaded)

### <h2 id="2.2">2.2 Sample level files</h2>

These four files are sample-specific.

* tumor.cram or bam (exist)
* tumor.crai or bai (along with tumor.carm or bam)
* normal.cram or bam (exist)
* normal.crai or bai (along with normal.carm or bam)

## <h2 id="3">3. Files Editing</h2>

Two input files, intervals.bed, and reference.dict, need editing. Here is the process of editing these two files.

### <h2 id="3.1">3.1 Bam files viewing</h2>

The first step of editing the files is to view the bam files you are analyzing to see the format people are using to represent different chromosomes.

By the way, bam files and cram files are the same in this pipeline, so I will always say bam files to stand for both bam and cram files in this documentation.

On CARC, you can load the module samtools for viewing the bam files.

```
module load samtools
```
Use this command to view the bam file.

```
samtools view filename.bam
```
Here are the two examples of the bam files.

![Image1](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/bam_example_1.png?raw=true)

![Image2](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/bam_example_2.png?raw=true)

As you can see in the third column of the two examples, there are two ways, "1" and "chr1", to represent the chromosome 1. It is important to know how to represent the chromosome number in the bam files.

### <h2 id="3.2">3.2 Intervals.bed</h2>

There are usually several bed files for a set of samples. I recommend you to use targets.bed or regions.bed.

Now, let's view a bed file.

![Bed](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/bed.png?raw=true)

Three places in the bed file may need to be edited. First, if the file has some header lines, you need to remove them. Second, the first column of the bed file is the chromosome number, so you need to unify the format of the chromosome number, "N" or "chrN" with the bam file of the set of samples. You can use the %s function in vim to edit the bed file. Third, you need to remove the lines, not on chromosomes 1-22 or x and y. Use this command to remove these lines as well as the header.

```
egrep "^chr[0-9XY]" filename.bed > filename_edit.bed
```
Or this command if the chromosome numbers are numbers only.

```
egrep "^[0-9XY]" filename.bed > filename_edit.bed
```

### <h2 id="3.3">3.3 Reference.dict</h2>

Let us view a dictionary file first.

![Dict](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/dict.png?raw=true)

The editing of the dictionary file is the same as the bed file. First, use the %s function in vim to unify the format of the chromosome number. Second, use this command to remove the header and lines, not on chromosomes 1-22 or x and y.

```
egrep "^@SQ SN:chr[0-9XY]" filename.dict > filename_edit.dict
```

Or this command if the chromosome numbers are numbers only.

```
egrep "^@SQ SN:[0-9XY]" filename.dict > filename_edit.dict
```

### <h2 id="3.4">3.4 Geneinfo.txt</h2>

Geneinfo.txt is the table of information on each gene. It includes genes' locations, names, and ids. Usually, we use GRCh38 as the geneinfo file, and you can directly download the GRCh38 geneinfo file from my GitHub.

![Geneinfo](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/geneinfo.png?raw=true)

You do not need to modify the chromosome format ("chrN" or "N") of the geneinfo file since I have a chromosome format judgment step in the CNTools R script.

If you want to use your geneinfo file, you need to remove the duplicate data and data, not on chromosomes 1-22 or X and Y.

### <h2 id="3.5">3.5 Prepare the bamIdsUniq.txt</h2>

Using this pipeline, you are required to create a bamIdsUniq.txt file. The file should have two columns. The first column is all the tumor bam file names, and the second column is all the normal bam file names. Use commas to separate two columns. Here is an example of what the file looks like.

![BamIdsUniq](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/bamidsuniq.png?raw=true)


## <h2 id="4">4. Environment Setting</h2>

### <h2 id="4.1">4.1 Installation of GATK</h2>

Download the GATK with the version of 4.2.6.1

```
wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
```

Unzip the GATK zip file.

```
unzip gatk-4.2.6.1.zip
```

Go into the `gatk-4.2.6.1` directory and test the GATK.

```
./gatk --help
```

Put the path of GATK into the .bashrc file.

```
vim ~/.bashrc
```

```
export PATH=$PATH:/add/your/path/here/gatk-4.2.6.1
```

```
source ~/.bashrc
```

### <h2 id="4.2">4.2 R environment setting</h2>

In this pipeline, I used R with the version 4.1.2. You need to install four packages, optparse, data.table, ggplot2, and CNTools, in the r/4.1.2.

```
install.packages("optparse")
install.packages("data.table")
install.packages("ggplot2")

if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("CNTools")
```

### <h2 id="4.3">4.3 HPC modules loading</h2>

In this pipeline, I load four modules on the HPC. They are avaliable on the USC CARC, so if you are using another HPC, please whether these modules are avaliable.

* module load gcc/11.2.0
* module load jdk/17.0.1
* module load openblas/0.3.18
* module load r/4.1.2

You don't need to load the modules yourself since I have written them inside the scripts. If some of the modules are not avaliable on the HPC you are using, try to install the programs yourself.

## <h2 id="5">5. Introduction of Scripts</h2>

This GATK CNV analysis pipeline includes two main steps put into two Slurm scripts, `gatk_cnv_prepare.slurm` and `gatk_cnv.slurm`. Both scripts have config files, `config_gatk_cnv.txt` and `config_gatk_cnv_prepare.txt`, for inputting the variables or files. In addition, there is an R script, `cntools_gatk.R`, that should be used in the second Slurm script. You can download these five files on GitHub. You should locate these files in the same directory when using the pipeline.

### <h2 id="5.1">5.1 Preparation step</h2>

The first step, also called the preparation step, is for the whole set of samples since they share the same reference and bed files.

This step aims to generate two intermedia files, `targets.preprocessed.interval_list` and `annotated_intervals.tsv` for the next step.

For inputting the variables and files, you do not need to modify the script, `gatk_cnv_prepare.slurm`, but write them into the config file `config_gatk_cnv_prepare.txt`. Now, let us have a look at the config file.

![Config_Prepare](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/config_prepare.png?raw=true)

There are three variables in this config file, `output_path`, `bed` and `refer`. `output_path` is the path of a directory to store the output files, `targets.preprocessed.interval_list` and `annotated_intervals.tsv`. `bed` is the absolute path of the edited .bed file. `refer` is the absolute path of the .fa file. One important thing is that you should put the .fai in the same directory as the .fa file.

After completing the config file, you can submit the Slurm script, `gatk_cnv_prepare.slurm`, to run the preparation step.

```
sbatch --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=4GB --time=1:00:00 gatk_cnv_prepare.slurm
```
You can modify the memory, time, and other parameters, but they are usually enough in this step.

Now, you should have the two new files, `targets.preprocessed.interval_list` and `annotated_intervals.tsv`, in the target directory and be able to go to the next step.

### <h2 id="5.2">5.2 Main step</h2>

You can do the main step with the two newly generated intermedia files. This step contains a Slurm script `gatk_cnv.slurm`, a config file `config_gatk_cnv.txt`, and a R script `cntools_gatk.R`.

For inputting the variables and files, you need to add them into the config file `config_gatk_cnv.txt`. Now, let us have a look at this one.

![Config](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/config.png?raw=true)

`BAMDIR` is the absolute path of the directory storing bam files. You need to put the set of bam files (bai files as well) in one directory.

`FILE` is the absolute path of the bamIdsUniq.txt file you create in chapter 3.5.

`intlist` and `inttsv` are the absolute path of the intermedia files generated in chapter 5.1.

`refer` is the absolute path of the .fa file. One important thing is that you should put the .fai in the same directory as the .fa file.

`output_path` is the absolute path of the directory for storing the GATK CNV pipeline outputs. Finally, there should be a set of directories named by the tumor bam file names in this directory. Inside individual directories, you can find the output files of each sample.

`dict` is the absolute path of the dictionary file you edit in chapter 3.3.

`geneinfo_path` is the absolute path of the geneinfo file mentioned in chapter 3.4. You can download the GRCh38 version of the geneinfo file through this GitHub.

`modelseg_mem` is the java memory requirement of the ModelSegments step, which usually requires larger memory. 70G of memory is usually enough for the general bam files, around 20G - 40G. However, larger java memory for this step is required for some large bam files, maybe larger than 50G. You can change the job partition into `largemem`, and set this parameter as 256G or larger. 

After completing the config file, you can submit the Slurm script, `gatk_cnv.slurm`, to run the main step.

```
sbatch --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=80GB --time=8:00:00 --array=[1-2]%2 gatk_cnv.slurm
```

The parameter of the array depends on the number of samples you want to do. In other word, the total number of the array is supposed to be the same as the number of rows in the `bamIdsUniq.txt` files. The details of `bamIdsUniq.txt` can be found in the chapter 3.5. For example, if you have 50 samples requiring CNV analysis, and you want to run five jobs each time, you are supposed to set the array parameter as `--array=[1-50]%5`.

I set the default value of the memory parameter as 80GB since the ModuleSegments step usually requires 70GB. If your sample's bam files are too large, you may face an error called `java out of memory`. You can change the parameter of `partition` into `largemem`, and set the larger value of both `mem` parameter and the `modelseg_mem` variable inside the `config_gatk_cnv.txt` file. Here is an example of the large bams' jobs.

```
sbatch --partition=largemem --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=500GB --time=8:00:00 --array=[1-2]%2 gatk_cnv.slurm
```

It usually takes 2-5 hours to run the pipeline once, meaning if you set the array as `--array=[1-50]%5`, it will take 20-50 hours to run all 50 samples. After running the main step, you will get a set of directories, which are the same number as the samples, inside the output directory you set, and the name of each directory is supposed to be the same as the tumor bam's filename.

![Output1](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/output1.png?raw=true)

Inside each individual directory, there should be 12 files (five `.seg` files, two `.tsv` files, four `.param` files and a `.png` image) and a directory called `cntools_result`.

![Output2](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/output2.png?raw=true)

This main step should have generated several intermedia files, but you cannot see them in the final results since I added a clean-up step in the script. Here I also list all intermedia files generated by this step in case you are interested in them.

```
tumor_bam_name_t.counts.hdf5
normal_bam_name_g.counts.hdf5
normal_bam_name_cnv.pon.hdf5
tumor_bam_name.standardizedCR.tsv
tumor_bam_name.denoisedCR.tsv
tumor_bam_name_t.allelicCounts.tsv
normal_bam_name_g.allelicCounts.tsv
```

## <h2 id="6">6. Introduction of Results</h2>

The results of the GATK can be divided into three categories, CNV plots, segment level CNV table, and gene level CNV table. First, let us have a look at the CNV plots.

### <h2 id="6.1">6.1 CNV plots</h2>

Here is an example of the GATK pipeline CNV plots.

![CNV_plots](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/Sample1.modeled.png?raw=true)

You can finally get two plots from running the GATK pipeline on a set of tumor normal samples. The x-axis of both plots shows the different segments on different chromosomes. The difference between the two plots is the y-axis. 

The y-axis of the first plot is `denoised copy ratio`, meaning the copy number of the tumor sample divided by the copy number of the normal sample. For example, 3/2 means a single copy gain, and 1/2 means a single copy gain. 

For the color, the official webpage explained this. "Different copy ratio segments are indicated by alternating blue and orange color groups
The denoised median is drawn in thick black"

The y-axis of the second plot is `alternate allele fraction`. Allele fraction can be used to infer whether a mutation at a specific locus is somatic or germline. It is calculated as the proportion of sequence reads mutated at a locus divided by the total coverage of the locus. 

An allele fraction of roughly 50% for diploid means that the mutation in this sample at the particular locus is a heterozygous germline mutation, i.e., recombination from the parents' germ cells. 

An allele fraction of 100% or 0 means that the allele is homozygous at the locus, but this cannot determine whether it is a somatic or germline mutation. The allele fraction of the reference genome is also 0. If the allele fraction is not 50%, 100%, or 0, then there are potential errors other than germline mutations, i.e., somatic mutations, at that locus for that sample. 

Allele fraction's representation of somatic variation requires consideration of three additional issues. First, tumor samples often contain contamination by normal cells with only germline mutations, so the correspondence between allele fraction and somatic mutations also needs to include contamination by the normal cells in the calculation. 

Second, changes in allele copy number or chromosome ploidy can intrinsically affect the allele fraction of heterozygous mutations. However, the allele fraction of homozygotes is independent of copy number or ploidy and is always 0 or 100%. 

Third, the allele fraction for germline mutations is not 50% because some large germline mutations, such as large insertions or deletions (INDEL), can alter allele copy numbers, even the ploidy. Therefore, the allele fraction should be used to determine tumors' somatic mutations and match normal samples.

### <h2 id="6.2">6.2 Segment level CNV table</h2>

As I mentioned in chapter 5.2, the GATK CNV pipeline has 12 output files. The `.seg` files are the result tables you are looking for. There are five `.seg` files in total.

```
tumor_bam_filename.modelBegin.seg
tumor_bam_filename.modelFinal.seg
tumor_bam_filename.cr.seg
tumor_bam_filename.cr.igv.seg
tumor_bam_filename.af.igv.seg
```

The first two tables are the initial and final log2 copy ratio results before and after the segmentation step. Both have the log2 copy ratio results at different positions, 10, 50, and 90, so it is hard to get final readable results for these two files.

`tumor_bam_filename.cr.seg` and `tumor_bam_filename.cr.igv.seg` essentially include the same data，segments' position information, number of the probs and the mean value of log2 copy ratio, but the latter is more arrangeable and perfect being the input file of the next step. `tumor_bam_filename.af.igv.seg` has same format as `tumor_bam_filename.cr.igv.seg`, and the only difference is that the log2 copy ratio in the last column is replaced with the allele fraction.

Therefore, `tumor_bam_filename.cr.igv.seg` and `tumor_bam_filename.af.igv.seg` are the more informative results. You can get the segment-level CNV information by having a look at these two tables. Here are the first ten rows of these two tables.

![Cr](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/cr.png?raw=true)

![Af](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/af.png?raw=true)

### <h2 id="6.3">6.3 Gene level CNV table</h2>

After running the GATK CNV pipeline, there is supposed to be a directory called cntools_result in the output directory, and you can find a file called `tumor_bam_filename.cntools.gatk.txt`.

![Cntools](https://github.com/DZBohan/GATK_CNV_Pipeline/blob/main/images/cntools.png?raw=true)

In this file, you can get the genenames, geneids, genes' position, and genes' log2 copy ratio.

