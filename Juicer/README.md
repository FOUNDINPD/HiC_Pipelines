# Juicer Scripts Overview
- **Author(s)**
	- credit to [aidenlab](https://github.com/aidenlab/juicer) for the juicer pipeline
	- jupyter notebooks written by Frank Grenn
- **Date Started** 
	- March 2020
- **Quick Description:** 
	- Scripts to run the samples through the juicer pipeline and analyze results.

## Scripts
listed in the general order they were ran. Everything coded to run on biowulf.

### 1) `Juicer_Pipeline_Setup.ipynb`
- sets up command to run juicer pipeline for a sample on biowulf. Uses an altered version of `juicer.sh` from the juicer pipeline that removes data on mitochondrial chromosomes. 

### 2) `Gather_Juicer_Stats.ipynb`
- code to collect output statistics from the juicer pipeline to compare

### 3) `HIC_overlap_examples.ipynb`
- example code for checking for overlap between sample `merged_loops.bedpe` files using bedtools.
- uses functions from `overlapBedpeAndBed.py`. Examples in the notebook run using test data in `example_test_files`.
- includes code to:
	- overlap sample bedpes with each other
	- overlap sample bedpes with one bedpe file
	- overlap sample bedpes with one bed file
	- randomly sample the sample bedpes and check for overlap with a bed file

### 4) `HIC_overlap_between_samples.ipynb`
- uses functions from `overlapBedpeAndBed.py` to check for overlap between samples using bedtools.
- output matrix files show loop count or percent overlap with respect to the column sample.
##### Output:
    - `sample_overlap_count.csv`
    - `sample_overlap_percent.csv`

### 5) `HIC_overlap_PsychENCODE.ipynb`
- uses functions from `overlapBedpeAndBed.py` to check for PsychENCODE loop overlap with samples
- includes code to process the PsychENCODE data for comparison (map their hg19 positions to hg38)
##### Output:
    - `psychencode_loop_overlap.csv`

### 6) `Juicer_Tools_Compare_Lists.ipynb`
- uses the [Compare](https://github.com/aidenlab/Juicer/wiki/Compare-Lists) tool from juicer to compare the sample bedpes. 
- Different from using bedtools for overlap, because this uses a radius threshold to check for euclidean distance between loops to check for overlap.
- output matrix files show loop count or percent overlap with respect to the column sample.
##### Output:
    - `compare_lists_25000`
    - `compare_lists_default`

### 7) `Juicer_Tools_Huccups_Diff.ipynb`
- uses the [hiccupsdiff](https://github.com/aidenlab/Juicer/wiki/HiCCUPSDiff) tool from juicer to look for differential loops between samples. 
- This reruns the Hiccups algorithm on the samples and looks for enrichment in annotated loop regions. Loops are only counted as being differential to one sample if the loop has no enrichment at all in the other sample. 
##### Output:
    - `hiccups_diff_progression.csv`