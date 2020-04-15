# QC and Preprocessing Scripts Overview
- **Author(s)**
	- jupyter notebooks written by Frank Grenn
- **Date Started** 
	- March 2020


## `Combine_Fastq_Lanes.ipynb`
- code to combine separate lane fastqs for each sample. 
- each sample started with 4 lanes per read. (read1 and read2, so 4x2=8 total fastq files per sample)
- combines fastqs to give one per read per sample. So each sample ends with a `samplename_R1.fastq` and `samplename_R2.fastq`

## `run_fastqc.ipynb`
- code to run fastqc tool on the sample fastqs

