#!/bin/bash

# Modified by Wolfgang mostly following Susan's modifications for the
# previous version of the pipeline.

##########
#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########
# Alignment script. Sets the reference genome and genome ID based on the input
# arguments (default human, MboI). Optional arguments are the queue for the 
# alignment (default short), description for stats file, 
# using the short read aligner, read end (to align one read end using short 
# read aligner), stage to relaunch at, paths to various files if needed,
# chunk size, path to scripts directory, and the top-level directory (default 
# current directory). In lieu of setting the genome ID, you can instead set the
# reference sequence and the chrom.sizes file path, but the directory 
# containing the reference sequence must also contain the BWA index files.
#
# Splits the fastq files, creates jobs to align them, creates merge jobs that
# wait for the alignment to finish, and creates a final merge job.
#
# Also creates "cleanup" jobs that at each stage, deletes jobs off the cluster
# if any one of them fails.
#
# If all is successful, takes the final merged file, removes name duplicates,
# removes PCR duplicates, and creates the hic job and stats job.  Final
# product will be hic file and stats file in the aligned directory.
#                                                                       
# [topDir]/fastq  - Should contain the fastq files. This code assumes that
#                   there is an "R" in the appropriate files, i.e. *R*.fastq
# From the top-level directory, the following two directories are created:
#                                                                              
# [topDir]/splits  - Where to write the scratch split files (fastq files and
#                    intermediate SAM files). This can be deleted after 
#                    execution.
# [topDir]/aligned - Where to write the final output files.
#
# The following globals should be set correctly before proceeding:
#
# splitsize - The number of lines that each split fastq should contain. Larger
#             means fewer files and longer overall, but too small means there
#             are so many jobs that the cluster won't run them. This can be
#             set with the -C command as well
# read1str  - portion of fastq filename that indicates this is the "read 1"
#             file; used to loop over only the read 1 and within that loop,
#             also align read 2 and merge.  If this is not set correctly,
#             script will not work. The error will often manifest itself
#             through a "*" in the name because the wildcard was not able to
#             match any files with the read1str.   
# Juicer version 1.5.6
shopt -s extglob
juicer_version="1.5.6"

###
### Biowulf settings
###

load_bwa="module load bwa/0.7.17" 
load_java="" 
load_gpu="module load CUDA/8.0" 
load_coreutils="module load coreutils/8.27"
# Juicer directory, contains scripts/, references/, and restriction_sites/
# can also be set in options via -D
juiceDir="/usr/local/apps/juicer/juicer-1.5.6"
# default queue, can also be set in options via -q
queue="norm"
queue_time="1440"
# default long queue, can also be set in options via -l
long_queue="norm"
long_queue_time="2880"

# set the locale
export LC_ALL=C

# default value for bwa threads
threads=12

###
### general settings
###

# size to split fastqs. adjust to match your needs. 4000000=1M reads per split
# can also be changed via the -C flag
splitsize=90000000

# fastq files should look like filename_R1.fastq and filename_R2.fastq 
# if your fastq files look different, change this value
read1str="_R1" 
read2str="_R2" 

# unique name for jobs in this run
groupname="a$(date +%s)"

## Default options, overridden by command line arguments

# top level directory, can also be set in options
topDir=$(pwd)
# restriction enzyme, can also be set in options
site="DpnII"
# genome ID, default to human, can also be set in options
genomeID="hg19"
# normally both read ends are aligned with long read aligner; 
# if one end is short, this is set                 
shortreadend=0
# description, default empty
about=""
nofrag=0

echo "Running juicer version ${juicer_version}"

## Read arguments                                                     
usageHelp="Usage: ${0##*/} [-g genomeID] [-d topDir] [-q queue] [-l long queue] [-s site]\n                 [-a about] [-R end] [-S stage] [-p chrom.sizes path]\n                 [-y restriction site file] [-z reference genome file]\n                 [-C chunk size] [-D Juicer scripts directory]\n                 [-Q queue time limit] [-L long queue time limit] [-b ligation] [-t threads]\n                 [-r] [-h] [-x]"
genomeHelp="* [genomeID] must be defined in the script, e.g. \"hg19\" or \"mm10\" (default \n  \"$genomeID\"); alternatively, it can be defined using the -z command"
dirHelp="* [topDir] is the top level directory (default\n  \"$topDir\")\n     [topDir]/fastq must contain the fastq files\n     [topDir]/splits will be created to contain the temporary split files\n     [topDir]/aligned will be created for the final alignment"
queueHelp="* [queue] is the queue for running alignments (default \"$queue\")"
longQueueHelp="* [long queue] is the queue for running longer jobs such as the hic file\n  creation (default \"$long_queue\")"
siteHelp="* [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" \n  (default \"$site\")"
aboutHelp="* [about]: enter description of experiment, enclosed in single quotes"
shortHelp="* -r: use the short read version of the aligner, bwa aln\n  (default: long read, bwa mem)"
shortHelp2="* [end]: use the short read aligner on read end, must be one of 1 or 2 "
stageHelp="* [stage]: must be one of \"merge\", \"dedup\", \"final\", \"postproc\", or \"early\".\n    -Use \"merge\" when alignment has finished but the merged_sort file has not\n     yet been created.\n    -Use \"dedup\" when the files have been merged into merged_sort but\n     merged_nodups has not yet been created.\n    -Use \"final\" when the reads have been deduped into merged_nodups but the\n     final stats and hic files have not yet been created.\n    -Use \"postproc\" when the hic files have been created and only\n     postprocessing feature annotation remains to be completed.\n    -Use \"early\" for an early exit, before the final creation of the stats and\n     hic files"
pathHelp="* [chrom.sizes path]: enter path for chrom.sizes file"
siteFileHelp="* [restriction site file]: enter path for restriction site file (locations of\n  restriction sites in genome; can be generated with the script\n  misc/generate_site_positions.py)"
chunkHelp="* [chunk size]: number of lines in split files, must be multiple of 4\n  (default ${splitsize}, which equals $(awk -v ss=${splitsize} 'BEGIN{print ss/4000000}') million reads)"
scriptDirHelp="* [Juicer scripts directory]: set the Juicer directory,\n  which should have scripts/ references/ and restriction_sites/ underneath it\n  (default ${juiceDir})"
refSeqHelp="* [reference genome file]: enter path for reference sequence file, BWA index\n  files must be in same directory"
queueTimeHelp="* [queue time limit]: time limit for queue, i.e. -W 12:00 is 12 hours\n  (default ${queue_time})"
longQueueTimeHelp="* [long queue time limit]: time limit for long queue, i.e. -W 168:00 is one week\n  (default ${long_queue_time})"
ligationHelp="* [ligation junction]: use this string when counting ligation junctions"
threadsHelp="* [threads]: number of threads when running BWA alignment. default: ${threads}"
excludeHelp="* -x: exclude fragment-delimited maps from hic file creation"
helpHelp="* -h: print this help and exit"
helpHelp2="* -H: print resource allocation table"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$genomeHelp"
    echo -e "$dirHelp"
    echo -e "$siteHelp"
    echo -e "$aboutHelp"
    echo -e "$shortHelp"
    echo -e "$shortHelp2"
    echo -e "$stageHelp"
    echo -e "$pathHelp"
    echo -e "$siteFileHelp"
    echo -e "$refSeqHelp"
    echo -e "$chunkHelp"
    echo -e "$scriptDirHelp"
    echo -e "$ligationHelp"
    echo -e "$threadsHelp"
    echo "$excludeHelp"
    echo -e "\n-- resource allocation options----------------------------------------"
    echo -e "$queueHelp"
    echo -e "$longQueueHelp"
    echo -e "$queueTimeHelp"
    echo -e "$longQueueTimeHelp"
    echo "----------------------------------------------------------------------"
    echo "$helpHelp"
    echo "$helpHelp2"
    echo "* -v: verbose"
    exit "$1"
}

###
### keep this in sync with any changes to jobs below !!!
###
printResourceTable () {
    echo -e "\nDefaults: \$queue_time      = $queue_time"
    echo -e "          \$queue           = $queue"
    echo -e "          \$long_queue_time = $long_queue_time"
    echo -e "          \$long_queue      = $long_queue"
    echo -e "          \$threads         = $threads (for bwa)"
    echo -e "          \$bwa_mem         = min(40G, \$threads * 5G)\n"
    echo -e "defaults can be changed with [-q queue] [-l long queue] and"
    echo -e "[-Q queue time limit] [-L long queue time limit]\n"
    cat <<END_OF_TABLE  | column -t -s'|'
Job|Partition|CPU|Mem[GiB]|GPU|lscratch|time[min]
--------------|---------|---|--------|---------|--------|---------------
_cmd|\$queue|2||||2
_split|\$queue|2||||\$queue_time / 3
_splitwait|\$queue|2||||5
_CountLigation|\$queue|2||||\$queue_time / 3
_align1|\$queue|\$threads|\$bwa_mem|||\$queue_time
_align2|\$queue|\$threads|\$bwa_mem|||\$queue_time
_merge|\$long_queue|8|40G||800|\$long_queue_time
_check|\$queue|||||\$queue_time
_fragmerge|\$long_queue|8|246G||800|\$long_queue_time
_dedup_guard|\$queue|2||||5
_dedup|\$queue|2|2G/CPU|||\$queue_time / 2
_post_dedup|\$queue|2||||100
_dupcheck|\$queue|2|1G/CPU|||\$queue_time / 3
_stats|\$queue|2|10G/CPU|||\$queue_time
_hic|\$long_queue|2|32G/CPU|||\$long_queue_time
_hic30|\$long_queue|2|32G/CPU|||\$long_queue_time
_hiccups_wrap|gpu|2|18G|gpu:k80:1||\$queue_time
_arrowhead_wrap|\$queue|2|18G|||\$queue_time
_prep_done|\$queue|2|2G/CPU|||\$queue_time
_prep_done|\$queue|2|2G/CPU|||\$queue_time
END_OF_TABLE
    exit
}

while getopts "d:g:R:a:hHrq:s:p:l:y:z:S:C:D:Q:L:b:t:xv" opt; do
    case $opt in
	g) genomeID=$OPTARG ;;
	h) printHelpAndExit 0;;
    H) printResourceTable ;;
	d) topDir=$OPTARG ;;
	l) long_queue=$OPTARG ;;
	q) queue=$OPTARG ;;
	s) site=$OPTARG ;;
	R) shortreadend=$OPTARG ;;
	r) shortread=1 ;;  #use short read aligner
	a) about=$OPTARG ;;
	p) genomePath=$OPTARG ;;  
	y) site_file=$OPTARG ;;
	z) refSeq=$OPTARG ;;
	S) stage=$OPTARG ;;
	C) splitsize=$OPTARG; splitme=1 ;;
	D) juiceDir=$OPTARG ;;
	Q) queue_time=$OPTARG ;;
	L) long_queue_time=$OPTARG ;;
	x) nofrag=1 ;;
	b) ligation=$OPTARG ;;
	t) threads=$OPTARG ;;
    v) echo "Running in verbose mode"; set -x ;;
	[?]) printHelpAndExit 1;;
    esac
done

if [[ ! "$threads" =~ [0-9]+ ]]; then
    echo "'$threads' is not a valid number of threads"
    exit 1
fi

if [[ ! "$queue_time" =~ ^[0-9]+$ ]]; then
    echo "-Q: queue_time needs to be specified in minutes"
    exit 1
fi

if [[ ! "$long_queue_time" =~ ^[0-9]+$ ]]; then
    echo "-L: long_queue_time needs to be specified in minutes"
    exit 1
fi


if [ ! -z "$stage" ]
then
    case $stage in
        merge) merge=1 ;;
        dedup) dedup=1 ;;
        early) earlyexit=1 ;;
        final) final=1 ;;
	postproc) postproc=1 ;; 
        *)  echo "$usageHelp"
	    echo "$stageHelp"
	    exit 1
    esac
fi

## Set reference sequence based on genome ID
if [ -z "$refSeq" ]
then 
    case $genomeID in
	mm9)	refSeq="${juiceDir}/references/mm9.fa";;
	mm10)	refSeq="${juiceDir}/references/mm10.fa";;
#	hg38)	refSeq="${juiceDir}/references/hg38/hg38.fa";;
	hg19)	refSeq="${juiceDir}/references/hg19.fa";;
#	hg18)	refSeq="${juiceDir}/references/hg18.fasta";;
	*)		echo "$usageHelp"
	    echo "$genomeHelp"
	    exit 1
    esac
else
    ## Reference sequence passed in, so genomePath must be set for the .hic 
    ## file to be properly created
    if [ -z "$genomePath" ]
    then
        echo "***! You must define a chrom.sizes file via the \"-p\" flag that delineates the lengths of the chromosomes in the genome at $refSeq";
        exit 1;
    fi
fi

## Check that refSeq exists 
if [ ! -e "$refSeq" ]; then
    echo "***! Reference sequence $refSeq does not exist";
    exit 1;
fi

## Check that index for refSeq exists
if [ ! -e "${refSeq}.bwt" ]; then
    echo "***! Reference sequence $refSeq does not appear to have been indexed. Please run bwa index on this file before running juicer.";
    exit 1;
fi

## Set ligation junction based on restriction enzyme
case $site in
    HindIII) ligation="AAGCTAGCTT";;
    DpnII) ligation="GATCGATC";;
    MboI) ligation="GATCGATC";;
    NcoI) ligation="CCATGCATGG";;
    none) ligation="XXXX";;
    *)  ligation="XXXX"
	echo "$site not listed as recognized enzyme. Using $site_file as site file"
	echo "Ligation junction is undefined"
esac

## If DNAse-type experiment, no fragment maps
if [ "$site" == "none" ]
then
    nofrag=1;
fi

## If short read end is set, make sure it is 1 or 2
case $shortreadend in
    0) ;;
    1) ;;
    2) ;;
    *)	echo "$usageHelp"
	echo "$shortHelp2"
	exit 1
esac

if [ -z "$site_file" ]
then
    site_file="${juiceDir}/restriction_sites/${genomeID}_${site}.txt"
fi

## Check that site file exists, needed for fragment number for merged_nodups
if [ ! -e "$site_file" ] && [ "$nofrag" -ne 1 ]
then
    echo "***! $site_file does not exist. It must be created before running this script."
    exit 1
fi

## Set threads for sending appropriate parameters to cluster and string for BWA call
threadstring="-t $threads"
bwa_mem=$((threads * 5000))

if [ $bwa_mem -gt 40000 ]
then
    bwa_mem=40000
fi

## Directories to be created and regex strings for listing files
splitdir=${topDir}"/splits"
donesplitdir=$topDir"/done_splits"
fastqdir=${topDir}"/fastq/*_R*.fastq*"
outputdir=${topDir}"/aligned"
tmpdir=${topDir}"/HIC_tmp"
debugdir=${topDir}"/debug"

## Check that fastq directory exists and has proper fastq files
if [ ! -d "$topDir/fastq" ]; then
    echo "Directory \"$topDir/fastq\" does not exist."
    echo "Create \"$topDir/fastq\" and put fastq files to be aligned there."
    echo "Type \"juicer.sh -h\" for help"
    exit 1
else 
    if stat -t ${fastqdir} >/dev/null 2>&1
    then
	echo "(-: Looking for fastq files...fastq files exist:"
    ls -lh ${fastqdir}
    else
	if [ ! -d "$splitdir" ]; then 
	    echo "***! Failed to find any files matching ${fastqdir}"
	    echo "***! Type \"juicer.sh -h \" for help"
	    exit 1		
	fi
    fi
fi

## Create output directory, only if not in postproc, dedup or final stages
if [[ -d "$outputdir" && -z "$final" && -z "$dedup" && -z "$postproc" ]] 
then
    echo "***! Move or remove directory \"$outputdir\" before proceeding."
    echo "***! Type \"juicer.sh -h \" for help"
    exit 1			
else
    if [[ -z "$final" && -z "$dedup" && -z "$postproc" ]]; then
        mkdir "$outputdir" || { echo "***! Unable to create ${outputdir}, check permissions." ; exit 1; } 
    fi
fi

## Create split directory
if [ -d "$splitdir" ]; then
    splitdirexists=1
else
    mkdir "$splitdir" || { echo "***! Unable to create ${splitdir}, check permissions." ; exit 1; }
fi

## Create temporary directory, used for sort later
if [ ! -d "$tmpdir" ] && [ -z "$final" ] && [ -z "$dedup" ] && [ -z "$postproc" ]; then
    mkdir "$tmpdir"
    chmod 777 "$tmpdir"
fi

## Create output directory, used for reporting commands output
if [ ! -d "$debugdir" ]; then
    mkdir "$debugdir"
    #chmod 777 "$debugdir"
fi

###
### Pipeline
###

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline
# If chunk size sent in, split. Otherwise check size before splitting
if [ -z $splitme ]
then
    fastqsize=$(ls -lL  ${fastqdir} | awk '{sum+=$5}END{print sum}')
    if [ "$fastqsize" -gt "2592410750" ]
    then
	splitme=1
    fi
fi

testname=$(ls -l ${fastqdir} | awk 'NR==1{print $9}')
if [ "${testname: -3}" == ".gz" ]
then
    read1=${splitdir}"/*${read1str}*.fastq.gz"
    gzipped=1
else
    read1=${splitdir}"/*${read1str}*.fastq"
fi





### JOB _cmd
# Add header containing command executed and timestamp:
jid=`sbatch <<- HEADER | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -t 2
	#SBATCH -c 2
	#SBATCH -o $debugdir/head-%j.out
	#SBATCH -e $debugdir/head-%j.err
	#SBATCH -J "${groupname}_cmd"
	date
	$load_bwa
	$load_java

	# Experiment description
	if [ -n "${about}" ]
	then
		echo -ne 'Experiment description: ${about}; '
	else
		echo -ne 'Experiment description: '
	fi

	# Get version numbers of all software
	echo -ne "Juicer version $juicer_version;" 
	bwa 2>&1 | awk '\\\$1=="Version:"{printf(" BWA %s; ", \\\$2)}'
	echo -ne "${threads} threads; "
	if [ -n "$splitme" ]
	then
		echo -ne "splitsize $splitsize; "
	fi  
	java -version 2>&1 | awk 'NR==1{printf("%s; ", \\\$0);}'
	${juiceDir}/scripts/juicer_tools -V 2>&1 | awk '\\\$1=="Juicer" && \\\$2=="Tools"{printf("%s; ", \\\$0);}'
	
	echo "$0 $@"
HEADER`
headfile="${debugdir}/head-${jid}.out"

## Record if we failed while aligning, so we don't waste time on other jobs
## Remove file if we're relaunching Juicer 
errorfile=${debugdir}/${groupname}_alignfail
if [ -f $errorfile ]
then
    rm $errorfile
fi

# Not in merge, dedup,  or final stage, i.e. need to split and align files.
if [ -z $merge ] && [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
    if [ "$nofrag" -eq 0 ]
    then
	echo -e "(-: Aligning files matching $fastqdir\n in queue $queue to genome $genomeID with site file $site_file"
    else
        echo -e "(-: Aligning files matching $fastqdir\n in queue $queue to genome $genomeID with no fragment delimited maps."
    fi
    
    ## Split fastq files into smaller portions for parallelizing alignment 
    ## Do this by creating a text script file for the job on STDIN and then 
    ## sending it to the cluster	
    ### JOB _split
    dependsplit="afterany"
    if [ ! $splitdirexists ]
    then
	echo "(-: Created $splitdir and $outputdir."
	if [ -n "$splitme" ]
        then
            for i in ${fastqdir}
            do
		filename=$(basename $i)
		filename=${filename%.*}      
                if [ -z "$gzipped" ]
                then	
		    jid=`sbatch <<- SPLITEND | egrep -o -e "\b[0-9]+$"
			#!/bin/bash -l
			#SBATCH -p $queue
			#SBATCH -t $((queue_time / 3))
			#SBATCH -c 2
			#SBATCH -o $debugdir/split-%j.out
			#SBATCH -e $debugdir/split-%j.err
			#SBATCH -J "${groupname}_split_${i}"
			date
			echo "Split file: $filename"
			${load_coreutils}
			split -a 3 -l $splitsize -d --additional-suffix=.fastq $i $splitdir/$filename
			date
SPLITEND`
		else
		    jid=`sbatch <<- SPLITEND | egrep -o -e "\b[0-9]+$"
			#!/bin/bash -l
			#SBATCH -p $queue
			#SBATCH -t $((queue_time / 3))
			#SBATCH -c 2
			#SBATCH -o $debugdir/split-%j.out
			#SBATCH -e $debugdir/split-%j.err
			#SBATCH -J "${groupname}_split_${i}"
			date
			echo "Split file: $filename"
			${load_coreutils}
			zcat $i | split -a 3 -l $splitsize -d --additional-suffix=.fastq - $splitdir/$filename
			date
SPLITEND`
		fi
		dependsplit="$dependsplit:$jid"
                # if we split files, the splits are named .fastq
                read1=${splitdir}"/*${read1str}*.fastq"
	    done
	    
		# dependencies don't work with srun from within an interactive session
	    #srun -c 2 -p "$queue" -t 1 -o $debugdir/wait-%j.out -e $debugdir/wait-%j.err -d $dependsplit -J "${groupname}_wait" sleep 1
		jid=`sbatch <<- SPLITWAITEND  | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $queue
		#SBATCH -t 5
		#SBATCH -c 2
		#SBATCH -o /dev/null
		#SBATCH -J "${groupname}_splitwait_${i}"
		#SBATCH -d ${dependsplit}
		touch ${splitdir}/SPLIT_DONE
SPLITWAITEND`
		echo "    Submitted job $jid which will wait for fastq file splitjobs to finish. START: $(date)"
		while true; do
			[[ -f ${splitdir}/SPLIT_DONE ]] && break
			sleep 1m
		done
		echo "              job $jid waiting for fastq file split to finish is done. END:   $(date)"
        echo "    Created read1 split files:"
        ls -lh ${read1} | sed 's/^/      /'
        else
            cp -rs ${fastqdir} ${splitdir}
            wait
        fi
    else
        ## No need to re-split fastqs if they already exist
        echo -e "---  Using already created files in $splitdir\n"
	# unzipped files will have .fastq extension, softlinked gz 
        testname=$(ls -l ${splitdir} | awk '$9~/fastq$/||$9~/gz$/{print $9; exit}')

        if [ ${testname: -3} == ".gz" ]
        then
            read1=${splitdir}"/*${read1str}*.fastq.gz"
        else
	    read1=${splitdir}"/*${read1str}*.fastq"
        fi
    fi
    
    ## Launch job. Once split/move is done, set the parameters for the launch. 
    echo "(-: Starting job to launch other jobs once splitting is complete"
    
    ## Loop over all read1 fastq files and create jobs for aligning read1,
    ## aligning read2, and merging the two. Keep track of merge names for final
    ## merge. When merge jobs successfully finish, can launch final merge job.
    ## ARRAY holds the names of the jobs as they are submitted
    ## Loop over all read1 fastq files and create jobs for aligning read1,
    ## aligning read2, and merging the two. Keep track of merge names for final
    ## merge. When merge jobs successfully finish, can launch final merge job. 
    countjobs=0
    declare -a ARRAY
    declare -a JIDS
    declare -a TOUCH

    dependmerge="afterok"

    for i in ${read1}
    do
	ext=${i#*$read1str}
	name=${i%$read1str*} 
	# these names have to be right or it'll break
	name1=${name}${read1str}
	name2=${name}${read2str}	
	jname=$(basename "$name")${ext}
    usegzip=0
    if [ "${ext: -3}" == ".gz" ]
    then
        usegzip=1
    fi

	# count ligations
    ### JOB _count_ligation
	jid=`sbatch <<- CNTLIG
		#!/bin/bash -l
		#SBATCH -p $queue
		#SBATCH -t $((queue_time / 3))
		#SBATCH -c 2
		#SBATCH -o $debugdir/count_ligation-%j.out
		#SBATCH -e $debugdir/count_ligation-%j.err
		#SBATCH -J "${groupname}${jname}_Count_Ligation"
		date
		export usegzip=${usegzip}; 
		export name=${name}; 
		export name1=${name1}; 
		export name2=${name2}; 
		export ext=${ext}; 
		export ligation=${ligation}; 
		${juiceDir}/scripts/countligations.sh
		date
CNTLIG`
    echo "Submitted ligation counting job ${groupname}${jname}_Count_Ligation"

	# align read1 fastq
	touchfile1=${tmpdir}/${jname}1

    ### JOB _align1
	jid=`sbatch <<- ALGNR1 | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $queue
		#SBATCH -t ${queue_time}
		#SBATCH -o $debugdir/align1-%j.out
		#SBATCH -e $debugdir/align1-%j.err
		#SBATCH -c ${threads}
		#SBATCH --mem=$bwa_mem
		#SBATCH -J "${groupname}_align1_${jname}"
		$load_bwa
		# Align read1
		date
		if [ -n "$shortread" ] || [ "$shortreadend" -eq 1 ]
		then
			echo 'Running command bwa aln $threadstring -q 15 $refSeq $name1$ext > $name1$ext.sai && bwa samse $refSeq $name1$ext.sai $name1$ext > $name1$ext.sam'
			bwa aln $threadstring -q 15 $refSeq $name1$ext > $name1$ext.sai \
                && bwa samse $refSeq $name1$ext.sai $name1$ext > $name1$ext.sam
			if [ \\\$? -ne 0 ]
			then
				touch $errorfile
				exit 1
			else
				touch $touchfile1
				echo "(-: Short align of $name1$ext.sam done successfully"
    			fi
		else
			echo 'Running command bwa mem $threadstring $refSeq $name1$ext > $name1$ext.sam '
			bwa mem $threadstring $refSeq $name1$ext > $name1$ext.sam
			if [ \\\$? -ne 0 ]
			then  
				touch $errorfile
				exit 1
			else
				touch $touchfile1
				echo "(-: Mem align of $name1$ext.sam done successfully"
			fi
		fi
		date
ALGNR1`
    echo "Submitted align1 job ${groupname}_align1_${jname}"

	dependalign="afterok:$jid"
	
    ### JOB _align2
	# align read2 fastq
	touchfile2=${tmpdir}/${jname}2
	jid=`sbatch <<- ALGNR2 | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $queue
		#SBATCH -t ${queue_time}
		#SBATCH -o $debugdir/align2-%j.out
		#SBATCH -e $debugdir/align2-%j.err
		#SBATCH -c ${threads}
		#SBATCH --mem=$bwa_mem
		#SBATCH -J "${groupname}_align2_${jname}"
		$load_bwa
		date
		# Align read2
		if [ -n "$shortread" ] || [ "$shortreadend" -eq 2 ]
		then		
			echo 'Running command bwa aln $threadstring -q 15 $refSeq $name2$ext > $name2$ext.sai && bwa samse $refSeq $name2$ext.sai $name2$ext > $name2$ext.sam '
			 $threadstribwa alnng -q 15 $refSeq $name2$ext > $name2$ext.sai \
                && bwa samse $refSeq $name2$ext.sai $name2$ext > $name2$ext.sam
			if [ \\\$? -ne 0 ]
			then 
				touch $errorfile
				exit 1
			else
				touch $touchfile2
				echo "(-: Short align of $name2$ext.sam done successfully"
			fi
		else	
			echo 'Running command bwa mem $threadstring $refSeq $name2$ext > $name2$ext.sam'
			bwa mem $threadstring $refSeq $name2$ext > $name2$ext.sam
			if [ \\\$? -ne 0 ]
			then 
				touch $errorfile
				exit 1
			else
				touch $touchfile2
				echo "(-: Mem align of $name2$ext.sam done successfully"
			fi		
		fi
		date		
ALGNR2`
    echo "Submitted align2 job ${groupname}_align2_${jname}"

	dependalign="$dependalign:$jid"

	touchfile3=${tmpdir}/${jname}3
    ### JOB _merge
	# wait for top two, merge
	jid=`sbatch <<- MRGALL | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $long_queue
		#SBATCH -t ${long_queue_time}
		#SBATCH -o $debugdir/merge-%j.out
		#SBATCH -e $debugdir/merge-%j.err
		#SBATCH --mem=40G
		#SBATCH -c 8 
		#SBATCH -d $dependalign
		#SBATCH -J "${groupname}_merge_${jname}"
		#SBATCH --gres=lscratch:800
		export LC_COLLATE=C
        ${load_coreutils}
		date
		if [ ! -f "${touchfile1}" ] || [ ! -f "${touchfile2}" ]
		then
			echo "***! Error, cluster did not finish aligning ${jname}"
			touch $errorfile 
			exit 1
		fi
		# sort read 1 aligned file by readname
		sort --parallel=8 -S 30G -T /lscratch/\\\${SLURM_JOB_ID} -k1,1 $name1$ext.sam \
            > $name1${ext}_sort.sam
		if [ \\\$? -ne 0 ]
		then 
			echo "***! Error while sorting $name1$ext.sam"
			touch $errorfile
			exit 1
		else
			echo "(-: Sort read 1 aligned file by readname completed."
		fi
		
		# sort read 2 aligned file by readname 
		sort --parallel=8 -S 30G -T /lscratch/\\\${SLURM_JOB_ID} -k1,1 $name2$ext.sam \
            > $name2${ext}_sort.sam
		if [ \\\$? -ne 0 ]
		then
			echo "***! Error while sorting $name2$ext.sam"
			touch $errorfile
			exit 1
		else
			echo "(-: Sort read 2 aligned file by readname completed."
		fi
		
		# remove header, add read end indicator toreadname
		awk 'NF >= 11{\\\$1 = \\\$1"/1";print}' ${name1}${ext}_sort.sam > ${name1}${ext}_sort1.sam
		awk 'NF >= 11{\\\$1 = \\\$1"/2";print}' ${name2}${ext}_sort.sam > ${name2}${ext}_sort1.sam

		# merge the two sorted read end files
		sort --parallel=8 -S 30G -T /lscratch/\\\${SLURM_JOB_ID} -k1,1 \
            -m $name1${ext}_sort1.sam $name2${ext}_sort1.sam > $name$ext.sam
		if [ \\\$? -ne 0 ]
		then
			echo "***! Failure during merge of read files"
			touch $errorfile
			exit 1
		else
			echo "$name$ext.sam created successfully."
		fi

		# call chimeric_blacklist.awk to deal with chimeric reads; sorted file is sorted by read name at this point
		touch ${name}${ext}_abnorm.sam ${name}${ext}_unmapped.sam
		awk -v "fname1"=${name}${ext}_norm.txt -v "fname2"=${name}${ext}_abnorm.sam -v "fname3"=${name}${ext}_unmapped.sam -f $juiceDir/scripts/chimeric_blacklist.awk ${name}${ext}.sam

		if [ \\\$? -ne 0 ] 
		then    
			echo "***! Failure during chimera handling of $name${ext}"
			touch $errorfile
			exit 1   
		fi  
		# if any normal reads were written, find what fragment they 
		# correspond to and store that
		# check if site file exists and if so write the fragment number
		# even if nofrag set
		# one is not obligated to provide a site file if nofrag set; 
		# but if one does, frag numbers will be calculated correctly
		if [ -e "$name${ext}_norm.txt" ] && [ "$site" != "none" ] && [ -e "$site_file" ]
		then
			${juiceDir}/scripts/fragment.pl ${name}${ext}_norm.txt ${name}${ext}.frag.txt $site_file
		elif [ "$site" == "none" ] || [ "$nofrag" -eq 1 ]
		then
			awk '{printf("%s %s %s %d %s %s %s %d", \\\$1, \\\$2, \\\$3, 0, \\\$4, \\\$5, \\\$6, 1); for (i=7; i<=NF; i++) {printf(" %s",\\\$i);}printf("\n");}' $name${ext}_norm.txt > $name${ext}.frag.txt
		else
			echo "***! No $name${ext}_norm.txt file created"
			touch $errorfile
			exit 1
		fi
		if [ \$? -ne 0 ]
		then
			echo "***! Failure during fragment assignment of $name${ext}"
			touch $errorfile
			exit 1 
		fi
		# sort by chromosome, fragment, strand, and position
		sort -S 2G -T $tmpdir -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $name${ext}.frag.txt > $name${ext}.sort.txt
		if [ \\\$? -ne 0 ]   
		then
			echo "***! Failure during sort of $name${ext}"
			touch $errorfile
			exit 1
		else
			rm $name${ext}_norm.txt $name${ext}.frag.txt
			rm $name1$ext.sa* $name2$ext.sa* $name1${ext}_sort*.sam $name2${ext}_sort*.sam
		fi
		touch $touchfile3
		date
MRGALL`
    echo "Submitted merge job ${groupname}_merge_${jname}"

	dependmerge="${dependmerge}:${jid}"
	ARRAY[countjobs]="${groupname}_merge_${jname}"
	JIDS[countjobs]="${jid}"
	TOUCH[countjobs]="$touchfile3"
        countjobs=$(( $countjobs + 1 ))
    done # done looping over all fastq split files
    # note that the dependmerge concatentation failed sometimes in UGER due to
    # the string becoming too long, so we had the below code instead.
    # list of all jobs.  hold next steps until they have finished
#    for (( i=0; i < countjobs; i++ ))
#    do
#        if [ $i -eq 0 ]; then
#            holdjobs="-hold_jid ${ARRAY[i]}"
#        else
#            holdjobs="${holdjobs},${ARRAY[i]}"
#        fi
#    done
    
    # list of all jobs. print errors if failed    
    for (( i=0; i < $countjobs; i++ ))
    do
	f=${TOUCH[$i]}
	msg="***! Error in job ${ARRAY[$i]}  Type squeue -j ${JIDS[$i]} to see what happened"
	
    ### JOB _check
	# check that alignment finished successfully
	jid=`sbatch <<- EOF
		#!/bin/bash
		#SBATCH -p $queue
		#SBATCH -t ${queue_time}
		#SBATCH -o $debugdir/aligncheck-%j.out
		#SBATCH -e $debugdir/aligncheck-%j.err
		#SBATCH -J "${groupname}_check"
		#SBATCH -d $dependmerge
		date
		echo "Checking $f"
		if [ ! -e $f ]
		then
			echo $msg
			touch $errorfile
		fi
		date
EOF`
    echo "Submitted check job ${groupname}_check"

	jid=$(echo $jid | egrep -o -e "\b[0-9]+$")
	dependmergecheck="${dependmerge}:${jid}"
    done
fi  # Not in merge, dedup,  or final stage, i.e. need to split and align files.

# Not in final, dedup, or postproc
if [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
    if [ -z $merge ]
    then
	sbatch_wait="#SBATCH -d $dependmergecheck"
    else
        sbatch_wait=""
    fi
    
    # merge the sorted files into one giant file that is also sorted.      jid=`sbatch <<- MRGSRT | egrep -o -e "\b[0-9]+$"
    
    ### JOB _fragmerge
    jid=`sbatch <<- EOF
		#!/bin/bash
		#SBATCH -o $debugdir/fragmerge-%j.out
		#SBATCH -e $debugdir/fragmerge-%j.err
		#SBATCH --mem 246G
		#SBATCH -p $long_queue
		#SBATCH -t ${long_queue_time}
		#SBATCH -c 8
		#SBATCH -J "${groupname}_fragmerge"
		#SBATCH --gres=lscratch:800
		${sbatch_wait}
		date
		if [ -f "${errorfile}" ]
		then
			echo "***! Found errorfile. Exiting." 
			exit 1 
		fi
		export LC_COLLATE=C
        ${load_coreutils}

        if ! sort --parallel=8 -S220G -T /lscratch/\\\${SLURM_JOB_ID} \
            -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n \
            $splitdir/*.sort.txt > $outputdir/merged_sort.txt
        then
            echo "***! Some problems occurred somewhere in creating sorted align files."
            touch $errorfile
            exit 1
        else
            echo "(-: Finished sorting all sorted files into a single merge."
            echo "now grep out chrM"
            mv $outputdir/merged_sort.txt $outputdir/merged_sort_chrM.txt
            grep -v "chrM" $outputdir/merged_sort_chrM.txt > $outputdir/merged_sort.txt
            wc -l $outputdir/merged_sort_chrM.txt
            wc -l $outputdir/merged_sort.txt
            echo "chrM removed from merged_sort.txt"
        fi
		date
EOF`
    echo "Submitted fragmerge job ${groupname}_fragmerge"

    jid=$(echo $jid | egrep -o -e "\b[0-9]+$")
    dependmrgsrt="afterok:$jid"
fi

# Remove the duplicates from the big sorted file
if [ -z $final ] && [ -z $postproc ]
then
    if [ -z $dedup ]
    then
        sbatch_wait="#SBATCH -d $dependmrgsrt"
    else
        sbatch_wait=""
    fi
    # Guard job for dedup. this job is a placeholder to hold any job submitted after dedup.
    # We keep the ID of this guard, so we can later alter dependencies of inner dedupping phase.
    # After dedup is done, this job will be released. 
    ### JOB _dedup_guard
    guardjid=`sbatch <<- DEDUPGUARD | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $debugdir/dedupguard-%j.out
	#SBATCH -e $debugdir/dedupguard-%j.err
	#SBATCH -t 5
	#SBATCH -c 2
	#SBATCH -H
	#SBATCH -J "${groupname}_dedup_guard"
	${sbatch_wait}
	date
DEDUPGUARD`
    echo "Submitted dedup_guard job in held state ${groupname}_dedup_guard"

    dependguard="afterok:$guardjid"

    ### JOB _dedup
    # if jobs succeeded, kill the cleanup job, remove the duplicates from the big sorted file
    jid=`sbatch <<- DEDUP | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -t $((queue_time / 2))
	#SBATCH --mem-per-cpu=2G
	#SBATCH -o $debugdir/dedup-%j.out
	#SBATCH -e $debugdir/dedup-%j.err
	#SBATCH -c 2
	#SBATCH -J "${groupname}_dedup"
	${sbatch_wait}
	date
    export LC_ALL=C
        if [ -f "${errorfile}" ]
        then 
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi 
	squeue -u $USER -o "%A %T %j %E %R" | column -t
	awk -v queue=$long_queue -v groupname=$groupname -v debugdir=$debugdir -v dir=$outputdir -v topDir=$topDir -v juicedir=$juiceDir -v site=$site -v genomeID=$genomeID -v genomePath=$genomePath -v user=$USER -v guardjid=$guardjid -f $juiceDir/scripts/split_rmdups.awk $outputdir/merged_sort.txt
	##Schedule new job to run after last dedup part:
	##Push guard to run after last dedup is completed:
	##srun --ntasks=1 -c 1 -p "$queue" -t 1 -o ${debugdir}/dedup_requeue-%j.out -e ${debugdir}/dedup-requeue-%j.err -J "$groupname_msplit0" -d singleton echo ID: $ echo "\${!SLURM_JOB_ID}"; scontrol update JobID=$guardjid dependency=afterok:\$SLURM_JOB_ID
	squeue -u $USER -o "%A %T %j %E %R" | column -t
	date
	scontrol release $guardjid
DEDUP`
    echo "Submitted dedup job ${groupname}_dedup"

    dependosplit="afterok:$jid"

    #Push dedup guard to run only after dedup is complete:
    scontrol update JobID=$guardjid dependency=afterok:$jid

    ### JOB _post_dedup
    #Wait for all parts of split_rmdups to complete:
    jid=`sbatch <<- MSPLITWAIT | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $debugdir/post_dedup-%j.out
	#SBATCH -e $debugdir/post_dedup-%j.err
	#SBATCH -t 100
	#SBATCH -c 2
	#SBATCH -J "${groupname}_post_dedup"
	#SBATCH -d ${dependguard}
	date
	rm -Rf $tmpdir ;
	find $debugdir -type f -size 0 | xargs rm
	squeue -u $USER -o "%A %T %j %E %R" | column -t
	date
MSPLITWAIT`
    echo "Submitted post_dedup job ${groupname}_post_dedup"

    dependmsplit="afterok:$jid"
    sbatch_wait="#SBATCH -d $dependmsplit"
else
    sbatch_wait=""
fi

if [ -z "$genomePath" ]
then
    #If no path to genome is give, use genome ID as default.
    genomePath=$genomeID
fi

# if early exit, we stop here, once the merged_nodups.txt file is created.
if [ -z "$earlyexit" ]
then
    ### JOB _dupcheck
    # Check that dedupping worked properly
    # in ideal world, we would check this in split_rmdups and not remove before we know they are correct
    jid=`sbatch <<- DUPCHECK | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -t $((queue_time / 3))
	#SBATCH -o $debugdir/dupcheck-%j.out
	#SBATCH -e $debugdir/dupcheck-%j.err
	#SBATCH -c 2
	#SBATCH --mem-per-cpu=1G
	#SBATCH -J "${groupname}_dupcheck"
	${sbatch_wait}

	date      
	ls -l ${outputdir}/merged_sort.txt | awk '{printf("%s ", \\\$5)}' > $debugdir/dupcheck-${groupname}
	ls -l ${outputdir}/merged_nodups.txt ${outputdir}/dups.txt ${outputdir}/opt_dups.txt | awk '{sum = sum + \\\$5}END{print sum}' >> $debugdir/dupcheck-${groupname}
	awk 'NR==1{if (NF == 2 && \\\$1 == \\\$2){print "Sorted and dups/no dups files add up"}else{print "Problem" >> ${errorfile}; print "***! Error! The sorted file and dups/no dups files do not add up, or were empty."}}' $debugdir/dupcheck-${groupname} 
        date                                                                                                           
DUPCHECK`
    echo "Submitted dupcheck job ${groupname}_dupcheck"
    sbatch_wait="#SBATCH -d afterok:$jid"

    #Skip if post-processing only is required
    if [ -z $postproc ]
    then
    ### JOB _stats
	jid=`sbatch <<- STATS | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $queue
		#SBATCH -t ${queue_time}
		#SBATCH -o $debugdir/stats-%j.out
		#SBATCH -e $debugdir/stats-%j.err
		#SBATCH -c 2
		#SBATCH --mem-per-cpu=10G
		#SBATCH -J "${groupname}_stats"
		${sbatch_wait}

		date
		if [ -f "${errorfile}" ]
		then 
			echo "***! Found errorfile. Exiting." 
			exit 1 
		fi 
		${load_java}
		export _JAVA_OPTIONS="-Xms2048m -Xmx16384m -XX:ParallelGCThreads=1"
		tail -n1 $headfile | awk '{printf"%-1000s\n", \\\$0}' > $outputdir/inter.txt 

		${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/stats_dups.txt $outputdir/dups.txt
		cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter.txt
		java -cp ${juiceDir}/scripts/ LibraryComplexity $outputdir inter.txt >> $outputdir/inter.txt
		${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter.txt -q 1 $outputdir/merged_nodups.txt

		tail -n1 $headfile | awk '{printf"%-1000s\n", \\\$0}' > $outputdir/inter_30.txt 

		cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter_30.txt
		java -cp ${juiceDir}/scripts/ LibraryComplexity $outputdir inter_30.txt >> $outputdir/inter_30.txt
		${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter_30.txt -q 30 $outputdir/merged_nodups.txt

		cat $splitdir/*_abnorm.sam > $outputdir/abnormal.sam
		cat $splitdir/*_unmapped.sam > $outputdir/unmapped.sam
		awk -f ${juiceDir}/scripts/collisions.awk $outputdir/abnormal.sam > $outputdir/collisions.txt
		date
STATS`
    echo "Submitted stats job ${groupname}_stats"

	dependstats="afterok:$jid"

    ### JOB _hic
	jid=`sbatch <<- HIC | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -t ${long_queue_time}
	#SBATCH -o $debugdir/hic-%j.out
	#SBATCH -e $debugdir/hic-%j.err	
	#SBATCH -c 2
	#SBATCH --mem-per-cpu=32G
	#SBATCH -J "${groupname}_hic"
	#SBATCH -d $dependstats
	${load_java}
	export _JAVA_OPTIONS="-Xms2048m -Xmx48192m -XX:ParallelGCThreads=1"
	date
	if [ -f "${errorfile}" ]
	then 
		echo "***! Found errorfile. Exiting." 
		exit 1 
	fi 
	if [ "$nofrag" -eq 1 ]
	then 
		${juiceDir}/scripts/juicer_tools48g pre -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath
	else
		${juiceDir}/scripts/juicer_tools48g pre -f $site_file -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath
fi
	date
HIC`
    echo "Submitted hic job ${groupname}_hic"

	dependhic="afterok:$jid"

    ### JOB _hic30
	jid=`sbatch <<- HIC30 | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -t ${long_queue_time}
	#SBATCH -o $debugdir/hic30-%j.out
	#SBATCH -e $debugdir/hic30-%j.err
	#SBATCH -c 2
	#SBATCH --mem-per-cpu=32G
	#SBATCH -J "${groupname}_hic30"
	#SBATCH -d ${dependstats}
	${load_java}
	export _JAVA_OPTIONS="-Xms2048m -Xmx48192m -XX:ParallelGCThreads=1"
	date
        if [ -f "${errorfile}" ]
        then 
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi 
        if [ "$nofrag" -eq 1 ]
        then 
	    ${juiceDir}/scripts/juicer_tools48g pre -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath
	else
	    ${juiceDir}/scripts/juicer_tools48g pre -f $site_file -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath
	fi
	date
HIC30`
    echo "Submitted hic job ${groupname}_hic30"

	dependhic30="${dependhic}:$jid"
	fi

	if [ -z $postproc ]
	then
		sbatch_wait="#SBATCH -d $dependhic30"
	else
		sbatch_wait=""
	fi
    ### JOB _hiccups_wrap
	jid=`sbatch <<- HICCUPS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p gpu
	#SBATCH -t ${queue_time}
	#SBATCH -c 2
	#SBATCH --mem=18g
	#SBATCH --gres=gpu:k80:1
	#SBATCH -o $debugdir/hiccups_wrap-%j.out
	#SBATCH -e $debugdir/hiccups_wrap-%j.err
	#SBATCH -J "${groupname}_hiccups_wrap"
	${sbatch_wait}
	${load_gpu}
	echo "load: $load_gpu"
	${load_java}
	date
	nvcc -V
        if [ -f "${errorfile}" ]
        then 
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi 
	${juiceDir}/scripts/juicer_hiccups.sh -j ${juiceDir}/scripts/juicer_tools -i $outputdir/inter_30.hic -m ${juiceDir}/references/motif -g $genomeID
	date
HICCUPS`
    echo "Submitted hiccups_wrap job ${groupname}_hiccups_wrap"
    dependhiccups="afterok:$jid"

    ### JOB _arrowhead_wrap
	jid=`sbatch <<- ARROWS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -t ${queue_time}
	#SBATCH -c 2
	#SBATCH --mem=18g
	#SBATCH -o $debugdir/arrowhead_wrap-%j.out
	#SBATCH -e $debugdir/arrowhead_wrap-%j.err
	#SBATCH -J "${groupname}_arrowhead_wrap"
	${sbatch_wait}
	${load_java}
	date
        if [ -f "${errorfile}" ]
        then 
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi 
	${juiceDir}/scripts/juicer_arrowhead.sh -j ${juiceDir}/scripts/juicer_tools -i $outputdir/inter_30.hic
	date;
ARROWS`
    echo "Submitted arrowhead job ${groupname}_arrowhead_wrap"
        dependarrows="${dependhiccups}:$jid"

    ### JOB _fincln
	jid=`sbatch <<- FINCLN1 | egrep -o -e "\b[0-9]+$"
	#!/bin/bash
	#SBATCH -p $queue
	#SBATCH -t ${queue_time}
	#SBATCH --mem-per-cpu=2G
	#SBATCH -o $debugdir/fincln-%j.out
	#SBATCH -e $debugdir/fincln-%j.err
	#SBATCH -c 2
	#SBATCH -J "${groupname}_prep_done"
	#SBATCH -d $dependarrows
	date
	export splitdir=${splitdir}; export outputdir=${outputdir}; ${juiceDir}/scripts/check.sh
	date
FINCLN1`
    echo "Submitted fincln job ${groupname}_prep_done"
else

    ### JOB _fincln
	jid=`sbatch <<- FINCLN1 | egrep -o -e "\b[0-9]+$" 
	#!/bin/bash
	#SBATCH -p $queue
	#SBATCH -t ${queue_time}
	#SBATCH --mem-per-cpu=2G
	#SBATCH -o $debugdir/fincln1-%j.out
	#SBATCH -e $debugdir/fincln1-%j.err
	#SBATCH -c 2
	#SBATCH -J "${groupname}_prep_done"     
	#SBATCH -d $dependarrows
	date
	export splitdir=${splitdir}; export outputdir=${outputdir}; export early=1; ${juiceDir}/scripts/check.sh
	date
FINCLN1`
    echo "Submitted fincln job ${groupname}_prep_done"
	dependprepdone="afterok:$jid"

fi
echo "(-: Finished adding all jobs... Now is a good time to get that cup of coffee.."
