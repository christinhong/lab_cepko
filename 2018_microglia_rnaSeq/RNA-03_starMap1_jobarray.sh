#!/bin/bash

#SBATCH -p short                # Required, partition/queue for submission
#SBATCH -t 0-10:00              # Required, time allowed for script to run in D-HH:MM
#SBATCH -c 4                    # number of cores/cpus per node (32 c/node)
#SBATCH --mem=48G               # total RAM requested per job (256 GB RAM/node)

#SBATCH -e /home/ch220/jobLogs/RNA-03_%A-%a.err     # standard err
#SBATCH -o /home/ch220/jobLogs/RNA-03_%A-%a.out     # standard out
#SBATCH --open-mode=append                          # append adds to outfile, truncate deletes old outfile first
#SBATCH --mail-type=END                             # email when job ends
#SBATCH --mail-user=christinhong@g.harvard.edu      # address for email
	

#### INTRO ####

# Christin Hong
# Last modified: 2018-10
# Harvard Medical School, Connie Cepko Lab

# Script for differential expression analysis of murine RNA-seq data.
    # Decided to keep flexibility of array command by leaving it outside this file. Then can choose each time which values to run, e.g. "sbatch --array=1-50 <script.sh>" for NextSeq samples, or "sbatch --array=51-60 <script.sh>" for HiSeq samples.


# Tasks
    # First pass STAR mapping of reads to m38 with STAR via job array (1 job per sample)


#### INFRASTRUCTURE ####

# Stop script if error occurs
set -Eeuo pipefail		# See https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/ and http://redsymbol.net/articles/unofficial-bash-strict-mode/


# GLOBAL VARIABLES

# Job-specific
export intCores=4
export pathLogs=/home/ch220/jobLogs


# Experiment-specific
export cepko=/n/data2/hms/genetics/cepko
export christin=${cepko}/christin

export pathProj=${christin}/2018_microglia_rnaSeq
export pathData=${pathProj}/data
export pathBash=${pathProj}/bash
export pathDoc=${pathProj}/doc

export pathOut=${pathProj}/output
export pathOutStar=${pathOut}/starMap
export pathOutStar2=${pathOut}/starMap/pass2
export fileSJ=${pathOutStar}/sj.txt                # File of STAR splice junctions

export pathData3=${pathOut}/bamsMerged


# References
export pathRef=${cepko}/resources/ref
export pathGen=${pathRef}/genome_m38
export fileGen=${pathGen}/Mus_musculus.GRCm38.dna.primary_assembly.fa
export pathStarInd=${pathGen}/STAR-m38


# Tools and general scripts with versions noted in their respective paths
export pathTools=${cepko}/resources/tools
export parallel=${pathTools}/parallel-20180722/src/parallel


# Modules loaded from O2
    # Need to load prereq modules first. Check for prereqs and command syntax with "module spider <tool name>".
module load gcc/6.2.0 python/2.7.12
module load fastqc/0.11.3
module load multiqc/1.5
module load star/2.5.4a
module load picard/2.8.0


# Other options
LC_COLLATE=C    # specifies sort order (numbers, uppercase, then lowercase)


# Notes
    # Export: Variables are inherited by child processes (e.g. subshells from GNU parallel).
    # Bash variables are untyped by default.  For more robust code, can declare data type with [declare] (see http://tldp.org/LDP/abs/html/declareref.html ), but I'm not sure how declare works with export.  May try later.
    # When possible, using full path to minimize confusion by shell, record tool versions, and increase clarity regarding dependencies.




#### START ####

echo "Starting first pass of STAR mapping on $(date '+%Y-%m-%d %H:%M:%S')"
echo "Using STAR index ${pathStarInd}"
echo


#### Aligning reads by STAR 2-pass mapping ####

# Move to folder FASTQ files
cd ${pathData}/S*."${SLURM_ARRAY_TASK_ID}"


# Check working directory
echo "Working directory is $PWD"
echo


# Get list of lanes, excluding hidden files/directories
echo "Finding *.fastq.gz"
echo

find . *.fastq.gz -name ".*" -prune -o -print > fq.txt

#

echo "Matching expected lane ID regex _L([0-9]{3})_, e.g. _L001_"
echo

while IFS='' read -r fq || [[ -n "$fq" ]]; do
    [[ $fq =~ _L([0-9]{3})_ ]] && echo ${BASH_REMATCH[0]}
done < fq.txt | sort | uniq > S"${SLURM_ARRAY_TASK_ID}"_lanes.txt



# Check output
echo "List of lane identifiers:"
cat S"${SLURM_ARRAY_TASK_ID}"_lanes.txt
echo


# Get read files per lane
     # Adapted from https://stackoverflow.com/questions/10929453/read-a-file-line-by-line-assigning-the-value-to-a-variable and https://stackoverflow.com/questions/25908070/how-to-get-the-output-of-find-as-a-space-separated-string

while IFS='' read -r lane || [[ -n "$lane" ]]; do
    find . *$lane*.fastq.gz  -name ".*" -prune -o -print | paste -sd " "
done <  S"${SLURM_ARRAY_TASK_ID}"_lanes.txt > S"${SLURM_ARRAY_TASK_ID}"_readPairs.txt


# Check output
echo "List of reads per lane:"
cat S"${SLURM_ARRAY_TASK_ID}"_readPairs.txt
echo



#### First pass mapping with STAR ####
	# Since STAR tends to make large intermediate files, will make output in scratch2

# Make directory for STAR output
mkdir ${pathOutStar}/S"${SLURM_ARRAY_TASK_ID}"


# Run STAR for each read pair
while IFS='' read -r pair || [[ -n "$pair" ]]; do

	# Check input
    echo
    echo "STAR input: $pair"
    echo

    # Name output
    outPre1=$(echo ${pair} | cut -d " " -f1)
    outPre2=${outPre1%_*}
    outPre3=${outPre2%_*}

    outPre="${outPre3}"_
    echo "Output bam is prefixed with ${outPre}"
    echo


    # Run STAR
    STAR \
        --runThreadN ${intCores} \
        --genomeDir ${pathStarInd} \
        --readFilesCommand zcat \
        --readFilesIn ${pair} \
        --genomeLoad NoSharedMemory \
        --outFileNamePrefix ${pathOutStar}/S"${SLURM_ARRAY_TASK_ID}"/${outPre} \
        --outSAMtype BAM Unsorted \
        --sjdbScore 2

done < S"${SLURM_ARRAY_TASK_ID}"_readPairs.txt

echo
echo "Done!"
echo

