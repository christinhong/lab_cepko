#!/bin/bash

#SBATCH -p short                # Required, partition/queue for submission
#SBATCH -t 0-10:00              # Required, time allowed for script to run in D-HH:MM
#SBATCH -c 4                    # number of cores/cpus per node (32 c/node)
#SBATCH --mem=48G               # total RAM requested per job (256 GB RAM/node)

#SBATCH -e /home/ch220/jobLogs/RNA-05_%A-%a.err     # standard err
#SBATCH -o /home/ch220/jobLogs/RNA-05_%A-%a.out     # standard out
#SBATCH --open-mode=append                          # append adds to outfile, truncate deletes old outfile first
#SBATCH --mail-type=END                             # email when job ends
#SBATCH --mail-user=christinhong@g.harvard.edu      # address for email
	

#### INTRO ####

# Christin Hong
# Last modified: 2018-08
# Harvard Medical School, Connie Cepko Lab

# Script for differential expression analysis of chick RNA-seq data. See project README.
    # Decided to keep flexibility of array command by leaving it outside this file. Then can choose each time which values to run, e.g. "sbatch --array=1-50 <script.sh>" for NextSeq samples, or "sbatch --array=51-60 <script.sh>" for HiSeq samples.


# Tasks
    # Second pass mapping of reads to Galgal5 with STAR via job array (1 job per sample)


#### INFRASTRUCTURE ####

# Stop script if error occurs
set -Eeuo pipefail		# See https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/ and http://redsymbol.net/articles/unofficial-bash-strict-mode/


# GLOBAL VARIABLES

export cepko=/n/data2/hms/genetics/cepko
export christin=${cepko}/christin


# Job-specific
export intCores=4
export pathLogs=/home/ch220/jobLogs


# Experiment-specific
export pathProj=${christin}/2018_chickRFZ_rnaSeq
export pathData=${pathProj}/data/*/*
export pathBash=${pathProj}/bash
export pathDoc=${pathProj}/doc

export pathOut=${pathProj}/output
export pathData2=${pathOut}/fq_trimmed
export pathOutStar=${pathOut}/starMap
export pathOutStar2=${pathOut}/starMap/pass2
export fileSJ=${pathOutStar}/sj_1-60.txt                # File of STAR splice junctions


# References
export pathRef=${cepko}/resources/ref
export pathGen=${pathRef}/genomes
export fileGen=${pathGen}/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa
export pathStarInd=${pathGen}/STAR-gg5-sjo49            # STAR indexed with --sjbdOverhang 49


# Tools and general scripts with versions in their respective paths 
export pathTools=${cepko}/resources/tools
export parallel=${pathTools}/parallel-20180722/src/parallel


# Modules loaded from O2
    # Need to load prereq modules first. Check for prereqs and command syntax with "module spider <tool name>".
module load gcc/6.2.0 python/2.7.12
module load fastqc/0.11.3
module load multiqc/1.5
module load trimmomatic/0.36    # Problematic syntax. Script will need to be manually updated if program is updated.
module load star/2.5.4a


# Other options
LC_COLLATE=C    # specifies sort order (numbers, uppercase, then lowercase)


# Notes
    # Export: Variables are inherited by child processes (e.g. subshells from GNU parallel).
    # Bash variables are untyped by default.  For more robust code, can declare data type with [declare] (see http://tldp.org/LDP/abs/html/declareref.html ), but I'm not sure how declare works with export.  May try later.
    # When possible, using full path to minimize confusion by shell, record tool versions, and increase clarity regarding dependencies.



#### START ####

echo "Starting second pass of multi-sample 2-pass STAR mapping on $(date '+%Y-%m-%d %H:%M:%S')"
echo "Using STAR index ${pathStarInd}"
echo


# Move to folder with trimmed FASTQ reads
cd ${pathData2}/S"${SLURM_ARRAY_TASK_ID}"


# Check working directory
echo "Working directory is $PWD"
echo



#### Multi-sample 2-pass mapping with STAR ####
    # Since STAR tends to make large intermediate files, will make output in scratch2

# Make directory for STAR output
mkdir ${pathOutStar2}/S"${SLURM_ARRAY_TASK_ID}"


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
    outPre4=${outPre3%_*}

    outPre="${outPre4}"_star2_
    echo "Output bam is prefixed with ${outPre}"
    echo

    # Run STAR
    STAR \
        --runThreadN ${intCores} \
        --genomeDir ${pathStarInd} \
        --readFilesCommand zcat \
        --readFilesIn ${pair} \
        --genomeLoad NoSharedMemory \
        --outFileNamePrefix ${pathOutStar2}/S"${SLURM_ARRAY_TASK_ID}"/${outPre} \
        --outSAMtype BAM Unsorted \
        --sjdbScore 2 \
        --sjdbFileChrStartEnd $(cat ${fileSJ}) \
        --quantMode TranscriptomeSAM GeneCounts

done < S"${SLURM_ARRAY_TASK_ID}"_readPairs.txt

# For --sjdbFileChrStartEnd $(cat ${pathOutStar}/sj_1-50.txt): Explanation of subshells vs. command substition (which acts like variables) to expand TXT of splice junctions from https://unix.stackexchange.com/questions/213530/difference-of-using-and-to-execute-a-series-of-commands


echo
echo "Done!"
echo

