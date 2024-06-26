#!/bin/bash

#SBATCH -p short                # Required, partition/queue for submission
#SBATCH -t 0-10:00              # Required, time allowed for script to run in D-HH:MM
#SBATCH -c 4                    # number of cores/cpus per node (32 c/node)
#SBATCH --mem=48G               # total RAM requested per job (256 GB RAM/node)

#SBATCH -e /home/ch220/jobLogs/RNA-05_%A-%a.err     # standard err
#SBATCH -o /home/ch220/jobLogs/RNA-05_%A-%a.out     # standard out
#SBATCH --open-mode=append                          # append adds to outfile, truncate deletes old outfile first

	

#### INTRO ####

# Christin M. Hong
# Last modified: 2018-10
# Harvard Medical School, Connie Cepko Lab

# Script for differential expression analysis of murine RNA-seq data.
    # Decided to keep flexibility of array command by leaving it outside this file. Then can choose each time which values to run, e.g. "sbatch --array=1-50 <script.sh>" for NextSeq samples, or "sbatch --array=51-60 <script.sh>" for HiSeq samples.


# Tasks
    # Second pass STAR mapping of reads to m38 with STAR via job array (1 job per sample)


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
export fileGTF=${pathGen}/Mus_musculus.GRCm38.94.gtf
export pathStarInd=${pathGen}/STAR-m38


# Tools and general scripts with versions noted in their respective paths
export pathTools=${cepko}/resources/tools
export qualimap=${pathTools}/qualimap_v2.2.1/qualimap

export parallel=${pathTools}/parallel-20180722/src/parallel
    # Note on GNU Parallel: --keep-order requires --joblog, and if the joblog file already exists, Parallel with terminate without an error notice. If using keep-order, it's best to tag the joblog with the job ID to generate a unique file per submission, e.g. with ${SLURM_ARRAY_JOB_ID}.



# Modules loaded from O2
    # Need to load prereq modules first. Check for prereqs and command syntax with "module spider <tool name>".
module load gcc/6.2.0 python/2.7.12
module load fastqc/0.11.3
module load multiqc/1.5
module load trimmomatic/0.36    # Problematic syntax. Version number will need to be manually updated in script if program is updated.

module load star/2.5.4a
module load picard/2.8.0        # Problematic syntax. Version number will need to be manually updated in script if program is updated.

module load samtools/1.9
module load R/3.5.1             # Used by Picard CollectMultipleMetrics and Qualimap




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
cd ${pathData}/S*."${SLURM_ARRAY_TASK_ID}"


# Check working directory
echo "Working directory is $PWD"
echo



#### Multi-sample 2-pass mapping with STAR ####

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

    outPre="${outPre3}"_star2_
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

# For $(cat ${fileSJ}): Explanation of subshells vs. command substition (which acts like variables) to expand TXT of splice junctions from https://unix.stackexchange.com/questions/213530/difference-of-using-and-to-execute-a-series-of-commands


echo
echo "Done!"
echo


# Due to the splice junction insertion step, for 96 samples, this pass takes ~3-4 hours to run.  It mainly depends on the queue - each job itself takes ~30 minutes.

