#!/bin/bash

#SBATCH -p short                # Required, partition/queue for submission
#SBATCH -t 0-10:00              # Required, time allowed for script to run in D-HH:MM
#SBATCH -c 4                    # number of cores/cpus per node
#SBATCH --mem=48G               # total RAM requested per job

#SBATCH -D /home/ch220/2018_chickRFZ_rnaSeq         # set working directory
#SBATCH --open-mode=append                          # append adds to outfile, truncate deletes old outfile first

#SBATCH -e /home/ch220/jobLogs/RNA-02_%A-%a.err     # standard err
#SBATCH -o /home/ch220/jobLogs/RNA-02_%A-%a.out     # standard out
#SBATCH --mail-type=END                             # email when job ends
#SBATCH --mail-user=christinhong@g.harvard.edu      # address for email
	

#####

# Christin Hong
# Last modified: 2018-08
# Harvard Medical School, Connie Cepko Lab

# Script for differential expression analysis of chick RNA-seq data. See project README.
    # Decided to keep flexibility of array command by leaving it outside this file. Then can choose each time which values to run, e.g. "sbatch --array=1-50 <script.sh>" for NextSeq samples, or "sbatch --array=51-60 <script.sh>" for HiSeq samples.


# Tasks
    # Prep paired input by lane
    # Gently clean reads with Trimmomatic


#### INFRASTRUCTURE ####

# Stop script if error occurs
set -Eeuo pipefail		# See https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/ and http://redsymbol.net/articles/unofficial-bash-strict-mode/


# GLOBAL VARIABLES

# Job-specific
export intCores=4
export pathLogs=/home/ch220/jobLogs


# Experiment-specific
export pathProj=/home/ch220/2018_chickRFZ_rnaSeq
export pathData=${pathProj}/data/*/*
export pathBash=${pathProj}/bash
export pathDoc=${pathProj}/doc


export pathOut=/n/scratch2/ch220
export pathData2=${pathOut}/fq_trimmed
export pathStarInd=${pathOut}/STAR-gg5-sjo49            # STAR indexed with --sjbdOverhang 49

export pathOutStar=${pathOut}/starMap                   
export fileSJ=${pathOutStar}/sj_1-60.txt                # File of STAR splice junctions


# References
export pathRef=/home/ch220/resources/ref
export fileGenome=${pathRef}/genomes/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa


# Tools and general scripts with versions in their respective paths 
export pathTools=/home/ch220/resources/tools
export parallel=${pathTools}/parallel-20180722/src/parallel


# Modules loaded from O2
    # Need to load prereq modules first. Check for prereqs and command syntax with "module spider <tool name>"
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

echo "Starting Trimmomatic on $(date '+%Y-%m-%d %H:%M:%S')"

# Add newline to output for readability
echo


# Move to data folder
cd ${pathData}/S*."${SLURM_ARRAY_TASK_ID}"


# Check working directory
echo "Working directory is $PWD"
echo


#### Prep paired end input and output ####
# Goal: Preserve paired end and lane information. To do that, want to match R1 and R2 on the lane identifier to feed them in as paired input.
    # It isn't elegant, but I think the code below is robust to missing partners (e.g. if an R2 is removed due to failing QC) and different lane numbers (though it assumes an "_L###_" format for lane IDs).


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
    find . *$lane*  -name ".*" -prune -o -print | paste -sd " "
done <  S"${SLURM_ARRAY_TASK_ID}"_lanes.txt > S"${SLURM_ARRAY_TASK_ID}"_readPairs.txt

     # Note: Could match on prefix with lane ID instead of lane id only, like I do for naming the output BAM files. But this way is less dependent on exactly where the underscores are and more specific for lanes.



# Check output
echo "List of reads per lane:"
cat S"${SLURM_ARRAY_TASK_ID}"_readPairs.txt
echo



#### Run Trimmomatic ####
# I haven't seen much justification for aggressive trimming, so I prefer to do gentle trimming to retain as much information as possible. Mainly for removing adapter sequences and low-quality bases (< Q5).

# In this case, from FastQC, T-7-G04_S34_L008_R1.fastq.gz looks like it's contaminated with "TruSeq Adapter, Index 9 (100% over 50bp)". Make reference file with contaminating adapter(s).
echo ">truseq9" > adapters.fa
echo "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGC" >> adapters.fa


# Large output files. Generate in scratch to avoid "No space left on disk" error.
mkdir ${pathOutTrimmed}/S"${SLURM_ARRAY_TASK_ID}"

echo "Called subscript written for trimmomatic/0.36"
echo
echo "Gently cleaning with Trimmomatic for paired end reads. Removing adapter sequences noted by FastQC. Trimming leading and trailing bases equal or below quality 3 (Q3), and trimming sequences once a 4 bp window falls below an average of Q5. "keepBothReads" set to TRUE to improve compatibility with STAR."

${parallel} -j ${intCores} --joblog "${pathLogs}/S"${SLURM_ARRAY_TASK_ID}"_parallel-trimmomatic.log" \
    --resume-failed --keep-order \
    'sh ${pathBash}/RNA_trimmomatic-pe.sh' \
    :::: S"${SLURM_ARRAY_TASK_ID}"_readPairs.txt


echo
echo "Done!"
echo
