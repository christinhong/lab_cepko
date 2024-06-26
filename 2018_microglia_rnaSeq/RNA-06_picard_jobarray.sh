#!/bin/bash

#SBATCH -p short                # Required, partition/queue for submission
#SBATCH -t 0-10:00              # Required, time allowed for script to run in D-HH:MM
#SBATCH -c 4                    # number of cores/cpus per node (32 c/node)
#SBATCH --mem=48G               # total RAM requested per job (256 GB RAM/node)

#SBATCH -e /home/ch220/jobLogs/RNA-06_%A-%a.err     # standard err
#SBATCH -o /home/ch220/jobLogs/RNA-06_%A-%a.out     # standard out
#SBATCH --open-mode=append                          # append adds to outfile, truncate deletes old outfile first

	

#### INTRO ####

# Christin M. Hong
# Last modified: 2018-10
# Harvard Medical School, Connie Cepko Lab

# Script for differential expression analysis of murine RNA-seq data.
    # Decided to keep flexibility of array command by leaving it outside this file. Then can choose each time which values to run, e.g. "sbatch --array=1-50 <script.sh>" for NextSeq samples, or "sbatch --array=51-60 <script.sh>" for HiSeq samples.


# Tasks
    # BAM annotation with Picard



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

echo 'Adding read groups to BAMs with Picard'
echo

cd ${pathOutStar2}/S"${SLURM_ARRAY_TASK_ID}"


echo "Working directory is ${PWD}"
echo


# Split filename on underscores and pass to array.
    # RGID: Unique read group identifier
    # RGSM: Sample/tissue that read group came from
    # RGPU: Platform, e.g. {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}
    # RGPL: Tech used to produce reads. Valid values: ILLUMINA, SOLID, LS454, HELICOS and PACBIO.
    # RGLB: DNA library prep identifier. I'm also using this as a batch identifier.
    
# Note that sorting by coordinate requires enough memory to load the reference genome (>30 GB).


echo "Cleaning BAMs and adding read groups from filename with Picard"
echo

find . *_Aligned.out.bam -name ".*" -prune -o -print > bams.txt

${parallel} -j ${intCores} \
    --resume-failed --keep-order \
    --joblog "${pathLogs}/S"${SLURM_ARRAY_TASK_ID}"-parallel_picardAnno_${SLURM_ARRAY_JOB_ID}.log" \
    'sh ${pathBash}/RNA_picard.sh' \
    :::: bams.txt


echo "Done annotating BAMs!"



#### Check for successful read group assignment in first file ####
# samtools view -H ${pathData3}/*_p2rg.bam | grep '@RG'

# Check for completion of script
# grep -rl "Done annotating BAMs!" ${pathLogs}/RNA-06*.out | wc -l # Should be the number of samples (95)

