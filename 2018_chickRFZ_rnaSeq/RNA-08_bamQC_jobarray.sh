#!/bin/bash

#SBATCH -p short                # Required, partition/queue for submission
#SBATCH -t 0-10:00              # Required, time allowed for script to run in D-HH:MM
#SBATCH -c 1                    # number of cores/cpus per node (32 c/node)
#SBATCH --mem=48G               # total RAM requested per job (256 GB RAM/node)

#SBATCH -e /home/ch220/jobLogs/RNA-08_%A-%a.err        # standard err
#SBATCH -o /home/ch220/jobLogs/RNA-08_%A-%a.out        # standard out
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



#### INFRASTRUCTURE ####

# Stop script if error occurs
set -Eeuo pipefail		# See https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/ and http://redsymbol.net/articles/unofficial-bash-strict-mode/


# Timing script
res1=$(date +%s)


# GLOBAL VARIABLES

export cepko=/n/data2/hms/genetics/cepko
export christin=${cepko}/christin


# Job-specific
export intCores=1
export intMem=48
export pathLogs=/home/ch220/jobLogs


# Experiment-specific
export pathProj=${christin}/2018_chickRFZ_rnaSeq
export pathData=${pathProj}/data/*/*
export pathBash=${pathProj}/bash

export pathDoc=${pathProj}/doc
export pathMulti=${pathDoc}/multiQC
export pathBamQC=${pathDoc}/bamQC

export pathOut=${pathProj}/output
export pathData2=${pathOut}/fq_trimmed

export pathOutStar=${pathOut}/starMap
export pathOutStar2=${pathOut}/starMap/pass2
export fileSJ=${pathOutStar}/sj_1-60.txt                # File of STAR splice junctions

export pathData3=${pathOut}/bamsMerged


# References
export pathRef=${cepko}/resources/ref
export pathGen=${pathRef}/genomes
export fileGen=${pathGen}/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa
export fileGTF=${pathGen}/Gallus_gallus.Gallus_gallus-5.0.93.gtf
export pathStarInd=${pathGen}/STAR-gg5-sjo49            # STAR indexed with --sjbdOverhang 49


# Tools and general scripts with versions in their respective paths 
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
module load samtools/1.3.1
module load R/3.5.1             # Used by Picard CollectMultipleMetrics and Qualimap


# Other options
LC_COLLATE=C    # specifies sort order (numbers, uppercase, then lowercase)


# Notes
    # Export: Variables are inherited by child processes (e.g. subshells from GNU parallel).
    # Bash variables are untyped by default.  For more robust code, can declare data type with [declare] (see http://tldp.org/LDP/abs/html/declareref.html ), but I'm not sure how declare works with export.  May try later.
    # When possible, using full path to minimize confusion by shell, record tool versions, and increase clarity regarding dependencies.


#### START ####

cd ${pathData3}/Sample_*."${SLURM_ARRAY_TASK_ID}"

echo "Working directory is ${PWD}"
echo

sBam=$(find . *_p3merged.bam -name ".*" -prune -o -print | cut -d "_" -f1)


echo "Sorting ${sBam}"
echo    
java -jar ${PICARD}/picard-2.8.0.jar SortSam \
    I=${sBam}_p3merged.bam \
    O=${sBam}_p4sorted.bam \
    SORT_ORDER=coordinate
# Note that sorting by coordinate requires a large amount of memory (>20 GB). Check exit values afterwards to confirm successful completion.


echo "Marking duplicates in ${sBam}"
echo
java -jar ${PICARD}/picard-2.8.0.jar MarkDuplicates \
    I=${sBam}_p4sorted.bam \
    O=${sBam}_p5mDups.bam \
    M=${pathBamQC}/${sBam}_p5mDups_MarkDuplicatesMetrics.txt


echo "Indexing ${sBam}"
echo
java -jar ${PICARD}/picard-2.8.0.jar BuildBamIndex \
    I=${sBam}_p5mDups.bam
        
        
echo "Estimating library complexity of ${sBam}"
echo
java -jar ${PICARD}/picard-2.8.0.jar EstimateLibraryComplexity \
    I=${sBam}_p5mDups.bam \
    O=${pathBamQC}/${sBam}_p5mDups_libComplexityMetrics.txt


echo "Collecting other Picard metrics on ${sBam}"
echo
java -jar ${PICARD}/picard-2.8.0.jar CollectMultipleMetrics \
    I=${sBam}_p5mDups.bam \
    O=${pathBamQC}/${sBam}_p5mDups_multipleMetrics \
    R=${fileGen} \
    PROGRAM=null \
    PROGRAM=CollectAlignmentSummaryMetrics \
    PROGRAM=CollectInsertSizeMetrics \
    PROGRAM=QualityScoreDistribution \
    PROGRAM=MeanQualityByCycle \
    PROGRAM=CollectBaseDistributionByCycle \
    PROGRAM=CollectGcBiasMetrics \
    PROGRAM=CollectQualityYieldMetrics
        

echo "Running Qualimap BAM QC"
echo
${qualimap} bamqc \
    -bam ${sBam}_p5mDups.bam \
    -outdir ${pathBamQC}/qualimapBamQC_${sBam} \
    -gff ${fileGTF} \
    --outside-stats \
    --collect-overlap-pairs \
    --paint-chromosome-limits \
    --java-mem-size=${intMem}G


echo "Running Qualimap RNA-seq QC"
echo
${qualimap} rnaseq \
    -bam ${sBam}_p5mDups.bam \
    -outdir ${pathBamQC}/qualimapRnaQC_${sBam} \
    -gtf ${fileGTF} \
    --paired \
    --java-mem-size=${intMem}G


# Timing script.  This takes ~30-60 minutes per job.
res2=$(date +%s)
echo "Start time: $res1"
echo "Stop time:  $res2"
timeSec=$(echo "$res2 - $res1" | bc )
echo "Elapsed minutes:  $(echo "scale=3; ${timeSec} / 60" | bc )"


# EOM
echo "BAM QC done!"
echo
