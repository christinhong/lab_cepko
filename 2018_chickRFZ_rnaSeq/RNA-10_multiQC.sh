# Christin M. Hong
# Last updated: 2018-08
# Harvard Medical School, Connie Cepko Lab


# Script 3 for differential expression analysis of chick RNA-seq data.
    # Running bash on HMS O2 cluster (Slurm).
    # Tasks
        # Check for successful completion of STAR mapping


#### INFRASTRUCTURE ####

# Request resources on interactive node
srun --pty -p interactive -c 8 --mem=64G -t 0-11:00 /bin/bash


# GLOBAL VARIABLES

export cepko=/n/data2/hms/genetics/cepko
export christin=${cepko}/christin


# Job-specific
export intCores=8
export pathLogs=/home/ch220/jobLogs


# Experiment-specific
export pathProj=/home/ch220/2018_chickRFZ_rnaSeq # Project's working directory, set for script with "SBATCH -D"
export pathData=${pathProj}/data/*/*
export pathBash=${pathProj}/bash
export pathDoc=${pathProj}/doc

export pathOut=/n/scratch2/ch220
export pathData2=${pathOut}/fq_trimmed


# References
export pathRef=/home/ch220/resources/ref
export fileGenome=${pathRef}/genomes/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa


# Tools and general scripts with versions in their respective paths 
export pathTools=/home/ch220/resources/tools
export parallel=${pathTools}/parallel-20180722/src/parallel


# Modules loaded from O2
    # Need to load prereq modules first
    # Check for prereqs and command syntax with "module spider <tool name>"

module load gcc/6.2.0 python/2.7.12
module load trimmomatic/0.36
module load fastqc/0.11.3
module load multiqc/1.5


# Other options
LC_COLLATE=C    # specifies sort order (numbers, uppercase, then lowercase)


# Notes
    # Export: Variables are inherited by child processes (e.g. subshells from GNU parallel).
    # Bash variables are untyped by default.  For more robust code, can declare data type with [declare] (see http://tldp.org/LDP/abs/html/declareref.html ), but I'm not sure how declare works with export.  May try later.
    # When possible, using full path to minimize confusion by shell, record tool versions, and increase clarity regarding dependencies.


#### START ####

# Check exit codes - should all be 0:0
sacct -j <jobid>

# Move Trimmomatic output (job error logs) for multiQC
	# Inelegant solution...would be better to output the error logs to doc directly, or find them within the jobLogs folder automatically.  But this'll do for now.
mv ~/jobLogs/*.err ~/2018_chickRFZ_rnaSeq/doc/trimmomatic/


# Run FastQC on trimmed FASTQ files
mkdir /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_trimmed

tmux

${parallel} -j ${intCores} --verbose --joblog "${pathLogs}/parallel-fastQC_trimmed.log" --resume-failed --keep-order "fastqc -o /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_trimmed {}" ::: /n/scratch2/ch220/fq_trimmed/*/*_paired_trimmed.fq.gz


# Compile reports with MultiQC and check output
cd /home/ch220/2018_chickRFZ_rnaSeq/doc

# Pre-trimmomatic
multiqc . --ignore fastqc_trimmed/ --ignore trimmomatic/ -o multiQC -n multiQC-01_origFastQC

# Move old output to archive folder
mkdir /home/ch220/2018_chickRFZ_rnaSeq/doc/archive

mv -n fastqc_nextSeq archive/

mv -n fastqc_hiSeq archive/

# Post-trimmomatic
multiqc \
    /home/ch220/2018_chickRFZ_rnaSeq/doc \
    --ignore /home/ch220/2018_chickRFZ_rnaSeq/doc/archive \
    --ignore /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_orig \
    -o /home/ch220/2018_chickRFZ_rnaSeq/doc/multiQC \
    -n multiQC-02_afterTrimmomatic_$(date '+%Y-%m-%d')



# Looks good!  Prep for STAR mapping
mkdir /n/scratch2/ch220/starMap


#### After 2-pass STAR mapping ####

# MultiQC after STAR takes a couple hours. Definitely run within tmux.

sacct -j 21496183            # Exit values should be 0:0

cat ~/jobLogs/RNA-05_*.err  # Should all be empty

# MultiQC reports: It's easier to see the details from each step if the reports aren't all concatenated together. Maybe break up as:

    # fastQC
    # trimmomatic-fastQC
    # starPass1
    # starPass2


multiqc \
    /home/ch220/2018_chickRFZ_rnaSeq/doc \
    /n/scratch2/ch220/starMap/ \
    --ignore /n/scratch2/ch220/starMap/pass2 \
    --ignore /home/ch220/2018_chickRFZ_rnaSeq/doc/archive \
    --ignore /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_orig \
    -o /home/ch220/2018_chickRFZ_rnaSeq/doc/multiQC \
    -n multiQC-02_STARpass1_$(date '+%Y-%m-%d')


multiqc \
    /n/scratch2/ch220/starMap/pass2 \
    -o /home/ch220/2018_chickRFZ_rnaSeq/doc/multiQC \
    -n multiQC-03_STARpass2_$(date '+%Y-%m-%d')



#### Downloading data from home terminal ####
rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/home/ch220/2018_chickRFZ_rnaSeq/doc/multiQC/*" "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_chickRFZ_rnaSeq/doc/"

