# Christin M. Hong
# Last updated: 2018-08
# Harvard Medical School, Connie Cepko Lab


# Script 3 for differential expression analysis of chick RNA-seq data.
    # Running bash on HMS O2 cluster (Slurm).
    # Tasks
        # Check for successful completion of STAR mapping


#### INFRASTRUCTURE ####

# Request resources on interactive node
srun --pty -p interactive -c 4 --mem=36G -t 0-11:00 /bin/bash


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

export pathData3=${pathOut}/bamsMerged


# References
export pathRef=${cepko}/resources/ref
export pathGen=${pathRef}/genomes
export fileGen=${pathGen}/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa
export fileGTF=${pathGen}/Gallus_gallus.Gallus_gallus-5.0.93.gtf
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
module load picard/2.8.0
module load samtools/1.3.1
module load R/3.5.1             # Used by Picard CollectMultipleMetrics and Qualimap


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

${parallel} -j ${intCores} --verbose --joblog "${pathLogs}/parallel-fastQC_trimmed.log" --resume-failed --keep-order "fastqc -o ${pathDoc}/fastqc_trimmed {}" ::: ${pathData2}/*/*_paired_trimmed.fq.gz


# Compile reports with MultiQC and check output
cd ${pathDoc}

# Pre-trimmomatic
multiqc \
    ${pathDoc}/fastqc_orig \
    -o ${pathDoc}/multiQC \
    -n multiQC-01_fastqc-orig_$(date '+%Y-%m-%d')

    #multiqc . --ignore fastqc_trimmed/ --ignore trimmomatic/ -o multiQC -n multiQC-01_origFastQC


# Move old output to archive folder
mkdir /home/ch220/2018_chickRFZ_rnaSeq/doc/archive

mv -n fastqc_nextSeq archive/

mv -n fastqc_hiSeq archive/

# Post-trimmomatic
multiqc \
    ${pathDoc}/trimmomatic \
    ${pathDoc}/fastqc_trimmed \
    -o ${pathDoc}/multiQC \
    -n multiQC-02_fastqc-trimmed_$(date '+%Y-%m-%d')


# Looks good!  Prep for STAR mapping
mkdir ${pathOut}/starMap


#### After STAR mapping ####

# MultiQC after STAR takes a couple hours. Definitely run within tmux.

sacct -j 22694194           # Exit values should be 0:0

cat ~/jobLogs/RNA-05_*.err  # Should all be empty

# MultiQC reports: It's easier to see the details from each step if the reports aren't all concatenated together. Maybe break up as:

    # fastQC-orig
    # trimmomatic-fastQC
    # starPass1
    # starPass2
    # picard-qualimap


multiqc \
    ${pathOutStar} \
    --ignore ${pathOutStar2} \
    -o ${pathDoc}/multiQC \
    -n multiQC-03_STARpass1_$(date '+%Y-%m-%d')


multiqc \
    ${pathOutStar2} \
    -o ${pathDoc}/multiQC \
    -n multiQC-04_STARpass2_$(date '+%Y-%m-%d')



#### Checking for rRNA ####
# ribokmers.fa.gz is from the silva rRNA database at https://www.arb-silva.de/. See https://www.biostars.org/p/159959/#175854

module load java

# For separating rRNA from FASTQ
/home/ch220/resources/tools/bbmap/bbduk.sh in=RFZ-1-A01_S1_L001_R1_001.fastq.gz outm=ribo.fa outu=nonribo.fa k=31 ref=/home/ch220/resources/ref/ribokmers.fa.gz

# For simply getting rRNA contamination stats
/home/ch220/resources/tools/bbmap/bbduk.sh in=RFZ-1-A01_S1_L001_R2_001.fastq.gz k=31 ref=/home/ch220/resources/ref/ribokmers.fa.gz


# Interesting! Doesn't seem like it's rRNA after all...in that case, I have no idea what happened with these libraries?

### Nevermind, I BLASTed the actual sequences from the "Overrepresented sequences" portion of the FastQC report, and these DO look like rRNA. "PREDICTED: Gallus gallus 18S ribosomal RNA."

# Well, good to know that FastQC + BLAST are always reliable.



#### BAM QC ####
sacct -j 22823785


multiqc \
    ${pathDoc}/bamQC \
    -o ${pathDoc}/multiQC \
    -n multiQC-05_bamQC_$(date '+%Y-%m-%d')




#### Downloading data from home terminal ####
rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/christin/2018_chickRFZ_rnaSeq/doc/multiQC/*" "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_chickRFZ_rnaSeq/doc/multiQC/"

rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/christin/2018_chickRFZ_rnaSeq/doc/bamQC" "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_chickRFZ_rnaSeq/doc/"


# rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/n/scratch2/ch220/bamsForJiho/sorted/*" "/home/christin/Desktop/bams/"

# rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_trimmed/*.html" "/home/christin/Desktop/bams/fastqc_trimmed_output/"



# Upload
rsync -avr --progress "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/resources/ref/*" "ch220@transfer.rc.hms.harvard.edu:/home/ch220/resources/ref/"

rsync -av --progress "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_chickRFZ_rnaSeq/" "ch220@transfer.rc.hms.harvard.edu:/home/ch220/2018_chickRFZ_rnaSeq"

#

rsync -av --progress "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_chickRFZ_rnaSeq/" "ch220@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/christin/2018_chickRFZ_rnaSeq"


# For files from a range of folders
rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/n/scratch2/ch220/starMap/pass2/S{51..60}/*_Aligned.out.bam" "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_chickRFZ_rnaSeq/doc/BAMs/"

or

rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/n/scratch2/ch220/starMap/pass2/S{51..60}/*_Aligned.out.bam" "/n/scratch2/ch220/bamsForJiho/"


rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/home/ch220/2018_chickRFZ_rnaSeq/data/*" "/n/data2/hms/genetics/cepko/christin/2018_chickRFZ_rnaSeq/data/"

