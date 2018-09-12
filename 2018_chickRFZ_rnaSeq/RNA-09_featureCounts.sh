# Christin M. Hong
# Last updated: 2018-08
# Harvard Medical School, Connie Cepko Lab


# Script 9 for differential expression analysis of chick RNA-seq data.
    # Running bash on HMS O2 cluster (Slurm).


#### INFRASTRUCTURE ####

# Request resources on interactive node
srun --pty -p interactive -c 4 --mem=48G -t 0-11:00 /bin/bash


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
module load subread/1.6.2


# Other options
LC_COLLATE=C    # specifies sort order (numbers, uppercase, then lowercase)


# Notes
    # Export: Variables are inherited by child processes (e.g. subshells from GNU parallel).
    # Bash variables are untyped by default.  For more robust code, can declare data type with [declare] (see http://tldp.org/LDP/abs/html/declareref.html ), but I'm not sure how declare works with export.  May try later.
    # When possible, using full path to minimize confusion by shell, record tool versions, and increase clarity regarding dependencies.



#### START ####

cd ${pathData3}

tmux

mkdir ${pathDoc}/featureCounts

find ./*/*_p5mDups.bam -name ".*" -prune -o -print > bamsFinal.txt

printf "Running featureCounts for summarizing multiple paired end datasets. 
Counting at gene (meta-feature) level. 
Counting reads mapping to EXONS only.
Counting duplicates. 
Allowing reads to overlap multiple metafeatures. 
Not counting chimeric fragments or multimapping reads."

featureCounts \
    -p \
    -T ${intCores} \
    -t exon \
    -g gene_id \
    -a ${fileGTF} \
    -G ${fileGen} \
    -C \
    -O \
    --byReadGroup \
    -o ${pathDoc}/featureCounts/featureCounts.txt \
    $(cat bamsFinal.txt) \
    > ${pathDoc}/featureCounts/featureCounts.out \
    2> ${pathDoc}/featureCounts/featureCounts.err



echo "Running MultiQC for BAM QC and featureCounts"
multiqc \
    ${pathDoc}/bamQC \
    ${pathDoc}/featureCounts \
    -o ${pathDoc}/multiQC \
    -n multiQC-05_bamQC-featureCounts_$(date '+%Y-%m-%d')


echo "Done with featureCounts!  Ready to move to R."



: << "Comment"
featureCount explanations:

* Counting the number of reads that align to a gene (meta-feature) rather than counting by alignment to an exon.

* -t exon: Decided to count only reads mapping to exons to minimize ambiguity. From the Qualimap BAM QC data, expect 60-80% of reads to be counted.

* -g gene_id: Identifier in the GTF being used for the genes/meta-feature.

* -C: Don't count chimeric fragments = fragments that span multiple chromosomes. I can see chimeric fragments being interesting in cancer research, but I don't see any reason they'd be present here.

* -O: Count fragments that overlap multiple features, e.g. fragments that map to more than one gene. Allowed since I can see thes potentially capturing genes that sit close to each other on the genome.

* --byReadGroup: Counting by read group annotation (see Picard's AddOrReplaceReadGroups above) so I can more easily analyze for differences between them later on. 

* To maximize accuracy, not counting multimapping reads = reads that map to more than one location.

* NOTE: If there are issues while running, can add a --verbose option for easier debugging.




On featureCounts (http://bioinf.wehi.edu.au/featureCounts/):
A read is said to overlap a feature if at least one read base is found to overlap the feature. For paired-end data, a fragment (or template) is said to overlap a feature if any of the two reads from that fragment is found to overlap the feature.

By default, featureCounts does not count reads overlapping with more than one feature (or more than one meta-feature when summarizing at meta-feature level). Users can use the -O option to instruct featureCounts to count such reads (they will be assigned to all their overlapping features or meta-features).

Note that, when counting at the meta-feature level, reads that overlap multiple features of the same meta-feature are always counted exactly once for that meta-feature, provided there is no overlap with any other meta-feature. For example, an exon-spanning read will be counted only once for the corresponding gene even if it overlaps with more than one exon. 



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

Comment
