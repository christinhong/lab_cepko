# Christin M. Hong
# Last updated: 2018-08
# Harvard Medical School, Connie Cepko Lab


# Script 9 for differential expression analysis of chick RNA-seq data.
    # Running bash on HMS O2 cluster (Slurm).


#### INFRASTRUCTURE ####

# Request resources on interactive node
srun --pty -p interactive -c 4 --mem=48G -t 0-11:00 /bin/bash


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
export pathBamQC=${pathDoc}/bamQC


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
module load trimmomatic/0.36    # Problematic syntax. Version number will need to be manually updated in scripts when program is updated.

module load star/2.5.4a
module load picard/2.8.0        # Problematic syntax. Version number will need to be manually updated in scripts when program is updated.

module load samtools/1.9
module load R/3.5.1             # Used by Picard CollectMultipleMetrics and Qualimap
module load subread/1.6.2



# Other options
LC_COLLATE=C    # specifies sort order (numbers, uppercase, then lowercase)


# Notes
    # Export: Variables are inherited by child processes (e.g. subshells from GNU parallel).
    # Bash variables are untyped by default.  For more robust code, can declare data type with [declare] (see http://tldp.org/LDP/abs/html/declareref.html ), but I'm not sure how declare works with export.  May try later.
    # When possible, using full path to minimize confusion by shell, record tool versions, and increase clarity regarding dependencies.



#### START ####

# Check exit values of RNA-08 (bam sorting, merging, and QC collection): sacct -j 27618595


# Run featureCounts (could job array/parallelize this, but it runs fast enough already)

printf "Running featureCounts for summarizing multiple paired end datasets. 
Counting at gene (meta-feature) level. 
Counting reads mapping to EXONS only.
Counting duplicates. 
Allowing reads to overlap multiple metafeatures. 
Not counting chimeric fragments or multimapping reads."


cd ${pathData3}

tmux


# Counts for Sean

find ./Sample_S{001..057}*/*_p5mDups.bam -name ".*" -prune -o -print > bamsFinal_microglia.txt


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
    -o ${pathDoc}/featureCounts/featureCounts_microglia.txt \
    $(cat bamsFinal_microglia.txt) \
    > ${pathDoc}/featureCounts/featureCounts_microglia.out \
    2> ${pathDoc}/featureCounts/featureCounts_microglia.err



# Counts for Sawyer

find ./Sample_S{058..095}*/*_p5mDups.bam -name ".*" -prune -o -print > bamsFinal_cones.txt

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
    -o ${pathDoc}/featureCounts/featureCounts_cones.txt \
    $(cat bamsFinal_cones.txt) \
    > ${pathDoc}/featureCounts/featureCounts_cones.out \
    2> ${pathDoc}/featureCounts/featureCounts_cones.err



echo
echo "Done with featureCounts!  Ready to move to R."
echo


: << "COMMENT"
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

COMMENT




#### Aggregate QC reports with MultiQC ####

# Note: By default, MultiQC's output is in whatever directory it runs in, and its input is its directory + all subdirectories.


echo
echo "Running MultiQC for STAR mapping"
echo 

# Microglia
multiqc \
    ${pathOutStar}/S{1..57} \
    --ignore ${pathOutStar2} \
    -o ${pathDoc}/multiQC \
    -n multiQC-microglia-03_STARpass1_$(date '+%Y-%m-%d')


multiqc \
    ${pathOutStar2}/S{1..57} \
    -o ${pathDoc}/multiQC \
    -n multiQC-microglia-04_STARpass2_$(date '+%Y-%m-%d')




# Cones
multiqc \
    ${pathOutStar}/S{58..95} \
    --ignore ${pathOutStar2} \
    -o ${pathDoc}/multiQC \
    -n multiQC-cones-03_STARpass1_$(date '+%Y-%m-%d')

multiqc \
    ${pathOutStar2}/S{58..95} \
    -o ${pathDoc}/multiQC \
    -n multiQC-cones-04_STARpass2_$(date '+%Y-%m-%d')



echo
echo "Running MultiQC for Picard, Qualimap, and featureCounts"
echo

# Microglia
multiqc \
    ${pathDoc}/bamQC/S{001..057}* \
    ${pathDoc}/featureCounts/*microglia* \
    -o ${pathDoc}/multiQC \
    -n multiQC-microglia-05_bamQC-featureCounts_$(date '+%Y-%m-%d')


# Cones
multiqc \
    ${pathDoc}/bamQC/S{058..095}* \
    ${pathDoc}/featureCounts/*cones* \
    -o ${pathDoc}/multiQC \
    -n multiQC-cones-05_bamQC-featureCounts_$(date '+%Y-%m-%d')





#### Download data from home terminal ####
# rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/christin/2018_microglia_rnaSeq/doc/multiQC/*" "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_microglia_rnaSeq/doc/multiQC/"


echo "Done with MultiQC!"
echo
