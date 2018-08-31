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

# Testing samtools mpileup for extracting regions within genes of interest that may be useful for RNA-FISH

# Checking to see which contigs are available in file
samtools view -H RFZ-6-F01_p5mDups.bam | grep "@SQ" | head

# Ah, chromosomes have only numbers, no "chr" prefix.  OK.

samtools tview RFZ-6-F01_p5mDups.bam ${fileGen}

g

6:22,703,402

# Hmm...well, it's kind of a mess to view alignment that way.  No wonder IGV is so popular.


# samtools mpileup output: "Each line consists of 
# 1. chromosome
# 2. 1-based coordinate
# 3. reference base
# 4. the number of reads covering the site
# 5. read bases and
# 6. base qualities." -http://samtools.sourceforge.net/pileup.shtml

# So maybe I can extract the genomic locations, then filter for >10 on column 4 (number of reads covering the site), and get the reference base sequence from there?


# FGF8
samtools mpileup -f ${fileGen} -r 6:22,703,402-22,715,085 RFZ-6-F01_p5mDups.bam -o test.txt

# From IGV, this is mainly between 6:22,711,000-22,713,100
samtools tview -p 6:22,711,000-22,713,100 RFZ-6-F01_p5mDups.bam ${fileGen}

samtools mpileup -f ${fileGen} -r 6:22,711,000-22,713,100 RFZ-6-F01_p5mDups.bam -o test.txt



# Make list of bams of interest
cd ${pathData3}

ls */*-6-*_p5mDups.bam > pileup_bams.txt

ls */*-7-*_p5mDups.bam >> pileup_bams.txt

#### Getting coordinates for each gene from IGV (easiest way)
# Genes of interest from Connie:
    # FGF8 (high in RFZ, T6, and D6): chr6:22,703,402-22,715,085
    # VAX1 (high in V and T6): chr6:29,104,869-29,112,316
    # EMB (no clear pattern): chrZ:14,734,701-14,763,276
    # CYP26A1 (high in V, some in RFZ6): chr6:20,553,756-20,560,848
    # CYP26C1 (high in RFZ): chr6:20,558,272-20,569,872
    # TBX2 (): chr19:7,639,044-7,651,850
    # TBX3 (): chr15:12,275,196-12,289,707
    # TBX5 (): chr15:12,378,768-12,423,152
    # RALDH1: Doesn't seem to be annotated in the Galgal5 genome.


# Make file of genes and their coordinates.  Remember to REMOVE the "chr" prefix.
nano pileup_genes.txt

# copy in below
fgf8 6:22,703,402-22,715,085
vax1 6:29,104,869-29,112,316
emb Z:14,734,701-14,763,276
cyp26a1 6:20,553,756-20,560,848
cyp26c1 6:20,558,272-20,569,872
tbx2 19:7,639,044-7,651,850
tbx3 15:12,275,196-12,289,707
tbx5 15:12,378,768-12,423,152


# Make output folder(s)

nReads=30

outdir=pileup${nReads}


mkdir pileup10

mkdir pileup30



# Run script below
while IFS=' ' read -ra gene || [[ -n "${gene}" ]]; do

varGene=$(echo "${gene[0]}")
geneRange=$(echo "${gene[1]}")

echo "Gene of interest is ${varGene}."
echo "Genomice coordinates from IGV are ${geneRange}."

echo "Minimum coverage per nucleotide is ${nReads}."


# Get alignments to genes of interest for each BAM of interest
while IFS='' read -r bam || [[ -n "$bam" ]]; do

    out=$(echo ${bam} | cut -d "/" -f2 | cut -d "_" -f1)
    
        echo "Filtering for nucleotide bases in gene with at least ${nReads} reads"

    samtools mpileup -f ${fileGen} -r ${geneRange} ${bam} | \
    awk '{$4 = sprintf("%05d", $4); print}' | \
    awk -v nReads="$nReads" '($4 >= nReads)' > ${outdir}/${varGene}_min${nReads}_${out}.txt


    # Checking if any bases in gene have aligned reads and acting accordingly
    if [[ -s ${outdir}/${varGene}_min${nReads}_${out}.txt ]] ; then

        echo "${varGene}_min${nReads}_${out} has aligned reads. Getting range and sequence."

        chrMin=$(cat ${outdir}/${varGene}_min${nReads}_${out}.txt | awk 'NR==1{min = $1 + 0; next} {if ($1 < min) min = $1;} END {print min}' )    
        pMin=$(cat ${outdir}/${varGene}_min${nReads}_${out}.txt | awk 'NR==1{min = $2 + 0; next} {if ($2 < min) min = $2;} END {print min}' )

        chrMax=$(cat ${outdir}/${varGene}_min${nReads}_${out}.txt | awk 'NR==1{max = $1 + 0; next} {if ($1 > max) max = $1;} END {print max}' )
        pMax=$(cat ${outdir}/${varGene}_min${nReads}_${out}.txt | awk 'NR==1{max = $2 + 0; next} {if ($2 > max) max = $2;} END {print max}')

        pRange=$(echo "$pMax - $pMin" | bc )

        printf "%s\n" "First position with at least ${nReads} aligned reads is ${chrMin}:${pMin}." \
        "Last position with at least ${nReads} aligned reads is ${chrMax}:${pMax}." \
        "Range is ${pRange} bases." \
        "Assuming range is on the first noted chromosome: ${chrMin}. Getting sequence from reference genome..." \
        " " \
        > ${outdir}/${varGene}_min${nReads}_${out}_sequence.txt

        samtools faidx ${fileGen} ${chrMin}:${pMin}-${pMax} >> ${outdir}/${varGene}_min${nReads}_${out}_sequence.txt
            # from https://davetang.org/wiki/tiki-index.php?page=SAMTools#Extracting_SAM_entries_mapping_to_a_specific_loci


    else

        echo "${varGene}_min${nReads}_${out}.txt is empty. No aligned reads to gene. Deleting file."
        rm -f ${outdir}/${varGene}_min${nReads}_${out}.txt

    fi

    echo "Done with ${bam}"
    echo

done < pileup_bams.txt

done < pileup_genes.txt



#### Downloading data from home terminal ####
rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/christin/2018_chickRFZ_rnaSeq/doc/multiQC/*" "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_chickRFZ_rnaSeq/doc/multiQC/"

rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/christin/2018_chickRFZ_rnaSeq/output/bamsMerged/pileup*" "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_chickRFZ_rnaSeq/doc/"


# Upload
rsync -avr --progress "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/resources/ref/*" "ch220@transfer.rc.hms.harvard.edu:/home/ch220/resources/ref/"


# For files from a range of folders
rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/n/scratch2/ch220/starMap/pass2/S{51..60}/*_Aligned.out.bam" "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_chickRFZ_rnaSeq/doc/BAMs/"

