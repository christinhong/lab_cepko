# Christin M. Hong
# Last updated: 2018-10
# Harvard Medical School, Connie Cepko Lab


# Script 1 for differential expression analysis of mouse RNA-seq data.
    # Running bash on HMS O2 cluster (Slurm).
    # Tasks
        # Setting up project directory on cluster
        # Importing Mus musculus reference genome assembly and annotation (m38)
        # Indexing m38 genome with STAR
        # Formatting sample names with annotations
        # Organizing samples for later job arrays
        # FastQC with MultiQC compilation of samples
        # Prep samples for STAR aligner job array



#### NOTE ####

# Hi! As RNA-seq has become ubiquitous in biology, I'm happy to provide these scripts for your personal use.  I hope you find them helpful.

# If you do, it'd mean a lot to me to be mentioned in the acknowledgments.  (No pressure, of course, it'd just be nice.)  

# Best of luck with your analysis!!  \^o^/  ---Christin



#### Intro ####

# I like to write scripts the same way I write protocols - so a complete novice to the task can pick it up and go. (Also so when I look back years later, I can actually understand what I did and why I did it.)

# This script is NOT meant for submitting as a job. This is a record of the initial commands I used to set up analysis on the cluster. These commands are all run through the interactive partition (see detail in script).

# This text is saved as .sh because I'm using an editor that does syntax highlighting (gedit), and the ".sh" tells the editor to syntax highlight for bash. I strongly recommend reading and writing scripts in editors that do syntax highlighting. 

# I've recorded full paths instead of shortcuts (e.g. /home/ch220/ instead of ~/) and excluded most interactive convenience tricks for transparency and readability. Feel free to make use of variables, tmux, etc. for greater convenience at will.

# On OS: Bash is designed for the Linux OS. I'm working from a Linux Mint OS (derived from Linux Ubuntu). If you're working on a PC through Windows Powershell, terminal commands that run on your own computer (not the cluster) may be different. Windows 10 comes with the option to partially enable an Ubuntu subsystem. It isn't perfect, but it might be easier to work from that Ubuntu subsystem instead of the default Windows command line. (Macs should work the same as Linux; they built in Linux OS support while developing.)



#### START ####

# Login with eCommons ID and password
    # Can automate connection later, see http://slowkow.com/notes/ssh-tutorial/
ssh ch220@o2.hms.harvard.edu


# Request resources on interactive node
srun --pty -p interactive -c 8 --mem=64G -t 0-11:00 /bin/bash


#### Set variables. Keep versions in their respective paths for future reference. ####

# Job-specific
export intCores=8
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


# References
export pathRef=${cepko}/resources/ref
export pathGen=${pathRef}/genome_m38
export fileGen=${pathGen}/Mus_musculus.GRCm38.dna.primary_assembly.fa
export pathStarInd=${pathGen}/STAR-m38


# Tools and general scripts with versions noted in their respective paths 
export pathTools=${cepko}/resources/tools
export parallel=${pathTools}/parallel-20180722/src/parallel


# Modules loaded from O2
    # Need to load prereq modules first. Check for prereqs and command syntax with "module spider <tool name>".
module load gcc/6.2.0 python/2.7.12
module load fastqc/0.11.3
module load multiqc/1.5
module load star/2.5.4a
module load picard/2.8.0


# Other options
LC_COLLATE=C    # specifies sort order (numbers, uppercase, then lowercase)



# Notes
    # Export: Variables are inherited by child processes (e.g. subshells from GNU parallel).
    # Bash variables are untyped by default.  For more robust code, can declare data type with [declare] (see http://tldp.org/LDP/abs/html/declareref.html ), but I'm not sure how declare works with export.  May try later.
    # When possible, using full path to minimize confusion by shell, record tool versions, and increase clarity regarding dependencies.




#### Start tmux ####
tmux

# Like screen, tmux allows a process to run in the background. Unlike screen, tmux keeps the process running even after logging out of the parent shell.
    # Will inherit loaded modules and exported variables from the parent shell (but I think they have to be loaded while no tmux child shells are open). 
    # Horizontal split with Ctrl+b, ".  Vertical split with Ctrl+b, %.
    # Switch pane with Ctrl+b, arrow key.
    # Can detach from tmux with Ctrl+b, then d. Reattach with tmux attach. 
    # List windows with tmux ls. 
    # Exit window with "exit."
    # For more, see https://lukaszwrobel.pl/blog/tmux-tutorial-split-terminal-windows-easily/




#### START ONE-TIME SETUP ####


#### Set up project directory ####

# Some bash basics
    # "cd <folder name>" to move into folder
    # "ls <folder name>" to list contents of folder
    # "mkdir <folder name>" to create a folder
    # "rm -rf <folder name>" to delete a folder and ALL contents. WARNING: Be careful with this command!!
        # "rm <file name>" to delete a file
    # "mv <old folder name> <new folder name>" to move and rename files/folders


# Christin's template RNA-seq directory structure
    # jobLogs/
    # resources/: General use analysis resources
        # ref/: Reference files
            # genomes/
        # tools/: General use software and scripts

    # projectDir/: Top level analysis scripts, including this one. Ideally includes a README with notes on project, collaborators, data collection, etc.
        # bash/: Bash subscripts
        # data/: Data that's safe to modify
            # orig/: Symlink to original data. Once data is in, can lock as READ ONLY. Best practice: Before any analysis, backup all data to an external location.
        # doc/: Reports, documents, Markdown, LaTeX, etc.
            # bamQC/
            # fastqc_orig/
            # fastqc_trimmed/
            # featureCounts/
            # graphs/
            # multiQC/
            # trimmomatic/
        # output/: Processed data and other output. May also write temp output to /n/scratch2/ch220 (10 TB of file space/user, files auto-purged after 30 days of no access).
            # bamsMerged/
            # fq_trimmed/
            # starMap/
        # R/: R subscripts



#### Install GNU Parallel ####
# GNU Parallel is great for parellelizing across cores in a node (a.k.a. across threads in a computer, semi-independent tasks).
    # Job arrays are great for parellelizing across nodes in a cluster (a.k.a. across computers in a network, completely independent tasks).

cd /n/data2/hms/genetics/cepko/resources/tools

# Install latest GNU Parallel per https://www.gnu.org/software/parallel/parallel_tutorial.html
(wget -O - pi.dk/3 || curl pi.dk/3/ || \
   fetch -o - http://pi.dk/3) | bash



#### Import reference genome assembly and annotation ####

    # Importing Ensembl Mus musculus unmasked primary assembly m38 with GTF annotation release 94.

    # For reasoning, see http://genomespot.blogspot.com/2015/06/mapping-ngs-data-which-genome-version.html
    # I'm not against using the toplevel assembly with the haplotypes as a read sink for multimappers, since they exist in the actual sample (though they may not be translated).  That said, standard practice is to use primary assembly to simplify mapping, so going with primary assembly for now.



# From http://useast.ensembl.org/Mus_musculus/Info/Index -> "Download DNA sequence (FASTA)" + additional directions for using rsync on the FTP site at https://useast.ensembl.org/info/data/ftp/rsync.html

rsync -av "rsync://ftp.ensembl.org/ensembl/pub/release-94/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz" "/n/data2/hms/genetics/cepko/resources/ref/genome_m38/"
 

# Checksum for m38: 23924 787097 Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
cd ${pathRef}/genome_m38/

sum Mus_musculus.GRCm38.dna.primary_assembly.fa.gz # Good


# For the annotation, from http://useast.ensembl.org/Mus_musculus/Info/Index -> "Download GTF or GFF3" + additional directions for using rsync on the FTP site at https://useast.ensembl.org/info/data/ftp/rsync.html

rsync -av "rsync://ftp.ensembl.org/ensembl/pub/release-94/gtf/mus_musculus/Mus_musculus.GRCm38.94.gtf.gz" "/n/data2/hms/genetics/cepko/resources/ref/genome_m38/"


# Checksum for release 94: 15062 28709 Mus_musculus.GRCm38.94.gtf.gz
sum Mus_musculus.GRCm38.94.gtf.gz # Good



#### Build STAR genome index ####
    # STAR writes large temp files that require >100 GB disk space.  That's greater than the size of my entire home folder, so build accordingly.


# Unzip files (STAR uses uncompressed files)
gzip -d Mus*.gz


# Load STAR and prep folders
cd ${pathRef}/genome_m38
mkdir STAR-m38

cd STAR-m38


# Index genome
STAR \
    --runThreadN ${intCores} \
    --runMode genomeGenerate \
    --genomeDir ${pathRef}/genome_m38/STAR-m38 \
    --genomeFastaFiles ${fileGen} \
    --sjdbGTFfile ${pathRef}/genome_m38/Mus_musculus.GRCm38.94.gtf \
    --sjdbOverhang 100

    # Do this with at least 36 GB RAM to ensure that entire genome can be loaded and read.

    # On --sjdbOverhang: STAR's manual notes that, "In case of reads of varying length, the ideal value is max(ReadLength)-1. In most cases, the default value of 100 will work as well as the ideal value."

        # I've tested changing the splice junction overhang to 31 for 32 bp reads and 49 for 50 bp reads, and I didn't see any difference, so I think leaving it at 100 is fine.  Hardcoding it here for future reference, just in case.



### Create sequence dictionary for Picard ####

echo 'Picard: Creating sequence dictionary'

java -jar ${PICARD}/picard-2.8.0.jar CreateSequenceDictionary R=${fileGen} O=${fileGen}.dict



#### DONE WITH ONE-TIME SETUP ####

# Note to self: Previous code should be in a separate "RNA-00" file.




#### Importing and annotating data files from Broad Genomics Platform ####

# Import
    # For importing the data, BGP provided a link to a 30-day get site.  The best option is to log into "transfer.rc.hms.harvard.edu" (instead of o2.hms.harvard.edu) and download it with BGP's link through a tmux window directly onto the lab server.


#### Check md5sums (Sean) ####
    # Sean Wang checked the md5sum hashes to confirm that they matched the hashes provided by BGP.
        # Note that the fastq.gz files needed to be unzipped before their hashes matched.
    
    # I've copied the code/records he sent me below:

    cd /n/data2/hms/genetics/cepko/Sean/2018_microglia_rnaSeq/data

    mkdir unzip # to store unzipped files

    cp /n/data2/hms/genetics/cepko/Sean/2018_microglia_rnaSeq/data/orig /n/data2/hms/genetics/cepko/Sean/2018_microglia_rnaSeq/data/unzip # copies zipped files

    gunzip -r /unzip # unzip files

    find -type f -exec md5sum "{}" +  >> md5sums.txt # performed in /unzip, creates txt file with md5sums

    scp -r skw18@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/Sean/2018_microglia_rnaSeq/data/unzip/orig/md5sums.txt /mnt/c/Users/seankwang/Desktop/HMS/Research/Cepko/RNAseq # from home terminal, downloads md5sums.txt


    # Compared md5sum values to MANIFEST file using Excel: for both md5sum.txt and MANIFEST, copied all contents into Excel and used Data > Text to Columns to delimit file names from md5sum values. Sorted lists alphabetically and then confirmed values equivalent using an IF formula

    rm -r /n/data2/hms/genetics/cepko/Sean/2018_microglia_rnaSeq/data/unzip # removes unzipped files



#### Data annotation ####
    # Looks to me that BGP provided both the multiplexed and demultiplexed data (and some mysterious "barcode_#.fastq.gz" files), but none of the data is annotated (the filenames are inscrutable without a key).  First things first: Annotating.

    # info_logs_logs.tar.gz file is empty! Would expect that to have the sample annotations in a machine-friendly format. Emailed Kristina to ask about it.

    # Otherwise, barcode key annotations were provided in Kristina's email. I can manually format them to make it shell-friendly and work off the existing filename plus the emailed key plus our submitted manifest.


# Final data structure
    # Paired end
    # Each sample has its own folder (so for paired-end data across 4 lanes, each folder would have 8 fastq.gz files).
    # File name syntax: group-replicate-well_sampleNumber(S##)_laneNumber(L###)_R1orR2_batchNumber.fastq.gz


# Manually created key to match barcode to sample metadata from Kristina's email (find and replace is great here.  Also manually padded with leading zero to sort to match manifest order.) and our submitted sample manifest. A little nerve-wracking to make this manually, but I checked five random rows and it looks good.


# Symlink to original data
ln -s /n/data2/hms/genetics/cepko/Sean/2018_microglia_rnaSeq/data/orig ${pathData}/


# Copy demultiplexed data from orig to data dir for renaming and organizing
rsync -av --progress ${pathData}/orig/*.unmapped.[1-2].fastq.gz ${pathData}/


# Download list of sample files from O2 cluster to local terminal
find ${pathData} -name '*.fastq.gz' -printf '%P\n' | sort > sampleFiles.txt

rsync -av --progress "ch220@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/christin/2018_microglia_rnaSeq/data/sampleFiles.txt" "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_seanWang_microglia/" 


# Made final key in R because it's better for formatting a data frame.  See "RNA-01_creatingFilenameKey.R".


# Upload to HMS O2
rsync -av --progress "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_seanWang_microglia/sampleNames.csv" "ch220@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/christin/2018_microglia_rnaSeq/data/"

rsync -av --progress "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_seanWang_microglia/sampleAnno.csv" "ch220@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/christin/2018_microglia_rnaSeq/data/"




#### Rename files ####

# From https://unix.stackexchange.com/questions/57754/how-to-rename-files-with-sed-and-csv

sed 's/^/mv -vi "/;s/,/" "/;s/$/";/' < sampleNames.csv | bash -

    # (Stackexchange.  <3 )



#### Sort into folders by sample number (zero padded for autosorting) ###

# Make folders for 96 samples (since I already know we only have 96)
for ((x=1;x<=96;x+=1)); do mkdir `printf "S%03d" $x`; done



## Make list of new files
find ${pathData} -name "*.fastq.gz" -prune -printf '%P\n' | sort > sampleFiles2.txt


echo "Making lists of files matching expected sample regex S([0-9]{3}), e.g. S001."
echo
echo "WARNING: Make sure no other fields in file name match this format!!!"
echo
# Can change 'S' prefix for this field when renaming files in earlier step.  Maybe I'll use Q next time for robustness...


while IFS='' read -r fq || [[ -n "$fq" ]]; do \
    [[ $fq =~ S([0-9]{3}) ]] && echo ${BASH_REMATCH[0]} ;
done < sampleFiles2.txt | sort | uniq > sampleRegex.txt


# Check output
echo "List of sample identifiers:"
cat sampleRegex.txt
echo


echo "\${parallel} has a value of $parallel"

# Move files by sample regex
${parallel} -j ${intCores} --verbose --joblog "${pathLogs}/parallel-sampleSorting.log" --resume-failed --keep-order "mv -vi ${pathData}/*_{}_*.fastq.gz ${pathData}/{}/" :::: sampleRegex.txt



#### Append identifier for job array ####

# Append job array index numbers (not 0 padded!!) to folder names
    # Beauty in one line. Adapted from https://stackoverflow.com/questions/3211595/renaming-files-in-a-folder-to-sequential-numbers
    
ls -d ${pathData}/S* | cat -n | while read n f; do mv -n "$f" "$f.$n"; done 


echo "Success!  Samples are now annotated and organized."



#### Run FastQC ####

${parallel} -j ${intCores} --verbose --joblog "${pathLogs}/parallel-fastqc.log" --resume-failed --keep-order "fastqc -o ${pathDoc}/fastqc_orig {}" ::: ${pathData}/S*/*.fastq.gz


# Check for non-zero exit values

    # For a 96-well plate with 8 files per well (paired end * 4 lanes), should have 768 fastq.gz files.

head ${pathLogs}/parallel-fastqc.log

cat ${pathLogs}/parallel-fastqc.log | awk -F '\t' '$7 == 0 {print $9 }' |  wc -l # 768

cat ${pathLogs}/parallel-fastqc.log | awk -F '\t' '$7 != 0 {print $9 }' |  wc -l # 1 (colnames)


# Check for FastQC files
cd ${pathDoc}/fastqc_orig

ls *fastqc.html | wc -l # 768


# Separate by project
mkdir microglia
mkdir cones

${parallel} -j ${intCores} --verbose --joblog "${pathLogs}/parallel-fastqcSort1.log" --resume-failed --keep-order "mv -vi ${pathDoc}/fastqc_orig/*_{}_* ${pathDoc}/fastqc_orig/microglia" ::: S{001..057}

${parallel} -j ${intCores} --verbose --joblog "${pathLogs}/parallel-fastqcSort2.log" --resume-failed --keep-order "mv -vi ${pathDoc}/fastqc_orig/*_{}_* ${pathDoc}/fastqc_orig/cones" ::: S{058..095}



#### Aggregate FastQC reports with MultiQC ####

# Note: By default, MultiQC's output is in whatever directory it runs in, and its input is its directory + all subdirectories.


# All
multiqc \
    ${pathDoc}/fastqc_orig \
    -o ${pathDoc}/multiQC \
    -n multiQC-01_fastqc-orig_$(date '+%Y-%m-%d')


# Microglia
multiqc \
    ${pathDoc}/fastqc_orig/microglia \
    -o ${pathDoc}/multiQC \
    -n multiQC-01_fastqc-orig-microglia_$(date '+%Y-%m-%d')


# Cones
multiqc \
    ${pathDoc}/fastqc_orig/cones \
    -o ${pathDoc}/multiQC \
    -n multiQC-01_fastqc-orig-cones_$(date '+%Y-%m-%d')



# From home terminal
rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/christin/2018_microglia_rnaSeq/doc/multiQC/*" "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_microglia_rnaSeq/doc/multiQC/"

rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/christin/2018_microglia_rnaSeq/doc/fastqc_orig/*" "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_microglia_rnaSeq/doc/fastqc_orig/"



: <<'THOUGHTS'
First off, no need for trimmomatic.  I'm guessing that the processing pipeline used by Broad Genomics Platform trims itself, because all the reads are 38 bp and have Q30+ scores (except R2 reads off the negative control).

GC content plot looks strange...  Going to rerun MultiQC on them separately to look more closely.  


The microglia data is S001-S057, while Sawyer's cone data is S058-S095.  (S096 is the negative control.)



Microglia:

%GC outside 45-55%:
    Mouse 64 (P20 fractalkine): S041 (microglia), S008 (retinal background): 39%
    Mouse 89 (P70 GFP): S054: 40%

%Dups for Mouse 64 and 89 is also abnormally high - ~50% instead of around 20-30% for other samples.

Overreprented Sequences, >3%:
    Mouse 64 (P20 fractalkine): S041 (microglia), S008 (retinal background).
    Mouse 89 (P70 GFP): S054.


Other samples that MultiQC suggests checking for GC content:
    S014 Mouse 3 A: 45%
    S035 Mouse 48 A: 47%
    S013 Mouse 2 A: 45-46%
    S039 Mouse 55 A: 46%
    S043 Mouse 70 A: 46% (slightly elevated dups, ~35%)
    S045 Mouse 72 A: 47% (slightly elevated dups, ~36%)
    S049 Mouse 84 A: 45-46%
    S044 Mouse 71 A: 46% (slightly elevated dups, ~34%)
    

Current thoughts: After reviewing the FastQC files, I may drop S041, S008, and S054 from the final analysis since we have the N, but I'll wait to see how the rest of the QC looks.  They're notably different from the rest of the samples, but I'm not sure if it's enough to matter in the final analysis.  Will keep the others - they have a little 2nd peak/shoulder, but they look essentially alright.



Cones:

%Dups looks good, no major outliers.

%GC also looks okay.  Lowest is 42% (S093, which also failed GC content, and S079, which passed GC content), while highest is 47%.  Quite a few failed the "Per Sequence GC Content" test, though.  It's actually a rather bizarre gradual shifted peak, not like the clear outliers in the microglia data...might be legit.

Overrepresented Sequences generally looks good.  S089, S090, S077, S078, and S093 are >3%, but the highest is ~4.5%, so it's pretty clean.

Looked into S093.  Its Per Sequence GC Content looks similar to S041, S008, and S054, but it looks like many of the cone samples have an imperfect distribution curve.  Again, not sure how much it matters for RNA-seq - may reflect true biological differences.  Nothing that screams "MUST REMOVE."



Upload next script(s):
rsync -av --progress "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_microglia_rnaSeq/" "ch220@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/christin/2018_microglia_rnaSeq/"

THOUGHTS

