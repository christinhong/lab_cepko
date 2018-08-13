# Christin M. Hong
# Last updated: 2018-08
# Harvard Medical School, Connie Cepko Lab


# Script 1 for differential expression analysis of chick RNA-seq data.
    # Running bash on HMS O2 cluster (Slurm).
    # Tasks
        # Setting up project directory on cluster
        # Importing Galgal5 reference genome assembly and annotation
        # Indexing Galgal5 genome with STAR
        # FastQC with MultiQC compilation of samples
        # Prep samples for STAR aligner job array


#### Intro notes ####
# Need to clean this log up!

# I like to write scripts the same way I write protocols - so a complete novice to the task can pick it up and go. (Also so when I look back years later, I can actually understand what I did and why I did it.)

# This script is NOT meant for submitting as a job. This is a record of the initial commands I used to set up analysis on the cluster. These commands are all run through the interactive partition (see detail in script).

# This text is saved as .sh because I'm using an editor that does syntax highlighting (gedit), and the ".sh" tells the editor to syntax highlight for bash. I strongly recommend reading and writing scripts in editors that do syntax highlighting. 

# I've recorded full paths instead of shortcuts (e.g. /home/ch220/ instead of ~/) and excluded most interactive convenience tricks for transparency and readability. Feel free to make variables, use tmux, etc. for greater convenience at will.

# On OS: Bash is designed for the Linux OS. I'm working from a Linux Mint OS (derived from Linux Ubuntu). If you're working on a PC through Windows Powershell, terminal commands that run on your own computer (not the cluster) may be different. Windows 10 comes with the option to partially enable an Ubuntu subsystem. It isn't perfect, but it might be easier to work from that Ubuntu subsystem instead of the default Windows command line. (Macs should work the same as Linux; they built in Linux OS support while developing.)



#### START ####
# Login with eCommons ID and password
    # Can automate connection later, see http://slowkow.com/notes/ssh-tutorial/
ssh ch220@o2.hms.harvard.edu


# Request resources on interactive node
srun --pty -p interactive -c 8 --mem=24G -t 0-11:00 /bin/bash
    # O2 has 7000 nodes. Each O2 node has 32 cores and 256 GB RAM, so this requests a quarter of a node with flexible memory distribution.


# Set job resource variables
export intCores=8



#### Set up project directory ####
# Bash basics
    # "cd <folder name>" to move into folder
    # "ls <folder name>" to list contents of folder
    # "mkdir <folder name>" to create a folder
    # "rm -rf <folder name>" to delete a folder and ALL contents. WARNING: Be careful with this command!!
        # "rm <file name>" to delete a file
    # "mv <old folder name> <new folder name>" to move and rename files/folders


# Christin's directory structure
    # jobLogs/
    # resources/: General use analysis resources
        # ref/: Reference files
            # genomes/
        # tools/: General use software and scripts

    # projectDir/: Top level analysis scripts, including this one. Ideally includes a README with notes on project, collaborators, data collection, etc.
        # bash/: Bash subscripts
        # data/: Original data and symlinks to original data. Once data is in, can lock as READ ONLY. Best practice: Before any analysis, backup all data to an external drive.
        # doc/: Reports, documents, Markdown, LaTeX, etc.
            # graphs/
        # output/: Processed data and other output. May also write temp output to /n/scratch2/ch220 (10 TB of file space/user, files auto-purged after 30 days of no access).
            # logs/
			# test/
        # R/: R subscripts



#### Install GNU Parallel ####
# GNU Parallel is great for parellelizing across cores in a node (a.k.a. across threads in a computer, semi-independent tasks).
    # Job arrays are great for parellelizing across nodes in a cluster (a.k.a. across computers in a network, completely independent tasks).

cd /home/ch220/resources/tools

# Install latest GNU Parallel per https://www.gnu.org/software/parallel/parallel_tutorial.html
(wget -O - pi.dk/3 || curl pi.dk/3/ || \
   fetch -o - http://pi.dk/3) | bash


# Set variables. Keep versions in their respective paths for future reference.
export pathRef=/home/ch220/resources/ref
export pathTools=/home/ch220/resources/tools
export parallel=${pathTools}/parallel-20180722/src/parallel



#### Import reference genome assembly and annotation ####
# Decided to go with Galgal5 because Galgal6 is still in the process of being annotated, as it came out 2018-03. 
    # For reference, the Galgal5 assembly came out in 2015-12 and its annotation was released in 2016-10, then patched until 2016-12.

# Transfer Ensembl's unmasked Galgal5 reference genome assembly to my cluster folder
    # Using unmasked because hard masking can decrease accuracy and soft masking doesn't make a difference. See http://genomespot.blogspot.com/2015/06/mapping-ngs-data-which-genome-version.html
    # From http://useast.ensembl.org/Gallus_gallus/Info/Index -> "Download DNA sequence (FASTA)"
rsync -av  "rsync://ftp.ensembl.org/ensembl/pub/release-93/fasta/gallus_gallus/dna/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa.gz" "/home/ch220/resources/ref/genomes/"

    # Note: This is a toplevel genome assembly. Toplevel isn't ideal because it includes N-padded haplotype/patch regions (areas that are highly variable or redundant, so hard to map to). For human and murine data, would choose the "primary assembly" option to exclude those regions, but primary assembly isn't available for chicken.


# Get Ensembl's Galgal5 GTF annotation file for splice junctions
    # Can't find this from the Galgal5 page. Have to go to http://www.ensembl.org/info/data/ftp/index.html to see all the available downloads and find the GTF file from there. (I actually found it by hunting through Google and working up from a mouse GTF question.)
rsync -av  "rsync://ftp.ensembl.org/ensembl/pub/release-93/gtf/gallus_gallus/Gallus_gallus.Gallus_gallus-5.0.93.gtf.gz" "/home/ch220/resources/ref/genomes/"

cd /home/ch220/resources/ref/genomes/

# Check checksums
sum Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa.gz    # Expect 58657 353241
sum Gallus_gallus.Gallus_gallus-5.0.93.gtf.gz             # Expect 29883 10614

# Unzip files (STAR uses uncompressed files)
gzip -d Gallus_gallus.Gallus_gallus-5.0*.gz



#### Build STAR genome index ####
    # This only needs to be done once per reference genome.
    # STAR writes large temp files that require >100 GB disk space, which is the size of my entire home folder. Decided write indices to scratch2 system, which allows up to 10 TB but deletes data after 30 days of non-use.

module load star/2.5.4a

tmux
# Like screen, tmux allows a process to run in the background. Unlike screen, tmux keeps the process running even after logging out of the parent shell. Will inherit loaded modules and exported variables from the parent shell (but I think they have to be loaded while no tmux child shells are open). Can detach from tmux with Ctrl+b, then d. Reattach with tmux attach. List windows with tmux ls. Kill the current window with Ctrl+b, x. For more, see https://lukaszwrobel.pl/blog/tmux-tutorial-split-terminal-windows-easily/

mkdir /n/scratch2/ch220
mkdir /n/scratch2/ch220/STAR-gg5

cd /n/scratch2/ch220/STAR-gg5

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir . --genomeFastaFiles /home/ch220/resources/ref/genomes/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa --sjdbGTFfile /home/ch220/resources/ref/genomes/Gallus_gallus.Gallus_gallus-5.0.93.gtf --sjdbOverhang 100

# On --sjdbOverhang: STAR's manual notes that, "In case of reads of varying length, the ideal value is max(ReadLength)-1. In most cases, the default value of 100 will work as well as the ideal value."
    # HiSeq reads are 50 bp while NextSeq are 32 bp. Not sure how much it matters, but I'll make indices for each vs. the default and compare the results.


mkdir /n/scratch2/ch220/STAR-gg5-sjdb49

cd /n/scratch2/ch220/STAR-gg5-sjdb49

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir . --genomeFastaFiles /home/ch220/resources/ref/genomes/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa --sjdbGTFfile /home/ch220/resources/ref/genomes/Gallus_gallus.Gallus_gallus-5.0.93.gtf --sjdbOverhang 49

#

mkdir /n/scratch2/ch220/STAR-gg5-sjo31

cd /n/scratch2/ch220/STAR-gg5-sjo31

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir . --genomeFastaFiles /home/ch220/resources/ref/genomes/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa --sjdbGTFfile /home/ch220/resources/ref/genomes/Gallus_gallus.Gallus_gallus-5.0.93.gtf --sjdbOverhang 31


# Note: While STAR was indexing in these tmux windows, I ran the sample files through FastQC in other tmux windows.


# When indexing is finished, make directory for index and copy final STAR index to it
    # WARNING: "If you are transferring files to /n/scratch2 using a tool and flag to preserve timestamps (e.g. rsync -a or -t), those files will also be subject to the deletion policy based on the original timestamp. . . . Please be very judicious about handling files when moving them to or generating them on /n/scratch2; as mentioned above, if you are affected by this behavior, the files are unrecoverable." -https://wiki.rc.hms.harvard.edu/display/O2/Filesystems#Filesystems-Scratchdirectory(/n/scratch2/ab123)

mkdir ${pathRef}/STAR-gg5

cp -R "/n/scratch2/ch220/STAR-gg5/." "${pathRef}/STAR-gg5/"

cp -R "/n/scratch2/ch220/STAR-gg5-sjdb49/." "${pathRef}/STAR-gg5/sjo49"



#### Format sample FASTQ file names ####
# Susana's files are fastq.bz2 while Nathan's are fastq.gz. Will convert all to fastq.gz - consistency always helps.
    # Note that I already manually renamed Susana's files to match the file name format of Nathan's files (for 20 files, manual renaming was faster for me than writing a script)
    # Keeping files in compressed format (*.fastq.gz) because most NGS tools work faster with compressed data

cd /home/ch220/2018_chickRFZ_rnaSeq/data

bzip2 -d /home/ch220/2018_chickRFZ_rnaSeq/data/susana/*/*/*.fastq.bz2
    # Decompressing one by one was pretty slow. Next time I'd parallelize.

${parallel} -j ${intCores} --verbose --joblog /home/ch220/2018_chickRFZ_rnaSeq/output/logs/parallel-gzip.log --resume-failed --keep-order "gzip {}" ::: /home/ch220/2018_chickRFZ_rnaSeq/data/susana/hiseq/*/*.fastq


# Appending batch number to hiSeq samples to match nextSeq sample syntax
rename .fastq.gz _001.fastq.gz /home/ch220/2018_chickRFZ_rnaSeq/data/susana/hiseq/*/*.fastq.gz


# Would be helpful if Nathan's two datasets had unique file names...
    # WARNING: Renaming is an easy way to accidentally erase data!! Be careful - make a copy of the data first to test renaming.
    # If data is accidentally overwritten, remember that data on group and home folders is backed up by snapshots.

rename _001.fastq.gz _002.fastq.gz /home/ch220/2018_chickRFZ_rnaSeq/data/nathan/nextSeq_b2/*/*.fastq.gz


    # Note: It's safer to prepend or append additional info since that doesn't risk changing the uniqueness of the names. This is a small number of files that are snapshotted, so I felt okay with this renaming strategy, but if this was a large number of files, automated (no manual checks), or unique data, I'd have simply added _b1 or _b2 before the .fastq.gz.
    # (That said, you should NEVER be working on unique data. Data is precious - back it up before any potential modifications.)



#### Run FastQC ####
# Can run FastQC on all samples simultaneously or separate by platform. The samples are QC'd independentely, but it can be easier to interpret the MultiQC reports when the platforms are seperate. This is particularly true because the NextSeq and HiSeq reads are different lengths (NextSeq: 32 bp, HiSeq: 50 bp), and the read lengths play a large role in calculating quality.

module load fastqc/0.11.3

mkdir /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_orig

mkdir /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_orig/ctrls
mkdir /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_orig/nextSeq
mkdir /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_orig/hiSeq


${parallel} -j ${intCores} --verbose --joblog "/home/ch220/jobLogs/parallel-fastqc-ctrls.log" --resume-failed --keep-order "fastqc -o /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_orig/ctrls {}" ::: /home/ch220/2018_chickRFZ_rnaSeq/dataCtrls/*/*/*.fastq.gz


${parallel} -j ${intCores} --verbose --joblog "/home/ch220/jobLogs/parallel-fastqc-nextSeq.log" --resume-failed --keep-order "fastqc -o /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_orig/nextSeq {}" ::: /home/ch220/2018_chickRFZ_rnaSeq/data/nathan/*/*/*.fastq.gz


${parallel} -j ${intCores} --verbose --joblog "/home/ch220/jobLogs/parallel-fastqc-hiSeq.log" --resume-failed --keep-order "fastqc -o /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_orig/hiSeq {}" ::: /home/ch220/2018_chickRFZ_rnaSeq/data/susana/*/*/*.fastq.gz


multiqc /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_orig -o /home/ch220/2018_chickRFZ_rnaSeq/doc/multiQC -n multiQC-01_fastqc-orig_$(date '+%Y-%m-%d')




#### Aggregate final FastQC reports with MultiQC ####

# Loading MultiQC
module spider multiqc/1.5   # Finding MultiQC's dependencies
module load gcc/6.2.0  python/2.7.12
module load multiqc/1.5

# By default, MultiQC's output is in whatever directory it runs in, and its input is its directory + all subdirectories. So for now:
cd /home/ch220/2018_chickRFZ_rnaSeq/doc
multiqc .

cd /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_nextSeq
multiqc .

cd /home/ch220/2018_chickRFZ_rnaSeq/doc/fastqc_hiSeq
multiqc .



#### Thoughts ####
# Considering the read lengths and that this data was collected 4 years ago, MultiQC reports look okay for all samples EXCEPT T-7-G04_S34_L008_R1 (failed overrepresented sequences check) and negative controls (as expected).

# The NextSeq samples have a pretty high level of duplicates, but that shouldn't affect the analysis - just means they were sequenced more deeply than necessary. The HiSeq data may be higher quality due to having a longer read length (50 bp vs. 32 bp), which will probably increase the accuracy when mapping. (~100 bp reads is thought to be ideal, as longer reads are more prone to sequencing bias.)

# I think I'll analyze the platforms separately first, especially since we aren't sure if they were collected at the same timepoints. I'll look into how to potentially combine them further after the separate analyses (for that, can start with https://www.biostars.org/p/251059/ and https://www.researchgate.net/post/How_can_you_combine_different_published_expression_datasets_and_analyze_them_in_R).

# Keep in mind that Susana's data analysis will now require adjusting for 1 unpaired read (T-7-G04_S34_L008_R2).



#### Rerun MultiQC after excluding samples and recheck output ####
# On cluster
cd /home/ch220/2018_chickRFZ_rnaSeq/doc
mkdir /home/ch220/2018_chickRFZ_rnaSeq/doc/failedQC
mv fastqc_nextSeq/Negative* failedQC/
mv fastqc_hiSeq/T-7-G04_S34_L008_R1* failedQC/

multiqc . --ignore failedQC/
# MultiQC automatically detects previous ouput and adjusts new output by appending a count to the filenames. If you want to overwrite old reports, use --force.


# From home terminal
scp -r ch220@transfer.rc.hms.harvard.edu:/home/ch220/2018_chickRFZ_rnaSeq/doc/multiqc* /home/christin/aaa_data_CepkoLab/2018_chickRFZ_rnaSeq_data/multiQC

# Looks good!


#### Separating controls ####
# Move samples and data with -n flag to prevent overwrites
mv -n /home/ch220/2018_chickRFZ_rnaSeq/data/nathan/nextSeq_b1/S0000 /home/ch220/2018_chickRFZ_rnaSeq/data/dataCtrls/nextSeq_b1

mv -n /home/ch220/2018_chickRFZ_rnaSeq/data/nathan/nextSeq_b2/S0000 /home/ch220/2018_chickRFZ_rnaSeq/data/dataCtrls/nextSeq_b2



#### Prep remaining samples for alignment job array ####
# List sample directories
ls -d /home/ch220/2018_chickRFZ_rnaSeq/data/*/*/S*


# Append job array index numbers (not 0 padded) to folder names
    # Beauty in one line. Adapted from https://stackoverflow.com/questions/3211595/renaming-files-in-a-folder-to-sequential-numbers
ls -d /home/ch220/2018_chickRFZ_rnaSeq/data/*/*/S* | cat -n | while read n f; do mv -n "$f" "$f.$n"; done 



#### Make output folder for Trimmomatic ####
# Many tools generate large output files. Generate in scratch to avoid "No space left on disk" error.
mkdir /n/scratch2/ch220/fq_trimmed





