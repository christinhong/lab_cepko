# 2018_chickRFZ_rnaSeq

Christin Hong

Harvard Medical School, Cepko Lab

Collaborators: Jiho Choi (main contact), Susana da Silva, Nathan Mundell

* Introduced to this project by Nico Lonfat's suggestion to Connie

*This README is written in GitHub Markdown. Outside GitHub, its general formatting can be viewed by copying and pasting into [http://markdownlivepreview.com/].*


---

## Changelog
Unreleased: Collected metrics, annotated, and merged BAMs with Picard

2018-08-12: README created. Completed FastQC, Trimmomatic, and STAR multi-sample 2-pass mapping after extensively testing parameters for mapping NextSeq data.

---

## TOC
1. [Intro](#intro)
1. [Data collection notes](#data-collection-notes)
1. [Controls](#controls)
1. [Cluster setup](#cluster-setup)
	1. [Directory structure](#directory-structure)
	1. [Reference files](#reference-files)
	1. [Raw data](#raw-data)
	1. [Sample script header](#sample-script-header)
1. [Acknowledgments](#acknowledgments)


## Intro


## Data collection notes
From Susana and Nathan



## Controls
From Connie, expect the following expression patterns:

* Fgf8: High in RFZ and maybe T
* Vax: High in V and low in D
* Tbx5: High in D and low in V
* Raldh1: High in D
* Emb: High in RFZ
* Cyp26c1: High in RFZ
* Cyp26a1: High in RFZ and N and T
* Tbx2, 3, and 5: More D than V, but on a gradient


## Cluster setup
The HMS O2 cluster has ~7000 nodes with 32 cores and 256 GB RAM each. It's managed by the Slurm job scheduler; for how to use Slurm, see [https://wiki.rc.hms.harvard.edu/display/O2/O2].

* Can request a max of 20 cores/node, but the more cores/memory requested, the longer the pend time.

* To run an interactive session for testing commands, I generally use `srun --pty -p interactive -c 8 --mem=64G -t 0-11:00 /bin/bash`. The interactive queue has a max of 12 hours.



### Directory structure
```
jobLogs/

resources/
resources/ref/
resources/ref/genomes/
resources/tools/

projectDir/
projectDir/bash/
projectDir/data/
projectDir/doc/
projectDir/doc/graphs/
projectDir/R/

/n/scratch2/ch220/
```

* `jobLogs/`: Where I write out sterr and stout from submitted jobs.

* `resources/`: General use analysis resources. Reference files, genomes, software and scripts.

* `projectDir/`: Contains top level analysis scripts and a README with project background, collaborators, data collection notes, changelog, etc.
	* `projectDir/bash/`: For bash subscripts
	* `projectDir/data/`: For raw data and symlinks to raw data. Once data is in, can lock as READ ONLY. Best practice: Before any analysis, backup all data to an external drive.
	* `projectDir/doc/`: Human-readable, things to share with other people. Reports, documents, Markdown, LaTeX, manuscripts, etc.
	* `projectDir/doc/graphs/`: Graphics, figures, charts, pdfs, etc.
	* `projectDir/R/`: For R subscripts
	
* `/n/scratch2/ch220/`: Due to limited personal space (100 GB), using scratch2 for output (10 TB of file space/user, files auto-purged after 30 days of no access). For processed data, logs, and other output. Can always be able to delete and regenerate the contents of this entire folder.


### Reference files
* Ensembl Galgal5 reference genome assembly
* Ensembl Galgal5 GTF annotation
* STAR genome index from genome assembly + annotation

## Raw data
* Genetics server at research.files.med.harvard.edu/genetics/??
* On HMS O2:
	* I've symlinked Nathan's NextSeq data with:

		`ln -s /n/data2/hms/genetics/cepko/Nathan/NextSeq/150227_NS500531_0022_AH5JJLBGXX/Data/Intensities/BaseCalls/Run_22/* /home/ch220/2018_chickRFZ_rnaSeq/data/nathan/nextSeq_b1`

		`ln -s /n/data2/hms/genetics/cepko/Nathan/NextSeq2/150326_NS500531_0031_AH2L3MBGXX/Data/Intensities/BaseCalls/Run_22/* /home/ch220/2018_chickRFZ_rnaSeq/data/nathan/nextSeq_b2`
	
		(Symlinked directories can be deleted with `rm <name of link>`, no / at the end. `rm <name of link>/` will be interpreted as the actual directory, which is precious and needs to be left alone.)
	
	* I've uploaded Susana's HiSeq data to `/home/ch220/2018_chickRFZ_rnaSeq/data/susana/hiSeq`
* Dropbox (Susana)
* eCommons (Nathan)


### Data format
* Paired end FASTQ files (R1 for forward read and R2 for reverse)
* FASTQ files are compressed as necessary by gunzip
* Each sample has its own folder
* Initial scripts format all filenames into the following structure:
	* `tissue-sample-well_sampleNumber_laneNumber_R1orR2_batchNumber.fastq.gz`
	* E.g. `Sample1.1/RFZ-1-A01_S1_L001_R1_001.fastq.gz`
* Samples are matched by lane assuming that lane regex is `_L###_`


## Sample script header
```bash
#!/bin/bash

#SBATCH -p short                # Required, partition/queue for submission
#SBATCH -t 0-10:00              # Required, time allowed for script to run in D-HH:MM
#SBATCH -c 4                    # number of cores/cpus per node (32 c/node)
#SBATCH --mem=48G               # total RAM requested per job (256 GB RAM/node)

#SBATCH -D /home/ch220/2018_chickRFZ_rnaSeq         # set working directory
#SBATCH --open-mode=append                          # append adds to outfile, truncate deletes old outfile first

#SBATCH -e /home/ch220/jobLogs/RNA-05_%A-%a.err     # standard err
#SBATCH -o /home/ch220/jobLogs/RNA-05_%A-%a.out     # standard out
#SBATCH --mail-type=END                             # email when job ends
#SBATCH --mail-user=christinhong@g.harvard.edu      # address for email
	

#### INTRO ####

# Christin M. Hong
# Last modified: 2018-08
# Harvard Medical School, Connie Cepko Lab

# Script for differential expression analysis of chick RNA-seq data.
    # Running bash on HMS O2 cluster (Slurm). Decided to keep flexibility of array command by leaving it outside this file. Then can choose each time which values to run.
        # Submit with "sbatch --array=1-50 <script.sh>" for NextSeq samples, or "sbatch --array=51-60 <script.sh>" for HiSeq samples.

    # Tasks
        # Second pass mapping of reads to Galgal5 with STAR via job array (1 job per sample)



#### INFRASTRUCTURE ####

# Stop script if error occurs
set -Eeuo pipefail		# See https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/ and http://redsymbol.net/articles/unofficial-bash-strict-mode/


# GLOBAL VARIABLES

# Job-specific
export intCores=4
export pathLogs=/home/ch220/jobLogs


# Experiment-specific
export pathProj=/home/ch220/2018_chickRFZ_rnaSeq # Project's working directory, set for script with "SBATCH -D"
export pathBash=${pathProj}/bash
export pathDoc=${pathProj}/doc

export pathData2=/n/scratch2/ch220/fq_trimmed

export pathOut=/n/scratch2/ch220
export pathOutStar=${pathOut}/starMap                   # S1-50 is from nextSeq. S51-60 is from hiSeq.
export fileSJ=${pathOutStar}/sj_51-60.txt               # File of splice junctions


export varReadL=49                                      # NextSeq reads are 32 bp. HiSeq reads are 50 bp.
export pathStarInd=${pathOut}/STAR-gg5-sjo${varReadL}   # Alt: STAR-gg5 (sjo100), STAR-gg5-sjo49


# References
export pathRef=/home/ch220/resources/ref
export fileGenome=${pathRef}/genomes/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa


# Tools and general scripts with versions in their respective paths 
export pathTools=/home/ch220/resources/tools
export parallel=${pathTools}/parallel-20180722/src/parallel


# Modules loaded from O2
    # "As Lmod is a hierarchical system, you may need to load a prerequisite module to be able to load you need (or see it in module avail). For example, if you'd like to load the Cython module, but don't know the prerequisite modules, run module spider cython. Then, you can load the prerequisites (listed in the module spider output), and will be able to load Cython." -https://wiki.rc.hms.harvard.edu/display/O2/Moving+from+Orchestra+to+O2

module load gcc/6.2.0 python/2.7.12 # Dependencies
module load fastqc/0.11.3
module load multiqc/1.5
module load trimmomatic/0.36
module load star/2.5.4a


# Other options
LC_COLLATE=C    # specifies sort order (numbers, uppercase, then lowercase)


# Notes
    # Export: Variables are inherited by child processes (e.g. subshells from GNU parallel).
    # Bash variables are untyped by default.  For more robust code, can declare data type with [declare] (see http://tldp.org/LDP/abs/html/declareref.html ), but I'm not sure how declare works with export.  May try later.
    # When possible, using full path to minimize confusion by shell, record tool versions, and increase clarity regarding dependencies.
```

I usually write bash in gedit and R in RStudio, but probably any editor with syntax highlighting will do.


## Acknowledgments
* Include citations for tools

