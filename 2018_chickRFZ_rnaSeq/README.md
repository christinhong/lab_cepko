# 2018_chickRFZ_rnaSeq

Christin M. Hong

HMS, Cepko Lab

Collaborators: Jiho Choi (main contact), Susana da Silva, Nathan Mundell

* Brought on this project by Nico Lonfat's suggestion to Connie

*This README is written in GitHub Markdown. Outside GitHub, its formatting can be viewed by copying and pasting into [http://markdownlivepreview.com/].*


---

## Changelog
Unreleased: Collected metrics, annotated, and merged BAMs with Picard

2018-08-12: README created. Done with FastQC/MultiQC, Trimmomatic, STAR multi-sample 2-pass mapping.

---

## TOC
1. [Intro](#intro)
1. [Raw data](#raw-data)
1. [Controls](#controls)
1. [Data collection notes](#data-collection-notes)
1. [Setup](#setup)
1. [TODOs](#todos)
1. [Acknowledgments](#acknowledgments)


---

## Intro

## Raw data

## Controls
From Connie


## Data collection notes
From Susana and Nathan

## Setup

## TODOs

## Acknowledgments
Include citations for tools


---


# 0) Directory structure
    # jobLogs
    # resources: General use analysis resources
        # ref: Reference files
            # genomes
        # tools: General use software and scripts

    # projectDir: Top level analysis scripts, including this one. Ideally includes a README with notes on project, collaborators, data collection, etc.
        # bash: Bash subscripts
        # data: Original data and symlinks to original data. Once data is in, can lock as READ ONLY. Best practice: Before any analysis, backup all data to an external drive.
        # doc: Reports, documents, Markdown, LaTeX, manuscripts, etc.
            # graphs: graphics, figures, charts, pdfs, etc.
        # output: For processed data, logs, and other output. Can always be able to delete and regenerate the contents of this entire folder. Storage in temp space is OK - may also write temp output to /n/scratch2/ch220 (10 TB of file space/user, files auto-purged after 30 days of no access).
        # logs
		# test
        # R: R subscripts


# 1) Reference files
    # Ensembl Galgal5 reference genome assembly
    # Ensembl GTF annotation
    # STAR genome index from genome assembly + annotation


# 2) FASTQ files in folders with numeric suffixes for job array, e.g. Sample1.1/RFZ-1-A01_S1_L001_R1_001.fastq.gz
    # For STAR, this pipeline assumes:
        # No spaces in file names
            # I adore bash, but it's terrible at handling filenames with spaces. If that's a major concern, I sincerely believe that it's easier to write a script to rename all the files at the very beginning or learn how to code in Python.
        # FASTQ files are compressed as "*.fastq.gz"
        # Lane IDs are formatted as "_L###_"
        # Prefix is before the final 2 underscores (so from previous example, prefix would be "RFZ-1-A01_S1_L001")


# 3) Resources in GLOBAL VARIABLES section



## Data collection ##


## Data structure ##
Nathan's data was symlinked from his group folder with:

```bash
ln -s /n/data2/hms/genetics/cepko/Nathan/NextSeq/150227_NS500531_0022_AH5JJLBGXX/Data/Intensities/BaseCalls/Run_22/* ~/2018_chickRFZ_rnaSeq/data/nathan/nextSeq_b1

ln -s /n/data2/hms/genetics/cepko/Nathan/NextSeq2/150326_NS500531_0031_AH2L3MBGXX/Data/Intensities/BaseCalls/Run_22/* ~/2018_chickRFZ_rnaSeq/data/nathan/nextSeq_b2
```

* Symlinked directories can be deleted with "rm [name of link]", no / at the end.  ("rm [name of link]/" will be interpreted as the actual directory, which is precious and needs to be left alone)
