# 2018_chickRFZ_rnaSeq

Christin Hong

HMS, Cepko Lab

Collaborators: Jiho Choi (main contact), Susana da Silva, Nathan Mundell

* Introduced to this project by Nico Lonfat's suggestion to Connie

*This README is written in GitHub Markdown. Outside GitHub, its general formatting can be viewed by copying and pasting into [http://markdownlivepreview.com/].*


---

## Changelog
Unreleased: Collected metrics, annotated, and merged BAMs with Picard

2018-08-12: README created. Completed FastQC, Trimmomatic, and STAR multi-sample 2-pass mapping.

---

## TOC
1. [Intro](#intro)
1. [Raw data](#raw-data)
1. [Controls](#controls)
1. [Data collection notes](#data-collection-notes)
1. [Setup](#setup)
1. [TODOs](#todos)
1. [Acknowledgments](#acknowledgments)


## Intro

## Raw data
Nathan's data was symlinked from his group folder with:

```bash
ln -s /n/data2/hms/genetics/cepko/Nathan/NextSeq/150227_NS500531_0022_AH5JJLBGXX/Data/Intensities/BaseCalls/Run_22/* ~/2018_chickRFZ_rnaSeq/data/nathan/nextSeq_b1

ln -s /n/data2/hms/genetics/cepko/Nathan/NextSeq2/150326_NS500531_0031_AH2L3MBGXX/Data/Intensities/BaseCalls/Run_22/* ~/2018_chickRFZ_rnaSeq/data/nathan/nextSeq_b2
```

* Symlinked directories can be deleted with "rm [name of link]", no / at the end.  ("rm [name of link]/" will be interpreted as the actual directory, which is precious and needs to be left alone)


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


## Data collection notes
From Susana and Nathan

## Setup

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
* `projectDir/bash/` and `projectDir/R/`: For bash and R subscripts
* `projectDir/data/`: For raw data and symlinks to raw data. Once data is in, can lock as READ ONLY. Best practice: Before any analysis, backup all data to an external drive.
* `projectDir/doc/`: Human-readable, things to share with other people. Reports, documents, Markdown, LaTeX, manuscripts, etc.
* `projectDir/doc/graphs/`: Graphics, figures, charts, pdfs, etc.

* `/n/scratch2/ch220/`: Due to limited personal space (100 GB), using scratch2 for output (10 TB of file space/user, files auto-purged after 30 days of no access). For processed data, logs, and other output. Can always be able to delete and regenerate the contents of this entire folder.

### Reference files
* Ensembl Galgal5 reference genome assembly
* Ensembl Galgal5 GTF annotation
* STAR genome index from genome assembly + annotation

### FASTQ files in folders with numeric suffixes for job array, e.g. Sample1.1/RFZ-1-A01_S1_L001_R1_001.fastq.gz
    # For STAR, this pipeline assumes:
        # No spaces in file names
            # I adore bash, but it's terrible at handling filenames with spaces. If that's a major concern, I sincerely believe that it's easier to write a script to rename all the files at the very beginning or learn how to code in Python.
        # FASTQ files are compressed as "*.fastq.gz"
        # Lane IDs are formatted as "_L###_"
        # Prefix is before the final 2 underscores (so from previous example, prefix would be "RFZ-1-A01_S1_L001")




## TODOs

## Acknowledgments
* Include citations for tools

