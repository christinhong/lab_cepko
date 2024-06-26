# 2018_chickRFZ_rnaSeq

Christin Hong

Harvard Medical School, Cepko Lab

Collaborators: Jiho Choi, Susana da Silva, Nathan Mundell

* Introduced to this project by Nico Lonfat's suggestion to Connie.

*README.md is written in GitHub Markdown. README.html is the .md in HTML (via RStudio) for easier viewing outside GitHub.*

---

## TOC
1. [Changelog](#changelog)
1. [Background](#background) - TODO
1. [Controls](#controls)
1. [Data collection](#data-collection)
1. [Cluster setup](#cluster-setup)
	1. [Directory structure](#directory-structure)
	1. [Reference files](#reference-files)
	1. [Raw data](#raw-data)
1. [Scripting tips](#scripting-tips)
	1. [Some bash commands](#some-bash-commands)
	1. [Parallelization](#parallelization)
	1. [Sample script header](#sample-script-header)
1. [Analysis](#analysis)
	1. [Pipelines](#pipelines)
	1. [FastQC and Trimmomatic](#fastqc-and-trimmomatic) - TODO
	1. [STAR mapping](#star-mapping)
	1. [On mapping NextSeq reads](#on-mapping-nextseq-reads)
	1. [Thoughts on mapped data](#thoughts-on-mapped-data)
	1. [BAM QC](#bam-qc)
	1. [featureCounts](#featurecounts)
	1. [PCA and dendrogram](#pca-and-dendrogram)
	1. [DESeq2](#deseq2)
	1. [Visualization](#visualization)
	1. [Complementary DE analyses](#complementary-de-analyses)
	1. [Probe design from transcriptome](#probe-design-from-transcriptome)
1. [Acknowledgments](#acknowledgments)

---

## Changelog
2018-09-12: Finished first DE analysis!

2018-09-05: Finished initial annotation and QC plots. Starting analysis with DESeq2.

2018-09-03: Done with featureCounts! On to analysis in R.

2018-08-30: Finished moving all work to Cepko group folder. Done annotating, merging BAMs, and collecting BAM QC metrics with Picard and Qualimap.
* Also wrote a small script for identifying regions with pileup of N+ aligned reads and extracting reference genome sequence for those regions, because it seemed easier than messing with IGV for probe development. Maybe it'll help?

*2018-08-17 to 2018-08-24: Away on family vacation.*

2018-08-13: HMS approved my request for increased group space! `/n/data2/hms/genetics/cepko` has gone from 1 to 10 TB. Will set up shop in group for easier resource sharing and avoiding the auto-purge in scratch2.

2018-08-12: Initialized project README. Tested parameters for mapping NextSeq data. Finalized FastQC, Trimmomatic, and STAR multi-sample 2-pass mapping. Code refactoring.

---

## Background
 :exclamation:


---

## Controls
From Connie, expect the following gene expression patterns:

* Fgf8: High in RFZ and maybe T
* Vax: High in V and low in D
* Tbx5: High in D and low in V
* Aldh1a1 (previously Raldh1): High in a specific region of D
* Emb: High in RFZ
* Cyp26c1: High in RFZ
* Cyp26a1: High in RFZ and N and T
* Tbx2: More D than V, but on a gradient
* Tbx3: More D than V, but on a gradient


---

## Data collection

**Organism**: Chick (Galgal5)

* **Groups for comparison**: 5 retinal tissue regions
	1. RFZ
	1. Dorsal
	1. Ventral
	1. Nasal
	1. Temporal

#### From Susana (retina 6-7)
* Identified RFZ progenitors at E5+1DIV by RARE reporter
* HMS Biopolymer Facility quantified RNA, diluted samples for SPIA amplification, and generated libraries
* Ran on HiSeq for 2x50 bp reads
* 2 retina * 5 tissues * 1 lane = 10 FASTQ files

#### From Nathan (retina 1-5)
* Identified RFZ progenitors at E5+1DIV by RARE reporter
* Input characteristics: 
	* RIN > 9 for all samples
	* ~125 ng RNA per sample
	* PolyA primer for selection
	* Libraries are not strand-specific
	* Generated 2 library sets following SmartSeq2 protocol: 1 set without RNA pre-amplification (failed due to reads only mapping to 3 loci) and 1 set with 12 cycles of pre-amplification with Kappa HF followed by Nextera XT reactions (OK)
* Ran on NextSeq for 2x32 bp reads
	* Expected ~685M reads/batch
* Batch 1: Mix of sets with and without RNA pre-amplification -> fewer reads than expected (~260M)
* Batch 2: Re-run of same pre-amplified set used in batch 1
* 5 retina * 5 tissues * 4 lanes * 2 runs = 200 FASTQ files


---

## Cluster setup
The HMS O2 cluster has ~7000 nodes with 32 cores and 256 GB RAM each. It's managed by the [Slurm job scheduler](https://wiki.rc.hms.harvard.edu/display/O2/O2).

On Linux and Macs, log in from the terminal with `ssh <eCommons username>@o2.hms.harvard.edu`

My home directory is at `/home/ch220`, as ch220 is my eCommons username. Quota is 100 GB.

The Cepko Lab group folder is at `/n/data2/hms/genetics/cepko` As of 2018-08, its quota is 10 TB.

* Can request a max of 20 cores/node, but the more cores/memory requested, the longer the pend time.

* To run an interactive session, I usually use `srun --pty -p interactive -c 8 --mem=64G -t 0-11:00 /bin/bash`. The interactive queue has a max of 12 hours.


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
* Ensembl Galgal5 toplevel reference genome assembly (no primary assembly available)
	* :star: **NOTE: On further thought, I think I also prefer using the toplevel assembly while viewing its repetitive regions as read sinks, and running STAR to exclude reads that map to too many loci.** It reduces the number of usable mapped reads, but it's also more conservative and reduces the risk of mapping reads incorrectly. Ideally, this would be with 2x100 (or 2x150) bp reads to maximize mapping accuracy. Then after mapping, reads that map primarily/strongly to problematic regions would be removed from downstream analysis. (This approach was also suggested by Babraham Bioinformatics, the devs of FastQC, at https://sequencing.qcfail.com/articles/genomic-sequence-not-in-the-genome-assembly-creates-mapping-artefacts/.)
* Ensembl Galgal5 GTF annotation
* STAR genome index from genome assembly + annotation


### Raw data
The raw data for this project can be found at:

* HMS Genetics server:  research.files.med.harvard.edu/genetics/Chick_Retina_Expression
* HMS O2:
	* I've symlinked Nathan's NextSeq data with:

		`ln -s /n/data2/hms/genetics/cepko/Nathan/NextSeq/150227_NS500531_0022_AH5JJLBGXX/Data/Intensities/BaseCalls/Run_22/* /home/ch220/2018_chickRFZ_rnaSeq/data/nathan/nextSeq_b1`

		`ln -s /n/data2/hms/genetics/cepko/Nathan/NextSeq2/150326_NS500531_0031_AH2L3MBGXX/Data/Intensities/BaseCalls/Run_22/* /home/ch220/2018_chickRFZ_rnaSeq/data/nathan/nextSeq_b2`
	
		(Symlinked directories can be deleted with `rm <name of link>`, no / at the end. `rm <name of link>/` will be interpreted as the actual directory, which is precious and needs to be left alone.)
	
	* I've uploaded Susana's HiSeq data to `/home/ch220/2018_chickRFZ_rnaSeq/data/susana/hiSeq`
* Dropbox (Susana)
* eCommons (Nathan)


#### Data format
* Paired end FASTQ.GZ files (R1 for forward read and R2 for reverse, GZ = compressed by gunzip)
* Each sample (unique combination of tissue and biological replicate) has its own folder
* Initial scripts format Hiseq filenames into the NextSeq naming syntax: 
	* `tissue-replicate-well_sampleNumber_laneNumber_R1orR2_batchNumber.fastq.gz`
	* E.g. `Sample1.1/RFZ-1-A01_S1_L001_R1_001.fastq.gz`
	* E.g. `Sample2.2/D-1-A02_S2_L001_R1_001.fastq.gz`
* Paired end reads are matched by lane assuming that lane regex is `_L###_`


---

## Scripting tips
* Be organized and consistent, especially with naming syntax. 
	* It's always worthwhile to make all your initial data files have a consistent, group-aware, and as-simple-as-possible filename syntax. Usually it's best to name from the most to least general characteristic. 
	* `tissue-replicate-well_sampleNumber_laneNumber_R1orR2_batchNumber` works well (as does `tissue-condition-replicate-well_sampleNumber_laneNumber_R1orR2_batchNumber`).

* Immediately strip whitespace from folder and filenames, and avoid having spaces in names as much as possible.
	* If this is infeasible, seriously consider coding in Python instead of bash. Bash wasn't developed for handling whitespace (but I still think it's a wonderful and under-appreciated language).

* Write in an editor with syntax highlighting. I usually write bash in gedit and R in RStudio, but probably any editor with syntax highlighting will do.

* Use [tmux](https://gist.github.com/MohamedAlaa/2961058) for interactive scripting.

* Use variables to minimize code duplication and set custom parameters. Keep them at the beginning of your scripts. Then later, when you need to change them, you won't have to hunt through your entire pipeline to find them all.

* Write scripts to be as modular and reusable as possible.

* Put `set -Eeuo pipefail` at the beginning of every bash script. See https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/ and http://redsymbol.net/articles/unofficial-bash-strict-mode/.

* Use `sacct -j <jobid>` to check whether jobs completed successfully. ExitValues should be 0:0.

* Remember version control with [GitHub](https://confluence.atlassian.com/bitbucketserver/basic-git-commands-776639767.html).


### Some bash commands
Common bash commands include:

* `cd <target dir>`: Change directory
	* Press tab to auto-complete a name, and press tab twice to see potential auto-complete options.
	* `cd ..`: Move up one directory
* `ls`: List contents of current directory
	* `ls *.bam > bams.txt`: List contents of current directory if it ends in ".bam" and write output to bams.txt
* `cat <file name>`: Display contents of file
* `head <file name>`: Display the first 10 lines of a file
* `mkdir <directory name>:` Make a new directory
* `mv -i <old file> <new file>`: Move a file or folder from one location to another, with a prompt if moving will overwrite a file. Also used for renaming.


### Parallelization
Job arrays parallelize across nodes. GNU Parallel parallelizes across cores.

I often use job arrays for between-sample analyses and GNU Parallel for within-sample analyses.

* Why?
	* A cluster is a collection of nodes. On HMS O2, each node has 32 cores and 256G memory.
	* This can also be said as, "A network is a collection of computers. On HMS O2, each computer has 32 threads and 256 GB RAM."
	* Each job is sent to 1 node (unless you're using mpi. That's all you).
	* Keep this structure in mind when parallelizing. You want to parallelize in a way that's analogous to the organization of your data.
	* Say you have ten donors, and you're analyzing all the chromosomes per donor. Each donor has their own folder named with a unique number, but their chromosomes have the same name (so both donor 1 and donor 2 have chr2). It's often a great idea to submit a job per donor, then within that job, simultaneously process the different chromosomes on different cores. On the other hand, submitting a job per chromosome while trying to analyze the different donors on different cores will generally make a mess.
		* Of course, submitting a job per chromosome while analyzing different donors across cores is possible. It all depends on how you've organized and annotated your data.

* Caveats
	* Jobs that require a large number of cores/memory (>8 cores and/or >48 GB RAM) can spend a long time languishing in the queue. Sometimes it's faster to reduce memory requirements by looping instead of parallelizing.
	* When a script generates large intermediate files (e.g. STAR), available disk space can become a limiting factor. Then it may be worthwhile to limit the number of jobs allowed to run from a job array at a time.


### Sample script header
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

# Christin Hong
# Last modified: 2018-08
# Harvard Medical School, Connie Cepko Lab

# Script for differential expression analysis of chick RNA-seq data. See project README.
    # Decided to keep flexibility of array command by leaving it outside this file. Then can choose each time which values to run, e.g. "sbatch --array=1-50 <script.sh>" for NextSeq samples, or "sbatch --array=51-60 <script.sh>" for HiSeq samples.

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
export pathProj=/home/ch220/2018_chickRFZ_rnaSeq
export pathData=${pathProj}/data/*/*
export pathBash=${pathProj}/bash
export pathDoc=${pathProj}/doc


export pathOut=/n/scratch2/ch220
export pathData2=${pathOut}/fq_trimmed
export pathStarInd=${pathOut}/STAR-gg5-sjo49            # STAR indexed with --sjbdOverhang 49

export pathOutStar=${pathOut}/starMap                   
export fileSJ=${pathOutStar}/sj_1-60.txt                # File of STAR splice junctions


# References
export pathRef=/home/ch220/resources/ref
export fileGenome=${pathRef}/genomes/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa


# Tools and general scripts with versions in their respective paths 
export pathTools=/home/ch220/resources/tools
export parallel=${pathTools}/parallel-20180722/src/parallel


# Modules loaded from O2
    # Need to load prereq modules first. Check for prereqs and command syntax with "module spider <tool name>"
module load gcc/6.2.0 python/2.7.12
module load fastqc/0.11.3
module load multiqc/1.5
module load trimmomatic/0.36    # Problematic syntax. Script will need to be manually updated if program is updated.
module load star/2.5.4a


# Other options
LC_COLLATE=C    # specifies sort order (numbers, uppercase, then lowercase)


# Notes
    # "${SLURM_ARRAY_TASK_ID}" sometimes needs to be quoted to be recognized by the shell. Just quote it all the time.
    # Export: Variables are inherited by child processes (e.g. subshells from GNU parallel).
    # Bash variables are untyped by default.  For more robust code, can declare data type with [declare] (see http://tldp.org/LDP/abs/html/declareref.html ), but I'm not sure how declare works with export.  May try later.
    # When possible, using full path to minimize confusion by shell, record tool versions, and increase clarity regarding dependencies.



#### START ####

echo "Starting second pass of multi-sample 2-pass STAR mapping on $(date '+%Y-%m-%d %H:%M:%S')"
echo "Using STAR index ${pathStarInd}"
echo
```

---

## Analysis

### Pipelines
- [ ] In development: Genome-based mapping: STAR multi-sample 2-pass mapping -> featureCounts -> DESeq2
- [ ] Clean up scripts with GNU Make and knitr
- [ ] *Hold: Genome-based mapping plus automated read count processing: (STAR) RSEM (used by Broad)*
- [ ] *Hold: Transcriptome-based mapping: (STAR vs. StringTie for identifying novel transcripts) -> Salmon -> tximport -> DESeq2*
	* https://combine-lab.github.io/salmon/faq/
	
For other options, see https://en.wikipedia.org/wiki/List_of_RNA-Seq_bioinformatics_tools


### FastQC and Trimmomatic
:exclamation:
 

### STAR mapping
Running STAR multi-sample 2-pass mapping, which increases sensitivity for novel splice junctions. From my understanding:

1. STAR is indexed with the published genome annotation, which identifies potential splice sites in the genome assembly.
1. STAR further detects and maps reads to novel splice junctions as it aligns, presumably based on a scoring system for the likelihood of novel splice junction vs. incorrect mapping. Due to this scoring system, reads with short overhangs (e.g. 5-10 bp) across the junctions are unlikely to be mapped to novel junctions.
1. Novel splice junctions that pass STAR's criteria for being legitimate (presumably based on the number and percentage of reads that mapped well across them) are outputted in a "SJ.out.tab" file.
1. Then it's possible to perform a second mapping pass with the SJ.out.tab files from the first mapping pass. In the second pass, STAR inserts the novel junctions into a private copy of the provided reference genome, then re-maps the reads against this personalized genome index.  This increases the sensitivity of read alignment to the novel junctions, presumably because the mapped reads are scored for "standard junction (reference genome) vs. incorrect mapping" rather than "novel junction vs. incorrect mapping."

From [STAR's publication](https://academic.oup.com/bioinformatics/article/29/1/15/272537#84630902): 

>One of the main inherent problems of all de novo RNA-seq aligners is the inability to accurately detect splicing events that involve short (<5–10 nt) sequence overhangs on the donor or acceptor sides of a junction. This causes a significant underdetection of splicing events, and also increases significantly the misalignment rate, as such reads are likely to be mapped with a few mismatches to a similar contiguous genomic region. In addition, this effect also biases the alignments toward processed pseudogenes, which are abundant in the human genome. Similarly to other RNA-seq aligners, to mitigate this problem, STAR has an option to obtain information about possible splice junction loci from annotation databases (Supplementary Section 4). It is also possible to run a second mapping pass, supplying it with splice junction loci found in the first mapping pass. In this case, STAR will not discover any new junctions but will align spliced reads with short overhangs across the previously detected junctions.

There is the question of how accurately STAR detects novel junctions. The developers have some data for that in their paper, but the real reason I'm using it is that STAR is the most popular aligner in the RNA-seq world. It's also used by Broad for their RNA-seq analysis (albeit in the more convenient single-sample 2-pass instead of multi-sample 2-pass). Part of this is simply because STAR is blazingly fast, but if there were any glaring errors in its mapping method, I'm fairly sure they would have been noticed by now. 

I also think STAR is the most flexible aligner. The trend is to move towards *de novo* genome/transcriptome-based analysis, e.g. Tophat searches the transcriptome first, then maps to the genome only if it can't find a good match in the transcriptome. Kallisto and Salmon map solely to the transcriptome. 

In comparison, STAR's default approach is to map to the annotated reference genome, then split the read across the reference genome if it can't find a good contiguous match to see if it can discover a novel splice junction. The transcriptome-based approaches make theorectical sense--in RNA-seq, the transcriptome is what matters--but STAR's approach is probably more forgiving of less-well-annotated genomes, and it's possibly more accurate for calling uniquely mapping reads. (See https://sequencing.qcfail.com/articles/mapping-to-a-transcriptome-can-incorrectly-report-reads-as-mapping-uniquely/, with the counterpoint being at https://cgatoxford.wordpress.com/2016/08/17/why-you-should-stop-using-featurecounts-htseq-or-cufflinks2-and-start-using-kallisto-salmon-or-sailfish/.)

But STAR *also* has the option of providing transcriptome-based counts, so it fits easily into the transcriptome-based approach. Hence its flexibility.

As for 2-pass: I think that if someone is using STAR, they've already bought into its novel splice junction detection method, because it'll do that on every run anyway. If we don't believe in its novel junction detection method, we have to use a different aligner. If we do believe that STAR's novel junction detection is accurate, then 2-pass mapping is part of using STAR well.

(As a genome-based, splice-aware aligner, I think STAR's default output makes it a hybrid between reference and *de novo* genome assembly. See https://biology.stackexchange.com/questions/56158/what-is-contigs-in-picards-reordersam.)


* Increasing mapping speed
	* To run this, I set STAR `--genomeLoad NoSharedMemory`. LoadAndRemove can be faster for the first pass, but it occasionally caused shared memory errors. For the second pass, junction insertion changes the genome, and I suspect that would also lead to shared memory errors.
	* With NoSharedMemory, every job needs ~30 GB of memory just to load the genome. I normally parallelize STAR mapping within each sample, but requesting 64+ GB of RAM significantly increases the amount of time a job spends in the queue. 
	* STAR itself is fast--mapping each read pair takes 2-20 minutes--so it's actually more efficient to loop through each read pair instead of parallelizing. This cuts the RAM requirement to what's necessary for mapping 1 read pair (I use 48 GB/job to be safe) instead of what would be necessary to map all the read pairs in parallel (e.g. 4 pairs per sample (from 4 flowcell lanes per pair) * >30 GB RAM = >120 GB RAM required). Less RAM/job required means faster job allocation from the queue.
	* Then I submit a job for each sample via job array. 
	* Average time for each pass across all samples is ~40-60 minutes, though it's taken >3 hours when the cluster is busy. (Like lab, equipment is most available during nights and weekends.)


### On mapping NextSeq reads
From STAR's output after the first pass, only ~20-25% of NextSeq reads are uniquely mapped vs. ~75-80% of HiSeq reads.

The NextSeq samples have ~40% of reads mapping "to too many loci" (>10 loci) vs 2-10% of HiSeq reads. NextSeq reads are also more likely to be unmapped due to being "too short" (e.g. 20% vs. 10%) or "other" (e.g. 15% vs. 2%). The second NextSeq batch has higher read counts, but the mapped proportions look the same.

Part of this is probably due to the NextSeq reads being 32 bp while the HiSeq are 50 bp + Galgal5 having only the toplevel genome assembly on Ensembl. Unlike primary assemblies, toplevel genome assemblies include haplotypes and repetitive patches. It's likely that the shorter reads map more easily to these repetitive regions.

I've tested the following:

1. Mapping reads to genomes indexed with `--sjdbOverhang 31` or `49` based on read length (default is 100, recommended is ((max read length)-1)).
	1. No effect.
1. Mapping only R1 reads to see if the sort order in paired FASTQ input was an issue (known STAR bug).
	1. No effect.
1. Relaxing the requirements for mapping length with `--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 25`. The average input length for the PE is ~60 (2x32). STAR defaults map if read aligns to at least 2/3 of its total input length, so this allows reads that match on at least 25 bp instead of at least 40 bp to be mapped.
	1. ~4% more reads map uniquely instead of falling in the "too short" category. Makes sense but doesn't solve the problem.
1. Setting `--outFilterMultimapNmax 20` (default is 10, meaning that once a read maps to >10 loci, it's marked as mapping to too many.  20 is seen at https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/).
	1. No effect.
1. Setting `--outFilterMultimapNmax 200`
	1. The "too many loci" reads are now "multiple loci" reads, which makes sense. Of course, this value is far too lenient, but good to know everything's working as expected.
1. Setting `--sjdbScore 1` (default is 2). This decreases mapping to splice junctions. From Alex Dobin (the developer of STAR) at https://groups.google.com/forum/#!msg/rna-star/O1oDItDltjY/0jSn0vy0ccgJ: "I think it is to be expected that some unique mappers become multi-mappers as you add more and more sjdb junctions, since this effectively adds more possibilities for the reads to align. Note, that by default the --sjdbScore = 2, which means that STAR will try to map aggressively to the sjdb junctions, preferring spliced alignment with 1 mismatch to an unspliced alignment without mismatches. You may want to try to reduce this parameter, though it will lead to yet another slight decrease in the % of unique mappers."
	1. No effect.

Unfortunately, I think these reads are genuinely mapping non-specifically. I'm tempted to relax the mapping length requirement, but quality > quantity. To be safe, I'll use the STAR defaults while using the genome index generated with `--sjdbOverhang 49` for processing all samples.

Hmm...this *could* be due to shorter read lengths + repetitive haplotypes/patches in the toplevel genome assembly, but the difference between the HiSeq and NextSeq data seems too dramatic. Will keep digging.


#### UPDATE
After looking around further, I discovered that another common reason for a high percentage of multimappers is rRNA "contamination" (incomplete rRNA depletion).

>2a. Contamination with exogenous sequences.
>You can try to BLAST a few of unmapped reads agains the full NCBI database to see if you get any good matches.
>
>2b. Contamination with ribosomal RNA.
>If your samples are "total RNA", depleted with Ribo-Zero or Ribo-Minus kits, it is possible that the depletion did not work well. rRNA are typically multi-mappers (and you get plenty of those), however, not all rRNA repeats make it into the main chromosomal assembly, and in this case they will not be mapped and will be reported as "alignment too short". We have recently had many cases like that in our lab for human tissues. I believe for the fly genome, the unplaced contigs are in chrU and chrUextra - please try to include them in the genome if you have not done so.
-https://groups.google.com/forum/#!msg/rna-star/7RwKkvNLmI4/REpWc1B4KDkJ

Checking for rRNA contamination in chick was less straightforward than I expected - seems like Galgal rRNA genes aren't really annotated. But the FastQC reports provide the top overrepresented sequences, so I copied the top 3-5 sequences from a sampling of NextSeq files into [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and searched the BLAST results for "gallus gallus." The first *Gallus gallus* result, if available, is listed below:

1. D-1-A02_S2_L001_S2_R1_paired_trimmed.fq.gz
	1. `GTACAAAGGGCAGGGACTTAATCAACGCGA`: PREDICTED: Gallus gallus 18S ribosomal RNA (LOC112533603), rRNA
	1. `ATCAGATACCGTCGTAGTTCCGACCATAAA`: PREDICTED: Gallus gallus 18S ribosomal RNA (LOC112533603), rRNA 
	1. `TCGTAGTTCCGACCATAAACGATGCCGA`: PREDICTED: Gallus gallus 18S ribosomal RNA (LOC112533603), rRNA
1. N-2-B05_S10_L001_S35_R1_paired_trimmed.fq.gz
	1. `GTACGGAAGCAGTGGTATCAACGCAGAGTA`: No hit
	1. `CTTCCGTACTCTGCGTTGATACCACTGCTT`: No hit
	1. `TCCGTACTCTGCGTTGATACCACTGCTTC`: No hit
	1. `GAGTACGGAAGCAGTGGTATCAACGCAGAG`: No hit
		* The sequences seem to be related, e.g. they all pull up multiple "common carp genome, scaffold" results. But nothing for chicken. Weird, but it is what it is.
	1. `GTACAAAGGGCAGGGACTTAATCAACGCGA`: PREDICTED: Gallus gallus 18S ribosomal RNA (LOC112533603), rRNA
1. RFZ-3-C01_S11_L003_S36_R2_paired_trimmed.fq.gz
	1. `GTACAAAGGGCAGGGACTTAATCAACGCGA`: PREDICTED: Gallus gallus 18S ribosomal RNA (LOC112533603), rRNA
	1. `GTACAGTGAAACTGCGAATGGCTCATTAAA`: PREDICTED: Gallus gallus 18S ribosomal RNA (LOC112533602), rRNA
	1. `GGGTAGACACAAGCTGAGCCAGTCAGTGTA`: PREDICTED: Gallus gallus 18S ribosomal RNA (LOC112533603), rRNA
1. T-4-D04_S19_L004_S19_R2_paired_trimmed.fq.gz
	1. `GTACAAAGGGCAGGGACTTAATCAACGCGA`: PREDICTED: Gallus gallus 18S ribosomal RNA (LOC112533603), rRNA
	1. `GTTAAGAGCATCGAGGGGGCGCCGAGAGA`: PREDICTED: Gallus gallus 18S ribosomal RNA (LOC112533603), rRNA
	1. `TCGTAGTTCCGACCATAAACGATGCCGA`: PREDICTED: Gallus gallus 18S ribosomal RNA (LOC112533603), rRNA 
1. V-5-E03_S23_L004_S48_R1_paired_trimmed.fq.gz
	1. `GTACAAAGGGCAGGGACTTAATCAACGCGA`: PREDICTED: Gallus gallus 18S ribosomal RNA (LOC112533603), rRNA
	1. `TCGTAGTTCCGACCATAAACGATGCCGA`: PREDICTED: Gallus gallus 18S ribosomal RNA (LOC112533603), rRNA
	1. `ATCAGATACCGTCGTAGTTCCGACCATAAA`: PREDICTED: Gallus gallus 18S ribosomal RNA (LOC112533603), rRNA 
	
HiSeq samples (retinas 6-7) had no notable overrepresented sequences.

Conclusion: NextSeq libraries probably have a high percentage of multimapping reads due to original input having high levels of rRNA.


### Thoughts on mapped data
NextSeq: The high percentage of reads that were "mapped to too many loci" in the first STAR pass seem to be split across "mapped to too many loci" and "reads unmapped: other" in the second STAR pass. But the percentage of uniquely mapped reads is stable (~20-25%), so that's probably fine.

HiSeq: There's a slight drop in uniquely mapping and rise in multimapping reads after the 2nd pass, which makes sense. It's a very minor shift (from ~75% to ~72%), so I would say these are also stable.

In terms of uniquely mapped reads, the NextSeq data consistently hovers around 0.5-1 M/BAM * 4 lanes = 2-4 M/sample.

From HiSeq, it's pretty variable, with the range being ~1-12 M/sample.

Ryoji's dendrogram seems to be based on log2+1 read counts, so I'd guess that's the main reason why the datasets looked so different? At any rate, all the samples have information. Will see how the analysis turns out.


### BAM QC

#### BAM QC with Picard (https://broadinstitute.github.io/picard/)

- [x] CreateSequenceDictionary (one-time generation of sequence dictionary for reference genome FASTA)
- [x] CleanSam (avoid errors from reads mapping over ends of chromosomes; happened in sample 15)
- [x] Annotate with AddOrReplaceReadGroups
- [x] Merge BAMs for each sample with samtools
- [x] SortSam to sort BAMs by coordinate 
	* Broad does ReorderSam -> MarkDuplicates instead of SortSam. SortSam sorts by coordinate, while ReorderSam sorts by the ordering in the provided reference genome file **and drops reads that don't map to the provided reference**. That seems counter-productive for STAR's novel junction discovery, so I've decided to use SortSam.
- [x] BuildBamIndex (faster with coordinate-sorted BAMs)
- [x] MarkDuplicates (output: txt file)
	* Removing duplicates is discouraged in RNA-seq analysis because RNA duplicates are most likely biologically informative.
	* "When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates. However, when the input is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads." -https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
- [x] EstimateLibraryComplexity (output: txt file)
- [x] CollectMultipleMetrics
	- [x] CollectAlignmentSummaryMetrics (main output: txt file)
	- [x] CollectInsertSizeMetrics (pdf)
	- [x] QualityScoreDistribution (not supported by MultiQC) (pdf)
	- [x] MeanQualityByCycle (not supported by MultiQC) (pdf)
	- [x] CollectBaseDistributionByCycle (pdf)
	- [x] CollectGcBiasMetrics (pdf)
		* (If necessary, using cqn with DESeq2: https://www.biostars.org/p/259378/ and in https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
		* From the GC% of 100 base windows, looks like chick has an AT-rich genome. From the normalized coverage, it also looks like there's some GC bias (preferential amplification of GC-rich fragmentS) due to how coverage generally increases with GC%. Most likely, the GC-bias will be consistent across samples, so there won't be a need to correct it out.
	- [x] CollectQualityYieldMetrics (txt)
		* May be redundant with FastQC since not doing GATK base recalibration, but technically it's measuring MapQ instead of Phred. See https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/
	* Not doing CollectRnaSeqMetrics because it requires a refFlat instead of GTF annotation file (https://github.com/broadinstitute/picard/issues/805)
	* Not doing CollectSequencingArtifactMetrics because we aren't interested in SNPs


#### Qualimap
Complementary BAM QC with Qualimap: http://qualimap.bioinfo.cipf.es/doc_html/analysis.html
- [x] BamQC
- [x] RNASeq


#### Processing QC output
- [x] Make high throughput
- [x] Run on all BAMs
- [x] Aggregate all BAM QC (Picard and Qualimap) with MultiQC
* Qualimap also has Counts QC post-featureCounts: http://qualimap.bioinfo.cipf.es/doc_html/samples.html


#### Further BAM QC
- [x] Test IGV with Oracle Java. Still laggy.
- [x] Test visualization with SeqMonk (which can also do read density vs. duplication, see https://bioinformatics.stackexchange.com/questions/722/visualisation-of-long-read-rna-seq-splicing). Definitely faster, but less user-friendly.
	* On technical duplicates: https://sequencing.qcfail.com/articles/libraries-can-contain-technical-duplication/ . Deduping isn't appropriate for RNA-seq (especially when these libraries are so saturated), but, "If your concern is with the inflated increase in power from duplication then a better solution might well be to quantitate the read counts as normal, but then try to estimate the overall level of duplication and divide all counts by this amount before moving on to doing statistical analyses.  This won’t change the magnitude of the changes seen, but will reduce the overall number of observations."
	* On "optical duplicates": https://sequencing.qcfail.com/articles/illumina-patterned-flow-cells-generate-duplicated-sequences/
	* An R version for prettier plots of BAM duplication from Babraham is dupRadar: https://sourceforge.net/projects/dupradar/ . But SeqMonk may be enough.
* For further analysis of Picard Metrics, can check out https://github.com/slowkow/picardmetrics


### featureCounts
The command I used for featureCounts was:

```bash
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
```

#### Explanations

* Counting the number of reads that align to a gene (meta-feature) rather than counting by alignment to an exon.
* -t exon: Decided to count only reads mapping to exons to minimize ambiguity.
* -g gene_id: Identifier in the GTF being used for the genes/meta-feature.
* -C: Don't count chimeric fragments = fragments that span multiple chromosomes. I can see chimeric fragments being interesting in cancer research, but I don't see any reason they'd be present here.
* -O: Count fragments that overlap multiple features, e.g. fragments that map to more than one gene. Allowed since I can see thes potentially capturing genes that sit close to each other on the genome.
* --byReadGroup: Counting by read group annotation (see Picard's AddOrReplaceReadGroups above) so I can more easily analyze for differences between them later on.
* To maximize accuracy, not counting multimapping reads = reads that map to more than one location.
* Including duplicates since for RNA-seq, duplicates may be biologically informative (they may come from actual duplicate mRNA transcripts).
* NOTE: If there are issues while running, can add a `--verbose` option for easier debugging.


#### Results

From the Qualimap BAM QC data, ~65-75% of total fragments mapped to exons. FeatureCounts counted ~60-70% of total fragments as mapping uniquely to exons.  Seems fine.  Moving to R!


### PCA and dendrogram
PCA on annotated, log2(x+1) transformed count data suggests that (as expected) batch is the major driver of PC1-3.  Oddly, RFZ-7 shows up with the Retina 6 samples, while the other Retina 7 samples cluster on their own.  Not sure what's going on there...**maybe a timepoint difference?**

I was hoping the major batch-related factor would be library size, but from the number of assigned fragments in featureCounts, it looks like that isn't the case.  Hmm.

Still, the tissues start clustering by PC4.  That's cool.

There might be a problem with D from Retina 7...?  It's further from the other D samples on PC3-7 than I'd expect.  But it isn't crazy far, so I'll leave it in for now.

If necessary, I can run ComBat and check that it doesn't remove the tissue-based differences, but it may not be necessary.  Since the tissues cluster so well on PC3-6, I have a hunch that these samples can actually be analyzed together.  That would be helpful, because if the batch-related differences are major enough to prevent that, the HiSeq samples can't be used for differential gene expression analysis.  It's 1 retina per batch for those, and we can't get statistics with N=1.

(If necessary, there's a package called GFOLD for analyzing single replicate data. Paper is here: https://academic.oup.com/bioinformatics/article/28/21/2782/235811.  X. Shirley Liu is an author, and it seems legit on Biostars.  It isn't really DE analysis, but it's an advanced method for ranking fold changes across samples for potential follow-up.)

- [x] Label Retina 7 into its own batch instead of keeping it in the same batch as Retina 6. I think they're from the same library but the samples were prepped with different protocols, as Susana mentioned that she thinks she submitted too much RNA for the Retina 7 samples.

- [x] Lanes all look consistent by PCA and dendrogam -> Sum lanes together for each sample for easier visualization downstream.


### DESeq2
See https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

edgeR, DESeq, and DESeq2 all have their own methods of correcting for library size, so they all want the raw count data.  (I could normalize the raw counts to TPM - see https://statquest.org/2015/07/09/rpkm-fpkm-and-tpm-clearly-explained/ and https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html - but DESeq2 counts are more sophisticated and preferred on Biostars, so no reason to write the extra code.)

* DESeq2 automatically filters out low count reads. DESeq offers a genefilter package for that. Not sure about edgeR.

- [x] GC bias is different between HiSeq and NextSeq datasets -> Add a new annotation for platform and use that as a covariate. 
	* Other options: cqn (https://support.bioconductor.org/p/95683/)? alpine (https://mikelove.wordpress.com/2016/09/26/rna-seq-fragment-sequence-bias/)?
	* For an explanation of bias, see https://mikelove.github.io/techbias/#/slide-1
- [x] Set RFZ as "control" and other tissues as other conditions to analyze data together? See https://www.biostars.org/p/110266/
- [x] Check DESeq2 results by PCA, heatmap with hierachical clustering, etc.
- [x] If the DESeq2 normalized counts still show a strong batch effect, research batch correction methods
- [x] Run svaseq. (Originally I was going to go with ComBat, but consensus by experienced bioinformaticians is that ComBat is too aggressive and over-fits, which creates beautiful plots, but significantly inflates false positives. I think it's considered only safe to use with spike-in controls. Otherwise, adjusting with SVA is preferred.
- [x] Figure out how to observe effects from svaseq in PCA.  :p  Eventually discovered its use with limma's removeBatchEffect
- [ ] On heatmaps: http://www.opiniomics.org/you-probably-dont-understand-heatmaps/
	* If the batches are still problematic, will probably have to analyze Retina 6 and Retina 7 separately. Can try GFOLD for getting gene rankings (GFOLD paper is at https://academic.oup.com/bioinformatics/article/28/21/2782/235811).
- [x] Add svaseq to DESeq2 to test correcting for NextSeq vs. HiSeq


### Visualization
See https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html

* NOISeq looks interesting for QC purposes? May throw that in to see if it'd be useful.

- [ ] Annotate gene names
- [ ] Get lists of DE genes for RFZ vs. other tissues plus intersect (Venn diagram Jiho showed me)
- [ ] Send to Jiho
- [ ] Heatmap
- [ ] PCA
- [ ] Volcano plot
- [ ] GO pathway analysis and/or KEGG?


### Polish
- [ ] Test removing R2 for D_7 (PCA outlier, 5'-3' bias is 2.36! See https://www.biostars.org/p/102812/ .  I think it also has the lowest GC content) and rerun.  Can check at BAM QC step to see if 5'-3' bias is still present)
- [ ] Test HiSeq alone
- [ ] Test NextSeq alone
- [ ] GNU Make
- [ ] knitr


### Complementary DE analyses
DESeq (more conservative) and edgeR (alternative library normalization method)



### Probe design from transcriptome
(isoform detection)

- [x] Detect pileups in genes of interest and reference genome sequence where they occur
- [ ] Read Sylvain's paper! Currently on bioRxiv. :D
- [ ] Check out OligoMiner
- [ ] Might actually fit into alternative RNA-seq pipelines (StringTie, Salmon)


---


## Acknowledgments
* Citing tools helps keep them funded. Keeping track below; citation format will vary.

Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, Thomas R. Gingeras; STAR: ultrafast universal RNA-seq aligner, Bioinformatics, Volume 29, Issue 1, 1 January 2013, Pages 15–21, https://doi.org/10.1093/bioinformatics/bts635

Tange, Ole. (2018). GNU Parallel 2018. GNU Parallel 2018 (p. 112). Ole Tange. http://doi.org/10.5281/zenodo.1146014

Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2015). "Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data." Bioinformatics, btv566

Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8

