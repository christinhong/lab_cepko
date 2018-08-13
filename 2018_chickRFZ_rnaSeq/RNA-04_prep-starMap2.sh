#### INTRO ####

# Christin M. Hong
# Last modified: 2018-08
# Harvard Medical School, Connie Cepko Lab


#### Request resources on interactive node ####
srun --pty -p interactive -c 4 --mem=24G -t 0-11:00 /bin/bash


#### INFRASTRUCTURE ####

# GLOBAL VARIABLES

# Job-specific
export intCores=4
export pathLogs=/home/ch220/jobLogs


# Experiment-specific
export pathProj=/home/ch220/2018_chickRFZ_rnaSeq # Project's working directory, can set for script with "SBATCH -D"
export pathBash=${pathProj}/bash
export pathDoc=${pathProj}/doc

export pathData2=/n/scratch2/ch220/fq_trimmed

export pathOut=/n/scratch2/ch220
export pathOutStar=${pathOut}/starMap                   # S1-50 is from nextSeq. S51-60 is from hiSeq.
export varReadL=31                                      # NextSeq reads are 32 bp. HiSeq reads are 50 bp.
export pathStarInd=${pathOut}/STAR-gg5-sjo${varReadL}   # Alt: STAR-gg5, STAR-gg5-sjo49


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



#### START ####

echo "Starting STAR alignment and mapping on $(date '+%Y-%m-%d %H:%M:%S')"

# Add newline to output for readability
echo


echo "Collecting 1st pass junctions from 32 bp samples (S1-50) for multi-sample 2-pass STAR mapping"

find ${pathOutStar}/S{1..50}/*SJ.out.tab -name ".*" -prune -o -print | paste -sd " " > ${pathOutStar}/sj_1-50.txt


echo "Collecting 1st pass junctions from 50 bp samples (S51-60) for multi-sample 2-pass STAR mapping"

find ${pathOutStar}/S{51..60}/*SJ.out.tab -name ".*" -prune -o -print | paste -sd " " > ${pathOutStar}/sj_51-60.txt


echo "Collecting 1st pass junctions from samples for multi-sample 2-pass STAR mapping"
find ${pathOutStar}/S{1..60}/*SJ.out.tab -name ".*" -prune -o -print | paste -sd " " > ${pathOutStar}/sj_1-60.txt


echo "Making directory for 2nd pass STAR samples"
mkdir ${pathOutStar}/pass2

