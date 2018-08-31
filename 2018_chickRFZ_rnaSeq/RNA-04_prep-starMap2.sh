#### INTRO ####

# Christin M. Hong
# Last modified: 2018-08
# Harvard Medical School, Connie Cepko Lab


#### Request resources on interactive node ####
srun --pty -p interactive -c 4 --mem=24G -t 0-11:00 /bin/bash


#### INFRASTRUCTURE ####

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


# References
export pathRef=${cepko}/resources/ref
export pathGen=${pathRef}/genomes
export fileGen=${pathGen}/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa
export pathStarInd=${pathGen}/STAR-gg5-sjo49            # STAR indexed with --sjbdOverhang 49


# Tools and general scripts with versions in their respective paths 
export pathTools=${cepko}/resources/tools
export parallel=${pathTools}/parallel-20180722/src/parallel


# Modules loaded from O2
    # Need to load prereq modules first. Check for prereqs and command syntax with "module spider <tool name>".
module load gcc/6.2.0 python/2.7.12
module load fastqc/0.11.3
module load multiqc/1.5
module load trimmomatic/0.36    # Problematic syntax. Version number will need to be manually updated in script if program is updated.
module load star/2.5.4a
module load picard/2.8.0        # Problematic syntax. Version number will need to be manually updated in script if program is updated.



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


#echo "Collecting 1st pass junctions from 32 bp samples (S1-50) for multi-sample 2-pass STAR mapping"
#find ${pathOutStar}/S{1..50}/*SJ.out.tab -name ".*" -prune -o -print | paste -sd " " > ${pathOutStar}/sj_1-50.txt

#echo "Collecting 1st pass junctions from 50 bp samples (S51-60) for multi-sample 2-pass STAR mapping"
#find ${pathOutStar}/S{51..60}/*SJ.out.tab -name ".*" -prune -o -print | paste -sd " " > ${pathOutStar}/sj_51-60.txt


echo "Collecting 1st pass junctions from samples for multi-sample 2-pass STAR mapping"
find ${pathOutStar}/S{1..60}/*SJ.out.tab -name ".*" -prune -o -print | paste -sd " " > ${pathOutStar}/sj_1-60.txt


echo "Making directory for 2nd pass STAR samples"
mkdir ${pathOutStar}/pass2

