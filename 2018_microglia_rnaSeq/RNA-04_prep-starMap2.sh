#### INTRO ####

# Christin M. Hong
# Last modified: 2018-10
# Harvard Medical School, Connie Cepko Lab



#### INFRASTRUCTURE ####

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




#### START ####


# Check for successful 1st pass mapping

cat ${pathLogs}/RNA-03*.err # Should return empty

grep -rl "ALL DONE!" ${pathOutStar}/S*/*_Log.out | wc -l # For 95 samples * 4 lanes, should return 380 (96 was our negative control -> not mapping it)


# For second pass: grep -rl "ALL DONE!" ${pathOutStar2}/S*/*_Log.out | wc -l

# To find files that do NOT have search term: grep -rL "ALL DONE!" ${pathOutStar2}/S*/*_Log.out



echo "Starting STAR alignment and mapping on $(date '+%Y-%m-%d %H:%M:%S')"

# Add newline to output for readability
echo


echo "Collecting 1st pass junctions from samples for multi-sample 2-pass STAR mapping"
find ${pathOutStar}/S{1..95}/*SJ.out.tab -name ".*" -prune -o -print | paste -sd " " > ${pathOutStar}/sj.txt


echo "Making directory for 2nd pass STAR samples"
mkdir ${pathOutStar}/pass2


# Upload next script(s) from local terminal:
# rsync -av --progress "/home/christin/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_microglia_rnaSeq/" "ch220@transfer.rc.hms.harvard.edu:/n/data2/hms/genetics/cepko/christin/2018_microglia_rnaSeq/"


