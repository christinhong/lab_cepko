#### For test code ####
sbatch --array=1-60 ~/2018_chickRFZ_rnaSeq/RNA-02_starMap1_jobarray.sh

sbatch --array=9 ~/2018_chickRFZ_rnaSeq/RNA-02_starMap1_jobarray.sh



#### Old but interesting code snippets ####

#### Commenting out below ####
: <<'COMMENT'


COMMENT

#

# String splitting
INPUT='someletters_12345_moreleters.ext'
SUBSTRING=$(echo $INPUT| cut -d'_' -f 2)
echo $SUBSTRING

#

# Remove first (empty) line of readPair file
    # Helps with debugging because it keeps the error logs empty during successful runs.
    # from https://stackoverflow.com/questions/339483/how-can-i-remove-the-first-line-of-a-text-file-using-bash-sed-script
tail -n +2 "S"${SLURM_ARRAY_TASK_ID}"_readPairs.txt" > "S"${SLURM_ARRAY_TASK_ID}"_readPairs.tmp" \
    && mv "S"${SLURM_ARRAY_TASK_ID}"_readPairs.tmp" "S"${SLURM_ARRAY_TASK_ID}"_readPairs.txt"


#

#### Sort and index BAM file ####
cd ${pathOut}/S"${SLURM_ARRAY_TASK_ID}"


# Make vars and output names
outBAMsorted="${outPrefix}"_Aligned.out.sorted.bam
printf "%s\n" "sorted bam: ${outBAMsorted}"


# Sort and index with samtools
samtools sort -o ${outBAMsorted} ${outBAM}

samtools index ${outBAMsorted}

#

# Inputting all R1 and R2 reads into STAR for a single BAM
ls *_R1_*.fastq.gz | tr -s '\r\n' ',' | sed -e 's/,$/\n/' > r1.csv

ls *_R2_*.fastq.gz | tr -s '\r\n' ',' | sed -e 's/,$/\n/' > r2.csv


${parallel} -j ${intCPUs} --verbose --joblog "${pathLogs}/parallel-starP1.log" --resume-failed --keep-order 'sh ${pathBash}/RNA-02_func01StarPass1.sh' :::: r1.csv :::: r2.csv

#

# Make CSV list of R1 reads in sample folder
ls *_R1*.fastq.gz | tr -s '\r\n' ',' | sed -e 's/,$/\n/' > S"${SLURM_ARRAY_TASK_ID}"_r1.csv

# Make CSV list of R2 reads in sample folder
ls *_R2*.fastq.gz | tr -s '\r\n' ',' | sed -e 's/,$/\n/' > S"${SLURM_ARRAY_TASK_ID}"_r2.csv

#

# Join R1 and R2 files line by line with space between the lines
paste -d " " S"${SLURM_ARRAY_TASK_ID}"_r1.sorted.txt S"${SLURM_ARRAY_TASK_ID}"_r2.sorted.txt > S"${SLURM_ARRAY_TASK_ID}"_readPairs.txt

#

    # In for this filename, when splitting on _, lane identifier is the 3rd substring.
find . *.fastq.gz -name ".*" -prune -o -print | awk -F '_' '{print $3}' | sort | uniq > S"${SLURM_ARRAY_TASK_ID}"_lanes.txt


# #### Transfer to home dir ####
# # Check size (will this fit in my home dir?)
# du -sh /n/scratch2/ch220/starMap # 27 G, cool. Though for a 100 GB home dir, this is still getting pretty tight...maybe I should make the scratch2 folder my temp output folder. Hmm.


# # Copy output to project output folder (remember to exclude timestamps!)
# rsync -rv /n/scratch2/ch220/starMap /home/ch220/2018_chickRFZ_rnaSeq/output/ # Probably should've logged onto interactive node first.



#### Download MultiQC results to PC to look them over by eye ####
# From shell on HOME COMPUTER, run:
rsync -avr --progress "ch220@transfer.rc.hms.harvard.edu:/home/ch220/2018_chickRFZ_rnaSeq/doc/multiqc_*_2*" "/home/christin/aaa_data_CepkoLab/2018_chickRFZ_rnaSeq_data/multiQC/"





