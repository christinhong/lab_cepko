#!/bin/bash

# Christin Hong
# Last updated: 2018-08
# Harvard Medical School, Connie Cepko Lab


#####

# Check input
echo
echo "Trimmomatic input: ${1}"
echo

# Naming output
outPre1=$(echo ${1} | cut -d " " -f1)
outPre2=${outPre1%.*}
outPre3=${outPre2%.*}

outPre="${outPre3}"_S"${SLURM_ARRAY_TASK_ID}_"
echo "Output FASTQ are prefixed with ${outPre}"
echo


# Running Trimmomatic
java -jar $TRIMMOMATIC/trimmomatic-0.36.jar PE \
	-threads ${intCores} \
	-phred33 \
	${1} \
	${pathData2}/S"${SLURM_ARRAY_TASK_ID}"/${outPre}R1_paired_trimmed.fq.gz \
	${pathData2}/S"${SLURM_ARRAY_TASK_ID}"/${outPre}R1_unpaired_trimmed.fq.gz \
	${pathData2}/S"${SLURM_ARRAY_TASK_ID}"/${outPre}R2_paired_trimmed.fq.gz \
	${pathData2}/S"${SLURM_ARRAY_TASK_ID}"/${outPre}R2_unpaired_trimmed.fq.gz \
	ILLUMINACLIP:adapters.fa:2:20:10:8:TRUE \
	LEADING:3 TRAILING:3 \
	SLIDINGWINDOW:4:5



