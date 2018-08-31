#!/bin/bash

# Christin Hong
# Last updated: 2018-08
# Harvard Medical School, Connie Cepko Lab


#####

# Check input
echo "Picard input: ${1}"
echo

while IFS='_' read -ra bam || [[ -n "${bam}" ]]; do

    # Naming output
    outPre=$(echo ${bam[0]}_${bam[5]}_B${bam[4]}_${bam[2]})
    echo "Output is prefixed with ${outPre}"
    echo

    java -jar ${PICARD}/picard-2.8.0.jar CleanSam \
        I=${1} \
        O=${outPre}_p1clean.bam
        
    java -jar ${PICARD}/picard-2.8.0.jar AddOrReplaceReadGroups \
        I=${outPre}_p1clean.bam \
        O=${pathData3}/${outPre}_p2rg.bam \
        SORT_ORDER=coordinate \
        RGID=${outPre} \
        RGSM=${bam[0]} \
        RGPU=${bam[2]} \
        RGPL=illumina \
        RGLB=Library${bam[4]}

done <<< ${1}


