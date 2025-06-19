#!/bin/bash

dataset="$1"

if [[ "$dataset" == "m5C-R1" ]]; then 
    classification_mapping_path="/data/fass5/projects/hv_rna_mod/data/classification_mapping/Pool2_R1"
elif [[ "$dataset" == "m5C-R2" ]]; then
    classification_mapping_path="/data/fass5/projects/hv_rna_mod/data/classification_mapping/Pool2_R2"
elif [[ "$dataset" == "psU-R1" ]]; then
    classification_mapping_path="/data/fass5/projects/hv_rna_mod/data/classification_mapping/Pool1_R1"
else 
    echo "Usage: m5C-R1 or m5C-R2 or psU-R1"
    exit 1
fi 

all_bams=($(ls "$classification_mapping_path" | grep 'ID.*\.bam')) 

function filter_bams(){
    # filter bams for rc read ids 
    for bam in "${all_bams[@]}"; do 
        echo $bam
        samtools view -b -f 16 "$classification_mapping_path/$bam" > "$classification_mapping_path/rcs/rc_${bam}" 
    done
}

filter_bams
# samtools view rc.bam | awk '{print $2}' | sort | uniq -c
