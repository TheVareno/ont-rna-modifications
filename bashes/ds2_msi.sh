#!/bin/bash
# map dataset 2 : R1, R2 with mod/can references 
# read -r -p "read set name ? " set 
# read -r -p "mod/can base ? " mode 

set="$1"
mode="$2"

if [[ -z "$set" || -z "$mode" ]]; then 
    echo "Usage: $0 <set> <mode>" 
    exit 1
fi 

refs_path="/data/fass5/projects/hv_rna_mod/data/reference/" 
raws_path="/home/hi68ren/Dokumente/MA/data/raw/" 
fastqs_path="/home/hi68ren/Dokumente/MA/data/called/fastqs/ds2/"

refs_mod=("ID1-GC_7k_1-5spacer_singleG_mod.fa" "ID2-GC_7k_1-5spacer_singleG_mod.fa" "ID3-GC_7k_1-5spacer_singleG_mod.fa" "ID4-GC_7k_1-5spacer_singleG_mod.fa")
refs_can=("ID1-GC_7k_1-5spacer_singleG_can.fa" "ID2-GC_7k_1-5spacer_singleG_can.fa" "ID3-GC_7k_1-5spacer_singleG_can.fa" "ID4-GC_7k_1-5spacer_singleG_can.fa") 

raw_R1="R1.blow5"
raw_R2="R2.blow5"

pool_R1=("$raw_R1" "R1.fastq")
pool_R2=("$raw_R2" "R2.fastq")

bam_files_all=($(ls | grep '\.bam'))

bam_files_R1=($(ls | grep '^R1.\.bam'))

bam_files_R2=($(ls | grep '^R2.\.bam'))

if [[ "$mode" == "mod" ]]; then 
    ref_array=("${refs_mod[@]}")  
elif [[ "$mode" == "can" ]]; then 
    ref_array=("${refs_mod[@]}") 
else 
    exit 1
fi 


if [[ "$set" == "R1" ]]; then 
    read_array=("${pool_R1[@]}")  
    bam_array=("${bam_files_R1[@]}")  
elif [[ "$set" == "R2" ]]; then 
    read_array=("${pool_R2[@]}") 
    bam_array=("${bam_files_R2[@]}")  
else 
    exit 1
fi 


function align_ds2() { 
    for ref in "${ref_array[@]}"; do 
        minimap2 -ax map-ont --rev-only "${refs_path}${ref}" "../../called/fastqs/ds2/${set}.fastq" > "./${set}_revcomp_${ref:0:3}_${mode}.bam"
    done 
} 

function sort_bams(){
    for bam in "${bam_files_all[@]}"; do 
        samtools sort -o "sorted-${bam}" "$bam" 
    done
}


bam_files_sorted=($(ls | grep 'sorted\.bam'))
mapfile -t bam_files_sorted < <(find . -maxdepth 1 -type f -name "sorted*.bam")

#echo "${#bams_files_sorted[@]}"

function index_bams(){
    for bam in "${bam_files_sorted[@]}"; do 
        echo "Indexing $bam..."
        samtools index "$bam"
    done
}

# R1 mod 
# R1 can 
# R2 mod 
# R2 can
function f5c_index(){
        ./../../../imp/tools/f5c-v1.5/f5c_x86_64_linux index -t 20 --slow5 "${raws_path}${read_array[0]}" "${fastqs_path}${read_array[1]}"   
}

function f5c_event_align(){

    for 
    ./../../../imp/tools/f5c-v1.5/f5c_x86_64_linux eventalign -t 40 --rna -r "${fastqs_path}${read_array[1]}" -b     

}

f5c_index    