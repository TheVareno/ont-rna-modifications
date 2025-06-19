#!/bin/bash 

raws_path='/data/fass5/projects/hv_rna_mod/data/raw/U_psU/'
can_mod_ref='/data/fass5/projects/hv_rna_mod/data/reference/all_5mers_single_T_3cycles.fa'
basecalled_data_path='/data/fass5/projects/hv_rna_mod/data/basecalled/psU/'
general_mapping='/data/fass5/projects/hv_rna_mod/data/general_mapping/'
basecaller_path="/data/fass5/projects/hv_rna_mod/data/basecaller/dorado-0.9.1-linux-x64/bin/"

#time ./dorado basecaller --emit-fastq ./rna002_70bps_hac@v3/ --device cuda:all  ../../../raw/psU-RNA_20201103_FAO12159.pod5 > ../../../basecalls/psU/psU-RNA_20201103_FAO12159.fastq 

function basecall(){

    model_mod="./rna004_130bps_sup@v5.1.0_pseU@v1"  
    model_can="./rna004_130bps_hac@v5.1.0"  
    
    raw_can='/data/fass5/projects/hv_rna_mod/data/raw/U_psU/U-RNA_20201103_FAO12153.pod5'
    raw_mod='/data/fass5/projects/hv_rna_mod/data/raw/U_psU/psU-RNA_20201103_FAO12159.pod5'
    file_name_can="${raw_can##*/}"
    file_name="${file_name_can%.*}"
    cd "$basecaller_path" && ./dorado basecaller --emit-fastq $model_can --device cuda:all $raw_can > "${basecalled_data_path}rna_sup_130bps_pseU_v1/${file_name}.fastq"  
    cd "$basecaller_path" && ./dorado basecaller $model_can --device cuda:all $raw_can > "${basecalled_data_path}rna_sup_130bps_pseU_v1/${file_name}.bam"
    
    file_name_mod="${raw_mod##*/}"
    file_name="${file_name_mod%.*}" 
    cd "$basecaller_path" && ./dorado basecaller  --emit-fastq  $model_can --device cuda:all $raw_mod > "${basecalled_data_path}rna_sup_130bps_pseU_v1/${file_name}.fastq"  
    cd "$basecaller_path" && ./dorado basecaller $model_can --device cuda:all $raw_mod > "${basecalled_data_path}rna_sup_130bps_pseU_v1/${file_name}.bam"
}

basecall

function convert_p2b(){

    raw_daten=($(ls "$path"*.pod5))  

    for raw_data in "${raw_daten[@]}"; do 
        output_file="${raw_data%.pod5}.blow5"
        blue-crab p2s "$raw_data" --output "$output_file"
    done
}


function map(){
    
    fastq_files=($(ls "$basecalled_data_path"*.fastq ))

    for fastq_file in "${fastq_files[@]}"; do
        minimap2 -ax map-ont -t 10 "$can_mod_ref" "$fastq_file" | samtools view -bo "${fastq_file%.fastq}.bam"           
    done

}


function sort() {
    echo "always sort then index!"
}


function index(){

    bam_files=($(ls "$general_mapping"*_sorted.bam)) 

    for bam_file in "${bam_files[@]}"; do
        samtools sort -o "$general_mapping$bam_file" "$general_mapping$bam_file"
        samtools index "$general_mapping$bam_file"
    done
}




