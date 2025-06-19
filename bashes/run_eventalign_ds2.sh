#!bin/bash 

# run f5c eventalign for all 


function run_eventalign_psU(){
    
    refs_path="/data/fass5/projects/hv_rna_mod/data/reference/"

    refs_AT_mod=($(ls $refs_path | grep -E ".AT.*mod\.fa$"))
    refs_AT_can=($(ls $refs_path | grep -E ".AT.*can\.fa$"))


    #./f5c_x86_64_linux index -t 20  --slow5 /data/fass5/projects/hv_rna_mod/data/raw/U_psU/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.blow5 /data/fass5/projects/hv_rna_mod/data/basecalled/psU/rna_hac_70bps/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.fastq 

    echo "eventalign for psU dataset 2..."

    ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna --min-mapq 0 \
            -r /data/fass5/projects/hv_rna_mod/data/basecalled/psU/rna_hac_70bps/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.fastq \
            -b /data/fass5/projects/hv_rna_mod/data/classification_mapping/Pool1_R1/rcs/rc_ID1-AT_7k_1-5spacer_singleA_can_sorted.bam \
            -g /data/fass5/projects/hv_rna_mod/data/reference/ID1-AT_7k_1-5spacer_singleA_can.fa \
            --slow5 /data/fass5/projects/hv_rna_mod/data/raw/U_psU/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.blow5 > /data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/ID1_can_new.tsv


    ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna --min-mapq 0 \
            -r /data/fass5/projects/hv_rna_mod/data/basecalled/psU/rna_hac_70bps/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.fastq \
            -b /data/fass5/projects/hv_rna_mod/data/classification_mapping/Pool1_R1/rcs/rc_ID1-AT_7k_1-5spacer_singleA_mod_sorted.bam \
            -g /data/fass5/projects/hv_rna_mod/data/reference/ID1-AT_7k_1-5spacer_singleA_mod.fa \
            --slow5 /data/fass5/projects/hv_rna_mod/data/raw/U_psU/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.blow5 > /data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/ID1_mod_new.tsv


    ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna --min-mapq 0 \
            -r /data/fass5/projects/hv_rna_mod/data/basecalled/psU/rna_hac_70bps/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.fastq \
            -b /data/fass5/projects/hv_rna_mod/data/classification_mapping/Pool1_R1/rcs/rc_ID2-AT_7k_1-5spacer_singleA_can_sorted.bam \
            -g /data/fass5/projects/hv_rna_mod/data/reference/ID2-AT_7k_1-5spacer_singleA_can.fa \
            --slow5 /data/fass5/projects/hv_rna_mod/data/raw/U_psU/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.blow5 > /data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/ID2_can_new.tsv


    ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna --min-mapq 0 \
            -r /data/fass5/projects/hv_rna_mod/data/basecalled/psU/rna_hac_70bps/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.fastq \
            -b /data/fass5/projects/hv_rna_mod/data/classification_mapping/Pool1_R1/rcs/rc_ID2-AT_7k_1-5spacer_singleA_mod_sorted.bam \
            -g /data/fass5/projects/hv_rna_mod/data/reference/ID2-AT_7k_1-5spacer_singleA_mod.fa \
            --slow5 /data/fass5/projects/hv_rna_mod/data/raw/U_psU/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.blow5 > /data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/ID2_mod_new.tsv


    ./f5c_x86_64_linux eventalign -t 40 \
        --signal-index --collapse-events --scale-events --print-read-names --rna --min-mapq 0 \
        -r /data/fass5/projects/hv_rna_mod/data/basecalled/psU/rna_hac_70bps/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.fastq \
        -b /data/fass5/projects/hv_rna_mod/data/classification_mapping/Pool1_R1/rcs/rc_ID3-AT_7k_1-5spacer_singleA_can_sorted.bam \
        -g /data/fass5/projects/hv_rna_mod/data/reference/ID3-AT_7k_1-5spacer_singleA_can.fa \
        --slow5 /data/fass5/projects/hv_rna_mod/data/raw/U_psU/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.blow5 > /data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/ID3_can_new.tsv


    ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna --min-mapq 0 \
            -r /data/fass5/projects/hv_rna_mod/data/basecalled/psU/rna_hac_70bps/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.fastq \
            -b /data/fass5/projects/hv_rna_mod/data/classification_mapping/Pool1_R1/rcs/rc_ID3-AT_7k_1-5spacer_singleA_mod_sorted.bam \
            -g /data/fass5/projects/hv_rna_mod/data/reference/ID3-AT_7k_1-5spacer_singleA_mod.fa \
            --slow5 /data/fass5/projects/hv_rna_mod/data/raw/U_psU/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.blow5 > /data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/ID3_mod_new.tsv

    ./f5c_x86_64_linux eventalign -t 40 \
        --signal-index --collapse-events --scale-events --print-read-names --rna --min-mapq 0 \
        -r /data/fass5/projects/hv_rna_mod/data/basecalled/psU/rna_hac_70bps/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.fastq \
        -b /data/fass5/projects/hv_rna_mod/data/classification_mapping/Pool1_R1/rcs/rc_ID4-AT_7k_1-5spacer_singleA_can_sorted.bam \
        -g /data/fass5/projects/hv_rna_mod/data/reference/ID4-AT_7k_1-5spacer_singleA_can.fa \
        --slow5 /data/fass5/projects/hv_rna_mod/data/raw/U_psU/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.blow5 > /data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/ID4_can_new.tsv

    ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna  --min-mapq 0 \
            -r /data/fass5/projects/hv_rna_mod/data/basecalled/psU/rna_hac_70bps/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.fastq \
            -b /data/fass5/projects/hv_rna_mod/data/classification_mapping/Pool1_R1/rcs/rc_ID4-AT_7k_1-5spacer_singleA_mod_sorted.bam \
            -g /data/fass5/projects/hv_rna_mod/data/reference/ID4-AT_7k_1-5spacer_singleA_mod.fa \
            --slow5 /data/fass5/projects/hv_rna_mod/data/raw/U_psU/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.blow5 > /data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/ID4_mod_new.tsv 

}

# psU ds1
run_eventalign_psU_ds1(){

    echo "eventalign for psU dataset 1..."
    ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna --min-mapq 0 \
            -r /data/fass5/projects/hv_rna_mod/data/basecalled/psU/rna_hac_70bps/psU-RNA_20201103_FAO12159.fastq \
            -b /data/fass5/projects/hv_rna_mod/data/general_mapping/psU-RNA_20201103_FAO12159_sort.bam\
            -g /data/fass5/projects/hv_rna_mod/data/reference/all_5mers_single_T_3cycles.fa \
            --slow5 /data/fass5/projects/hv_rna_mod/data/raw/U_psU/psU-RNA_20201103_FAO12159.blow5 > /data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/psU-RNA_20201103_FAO12159_new.tsv
}

# m5C ds2 
function run_eventalign(){

    refs_set=("ID1" "ID2" "ID3" "ID4")
    
    js_bams_path_R1="/data/fass5/projects/hv_rna_mod/data/classification_mapping/Pool2_R1/rcs/"
    js_bams_path_R2="/data/fass5/projects/hv_rna_mod/data/classification_mapping/Pool2_R2/rcs/"
    
    #cnt_list_bams_R1=($(ls "$js_bams_path_R1" | grep -E "ID1-GC_7k_1-5spacer_singleG_(can|mod)\.bam"))
    #cnt_list_bams_R2=($(ls "$js_bams_path_R1" | grep -E "ID1-GC_7k_1-5spacer_singleG_(can|mod)\.bam"))
    
    ID1_ref_mod="ID1-GC_7k_1-5spacer_singleG_mod.fa"
    ID2_ref_mod="ID2-GC_7k_1-5spacer_singleG_mod.fa"
    ID3_ref_mod="ID3-GC_7k_1-5spacer_singleG_mod.fa"
    ID4_ref_mod="ID4-GC_7k_1-5spacer_singleG_mod.fa"
    ID1_ref_can="ID1-GC_7k_1-5spacer_singleG_can.fa"
    ID2_ref_can="ID2-GC_7k_1-5spacer_singleG_can.fa"
    ID3_ref_can="ID3-GC_7k_1-5spacer_singleG_can.fa"
    ID4_ref_can="ID4-GC_7k_1-5spacer_singleG_can.fa"
    
    ID1_bam_mod="rc_ID1-GC_7k_1-5spacer_singleG_mod.bam"
    ID2_bam_mod="rc_ID2-GC_7k_1-5spacer_singleG_mod.bam"
    ID3_bam_mod="rc_ID3-GC_7k_1-5spacer_singleG_mod.bam"
    ID4_bam_mod="rc_ID4-GC_7k_1-5spacer_singleG_mod.bam"
    ID1_bam_can="rc_ID1-GC_7k_1-5spacer_singleG_can.bam"
    ID2_bam_can="rc_ID2-GC_7k_1-5spacer_singleG_can.bam"
    ID3_bam_can="rc_ID3-GC_7k_1-5spacer_singleG_can.bam"
    ID4_bam_can="rc_ID4-GC_7k_1-5spacer_singleG_can.bam"

    for data in "${daten[@]}"; do 
        
        if [[ "$data" == "R1" ]]; then

            ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna \
            -r "../../../data/called/fastqs/ds2/${data}.fastq" \
            -b "${js_bams_path_R1}${ID1_bam_mod}" \
            -g "${refs_path}${ID1_ref_mod}" \
            --slow5 "../../../data/raw/${data}.blow5" > "/home/hi68ren/Dokumente/MA/data/f5c/eventaligns/ds2/${data}_${ID1_ref_mod:0:3}_mod.tsv"

            ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna \
            -r "../../../data/called/fastqs/ds2/${data}.fastq" \
            -b "${js_bams_path_R1}${ID2_bam_mod}" \
            -g "${refs_path}${ID2_ref_mod}" \
            --slow5 "../../../data/raw/${data}.blow5" > "/home/hi68ren/Dokumente/MA/data/f5c/eventaligns/ds2/${data}_${ID2_ref_mod:0:3}_mod.tsv"

            ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna \
            -r "../../../data/called/fastqs/ds2/${data}.fastq" \
            -b "${js_bams_path_R1}${ID3_bam_mod}" \
            -g "${refs_path}${ID3_ref_mod}" \
            --slow5 "../../../data/raw/${data}.blow5" > "/home/hi68ren/Dokumente/MA/data/f5c/eventaligns/ds2/${data}_${ID3_ref_mod:0:3}_mod.tsv"

            ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna \
            -r "../../../data/called/fastqs/ds2/${data}.fastq" \
            -b "${js_bams_path_R1}${ID4_bam_mod}" \
            -g "${refs_path}${ID4_ref_mod}" \
            --slow5 "../../../data/raw/${data}.blow5" > "/home/hi68ren/Dokumente/MA/data/f5c/eventaligns/ds2/${data}_${ID4_ref_mod:0:3}_mod.tsv"

        elif [[ "$data" == "R2" ]]; then
             
            ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna \
            -r "../../../data/called/fastqs/ds2/${data}.fastq" \
            -b "${js_bams_path_R2}${ID1_bam_mod}" \
            -g "${refs_path}${ID1_ref_mod}" \
            --slow5 "../../../data/raw/${data}.blow5" > "/home/hi68ren/Dokumente/MA/data/f5c/eventaligns/ds2/${data}_${ID1_ref_mod:0:3}_mod.tsv"    
            
            ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna \
            -r "../../../data/called/fastqs/ds2/${data}.fastq" \
            -b "${js_bams_path_R2}${ID2_bam_mod}" \
            -g "${refs_path}${ID2_ref_mod}" \
            --slow5 "../../../data/raw/${data}.blow5" > "/home/hi68ren/Dokumente/MA/data/f5c/eventaligns/ds2/${data}_${ID2_ref_mod:0:3}_mod.tsv"    
            
            ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna \
            -r "../../../data/called/fastqs/ds2/${data}.fastq" \
            -b "${js_bams_path_R2}${ID3_bam_mod}" \
            -g "${refs_path}${ID3_ref_mod}" \
            --slow5 "../../../data/raw/${data}.blow5" > "/home/hi68ren/Dokumente/MA/data/f5c/eventaligns/ds2/${data}_${ID3_ref_mod:0:3}_mod.tsv"    
            
            ./f5c_x86_64_linux eventalign -t 40 \
            --signal-index --collapse-events --scale-events --print-read-names --rna \
            -r "../../../data/called/fastqs/ds2/${data}.fastq" \
            -b "${js_bams_path_R2}${ID4_bam_mod}" \
            -g "${refs_path}${ID4_ref_mod}" \
            --slow5 "../../../data/raw/${data}.blow5" > "/home/hi68ren/Dokumente/MA/data/f5c/eventaligns/ds2/${data}_${ID4_ref_mod:0:3}_mod.tsv"    
        
        fi
    done     
}



run_eventalign_psU_ds1

