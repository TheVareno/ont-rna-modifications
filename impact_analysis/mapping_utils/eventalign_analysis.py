"""
author: Hadi Vareno
e-mail: mohammad.noori.vareno@uni-jena.de
github: https://github.com/TheVareno
"""

import numpy as np 
import pandas as pd # type: ignore   
from read5.Reader import read  # type: ignore
import scipy.stats as sp # type: ignore
from scipy.optimize import curve_fit # type: ignore # optimized gauss parameters finder!
from scipy.stats import median_abs_deviation as mad # type: ignore
from mapping_utils import bc_bams_analysis as bcbams # type: ignore 
from mapping_utils import guassian_utils as guass # type: ignore
import pysam # type: ignore
import argparse
import os

def setup_working_dir()-> None:
    os.chdir('/home/hi68ren/Dokumente/MA/imp/scripts/annotation')
    #os.chdir('/Users/hadivareno/Documents/Projects/MA/imp/scripts/annotation')

def find_motif_sliding_C(kmer: str)-> bool:
    """
    Recieves:
        - the pore model (kmer) 
    
    Performs:
        - filtering after the motif of C*DDDD / DDC*DD / DC*DDD / DDDC*D / DDDDC* 
    
    Return:
        - if the kmer has motif pattern C*DDDD / DDC*DD / DC*DDD / DDDC*D / DDDDC* (405 motifs)
    """
    try:    
        
        if kmer[0] == 'C' and kmer[1] != 'C' and kmer[2] != 'C' and kmer[3] != 'C' and kmer[4] != 'C': 
            return True 
        
        elif kmer[0] != 'C' and kmer[1] == 'C' and kmer[2] != 'C' and kmer[3] != 'C' and kmer[4] != 'C':
            return True
        
        elif kmer[0] != 'C' and kmer[1] != 'C' and kmer[2] == 'C' and kmer[3] != 'C' and kmer[4] != 'C':
            return True
        
        elif kmer[0] != 'C' and kmer[1] != 'C' and kmer[2] != 'C' and kmer[3] == 'C' and kmer[4] != 'C':
            return True
        
        elif kmer[0] != 'C' and kmer[1] != 'C' and kmer[2] != 'C' and kmer[3] != 'C' and kmer[4] == 'C':
            return True
        
        else: 
            return False
         
    except TypeError: 
        return False 


def find_motif_sliding_U(kmer: str)-> bool:
    try: 
        
        if kmer[0] == 'T' and kmer[1] != 'T' and kmer[2] != 'T' and kmer[3] != 'T' and kmer[4] != 'T': 
            return True 

        elif kmer[0] != 'T' and kmer[1] == 'T' and kmer[2] != 'T' and kmer[3] != 'T' and kmer[4] != 'T':
            return True

        elif kmer[0] != 'T' and kmer[1] != 'T' and kmer[2] == 'T' and kmer[3] != 'T' and kmer[4] != 'T':
            return True

        elif kmer[0] != 'T' and kmer[1] != 'T' and kmer[2] != 'T' and kmer[3] == 'T' and kmer[4] != 'T':
            return True

        elif kmer[0] != 'T' and kmer[1] != 'T' and kmer[2] != 'T' and kmer[3] != 'T' and kmer[4] == 'T':
            return True

        else: 
            return False
        
    except TypeError: 
        return False 
 

def find_GC_rich_motifs(kmer: str, threshold: int = 3) -> bool: 
    """
    Receives:
        - kmer: A 5-mer sequence (str)
        - threshold: Minimum number of G or C bases to qualify as GC-rich which is 3 
    
    Performs:
        - Counts occurrences of 'G' and 'C' in the kmer
        - Checks if the count meets or exceeds the threshold
    
    Returns:
        - True if the kmer is GC-rich, False otherwise
    """
    return len(kmer) == 5 and sum(base in {'G', 'C'} for base in kmer) >= threshold



def find_DDCDD_motifs(kmer: str)-> bool:
    """
    Recieves:
        - the pore model (kmer) 
    
    Performs:
        - filtering after the motif of DDCDD / DDC*DD
    
    Return:
        - if the kmer has motif pattern DDCDD / DDC*DD
    """
    return (
        len(kmer) == 5 and
        kmer[0] != 'C' and
        kmer[1] != 'C' and 
        kmer[2] == 'C' and      
        kmer[3] != 'C' and 
        kmer[4] != 'C' 
    ) 


def find_UUTUU_motifs(kmer: str)-> bool:
    """
    Recieves:
        - the pore model (kmer) 
    
    Performs:
        - filtering after the motif of DDCDD / DDC*DD
    
    Return:
        - if the kmer has motif pattern DDCDD / DDC*DD
    """
    return (
        len(kmer) == 5 and
        kmer[0] != 'T' and
        kmer[1] != 'T' and 
        kmer[2] == 'T' and 
        kmer[3] != 'T' and 
        kmer[4] != 'T' 
    ) 


def find_psU_pockets_comp_motifs(dorado_bam_path)-> bool:
    """
    - Conceptual idea to find the motifs that biologically are the target of riboprotein complex responsible for modification of U. 
    
    - The target is known as Pseudouridine pocket. 
    
    """
    samfile = pysam.AlignmentFile(dorado_bam_path, 'rb', check_sq=False)
    base_pairs = ['AT', 'TA', 'GC', 'CG']
    readnames_pocketmotifs = {'read_name': [], 'read_position': [], 'pocket_motifs': [] }
    
    for read in samfile.fetch(until_eof=True):
        read_name = read.query_name   
        read_position = read.reference_start
        read_sequece = str(read.query_sequence)  
        current_read_pocket_motifs = []
        for i, base in enumerate(read_sequece): 
            if base == 'T' and i >= 6 :
                base_1_before = read_sequece[i - 1]
                base_1_after = read_sequece[i + 1] 
                bp_1_around = base_1_before + base_1_after 
                if bp_1_around not in base_pairs:
                    base_2_before = read_sequece[i - 2]
                    base_2_after = read_sequece[i + 2] 
                    bp_2_around = base_2_before + base_2_after 
                    if bp_2_around not in base_pairs: 
                        base_3_before = read_sequece[i - 3]
                        base_3_after = read_sequece[i + 3] 
                        bp_3_around = base_3_before + base_3_after 
                        if bp_3_around not in base_pairs:
                            base_4_before = read_sequece[i - 4]
                            base_4_after = read_sequece[i + 4] 
                            bp_4_around = base_4_before + base_4_after
                            if bp_4_around not in base_pairs:
                                base_5_before = read_sequece[i - 5]
                                base_5_after = read_sequece[i + 5] 
                                bp_5_around = base_5_before + base_5_after
                                if bp_5_around not in base_pairs:
                                    pocket_motif = base_2_before + base_1_before + 'T' + base_1_after + base_2_after 
                                    current_read_pocket_motifs.append(pocket_motif) 
                                    
            readnames_pocketmotifs['read_name'].append(read_name)     
            readnames_pocketmotifs['read_position'].append(read_position)
            readnames_pocketmotifs['pocket_motifs'].append(current_read_pocket_motifs)                                                        

    return pd.DataFrame(readnames_pocketmotifs)
    

def filter_eventalign_psU_pockets_comp_motifs(eventalign_dataset_path: str, dorado_bam_path: str):
    """
    - Implementation
    """
    psU_pockets_comp_motifs = find_psU_pockets_comp_motifs(dorado_bam_path)  
    
    for i, chunk in enumerate(pd.read_csv(eventalign_dataset_path, sep='\t', chunksize=5000)):
        
        print(f'Proccessing the {i}. chunk.') if i % 1000 == 0 else None
                        
        eventalign_grouped_df = chunk.groupby(['read_name']) 
        
        for read_tuple, sub_df in eventalign_grouped_df: 

            for read_name in psU_pockets_comp_motifs['read_name']:
                
                if read_name == read_tuple[0]: 
                   
                   for kmer in sub_df['model_kmer']:
                       
                       ls_motifs_crnt_read_name = psU_pockets_comp_motifs['pocket_motifs'][psU_pockets_comp_motifs['read_name'] == read_name].tolist()
                       ref_pos = sub_df['position'][sub_df['read_name'] == read_name]
                       pocket_position = psU_pockets_comp_motifs['read_position'][psU_pockets_comp_motifs['read_name'] == read_name ]
                       if kmer in ls_motifs_crnt_read_name and ref_pos == pocket_position:                
                            raw_signal = read.getpASignal(read_tuple[0])
                            norm_signal = bcbams.normalize_signal(raw_signal, float(sub_df['mean'].iloc[0]), float(sub_df['std_dev'].iloc[0]))
            
            for _, row in sub_df.iterrows():  
                start_index, end_index = int(row['start_idx']), int(row['end_idx']) 
                sub_df['kmer'].append(row['model_kmer'])
                sub_df['signal_values'].append(norm_signal[start_index:end_index])



def extract_dwell_time_kmerwise(eventalign_dataset_path: str): 
    """
    Recieves: eventalign output dataset. 

    Performs: extracts dwell time based on event length from eventalign analysis

    Returns: pandas dataframe kmers with their correspoding dwell time  
    """
    kmer_dwell_time = {'kmer': [], 'dwell_times': []}

    for i, chunk in enumerate(pd.read_csv(eventalign_dataset_path, sep='\t', chunksize=5000)):
        print(f'Proccessing the {i}. chunk.') if i % 1000 == 0 else None
        
        motif_df = chunk[chunk['model_kmer'].apply(find_motif_sliding_U)]    

        for kmer, gp in motif_df.groupby('model_kmer'):
            kmer_dwell_time['kmer'].append(kmer) 
            dwell_times = gp['event_length'].tolist()
            kmer_dwell_time['dwell_times'].append(dwell_times)

    return pd.DataFrame(kmer_dwell_time)



#! main eventalign dataset analysis in chunks
def extract_normalised_signalvalues_kmerwise(eventalign_dataset_path: str,
                  raw_signal_path: str, 
                  bam_file_path: str)-> pd.DataFrame: 
    
    """
    Recieves: 
        - eventaling output dataset, 
        - raw signal data path, 
        - dorado output bam file,
        - mode either canonical or modified.
    
    Performs: 
        - making read object from Read5
        - calling analyse_dorodo_bam
        - calling quality filtering based on basecaller bam data
        - processing of eventalign dataset in chunkwise, for each chunck: 
            - filtering after motif 
            - inner merge of both filtered dataframes (once after bc quality, once after motif)
            - getting the corresponding raw signal value of each motif in eventalign datagframe -> calls read from Read5
            - normalization of signal value of each motif.
    
    Returns: 
        - pandas dataframe with DDCDD / DDC*DD kmer pattern with its corresponding normalized signal values. 
    """
    
    read_obj = read(raw_signal_path)
    bam_df = bcbams.analyse_dorado_bam(bam_file_path) 
    filtered_reads = bcbams.filter_low_quality_reads_C(bam_df)  

    kmer_signal = {'kmer': [], 'signal_values': []} 
    
    for i, chunk in enumerate(pd.read_csv(eventalign_dataset_path, sep='\t', chunksize=5000)):
        
        print(f'Proccessing the {i}. chunk.') if i % 2000 == 0 else None
    
        motifs_df = chunk[chunk['model_kmer'].apply(find_GC_rich_motifs)]
        
        merged_df = motifs_df.merge(filtered_reads, on='read_name', how='inner')                
            
        grouped_df = merged_df.groupby(['read_name']) 
        
        for read_tuple, sub_df in grouped_df: 
           
            raw_signal = read_obj.getpASignal(read_tuple[0])
            norm_signal = bcbams.normalize_signal(raw_signal, float(sub_df['mean'].iloc[0]), float(sub_df['std_dev'].iloc[0]))
            
            for _, row in sub_df.iterrows():  
                start_index, end_index = int(row['start_idx']), int(row['end_idx']) 
                kmer_signal['kmer'].append(row['model_kmer'])
                #ls_norm_signal = norm_signal[start_index:end_index]
                #str_norm_signal = ','.join(str(x) for x in ls_norm_signal)
                kmer_signal['signal_values'].append(norm_signal[start_index:end_index])
                #string = ','.join(str(x) for x in arr)

    return pd.DataFrame(kmer_signal)
    #kmer_signal_df = pd.DataFrame(kmer_signal)
    #kmer_signal_df.to_csv(f'../../data/annotate/ds1/{mode}_kmer_normalized_signal_values.csv')
    

#! psU only func
def extract_normalised_signalvalues_kmerwise_NO_CHUNK(eventalign_dataset_path: str,
                  raw_signal_path: str, 
                  bam_file_path: str)-> pd.DataFrame: 

    """
    Same as above but for light eventalign dataset
    """
    
    read_obj = read(raw_signal_path)
    bam_df = bcbams.analyse_dorado_bam(bam_file_path) 
    filtered_reads = bcbams.filter_low_quality_reads_psU(bam_df)  
        
    kmer_signal = {'kmer': [], 'signal_values': []} 
    ea_df = pd.read_csv(eventalign_dataset_path, sep='\t')
    
    ea_df = ea_df[ea_df['model_kmer'].apply(find_motif_sliding_U)]
        
    merged_df = ea_df.merge(filtered_reads, on='read_name', how='inner')                

    grouped_df = merged_df.groupby(['read_name']) 
    
    for read_name, sub_df in grouped_df: 
        raw_signal = read_obj.getpASignal(read_name[0])
        norm_signal = bcbams.normalize_signal(raw_signal, float(sub_df['mean'].iloc[0]), float(sub_df['std_dev'].iloc[0]))
        
        for _, row in sub_df.iterrows():  
            start_index, end_index = int(row['start_idx']), int(row['end_idx']) 
            kmer_signal['kmer'].append(row['model_kmer'])
            kmer_signal['signal_values'].append(norm_signal[start_index:end_index])
                
    return pd.DataFrame(kmer_signal)


def fit_guass(kmer_motif: str, bunch_can: np.ndarray, bunch_mod: np.ndarray)-> dict:
    """
    Recieves: 
        - motif of kmer on from DDCDD 
        - concatenated normalized signal values of DDCDD kmers as np array  
        - concatenated normalized signal values of DDC*DD kmers as np array 
        
    Performs: 
        - fit gaussian distribution to each recieved bunch of signal values : DDCDD / DDC*DD
        
    Returns: 
        - parameters of fitted gauss : mean and std with motif bases plus postion on reference sequence 
    """
    
    fit_info = {'kmer': kmer_motif, 
              'can_params': [], 
              'mod_params': []}
    
    try:
        #! can    
        hist_can, bin_edges_can = np.histogram(bunch_can, bins=50, density=True)
        bin_centers_can = (bin_edges_can[:-1] + bin_edges_can[1:]) / 2  
        mu_init_can, sd_init_can, amplitude_init_can = np.mean(bunch_can), np.std(bunch_can), np.max(hist_can)
        popt_can, _ = curve_fit(guass.gaussian, bin_centers_can, hist_can, p0=[mu_init_can, sd_init_can, amplitude_init_can]) #   bounds=([0, 0, -np.inf], [10, 5, np.inf]) 
        mu_can, sd_can, _ = popt_can
        sk_can = guass.cal_skewness(bunch_can, mu_can, sd_can)
        ku_can = guass.cal_kurtosis(bunch_can, mu_can, sd_can)
        
        fit_info['can_params'] = [mu_can, sd_can, sk_can, ku_can]
    
        #! mod
        hist_mod, bin_edges_mod = np.histogram(bunch_mod, bins=50, density=True)
        bin_centers_mod = (bin_edges_mod[:-1] + bin_edges_mod[1:]) / 2  
        mu_init_mod, sd_init_mod, amplitude_init_mod = np.mean(bunch_mod), np.std(bunch_mod), np.max(hist_mod)
        popt_mod, _ = curve_fit(guass.gaussian, bin_centers_mod, hist_mod, p0=[mu_init_mod, sd_init_mod, amplitude_init_mod])
        mu_mod, sd_mod, _ = popt_mod
        sk_mod = guass.cal_skewness(bunch_mod, mu_mod, sd_mod)
        ku_mod = guass.cal_kurtosis(bunch_mod, mu_mod, sd_mod)
        fit_info['mod_params'] = [mu_mod, sd_mod, sk_mod, ku_mod]
        
        return fit_info
    
    except RuntimeError: 
        
        fit_info['can_params'] = [0, 0, 0, 0]
        fit_info['mod_params'] = [0, 0, 0, 0]
        print('no fit info! runtime error')
        return fit_info 
    

    
def make_context_table(fit_info: dict, context_tb: dict, mode: str)-> dict: 
    """
    Recieves: 
        - fitted gauss parametes as dict object 
        - empty context table
    
    Performs: 
        - filling the context table : placing each base on entry from 1 to 5 except 3 with corresponding gauss parameters for Can and Mod  
    
    Returns: 
        - filled context-wise gaussian parameters
    """    

    context_tb['base_1'].append(fit_info['kmer'][0])
    context_tb['base_2'].append(fit_info['kmer'][1])
    context_tb['base_3'].append(fit_info['kmer'][2])    
    context_tb['base_4'].append(fit_info['kmer'][3])
    context_tb['base_5'].append(fit_info['kmer'][4])
    context_tb['mu_can'].append(fit_info['can_params'][0])
    context_tb['mu_mod'].append(fit_info['mod_params'][0])
    context_tb['sd_can'].append(fit_info['can_params'][1])
    context_tb['sd_mod'].append(fit_info['mod_params'][1])
    context_tb['sk_can'].append(fit_info['can_params'][2])
    context_tb['sk_mod'].append(fit_info['mod_params'][2])
    context_tb['ku_can'].append(fit_info['can_params'][3])
    context_tb['ku_mod'].append(fit_info['mod_params'][3])
    context_tb['med_can'].append(fit_info['can_params'][4])
    context_tb['med_mod'].append(fit_info['mod_params'][4])
    context_tb['mad_can'].append(fit_info['can_params'][5])
    context_tb['mad_mod'].append(fit_info['mod_params'][5])
    if mode == 'dwell_time':
        context_tb['mean_can'].append(fit_info['can_params'][6])
        context_tb['mean_mod'].append(fit_info['mod_params'][6])
        context_tb['std_can'].append(fit_info['can_params'][7])
        context_tb['std_mod'].append(fit_info['mod_params'][7])
        context_tb['min_can'].append(fit_info['can_params'][8])
        context_tb['min_mod'].append(fit_info['mod_params'][8])
        context_tb['max_can'].append(fit_info['can_params'][9])
        context_tb['max_mod'].append(fit_info['mod_params'][9])
    

    
def concatenate_fit_dwell_times(kmer_dwell_time_can: pd.DataFrame, kmer_dwell_time_mod: pd.DataFrame, modi) -> pd.DataFrame: 
    
    filtered_df_can = kmer_dwell_time_can[kmer_dwell_time_can['kmer'].map(kmer_dwell_time_can['kmer'].value_counts()) > 30] 
    sampled_df_can = (filtered_df_can.groupby('kmer').sample(n=30, replace=False).reset_index(drop=True))
    
    if modi == 'm5C':
        filtered_df_mod = kmer_dwell_time_mod[kmer_dwell_time_mod['kmer'].map(kmer_dwell_time_mod['kmer'].value_counts()) > 30]
        sampled_df_mod = (filtered_df_mod.groupby('kmer').sample(n=30, replace=False).reset_index(drop=True))
    
    if modi == 'psU': 
        filtered_df_mod = kmer_dwell_time_mod[kmer_dwell_time_mod['kmer'].map(kmer_dwell_time_mod['kmer'].value_counts()) > 0]
        sampled_df_mod = (filtered_df_mod.groupby('kmer').sample(n=1, replace=False).reset_index(drop=True))
        
    
    context_tb = {'base_1': [], 'base_2': [], 'base_3': [], 'base_4': [], 'base_5': [], 'mu_can': [], 'mu_mod': [], 'sd_can': [], 'sd_mod': [], 
                'sk_can': [], 'sk_mod': [], 'ku_can': [], 'ku_mod': [], 'med_can': [], 'med_mod': [], 'mad_can': [], 'mad_mod':[], 
                'mean_can': [], 'mean_mod': [], 'std_can': [], 'std_mod': [], 'min_can': [], 'min_mod': [], 'max_can': [], 'max_mod': [] }    
    
    for kmer_mod, sub_df_mod in sampled_df_mod.groupby('kmer'): 
        
        for kmer_can, sub_df_can in sampled_df_can.groupby('kmer'): 
                    
            if kmer_can == kmer_mod : 
                
                print(f'5mer : {kmer_can.upper()}...')
                
                bunch_mod = np.concatenate(sub_df_mod['dwell_times'].tolist())  
                med_mod = np.median(bunch_mod) 
                mean_mod = np.mean(bunch_mod) 
                std_mod = np.std(bunch_mod) 
                std_mod = np.std(bunch_mod) 
                min_mod = np.min(bunch_mod) 
                max_mod = np.max(bunch_mod) 
                mad_mod = mad(bunch_mod)

                bunch_can = np.concatenate(sub_df_can['dwell_times'].tolist()) 
                med_can = np.median(bunch_can) 
                mean_can = np.mean(bunch_can) 
                mad_can = mad(bunch_can)
                std_can = np.std(bunch_can) 
                std_can = np.std(bunch_can) 
                min_can = np.min(bunch_can) 
                max_can = np.max(bunch_can) 
                
                gauss_fit_info = fit_guass(kmer_motif=kmer_can, bunch_mod=bunch_mod, bunch_can=bunch_can)
                                      
                gauss_fit_info['can_params'].append(med_can)
                gauss_fit_info['can_params'].append(mad_can)
                gauss_fit_info['can_params'].append(mean_can)
                gauss_fit_info['can_params'].append(std_can)
                gauss_fit_info['can_params'].append(min_can)
                gauss_fit_info['can_params'].append(max_can)
                
                gauss_fit_info['mod_params'].append(med_mod)
                gauss_fit_info['mod_params'].append(mad_mod)
                gauss_fit_info['mod_params'].append(mean_mod)
                gauss_fit_info['mod_params'].append(std_mod)
                gauss_fit_info['mod_params'].append(min_mod)
                gauss_fit_info['mod_params'].append(max_mod)
                
                make_context_table(gauss_fit_info, context_tb, 'dwell_time')
                
    return pd.DataFrame(context_tb, columns=context_tb.keys())            


def concatenate_fit_signal_values(kmer_sigval_can: pd.DataFrame, kmer_sigval_mod: pd.DataFrame, modi: str) -> pd.DataFrame:                      
    """
    Recieves: 
        - the dataframe of kmers and corresponding normalized signal values for Can
        - the dataframe of kmers and corresponding normalized signal values for Mod
    
    Performs: 
        - filtering each dataframe after kmers where value counts are above 100 
        - taking samples for of 100 kmer-signalvalue pairs 
        - concatenation of samples in one array 
        - fitting the concatenated array to gauss 
        - reporting the gauss parameters (1. and 2. moments) 
        
    Returns: 
        - pandas dataframe for pattern discovery 
    """
    
    filtered_df_can = kmer_sigval_can[kmer_sigval_can['kmer'].map(kmer_sigval_can['kmer'].value_counts()) > 30] 
    sampled_df_can = (filtered_df_can.groupby('kmer').sample(n=30, replace=False).reset_index(drop=True))
    
    if modi == 'm5C': 
        filtered_df_mod = kmer_sigval_mod[kmer_sigval_mod['kmer'].map(kmer_sigval_mod['kmer'].value_counts()) > 30]
        sampled_df_mod = (filtered_df_mod.groupby('kmer').sample(n=30, replace=False).reset_index(drop=True))
    
    if modi == 'psU': 
        filtered_df_mod = kmer_sigval_mod[kmer_sigval_mod['kmer'].map(kmer_sigval_mod['kmer'].value_counts()) > 0] 
        sampled_df_mod = (filtered_df_mod.groupby('kmer').sample(n=1, replace=False).reset_index(drop=True))

    
    context_tb = {'base_1': [], 'base_2': [], 'base_3': [], 'base_4': [], 'base_5': [], 
                'mu_can': [], 'mu_mod': [], 'sd_can': [], 'sd_mod': [], 
                'sk_can': [], 'sk_mod': [], 'ku_can': [], 'ku_mod': [], 
                'med_can': [], 'med_mod': [], 'mad_can': [], 'mad_mod':[]}    
    
    #sampled_df_m5c['signal_values'] = sampled_df_m5c['signal_values'].apply(ast.literal_eval)
    #sampled_df_c['signal_values'] = sampled_df_c['signal_values'].apply(ast.literal_eval)

    for kmer_mod, sub_df_mod in sampled_df_mod.groupby('kmer'): 
        
        for kmer_can, sub_df_can in sampled_df_can.groupby('kmer'): 
                    
            if kmer_can == kmer_mod : 
                
                print(f'5mer : {kmer_can.upper()}...')
                
                bunch_mod = np.concatenate(sub_df_mod['signal_values'].tolist())  
                print(len(bunch_mod))
                med_mod = np.median(bunch_mod) 
                mad_mod = mad(bunch_mod)

                bunch_can = np.concatenate(sub_df_can['signal_values'].tolist()) 
                med_can = np.median(bunch_can) 
                mad_can = mad(bunch_can)

                gauss_fit_info = fit_guass(kmer_motif=kmer_can, bunch_mod=bunch_mod, bunch_can=bunch_can)
                
                if len(gauss_fit_info) != 0:                      
                    gauss_fit_info['can_params'].append(med_can)
                    gauss_fit_info['can_params'].append(mad_can)
                    gauss_fit_info['mod_params'].append(med_mod)
                    gauss_fit_info['mod_params'].append(mad_mod)
                    make_context_table(gauss_fit_info, context_tb, 'signal_value')
                else:
                    pass
                
    return pd.DataFrame(context_tb, columns=context_tb.keys())            


def make_gaussian_parameters_kmerwise_table(eventalign_dataset_path: str, raw: str, dorado_bam: str, modification: str)-> pd.DataFrame:
    
    kmer_gauss_params = {'kmer': [], 'mu': [], 'sigma':[]}
    
    if modification == 'm5C': 
        m5C_kmer_signal_values_df = extract_normalised_signalvalues_kmerwise(eventalign_dataset_path, raw, dorado_bam)  
        m5C_kmer_signal_values_df_grouped_kmerwise = m5C_kmer_signal_values_df.groupby(['kmer'])
        for kmer_index, sub_gp in m5C_kmer_signal_values_df_grouped_kmerwise:
            try: 
                current_kmer_all_signal_values = np.concatenate(sub_gp['signal_values'].tolist())
                hist, bin_edges = np.histogram(current_kmer_all_signal_values, bins=50, density=True)
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  
                mu_init, sd_init, amplitude_init = np.mean(current_kmer_all_signal_values), np.std(current_kmer_all_signal_values), np.max(hist)
                popt_can, _ = curve_fit(guass.gaussian, bin_centers, hist, p0=[mu_init, sd_init, amplitude_init]) # bounds=([0, 0, -np.inf], [10, 5, np.inf]) 
                mu, sd, _ = popt_can
                kmer_gauss_params['kmer'].append(kmer_index[0])        
                kmer_gauss_params['mu'].append(mu)
                kmer_gauss_params['sigma'].append(sd)
            
            except RuntimeError: 
                print(f"Error - curve_fit failed for kmer: {kmer_index[0]}")
            
    if modification == 'psU': 
        psU_kmer_sigval = extract_normalised_signalvalues_kmerwise_NO_CHUNK(eventalign_dataset_path, raw, dorado_bam)
        psU_kmer_signal_values_df_grouped_kmerwise = psU_kmer_sigval.groupby(['kmer'])
        for kmer_index, sub_gp in psU_kmer_signal_values_df_grouped_kmerwise:
            try: 
                current_kmer_all_signal_values = np.concatenate(sub_gp['signal_values'].tolist())
                hist, bin_edges = np.histogram(current_kmer_all_signal_values, bins=50, density=True)
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  
                mu_init, sd_init, amplitude_init = np.mean(current_kmer_all_signal_values), np.std(current_kmer_all_signal_values), np.max(hist)
                popt_can, _ = curve_fit(guass.gaussian, bin_centers, hist, p0=[mu_init, sd_init, amplitude_init]) # bounds=([0, 0, -np.inf], [10, 5, np.inf]) 
                mu, sd, _ = popt_can
                kmer_gauss_params['kmer'].append(kmer_index[0])        
                kmer_gauss_params['mu'].append(mu)
                kmer_gauss_params['sigma'].append(sd)
            
            except RuntimeError: 
                print(f"Error - curve_fit failed for kmer: {kmer_index[0]}")
                
    return pd.DataFrame(kmer_gauss_params)    


def main()-> None: 
    
    parser = argparse.ArgumentParser(description="Process and Save output file.")
    parser.add_argument("--modi", type=str, choices=['m5C', 'psU'], help="Modification type")
    parser.add_argument("--ds", type=str, choices=['ds1', 'ds2'], help="Dataset type")
    parser.add_argument("--feat", type=str, choices=['signal-value', 'dwell-time'], help="Analysis type")

    args = parser.parse_args()
    
    setup_working_dir()
    
    #! m⁵C - ds1 
    
    c_ea = '../../../data/m5C/f5c/eventaligns/ds1/c_idxs_readids.tsv' 
    c_raw = '../../../data/m5C/raw/C-RNA_20200824_ACT688.pod5' 
    c_called = '../../../data/m5C/called/bams/ds1/c/c-rna-5mer.bam'
    
    m5c_ea = '../../../data/m5C/f5c/eventaligns/ds1/m5c_idxs_readids.tsv' 
    m5c_raw = '../../../data/m5C/raw/m5C-RNA_20200803_ACT409.pod5' 
    m5c_called = '../../../data/m5C/called/bams/ds1/m5c/m5c-rna-5mer.bam'

    #! m⁵C - ds2 
    ea_R1_mod = '../../../data/m5C/f5c/eventaligns/ds2/R1_events_mod.tsv'
    ea_R1_can = '../../../data/m5C/f5c/eventaligns/ds2/R1_events_can.tsv'
    ea_R2_mod = '../../../data/m5C/f5c/eventaligns/ds2/R2_events_mod.tsv'
    ea_R2_can = '../../../data/m5C/f5c/eventaligns/ds2/R2_events_can.tsv'

    raw_R1 = '../../../data/m5C/raw/RNAmod_20230826_Pool2_R1_FAW83677_LGWG.pod5'
    raw_R2 = '../../../data/m5C/raw/RNAmod_20230827_Pool2_R2_FAW83677_LGWG.pod5'

    bam_R1 = '../../../data/m5C/called/bams/ds2/R1.bam'
    bam_R2 = '../../../data/m5C/called/bams/ds2/R2.bam'

    
    if args.modi == 'm5C' and args.ds == 'ds1' and args.feat == 'signal-value': 
        
        c_kmer_signal_values_df = extract_normalised_signalvalues_kmerwise(c_ea, c_raw, c_called)  
        m5c_kmer_signal_values_df = extract_normalised_signalvalues_kmerwise(m5c_ea, m5c_raw, m5c_called) 
        context_tb = concatenate_fit_signal_values(c_kmer_signal_values_df, m5c_kmer_signal_values_df, args.modi)                    
        context_tb.to_csv('/data/fass5/projects/hv_rna_mod/data/ea_outputs/m5C/ds1_CG_rich_motifs_gauss_info_4_moments_signal_value.csv', header=True, index=True)
    
    if args.modi == 'm5C' and args.ds == 'ds1' and args.feat == 'dwell-time': 
        
        c_kmer_dwell_time_df = extract_dwell_time_kmerwise(c_ea)
        m5c_kmer_dwell_time_df = extract_dwell_time_kmerwise(m5c_ea)
        context_tb_dwell_time = concatenate_fit_dwell_times(c_kmer_dwell_time_df, m5c_kmer_dwell_time_df, args.modi)
        context_tb_dwell_time.to_csv('/data/fass5/projects/hv_rna_mod/data/ea_outputs/m5C/ds1_sliding_C_gauss_info_4_moments_dwell_time.csv', header=True, index=True)
    
    if args.modi == 'm5C' and args.ds == 'ds2' and args.feat == 'signal-value': 
    
        m5c_R1_kmer_sigval_df = extract_normalised_signalvalues_kmerwise(ea_R1_mod, raw_R1, bam_R1) 
        m5c_R2_kmer_sigval_df = extract_normalised_signalvalues_kmerwise(ea_R2_mod, raw_R2, bam_R2)
        c_R1_kmer_sigval_df = extract_normalised_signalvalues_kmerwise(ea_R1_can, raw_R1, bam_R1) 
        c_R2_kmer_sigval_df = extract_normalised_signalvalues_kmerwise(ea_R2_can, raw_R2, bam_R2) 
        c_kmer_sigval_df = pd.concat([c_R1_kmer_sigval_df, c_R2_kmer_sigval_df])
        m5c_kmer_sigval_df = pd.concat([m5c_R1_kmer_sigval_df, m5c_R2_kmer_sigval_df])
        context_df = concatenate_fit_signal_values(c_kmer_sigval_df, m5c_kmer_sigval_df, args.modi) 
        context_df.to_csv('/data/fass5/projects/hv_rna_mod/data/ea_outputs/m5C/ds2_CG_rich_motifs_gauss_info_4_moments_signal_value.csv', header=True, index=True)

        
    if args.modi == 'm5C' and args.ds == 'ds2' and args.feat == 'dwell-time':     
        
        ea_mod = pd.concat([pd.read_csv(ea_R1_mod, sep='\t'), pd.read_csv(ea_R2_mod, sep='\t')])
        ea_can = pd.concat([pd.read_csv(ea_R1_can, sep='\t'), pd.read_csv(ea_R2_can, sep='\t')])
        ea_can.to_csv('/data/fass5/projects/hv_rna_mod/data/ea_outputs/m5C/concatenated_can_R1_R2_eventaligns.csv')
        ea_mod.to_csv('/data/fass5/projects/hv_rna_mod/data/ea_outputs/m5C/concatenated_mod_R1_R2_eventaligns.csv')
        c_kmer_dwell_time_df = extract_dwell_time_kmerwise('/data/fass5/projects/hv_rna_mod/data/ea_outputs/m5C/concatenated_can_R1_R2_eventaligns.csv')
        m5c_kmer_dwell_time_df = extract_dwell_time_kmerwise('/data/fass5/projects/hv_rna_mod/data/ea_outputs/m5C/concatenated_mod_R1_R2_eventaligns.csv')
        context_tb_dwell_time = concatenate_fit_dwell_times(c_kmer_dwell_time_df, m5c_kmer_dwell_time_df, args.modi)
        context_tb_dwell_time.to_csv('/data/fass5/projects/hv_rna_mod/data/ea_outputs/m5C/ds2_sliding_C_gauss_info_4_moments_dwell_time.csv', header=True, index=True)

    
    #! psU - ds1 
    U_map_bam = '/data/fass5/projects/hv_rna_mod/data/general_mapping/U-RNA_20201103_FAO12153.bam'
    psU_map_bam = '/data/fass5/projects/hv_rna_mod/data/general_mapping/psU-RNA_20201103_FAO12159.bam'
    
    U_dorado_bam = '/data/fass5/projects/hv_rna_mod/data/basecalled/psU/rna_hac_70bps/U-RNA_20201103_FAO12153.bam'
    psU_dorado_bam = '/data/fass5/projects/hv_rna_mod/data/basecalled/psU/rna_hac_70bps/psU-RNA_20201103_FAO12159.bam'
    
    U_raw = '/data/fass5/projects/hv_rna_mod/data/raw/U_psU/U-RNA_20201103_FAO12153.pod5'
    psU_raw = '/data/fass5/projects/hv_rna_mod/data/raw/U_psU/psU-RNA_20201103_FAO12159.pod5'
    
    U_ea = '/data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/U-RNA_20201103_FAO12153.tsv'
    psU_ea = '/data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/psU-RNA_20201103_FAO12159.tsv'
    
    #! psU - ds2 
    ds2_raw = '/data/fass5/projects/hv_rna_mod/data/raw/U_psU/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.blow5'
    ds2_dorado_bam = '/data/fass5/projects/hv_rna_mod/data/basecalled/psU/rna_hac_70bps/RNAmod_20230826_Pool1_R1_FAW77078_LGWG.bam'
    
    """
    mod_tsvs = []
    can_tsvs = []
    for ref in 'ID1 ID2 ID3 ID4'.split(' '):
        mod_df = pd.read_csv(f'/data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/{ref}_mod.tsv', sep='\t') 
        can_df = pd.read_csv(f'/data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/{ref}_can.tsv', sep='\t') 
        mod_tsvs.append(mod_df)
        can_tsvs.append(can_df)
    
    concant_mod_df = pd.concat(mod_tsvs)
    concant_mod_df.to_csv('/data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/concat_mod_ds2.tsv', sep='\t', header=True)
    concant_can_df = pd.concat(can_tsvs)
    concant_can_df.to_csv('/data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/concat_can_ds2.tsv', sep='\t', header=True)
    """
    concat_can_df = '/data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/concat_can_ds2.tsv'
    concat_mod_df = '/data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/concat_mod_ds2.tsv'
    
    if args.modi == 'psU' and args.ds == 'ds1' and args.feat == 'signal-value':     
        U_kmer_sigval = extract_normalised_signalvalues_kmerwise(U_ea, U_raw, U_dorado_bam)
        psU_kmer_sigval = extract_normalised_signalvalues_kmerwise_NO_CHUNK(psU_ea, psU_raw, psU_dorado_bam)
        cntx_df = concatenate_fit_signal_values(U_kmer_sigval, psU_kmer_sigval, args.modi) 
        cntx_df.to_csv('/data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/cntx_dfs/ds1_sliding_U_gauss_info_4_moments_signal_value.csv', header=True, index=False)                               
    
    if args.modi == 'psU' and args.ds == 'ds1' and args.feat == 'dwell-time':         
        U_kmer_dwell_time_df = extract_dwell_time_kmerwise(U_ea)
        psU_kmer_dwell_time_df = extract_dwell_time_kmerwise(psU_ea)
        context_tb_dwell_time = concatenate_fit_dwell_times(U_kmer_dwell_time_df, psU_kmer_dwell_time_df, args.modi)
        context_tb_dwell_time.to_csv('/data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/cntx_dfs/ds1_sliding_U_gauss_info_4_moments_dwell_time.csv', header=True, index=True)
    
    if args.modi == 'psU' and args.ds == 'ds2' and args.feat == 'dwell-time':         
        U_kmer_dwell_time_df = extract_dwell_time_kmerwise(concat_can_df)
        psU_kmer_dwell_time_df = extract_dwell_time_kmerwise(concat_mod_df)
        context_tb_dwell_time = concatenate_fit_dwell_times(U_kmer_dwell_time_df, psU_kmer_dwell_time_df, args.modi)
        context_tb_dwell_time.to_csv('/data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/cntx_dfs/ds2_sliding_U_gauss_info_4_moments_dwell_time.csv', header=True, index=True)

    if args.modi == 'psU' and args.ds == 'ds2' and args.feat == 'signal-value':         
        U_kmer_sigval_ds2 = extract_normalised_signalvalues_kmerwise(concat_can_df, ds2_raw, ds2_dorado_bam) 
        psU_kmer_sigval_ds2 = extract_normalised_signalvalues_kmerwise_NO_CHUNK(concat_mod_df, ds2_raw, ds2_dorado_bam)
        cntx_df = concatenate_fit_signal_values(U_kmer_sigval_ds2, psU_kmer_sigval_ds2, args.modi)
        cntx_df.to_csv('/data/fass5/projects/hv_rna_mod/data/ea_outputs/psU/cntx_dfs/ds2_sliding_U_gauss_info_4_moments_signal_value.csv', header=True, index=False)                               

    
    """
    print('making gauss table m⁵C...')
    
    kmer_gauss_params_ds1 = make_gaussian_parameters_kmerwise_table(m5c_ea, m5c_raw, m5c_called, 'm5C')
    kmer_gauss_params_ds1.to_csv('/data/fass5/projects/hv_rna_mod/data/mod_gausses/m5C_ds1.csv')
    
    kmer_gauss_params_ds2_R1 = make_gaussian_parameters_kmerwise_table(ea_R1_mod, raw_R1, bam_R1, 'm5C')
    kmer_gauss_params_ds2_R2 = make_gaussian_parameters_kmerwise_table(ea_R2_mod, raw_R2, bam_R2, 'm5C')
    kmer_gauss_params_ds2_R1.to_csv('/data/fass5/projects/hv_rna_mod/data/mod_gausses/m5C_ds2_R1.csv')
    kmer_gauss_params_ds2_R2.to_csv('/data/fass5/projects/hv_rna_mod/data/mod_gausses/m5C_ds2_R2.csv')  
    """
    
    """
    print('making gauss table for psU ...')
    
    kmer_gauss_params_ds1 = make_gaussian_parameters_kmerwise_table(psU_ea, psU_raw, psU_dorado_bam, 'psU')
    kmer_gauss_params_ds1.to_csv('/data/fass5/projects/hv_rna_mod/data/mod_gausses/psU_ds1.csv')
    
    kmer_gauss_params_ds2 = make_gaussian_parameters_kmerwise_table(concat_mod_df, ds2_raw, ds2_dorado_bam, 'psU')
    kmer_gauss_params_ds2.to_csv('/data/fass5/projects/hv_rna_mod/data/mod_gausses/psU_ds2.csv')
    """
    

if __name__ == '__main__':
    main()
    



