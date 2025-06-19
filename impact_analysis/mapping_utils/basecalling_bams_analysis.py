
"""
author: Hadi Vareno
e-mail: mohammad.noori.vareno@uni-jena.de
github: https://github.com/TheVareno
"""


import numpy as np 
import pandas as pd # type: ignore   
import pysam # type: ignore 
import matplotlib.pyplot as plt 
from matplotlib.gridspec import GridSpec
import seaborn as sns # type: ignore
import os 
sns.set_style('darkgrid')

def setup_working_dir()-> None:
    os.chdir('/home/hi68ren/Dokumente/MA/imp/scripts')


def analyse_dorado_bam(path_to_bam: str, mode='rb')-> pd.DataFrame:
    """ 
    Recieves: 
        - path to bam files result of basecalling
        
    Performs: 
        - basecalling info extraction: mean std and quality score using paysam lib! 
    
    Returns: 
        - pandas dataframe consisting of read id, mean , std and phred score of quality of basecalling.  
    """
    
    samfile = pysam.AlignmentFile(path_to_bam, mode, check_sq=False)
    data = [] 
    
    for read in samfile.fetch(until_eof=True): # read type:  <class 'pysam.libcalignedsegment.AlignedSegment'>

        read_name = read.query_name 
        tags = dict(read.tags) 
        
        mean = tags.get('sm', np.nan)
        std_dev = tags.get('sd', np.nan) 
        qual_score = tags.get('qs', np.nan) # quality score of base calling -> lower in case of modifications
        
        data.append([read_name, mean, std_dev, qual_score])

    return pd.DataFrame(data, columns=['read_name', 'mean', 'std_dev', 'qual_score'])


def visualise_basecaling_quality_scores(path_to_can_bam, path_to_mod_bam)-> None:
    
    """
    Plot of ditribution of quality score of basecalling using among can / mod datasets.  
    
    """
    can_tags_df = analyse_dorado_bam(path_to_can_bam)
    mod_tags_df = analyse_dorado_bam(path_to_mod_bam)
    
    fig = plt.figure(figsize=(30, 20))
    sns.histplot(can_tags_df['qual_score'], kde= True, fill=True, color='#6395ee', legend=False, label='canonical', binwidth=1, alpha=0.5) 
    plt.axvline(can_tags_df['qual_score'].mean(), ls='dashed', lw='1.1' , color='#fa003f', label='can qs mean' )
    sns.histplot(mod_tags_df['qual_score'], kde= True, fill=True ,color='#ffa800', legend=False, label='modified', binwidth=1, alpha=0.5) 
    plt.axvline(mod_tags_df['qual_score'].mean(), ls='dashed', lw='1.1', color='#4cbb17', label='mod qs mean')
    plt.title('Distribution of Mean of Base-calling Quality Score for Reads with Canonical and Modified Bases', fontweight='bold', fontsize=22)
    plt.xlabel('Mean of Quality Score', fontweight='bold', fontsize=16)
    plt.ylabel('Frequency / Density', fontweight='bold', fontsize=16)
    can_qs_mean = can_tags_df['qual_score'].mean()
    can_qs_median = can_tags_df['qual_score'].median()
    can_qs_min = can_tags_df['qual_score'].min()
    can_qs_max = can_tags_df['qual_score'].max()
    mod_qs_mean = mod_tags_df['qual_score'].mean()
    mod_qs_median = mod_tags_df['qual_score'].median()
    mod_qs_min = mod_tags_df['qual_score'].min()
    mod_qs_max = mod_tags_df['qual_score'].max()
    note1 = f"Can:\n\nMean = {can_qs_mean:.2f}\n\nMedian = {can_qs_median:.2f}\n\nMin = {can_qs_min:.2f}\n\nMax = {can_qs_max:.2f}"
    note2 = f"Mod:\n\nMean = {mod_qs_mean:.2f}\n\nMedian = {mod_qs_median:.2f}\n\nMin = {mod_qs_min:.2f}\n\nMax = {mod_qs_max:.2f}"
    fig.text(0.53, 0.40, note1, fontsize=20, fontweight='bold', color='#111184')
    fig.text(0.68, 0.40, note2, fontsize=20, fontweight='bold', color='#c76e00')
    plt.legend(fontsize=18)
    plt.tight_layout(pad=2.0)
    plt.savefig(f'../../data/psU/qs_dist_set1.png', dpi=300)
    plt.show()
    plt.close() 



def filter_low_quality_reads_C(bam_df: pd.DataFrame, column = 'qual_score')-> pd.DataFrame:
    """
    Recieves: 
        - pandas dataframe with extracted basecalling info 
    
    Performs: 
        - Fitering of read ids based on quality score. (qualities over mean)
    
    Returns: 
        - filtered version of dataframe. 
    """
    qs_mean = bam_df[column].mean() 
    return bam_df[bam_df['qual_score'] > qs_mean]


def filter_low_quality_reads_U(bam_df: pd.DataFrame, column = 'qual_score')-> pd.DataFrame:
    """
    Recieves: 
        - pandas dataframe with extracted basecalling info 
    
    Performs: 
        - Fitering of read ids based on quality score. (qualities over mean)
    
    Returns: 
        - filtered version of dataframe. 
    """
    qs_mean = bam_df[column].mean() 
    theta_U = qs_mean + 5 # only for U 
    return bam_df[bam_df['qual_score'] > theta_U]



def filter_low_quality_reads_psU(bam_df: pd.DataFrame)-> pd.DataFrame:
    """
    Recieves: 
        - pandas dataframe with extracted basecalling info 
    
    Performs: 
        - Fitering of read ids based on quality score. (qualities over mean)
    
    Returns: 
        - filtered version of dataframe. 
    """
    return bam_df




def normalize_signal(raw_signal: np.array, mean: float, std: float)-> np.array: 
    """
    Recieves:
        - raw signal data as numpy array.
        - mean and std of signal values of each read  
    
    Perfomes:
        - calculation of z-score of each value of the array. 
    
    Returns:
        - Z normalized raw signal as numpy array. 
    """
    return ( raw_signal - mean ) / std



