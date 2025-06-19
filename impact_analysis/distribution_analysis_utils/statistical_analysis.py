"""
author: Hadi Vareno
e-mail: mohammad.noori.vareno@uni-jena.de
github: https://github.com/TheVareno
"""


import numpy as np 
from numpy import triu 
import pandas as pd  # type: ignore
from scipy import stats
import scikit_posthocs as sp # type: ignore
import matplotlib.pyplot as plt
import seaborn as sns # type: ignore
import os
import string
import warnings

sns.set_style('darkgrid') 
os.chdir('/home/hi68ren/Dokumente/MA/imp/scripts')
warnings.filterwarnings('ignore')


##################################################################
##################################################################
############################# EDA ################################
##################################################################
##################################################################


def can_mod_trends_plot_datasetwise(main_df: pd.DataFrame, modification: str, dataset: str) -> None:
    
    sns.set_theme(style="whitegrid", font_scale=1.5, rc={
        "axes.labelsize": 14,
        "axes.titlesize": 16,
        "legend.fontsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "figure.dpi": 300
    })

    settings = dict(kde=True, alpha=0.4, stat='density', bins=50)

    if modification == 'psU' and dataset == 'ds1': 
        
        fig, axs = plt.subplots(3, 2, figsize=(12, 10))  
        
        sns.histplot(main_df['mu_can'], label='Canonical', color='#6395EE', ax=axs[0, 0], **settings)
        axs[0, 0].axvline(main_df['mu_can'].mean(), color='red', linestyle='--', linewidth=1.3, label='Mean')
        axs[0, 0].axvspan(main_df['mu_can'].mean() - main_df['mu_can'].std(), main_df['mu_can'].mean() + main_df['mu_can'].std(), color='green',
                          linestyle='--', linewidth=1.3, alpha=0.2, label='Std. Dev.')

        sns.histplot(main_df['sd_can'], label='Canonical', color='#6395EE', ax=axs[0, 1], **settings)
        axs[0, 1].axvline(main_df['sd_can'].mean(), color='red', linestyle='--', linewidth=1.3, label='Mean')
        axs[0, 1].axvspan(main_df['sd_can'].mean() - main_df['sd_can'].std(), main_df['sd_can'].mean() + main_df['sd_can'].std(), color='green',
                          linestyle='--', linewidth=1.3, alpha=0.2, label='Std. Dev.')

        sns.histplot(main_df['mu_mod'], label='Modified', color='#FFA800', ax=axs[1, 0], **settings)
        axs[1, 0].axvline(main_df['mu_mod'].mean(), color='red', linestyle='--', linewidth=1.3, label='Mean')
        axs[1, 0].axvspan(main_df['mu_mod'].mean() - main_df['mu_mod'].std(), main_df['mu_mod'].mean() + main_df['mu_mod'].std(), color='green',
                          linestyle='--', linewidth=1.3, alpha=0.2, label='Std. Dev.')

        sns.histplot(main_df['sd_mod'], label='Modified', color='#FFA800', ax=axs[1, 1], **settings)
        axs[1, 1].axvline(main_df['sd_mod'].mean(), color='red', linestyle='--', linewidth=1.3, label='Mean')
        axs[1, 1].axvspan(main_df['sd_mod'].mean() - main_df['sd_mod'].std(), main_df['sd_mod'].mean() + main_df['sd_mod'].std(), color='green',
                          linestyle='--', linewidth=1.3, alpha=0.2, label='Std. Dev.')
        
        sns.histplot(main_df['med_can'], label='Canonical', color='#6395EE', ax=axs[2, 0], **settings)
        sns.histplot(main_df['med_mod'], label='Modified', color='#FFA800', ax=axs[2, 0], **settings)

        sns.histplot(main_df['mad_can'], label='Canonical', color='#6395EE', ax=axs[2, 1], **settings)
        sns.histplot(main_df['mad_mod'], label='Modified', color='#FFA800', ax=axs[2, 1], **settings)
        
        titles = ['Canonical Mean', 'Canonical Standard Deviation', 'Modified Mean', 'Modified Standard Deviation', 'Median', 'Median Absolute Deviation']
                
        for ax, title in zip(axs.flat, titles):
            title_print = f"Parameter: {title}"
            ax.set_title(title_print, fontweight='bold', fontsize=16)
            ax.set_xlabel(title, fontsize=14, fontweight='bold')
            ax.set_ylabel('Density', fontsize=14, fontweight='bold')  
            
            ax.legend()
            
            if title == 'Modified Mean': 
                ax.set_xlim(-13, 15)
            
            if title == 'Modified Standard Deviation': 
                ax.set_xlim(0, 10)
                
                        
    else: 
        
        fig, axs = plt.subplots(2, 2, figsize=(10, 8)) 
        
        sns.histplot(main_df['mu_can'], label='Canonical', color='#6395EE', ax=axs[0, 0], **settings)
        sns.histplot(main_df['mu_mod'], label='Modified', color='#FFA800', ax=axs[0, 0], **settings)

        sns.histplot(main_df['sd_can'], label='Canonical', color='#6395EE', ax=axs[0, 1], **settings)
        sns.histplot(main_df['sd_mod'], label='Modified', color='#FFA800', ax=axs[0, 1], **settings)

        sns.histplot(main_df['med_can'], label='Canonical', color='#6395EE', ax=axs[1, 0], **settings)
        sns.histplot(main_df['med_mod'], label='Modified', color='#FFA800', ax=axs[1, 0], **settings)

        sns.histplot(main_df['mad_can'], label='Canonical', color='#6395EE', ax=axs[1, 1], **settings)
        sns.histplot(main_df['mad_mod'], label='Modified', color='#FFA800', ax=axs[1, 1], **settings)
        
        titles = ['Mean', 'Standard Deviation', 'Median', 'Median Absolute Deviation']
    
        for ax, title in zip(axs.flat, titles):
            title_print = f'Parameter: {title}'
            ax.set_title(title_print, fontweight='bold', fontsize=16)
            ax.set_xlabel(title, fontsize=14, fontweight='bold')
            ax.set_ylabel('Density', fontsize=14, fontweight='bold')  
            ax.legend()

    if modification == 'psU' and dataset == 'ds1': subplot_labels = ['A', 'B', 'C', 'D', 'E', 'F'] 
    else: subplot_labels = ['A', 'B', 'C', 'D'] 
    
    for i, ax in enumerate(axs.flat): 
        ax.text(-0.13, 1.18, f'({subplot_labels[i]})',
                transform=ax.transAxes, 
                fontsize=20, 
                color='black',            
                fontweight='bold',       
                va='top',                
                ha='left'                
               )
    
    fig.suptitle(f'Distribution Comparison of Statistical Moments for {modification}',
                 fontweight='bold', fontsize=20, y=1.03)

    plt.tight_layout(pad=0.1)
    plt.subplots_adjust(top=0.90)
    fig.savefig(f"../../data/all_figs/moments_distribution_{modification}_{dataset}.png", dpi=600, bbox_inches='tight')



def can_mod_outlier_plot_datasetwise(main_df: pd.DataFrame, modification: str, dataset: str) -> None:
    
    sns.set_theme(style="whitegrid", font_scale=1.5, rc={
        "axes.labelsize": 14,
        "axes.titlesize": 16,
        "legend.fontsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "figure.dpi": 600
    })

    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    features = [
        ("mu", "Mean"),
        ("sd", "Standard Deviation"),
        ("med", "Median"),
        ("mad", "Median Absolute Deviation")
    ]

    flier_props = dict(marker='o',
                       markerfacecolor='red', 
                       markeredgecolor='red', 
                       markersize=6, 
                       linestyle='none'
                      )
    
    subplot_labels = list(string.ascii_uppercase[:len(features)]) 
    
    for i, (ax, (prefix, title)) in enumerate(zip(axs.flat, features)):
        data_to_plot = pd.DataFrame({
            "Canonical": main_df[f"{prefix}_can"],
            "Modified": main_df[f"{prefix}_mod"]
        })

        sns.boxplot(data=data_to_plot,
                    ax=ax,
                    palette=["#6395EE", "#FFA800"],
                    linewidth=1.3,
                    fliersize=2,
                    width=0.6,
                    fill=True, 
                    flierprops=flier_props
                    )

        ax.set_title(title, fontweight='bold', fontsize=20)
        ax.set_ylabel('Signal Value', fontsize=18, fontweight='bold')
        
        if modification == 'psU': 
            ax.set_ylim(-5, 5)
        
            if prefix == 'sd' or prefix == 'mad': 
                ax.set_ylim(-2, 2)
        
        
        ax.set_xlabel('', fontsize=15)
        ax.set_xticklabels(['Canonical', 'Modified'], fontsize=18, fontweight='bold')
        ax.grid(True) 
        
        ax.text(-0.13, 1.18, f'({subplot_labels[i]})',
                transform=ax.transAxes, 
                fontsize=20,             
                fontweight='bold',       
                va='top',                
                ha='left'                
               )

    fig.suptitle(f'Distribution of Statistical Features for {modification}',
                 fontweight='bold', fontsize=18, y=1.03)

    plt.tight_layout()
    plt.subplots_adjust(top=0.90)

    fig.savefig(f"../../data/all_figs/boxplot_{modification}_{dataset}.png", dpi=600, bbox_inches='tight')


##################################################################
##################################################################
######################## EUC. DISTANCE ###########################
##################################################################
##################################################################


# calculate distance column, default : euclidean
def make_distance_df(main_df: pd.DataFrame, mode: str, metric: str = 'euclidean') -> pd.DataFrame:
    
    if mode == 'contextwise':
        context_cols = ['base_1', 'base_2', 'base_4', 'base_5', 'modified_position']
    elif mode in ['sliding', 'CG_rich', 'psU_pocket']:
        context_cols = ['base_1', 'base_2', 'base_3', 'base_4', 'base_5', 'modified_position']
    else:
        raise ValueError(f"Unknown mode: {mode}")

    can_cols = ['mu_can', 'sd_can', 'med_can', 'mad_can']
    mod_cols = ['mu_mod', 'sd_mod', 'med_mod', 'mad_mod']

    X_can = main_df[can_cols].to_numpy()
    X_mod = main_df[mod_cols].to_numpy()

    if metric == 'euclidean':
        distances = np.linalg.norm(X_mod - X_can, axis=1)
    elif metric == 'manhattan':
        distances = np.sum(np.abs(X_mod - X_can), axis=1)
    elif metric == 'chebyshev':
        distances = np.max(np.abs(X_mod - X_can), axis=1)
    else:
        distances = np.linalg.norm(X_mod - X_can, axis=1)
        

    distance_df = main_df[context_cols].copy()
    distance_df[f'{metric}_distance'] = distances

    return distance_df



##################################################################
##################################################################
######################### STAT. TESTS ############################
##################################################################
##################################################################


def add_modification_position(df: pd.DataFrame, modification: str)-> pd.DataFrame:
    
    positions = []
    if modification == 'm5C': 
        # find where is the 'C'
        for index, row in df.iterrows(): # check each 5-mer
            if row['base_1'] == 'C':
                positions.append(1)
            elif row['base_2'] == 'C':
                positions.append(2)
            elif row['base_3'] == 'C':
                positions.append(3)
            elif row['base_4'] == 'C':
                positions.append(4)
            elif row['base_5'] == 'C':
                positions.append(5)
            else:
                positions.append(None) 
    
    if modification == 'psU': 
        for index, row in df.iterrows(): # check each 5-mer
            if row['base_1'] == 'U':
                positions.append(1)
            elif row['base_2'] == 'U':
                positions.append(2)
            elif row['base_3'] == 'U':
                positions.append(3)
            elif row['base_4'] == 'U':
                positions.append(4)
            elif row['base_5'] == 'U':
                positions.append(5)
            else:
                positions.append(None) 
        
    df['modified_position'] = positions 
    return df
    

def perform_kruskal(df: pd.DataFrame, modification: str):
    
    df = add_modification_position(df, modification)
    modified_positions = sorted(df['modified_position'].unique())
    distance_cols = ['mu_diff', 'sd_diff', 'med_diff', 'mad_diff']
    
    results = {}
    for col in distance_cols:
        groups = [df[df['modified_position'] == pos][col] for pos in modified_positions]
        h_stat, p_kruskal = stats.kruskal(*groups)
        results[col] = {
            'Kruskal H-statistic': h_stat, 'Kruskal p-value': p_kruskal
        }
    return pd.DataFrame(results)


def perform_mann_whitney(df: pd.DataFrame):
    distance_cols = ['mu_diff', 'sd_diff', 'med_diff', 'mad_diff']
    for i in [1, 2, 3, 4, 5]:
        print(f"\n--- Testing for modification at position {i} ---")
        for col in distance_cols:
            group_mod = df[df["modified_position"] == i][col]
            group_unmod = df[df["modified_position"] != i][col]
            stat, p_value = stats.mannwhitneyu(group_mod, group_unmod, alternative='two-sided')
            print(f"{col}: U = {stat:.3f}, p = {p_value:.6f}")
    
        
def check_assumptions(df: pd.DataFrame, modification: str):
    
    df = add_modification_position(df, modification)
    modification_positions = sorted(df['modified_position'].unique())
    distance_cols = ['mu_diff', 'sd_diff', 'med_diff', 'mad_diff']
    
    print("Shapiro-Wilk Normality Test (for each feature in each group):")
    for col in distance_cols:
        print(f"\nFeature: {col}")
        normality_met = True
        for pos in modification_positions:
            values = df[df['modified_position'] == pos][col]
            if len(values) < 3:  # Shapiro-Wilk requires at least 3 observations
                print(f"  Position {pos}: Not enough data for Shapiro-Wilk.")
                continue
            stat, p = stats.shapiro(values)
            print(f"  Position {pos}: Statistic={stat:.4f}, p={p:.4f}")
            if p < 0.05:
                normality_met = False
        if normality_met:
            print(f"Normality appears to be met for {col}.")
        else:
            print(f"Normality might be violated for {col}.")

    print("\nLevene's Test for Equal Variances (per feature):")
    for col in distance_cols:
        groups = [df[df['modified_position'] == pos][col] for pos in modification_positions]
        stat, p = stats.levene(*groups)
        print(f"Feature: {col} — Statistic={stat:.4f}, p={p:.4f}")
        if p >= 0.05:
            print(f"Homogeneity of variances appears to be met for {col}.")
        else:
            print(f"Homogeneity might be violated for {col}.")


# post-hoc test with Bonferroni correction
def perform_dunn_posthoc(df: pd.DataFrame): 
    
    posthoc = sp.posthoc_dunn(df, val_col='mu_diff', group_col='modified_position', p_adjust='bonferroni')
    print(posthoc)


# KW & MW test result
def get_result_table(modification: str, dataset: str, test: str): 
    
    m5c_ds1 = {
        "Mean":                 [0.000325, 0.2373, 0.0083, 0.3572, 0.0033, 0.9746],
        "Std. Dev.":            [0.001093, 0.0992, 0.0521, 0.0085, 0.3874, 0.0393],
        "Median":               [0.000414, 0.3000, 0.0091, 0.4416, 0.0024, 0.8697],
        "Median Abs. Dev.":     [0.000204, 0.0849, 0.0201, 0.0091, 0.1279, 0.0949]
    }
    
    m5c_ds2 = {
        "Mean":                 [ 0.118532, 0.153934, 0.152104, 0.320296, 0.040490, 0.853916],
        "Std. Dev.":            [ 0.658567, 0.824884, 0.985183, 0.462396, 0.712292, 0.367311],
        "Median":               [ 0.113127, 0.174174, 0.083379, 0.173167, 0.060400, 0.881472],
        "Median Abs. Dev.":     [ 0.937966, 0.947962, 0.838958, 0.925173, 0.685579, 0.441973]
    } 
    
    psu_ds1 = {
        "Mean":                 [ 0.025158, 0.713905, 0.002475, 0.194335, 0.090945, 0.776342],
        "Std. Dev.":            [ 0.821951, 0.955730, 0.467331, 0.761524, 0.548319, 0.335355],
        "Median":               [ 0.001779, 0.248965, 0.001025, 0.495672, 0.006430, 0.220600],
        "Median Abs. Dev.":     [ 0.002005, 0.064498, 0.495459, 0.033155, 0.000390, 0.791883]
    }
    
    psu_ds2 = {
        "Mean":                 [0.220092, 0.821580, 0.197074, 0.409911, 0.064275, 0.246725],
        "Std. Dev.":            [0.925804, 0.585791, 0.783831, 0.845599, 0.579965, 0.529491],
        "Median":               [0.229845, 0.556226, 0.122179, 0.176880, 0.153322, 0.519809],
        "Median Abs. Dev.":     [0.831694, 0.830672, 0.756248, 0.644721, 0.407503, 0.372409]
    } 
    
    m5c_ds1_dunn = {
        "Pos. 1": [1.000000, 1.000000, 1.000000, 0.003834, 1.000000],
        "Pos. 2": [1.000000, 1.000000, 0.270903, 0.000357, 0.299096],
        "Pos. 3": [1.000000, 0.270903, 1.000000, 0.544025, 1.000000],
        "Pos. 4": [0.003834, 0.000357, 0.544025, 1.000000, 0.496989],
        "Pos. 5": [1.000000, 0.299096, 1.000000, 0.496989, 1.000000],
    }
    
    
    psu_ds1_dunn = {
        "Pos. 1": [1.000000, 0.325621, 1.000000, 1.000000, 1.000000],
        "Pos. 2": [0.325621, 1.000000, 0.065059, 0.029196, 0.822167],
        "Pos. 3": [1.000000, 0.065059, 1.000000, 1.000000, 1.000000],
        "Pos. 4": [1.000000, 0.029196, 1.000000, 1.000000, 1.000000],
        "Pos. 5": [1.000000, 0.822167, 1.000000, 1.000000, 1.000000],
    }
    
    if modification == 'm5C' and dataset == 'ds1' and test == 'kw-mw': 
        return m5c_ds1

    elif modification == 'm5C' and dataset == 'ds2' and test == 'kw-mw': 
        return m5c_ds2 
    
    elif modification == 'psU' and dataset == 'ds1' and test == 'kw-mw': 
        return psu_ds1 
    
    elif modification == 'psU' and dataset == 'ds2' and test == 'kw-mw': 
        return psu_ds2 
    
    elif modification == 'm5C' and dataset == 'ds1' and test == 'dunn': 
        return m5c_ds1_dunn
    
    elif modification == 'psU' and dataset == 'ds1' and test == 'dunn': 
        return psu_ds1_dunn
    

def visualize_test_result(modification: str, dataset: str, test: str):

    results = get_result_table(modification, dataset, test)
    
    if test == 'kw-mw':
        columns = ["KW", "MW-Pos. 1", "MW-Pos. 2", "MW-Pos. 3", "MW-Pos. 4", "MW-Pos. 5"]
        pval_df = pd.DataFrame(results, index=columns).T

        def p_to_stars(p):
            if p < 0.001:
                return "***"
            elif p < 0.01:
                return "**"
            elif p < 0.05:
                return "*"
            else:
                return ""

        annotations = pval_df.copy()
        for row in pval_df.index:
            for col in pval_df.columns:
                p = pval_df.loc[row, col]
                annotations.loc[row, col] = f"{p:.5f}{p_to_stars(p)}"

        plt.figure(figsize=(10, 6))

        ax = sns.heatmap(
            pval_df.astype(float),
            annot=annotations,
            fmt="",
            cmap="coolwarm_r",
            cbar_kws={'label': 'p-value'},
            linewidths=1.0,
            linecolor='gray',
            vmin=0,
            vmax=1,
            annot_kws={"fontsize": 16, "fontweight": 'bold'}
        )

        plt.title(f"Statistical Test p-values for Signal-Based Features - {modification}", fontsize=18, fontweight='bold')
        plt.xlabel("Kruskal–Wallis(KW) and Mann–Whitney (MW) Statistical Test", fontsize=14, fontweight='bold')
        plt.ylabel("Shifts in Signal Feature", fontsize=14, fontweight='bold')
        plt.xticks(rotation=45, ha='right', fontsize=12, fontweight='bold', fontfamily='cursive')
        plt.yticks(fontsize=12, fontweight='bold', fontfamily='cursive')
        ax.figure.axes[-1].set_ylabel('p-value', size=14, weight='bold')
        plt.tight_layout()
        plt.savefig(f'../../data/all_figs/statistical_test_pvalues_heatmap_{modification}_{dataset}.png', dpi=600) 
    
    elif test == 'dunn':  
        
        columns = ["Pos. 1", "Pos. 2", "Pos. 3", "Pos. 4", "Pos. 5"]
        pval_df = pd.DataFrame(results, index=columns).T

        def p_to_stars(p):
            if p < 0.001:
                return "***"
            elif p < 0.01:
                return "**"
            elif p < 0.05:
                return "*"
            else:
                return ""

        annotations = pval_df.copy()
        for row in pval_df.index:
            for col in pval_df.columns:
                p = pval_df.loc[row, col]
                annotations.loc[row, col] = f"{p:.5f}{p_to_stars(p)}"


        plt.figure(figsize=(10, 7))

        ax = sns.heatmap(
            pval_df.astype(float),
            annot=annotations,
            fmt="",
            cmap="coolwarm_r",
            cbar_kws={'label': 'p-value'},
            linewidths=1.0,
            linecolor='gray',
            vmin=0,
            vmax=1,
            annot_kws={"fontsize": 16, "fontweight": 'bold'}
        )

        plt.title(f"Pairwise p-Values for Dunn's Post-hoc Test - {modification}", fontsize=18, fontweight='bold')
        plt.xlabel("Modification Position", fontsize=14, fontweight='bold')
        plt.ylabel("Modification Position", fontsize=14, fontweight='bold')
        plt.xticks(rotation=0, ha='right', fontsize=12, fontweight='bold', fontfamily='cursive')
        plt.yticks(fontsize=12, fontweight='bold', fontfamily='cursive')
        ax.figure.axes[-1].set_ylabel('p-value', size=14, weight='bold')
        plt.tight_layout()
        plt.savefig(f'../../data/all_figs/posthoc_heatmap_{modification}_{dataset}.png', dpi=600)




    