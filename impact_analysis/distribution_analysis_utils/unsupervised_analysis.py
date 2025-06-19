"""
author: Hadi Vareno
e-mail: mohammad.noori.vareno@uni-jena.de
github: https://github.com/TheVareno
"""

import numpy as np 
import pandas as pd  # type: ignore
from sklearn.cluster import KMeans # type: ignore
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from collections import Counter
from sklearn.metrics import silhouette_score # type: ignore
import plotly.express as px # type: ignore
import matplotlib.pyplot as plt
import seaborn as sns # type: ignore
import os
import warnings

warnings.filterwarnings('ignore')
os.chdir('/home/hi68ren/Dokumente/MA/imp/scripts')
sns.set_style('darkgrid') 

# encoding the nucleotide sequence
def perform_OHE(df: pd.DataFrame)-> pd.DataFrame: 
    return pd.get_dummies(df, dtype=int)


def find_k(
    df: pd.DataFrame, 
    features: list, 
    cluster_range: range, 
    modification: str, 
    dataset: str, 
    position: str,
    save_path: str = "../../data/all_figs"
) -> None:
    
    wcss = []                # Within-Cluster Sum of Squares
    silhouette_scores = []   # Silhouette Scores

    for k in cluster_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(df[features])
        wcss.append(kmeans.inertia_)
        silhouette_scores.append(silhouette_score(df[features], cluster_labels))

    plt.figure(figsize=(10, 4))
    sns.set_style('darkgrid') 

    # Elbow Method
    plt.subplot(1, 2, 1)
    plt.plot(cluster_range, wcss, marker='o', linestyle='--', linewidth=2.5, color='#4CBB17')
    plt.title('Elbow Method', fontsize=16, weight='bold')
    plt.xlabel('Number of Clusters (k)', fontsize=14, weight='bold')
    plt.xticks(fontsize=12, fontweight='bold')
    plt.yticks(fontsize=12, fontweight='bold')
    plt.ylabel('WCSS', fontsize=14, weight='bold')
    plt.grid(True)

    # Silhouette Score 
    plt.subplot(1, 2, 2)
    plt.plot(cluster_range, silhouette_scores, marker='o', linestyle='--', linewidth=2.5, color='#4CBB17')
    plt.title('Silhouette Scores', fontsize=16, weight='bold')
    plt.xlabel('Number of Clusters (k)', fontsize=14, weight='bold')
    plt.xticks(fontsize=12, fontweight='bold')
    plt.yticks(fontsize=12, fontweight='bold')
    plt.ylabel('Silhouette Score', fontsize=14, weight='bold')
    plt.grid(True)

    plt.suptitle(f'Finding Optimal Number of Clusters (k) for {modification} at Position {position}', fontsize=18, weight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    filename = f"{save_path}/finding_k_{modification}_{dataset}.png"
    plt.savefig(filename, dpi=600, bbox_inches='tight')



def perform_kmeans(df: pd.DataFrame, features, n_cluster: int): 
    kmeans = KMeans(n_clusters=n_cluster, random_state=42, n_init=10)
    df['kmeans_cluster'] = kmeans.fit_predict(df[features]) 
    df['kmeans_cluster'].value_counts() 
    
    
def visualize_kmeans_results(df: pd.DataFrame, modification: str, dataset : str, position: str):

    cluster_nuc_counts = df.groupby('kmeans_cluster')[
        [col for col in df.columns if col.startswith('base_')]
    ].sum()

    plt.figure(figsize=(10, 4)) 

    ax = sns.heatmap(cluster_nuc_counts,
                     annot=True,
                     cmap='Oranges',
                     fmt='.0f',
                     annot_kws={"fontsize": 16, "fontweight": 'bold'}) 

    plt.title(f'Nucleotide Frequencies by Cluster and {modification} at Position {position}', fontsize=18, fontweight='bold')
    plt.xlabel('Encoded Nucleotide Positions', fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=12, fontweight='bold', fontfamily='cursive')
    plt.ylabel('Cluster', fontsize=14, fontweight='bold')
    plt.yticks(fontsize=12, fontweight='bold')
    ax.figure.axes[-1].set_ylabel('Number of Nucleotides', size=14, weight='bold')
    plt.tight_layout()
    plt.savefig(f"../../data/all_figs/heatmap_kmeans_{modification}_{dataset}_pos{position}.png", dpi=600, bbox_inches='tight') 


def hierarchy_cluster(df: pd.DataFrame, modification_position: int, modification: str, dataset: str,
                     dist_col: str = 'euclidean_distance', save_path: str = '../../data/all_figs') -> None:
    
    X = df[[dist_col]].values
    linkage_matrix = linkage(X, method='ward')
    
    df['main_branch'] = fcluster(linkage_matrix, t=2, criterion='maxclust')
    
    branch_counts = df['main_branch'].value_counts()
    larger_branch = branch_counts.idxmax()
    
    larger_branch_mask = df['main_branch'] == larger_branch
    X_larger = df[larger_branch_mask][[dist_col]].values
    
    if len(X_larger) > 1:  
        linkage_matrix_larger = linkage(X_larger, method='ward')
        df.loc[larger_branch_mask, 'sub_branch'] = fcluster(linkage_matrix_larger, t=2, criterion='maxclust')
    else:
        df['sub_branch'] = 1 
    
    motif_labels = df[['base_1', 'base_2', 'base_3', 'base_4', 'base_5']].agg(''.join, axis=1)
    subgroup_motif = ''

    if modification == 'm5c':
        ls_motif_labels = list(motif_labels)
        mod_pos = ls_motif_labels[0].index('C')

        if modification_position == 1:
            subgroup_motif = 'C*DDDD'
        elif modification_position == 2:
            subgroup_motif = 'DC*DDD'
        elif modification_position == 3:
            subgroup_motif = 'DDC*DD'
        elif modification_position == 4:
            subgroup_motif = 'DDDC*D'
        elif modification_position == 5:
            subgroup_motif = 'DDDDC*'

        def get_context(row):
            bases = [row[f'base_{i}'] for i in range(1, 6)]
            mod_pos = bases.index('C')
            left = bases[mod_pos - 1] if mod_pos > 0 else '-'
            right = bases[mod_pos + 1] if mod_pos < 4 else '-'
            return f"{left}C{right}"

    elif modification == 'psu':
        ls_motif_labels = list(motif_labels)
        mod_pos = ls_motif_labels[0].index('U')

        if modification_position == 1:
            subgroup_motif = 'U*VVVV'
        elif modification_position == 2:
            subgroup_motif = 'VU*VVV'
        elif modification_position == 3:
            subgroup_motif = 'VVU*VV'
        elif modification_position == 4:
            subgroup_motif = 'VVVU*V'
        elif modification_position == 5:
            subgroup_motif = 'VVVVU*'

        def get_context(row):
            bases = [row[f'base_{i}'] for i in range(1, 6)]
            mod_pos = bases.index('U')
            left = bases[mod_pos - 1] if mod_pos > 0 else '-'
            right = bases[mod_pos + 1] if mod_pos < 4 else '-'
            return f"{left}U{right}"

    df['context'] = df.apply(get_context, axis=1)
    
    df['label'] = motif_labels + ' [' + df['context'] + '] '
    df['label'] += np.where(df['main_branch'] == larger_branch, 
                           '(MainB1-SubB' + df['sub_branch'].astype(str) + ')', 
                           '(MainB2)')

    plt.figure(figsize=(12, 5))
    dendrogram(linkage_matrix, labels=df['label'].values, leaf_rotation=90, leaf_font_size=8)
    
    top_contexts = {}
    contexts_main2 = df[df['main_branch'] != larger_branch]['context']
    count_main2 = Counter(contexts_main2)
    top_main2 = count_main2.most_common(5)
    top_contexts['main2'] = "Branch 2\n" + "\n".join([f"{ctx} ({cnt})" for ctx, cnt in top_main2])
    
    for sub_b in [1, 2]:
        contexts_sub = df[(df['main_branch'] == larger_branch) & (df['sub_branch'] == sub_b)]['context']
        count_sub = Counter(contexts_sub)
        top_sub = count_sub.most_common(5)
        top_contexts[f'main1_sub{sub_b}'] = f"Branch 1, Sub-branch {sub_b}\n" + "\n".join([f"{ctx} ({cnt})" for ctx, cnt in top_sub])
    
    annotation_styles = {
        'main2': {'color': 'darkred', 'fontsize': 16, 'fontweight': 'bold', 'fontfamily': 'cursive', 'alpha': 0.7},
        'main1_sub1': {'color': 'navy', 'fontsize': 16, 'fontweight': 'bold', 'fontfamily': 'cursive', 'alpha': 0.7},
        'main1_sub2': {'color': 'darkgreen', 'fontsize': 16, 'fontweight': 'bold', 'fontfamily': 'cursive', 'alpha': 0.7}
    }
    
    plt.text(0.65, 0.5, top_contexts['main2'], transform=plt.gca().transAxes,
             bbox=dict(boxstyle='round,pad=0.5', facecolor='mistyrose', edgecolor='darkred', alpha=0.7),
             **annotation_styles['main2'])
    
    plt.text(0.35, 0.5, top_contexts['main1_sub1'], transform=plt.gca().transAxes,
             bbox=dict(boxstyle='round,pad=0.5', facecolor='aliceblue', edgecolor='navy', alpha=0.7),
             **annotation_styles['main1_sub1'])
    
    plt.text(0.05, 0.5, top_contexts['main1_sub2'], transform=plt.gca().transAxes,
             bbox=dict(boxstyle='round,pad=0.5', facecolor='honeydew', edgecolor='darkgreen', alpha=0.7),
             **annotation_styles['main1_sub2'])
    
    annote_dic = { 1 : '(A)', 2 : '(B)', 3 : '(C)', 4 : '(D)', 5 : '(E)'}
    
    plt.text(-0.06, 1.06, f'{annote_dic[modification_position]}',
             transform=plt.gca().transAxes,  
             fontsize=16, fontweight='bold', color='black', ha='left', va='top')
    
    plt.title(f'Hierarchical Clustering of 5-mers for Subgroup of Motif {subgroup_motif}', fontsize=20, fontweight='bold')
    plt.xlabel('5-mers (labels hidden)', fontsize=14, fontweight='bold')
    # plt.xticks(fontsize=18, fontweight='bold', fontfamily='cursive', rotation=90)
    plt.xticks([])
    plt.ylabel('Distance', fontsize=14, fontweight='bold')
    plt.yticks(fontweight='bold', fontsize=12)
    # plt.ylim([0, 30])
    plt.tight_layout()
    
    os.makedirs(save_path, exist_ok=True)
    plt.savefig(f"{save_path}/hierarchical_{modification}_{dataset}_pos{modification_position}.png", 
                dpi=600, bbox_inches='tight')
  

