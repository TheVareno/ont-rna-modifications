"""
author: Hadi Vareno
e-mail: mohammad.noori.vareno@uni-jena.de
github: https://github.com/TheVareno
"""

import numpy as np 
import pandas as pd  # type: ignore
import statsmodels.api as sm # type: ignore
from sklearn.linear_model import LinearRegression # type: ignore
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score # type: ignore
from sklearn.model_selection import cross_val_score # type: ignore
import matplotlib.pyplot as plt
import seaborn as sns # type: ignore
import string
import os
import warnings

sns.set_style('darkgrid') 
os.chdir('/home/hi68ren/Dokumente/MA/imp/scripts')
warnings.filterwarnings('ignore')


def fit_eval_lm(X_train, y_train, X_test, y_test):
    
    model = LinearRegression()
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    cv_scores = cross_val_score(model, X_train, y_train, cv=3, scoring='r2')

    print("Model Evaluation Metrics:")
    print("-" * 30)
    print(f"  MSE: {mean_squared_error(y_test, y_pred):.5f}")
    print(f"  R²: {r2_score(y_test, y_pred):.5f}")
    print(f"  MAE: {mean_absolute_error(y_test, y_pred):.5f}")
    print("\nCross-Validation (R²):")
    print("-" * 30)
    print(f"  Mean CV Score: {np.mean(cv_scores):.5f}")
    print(f"  CV Std Dev: {np.std(cv_scores):.5f}")
    print(f"  CV Scores: {cv_scores}")

    return model.coef_[0][0], model.intercept_[0] 



def fit_lm_momentwise(main_df_ds1: pd.DataFrame, main_df_ds2: pd.DataFrame, moment: str, modification: str)-> None:

    main_df_ds1['source'] = 'ds1'
    main_df_ds2['source'] = 'ds2'
    combined_df = pd.concat([main_df_ds1, main_df_ds2], ignore_index=True)

    X_train = combined_df.loc[combined_df['source'] == 'ds1', f'{moment}_can'].values.reshape(-1, 1)
    X_test = combined_df.loc[combined_df['source'] == 'ds2', f'{moment}_can'].values.reshape(-1, 1)
    y_train = combined_df.loc[combined_df['source'] == 'ds1', f'{moment}_mod'].values.reshape(-1, 1)
    y_test = combined_df.loc[combined_df['source'] == 'ds2', f'{moment}_mod'].values.reshape(-1, 1)

    model_prop = fit_eval_lm(X_train, y_train, X_test, y_test)
    y_vals = model_prop[1] + model_prop[0] * combined_df[[f'{moment}_can']]

    fig = plt.figure(figsize=(12, 10))
    ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
    ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
    ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
    ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
    ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)

    axes = [ax1, ax2, ax3, ax4, ax5]
    hue_vars = ['base_1', 'base_2', 'base_3', 'base_4', 'base_5']
    subplot_labels = list(string.ascii_uppercase[:len(axes)]) 
    
    for i, ax in enumerate(axes):
        hue_var = hue_vars[i]
        sns.scatterplot(data=combined_df, y=f'{moment}_mod', x=f'{moment}_can', hue=hue_var, palette='viridis', style='source', ax=ax, legend=False) 
        ax.plot(combined_df[f'{moment}_can'], y_vals, color='red', linestyle='-', linewidth=0.8)
        ax.set_xlabel(f'canonical mean', fontsize=14, fontweight='bold')
        ax.set_ylabel('modified mean' if i == 0 or i == 3 else '', fontsize=14, fontweight='bold') # Set Y label only for first in each row

        

        ax.text(-0.10, 1.11, f'({subplot_labels[i]})', transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='left', zorder=3)
        ax.set_title(f'Hue on {hue_var.replace("_", " ").title()}', fontweight='bold', fontsize=16)
        ax.grid(True, linestyle='--', alpha=0.6) 

    source_handles = [
        plt.Line2D([0], [0], marker='o', color='gray', linestyle='None', markersize=8, label='Dataset 1'),
        plt.Line2D([0], [0], marker='X', color='gray', linestyle='None', markersize=8, label='Dataset 2')
    ]

    all_nucleotides = sorted(combined_df[['base_1', 'base_2', 'base_3', 'base_4', 'base_5']].stack().unique())
    nucleotide_colors = sns.color_palette('viridis', n_colors=len(all_nucleotides))
    nucleotide_handles = [] 
    
    for j, nucleotide in enumerate(all_nucleotides):
        nucleotide_handles.append(plt.Line2D([0], [0], marker='o', color=nucleotide_colors[j], linestyle='None', markersize=8, label=f'{nucleotide}'))

    regression_handle = plt.Line2D([0], [0], color='red', linestyle='-', linewidth=1.5, label=f'Fitted Line')

    all_handles = source_handles + nucleotide_handles + [regression_handle]
    all_labels = [h.get_label() for h in all_handles]

    fig.legend(handles=all_handles, labels=all_labels,
               loc='lower right',
               bbox_to_anchor=(0.98, 0.03),
               fontsize=12, 
               frameon=True,
               framealpha=0.9,
               borderpad=0.8, 
               labelspacing=0.8, 
               handletextpad=0.5)

    fig.suptitle(f'Modified Mean vs. Canonical Mean of Multiple Datasets of {modification}\nFitted Line: µ* = {model_prop[0]:.5f} × µ + {model_prop[1]:.5f}',
                 fontsize=20, fontweight='bold', y=1.02) 

    plt.tight_layout(rect=[0, 0.01, 1, 0.99]) 
    fig.savefig(f"../../data/all_figs/linear_model_main_{modification}.png", dpi=600, bbox_inches='tight')


    
def plot_flanking_context_for_base2(main_df_ds1: pd.DataFrame, main_df_ds2: pd.DataFrame, moment: str, modification: str) -> None:
    
    main_df_ds1['source'] = 'ds1'
    main_df_ds2['source'] = 'ds2' 
    combined_df = pd.concat([main_df_ds1, main_df_ds2], ignore_index=True)

    X_train = combined_df.loc[combined_df['source'] == 'ds1', f'{moment}_can'].values.reshape(-1, 1)
    X_test = combined_df.loc[combined_df['source'] == 'ds2', f'{moment}_can'].values.reshape(-1, 1)
    y_train = combined_df.loc[combined_df['source'] == 'ds1', f'{moment}_mod'].values.reshape(-1, 1)
    y_test = combined_df.loc[combined_df['source'] == 'ds2', f'{moment}_mod'].values.reshape(-1, 1)
    model_prop = fit_eval_lm(X_train, y_train, X_test, y_test)    
    
    y_vals = model_prop[1] + model_prop[0] * combined_df[f'{moment}_can']

    # m5C
    filtered_df = combined_df[combined_df['base_2'].isin(['G', 'U'])].copy() 
    filtered_df['highlight'] = 'none'
    filtered_df.loc[filtered_df['base_1'] == 'C', 'highlight'] = 'C at base_1'
    filtered_df.loc[filtered_df['base_3'] == 'C', 'highlight'] = 'C at base_3'

    dic_G = {'none': 'other', 'C at base_1': 'CG', 'C at base_3': 'GC'}
    dic_U = {'none': 'other', 'C at base_1': 'CU', 'C at base_3': 'UC'}
    
    filtered_df['context'] = filtered_df['base_1'] + filtered_df['base_3']
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    
    marker_dict = {
        'none': 'o',
        'C at base_1': 's',    
        'C at base_3': 'D',   
    }
    
    palette = {
        'none': '#999999',
        'C at base_1': '#6395EE',
        'C at base_3': '#FFA800',
    }

    for i, base in enumerate(['G', 'U']):
        ax = axes[i]
        subset_for_main_plot = filtered_df[filtered_df['base_2'] == base].copy() 

        if base == 'U':
            current_edgecolor = '#7CFC00' 
        elif base == 'G':
            current_edgecolor = '#0BDA51' 
        else:
            current_edgecolor = 'none'
        
            
        for i, label in enumerate(['none', 'C at base_1', 'C at base_3']):
            sub = subset_for_main_plot[subset_for_main_plot['highlight'] == label]
            
            if base == 'G':
                sns.scatterplot(
                    data=sub,
                    x=f'{moment}_can',
                    y=f'{moment}_mod',
                    ax=ax,
                    s=70 if label != 'none' else 30,
                    marker=marker_dict[label],
                    color=palette[label],
                    edgecolor=current_edgecolor if label != 'none' else 'none',
                    label=dic_G[label], 
                )
                
                ax.text(-0.10, 1.13, '(A)',
                transform=ax.transAxes, 
                fontsize=20, 
                color='black',            
                fontweight='bold',       
                va='top',                
                ha='left'                
               )            
            
            if base == 'U':
                sns.scatterplot(
                    data=sub,
                    x=f'{moment}_can',
                    y=f'{moment}_mod',
                    ax=ax,
                    s=70 if label != 'none' else 30,
                    marker=marker_dict[label],
                    color=palette[label],
                    edgecolor=current_edgecolor if label != 'none' else 'none',
                    label=dic_U[label], 
                )

                ax.text(-0.10, 1.13, '(B)',
                transform=ax.transAxes, 
                fontsize=20, 
                color='black',            
                fontweight='bold',       
                va='top',                
                ha='left'                
               )
            
        
        ax.plot(combined_df[f'{moment}_can'], y_vals, color='#FA5053', linestyle='-', alpha=0.7, linewidth=0.7)
        ax.set_title(f'2nd base = {base}', fontsize=16, fontweight='bold')
        ax.set_xlabel('canonical mean', fontsize=14, fontweight='bold')
        ax.set_ylabel('modified mean', fontsize=14, fontweight='bold')
        ax.legend(title='', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=14)
       
    plt.suptitle('Modified Mean vs. Canonical Mean for 2nd. Base = G or U\nHighlighting m5C/C Contexts at Positions 1 and 3', fontsize=20, fontweight='bold')
    plt.tight_layout(pad=1.0)
    fig.savefig(f"../../data/all_figs/lin_model_base_2_GU_{modification}.png", dpi=600, bbox_inches='tight')
    
    

def arrow_plot_flanking_C(main_df_ds1: pd.DataFrame, main_df_ds2: pd.DataFrame): 

    main_df_ds1['source'] = 'ds1'
    main_df_ds2['source'] = 'ds2'
    df = main_df_ds1.copy()
    
    df['mu_can_ds2'] = main_df_ds2['mu_can']
    df['mu_mod_ds2'] = main_df_ds2['mu_mod']

    df = df[df['base_2'].isin(['G', 'U'])].copy()

    df['combined_context'] = 'Other'
    df.loc[(df['base_2'] == 'G') & (df['base_1'] == 'C'), 'combined_context'] = 'CG'  
    df.loc[(df['base_2'] == 'G') & (df['base_3'] == 'C'), 'combined_context'] = 'GC' 
    df.loc[(df['base_2'] == 'U') & (df['base_1'] == 'C'), 'combined_context'] = 'CU' 
    df.loc[(df['base_2'] == 'U') & (df['base_3'] == 'C'), 'combined_context'] = 'UC' 

    df = df[df['combined_context'] != 'Other'].copy()

    combined_context_colors = {
        'CG': '#1f77b4',  
        'GC': '#87CEEB',  
        'CU': '#C76E00',  
        'UC': '#FFA800', 
    }

    fig = plt.figure(figsize=(8, 6))

    for _, row in df.iterrows():
        x_start = row['mu_can']
        y_start = row['mu_can_ds2']
        x_end = row['mu_mod']
        y_end = row['mu_mod_ds2']
        color = combined_context_colors[row['combined_context']]

        plt.arrow(x_start, y_start, x_end - x_start, y_end - y_start,
                  color=color, alpha=0.8, head_width=0.05, linewidth=1.8,
                  length_includes_head=True, zorder=2) 

    
    min_val = min(df['mu_can'].min(), df['mu_can_ds2'].min())
    max_val = max(df['mu_can'].max(), df['mu_can_ds2'].max())
    
    plt.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.6,
             label='Canonical Reference (ds1 = ds2)', zorder=1) 
   
    
    plt.text(0.85, 0.25, 'Arrow $\u2192$: µ* - µ', transform=plt.gca().transAxes,  
             fontsize=17, fontweight='bold', fontfamily = 'cursive', color='#6D8196', ha='right', va='bottom')
    

    plt.xlabel('Canonical Mean (Dataset 1)', fontsize=14, fontweight='bold') 
    plt.ylabel('Canonical Mean (Dataset 2)', fontsize=14, fontweight='bold') 
    plt.title('Paired Arrow Plot: Shift from Canonical to m5C-Modified Mean', fontsize=15, fontweight='bold')

    legend_elements = []
    for label in sorted(combined_context_colors.keys()):
        color = combined_context_colors[label]
        legend_elements.append(plt.Line2D([0], [0], color=color, lw=3, label=label,
                                          marker='>', markersize=11, markeredgecolor=color,
                                          markerfacecolor=color, linestyle='-')) 

    legend_elements.append(plt.Line2D([0], [0], color='k', lw=1.5, linestyle='--', label='Canonical Reference'))

    plt.legend(handles=legend_elements, title='Base at Position 2 & Flanking C/m5C',
               loc='best', frameon=True, fontsize=14, title_fontsize=14)

    plt.grid(True, linestyle='--', alpha=0.4, zorder=0) 
    plt.tight_layout(pad=0.5)

    fig.savefig(f"../../data/all_figs/arrow_plot.png", dpi=600, bbox_inches='tight')
    
    
    

