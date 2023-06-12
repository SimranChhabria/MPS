#!/usr/bin/env python3

import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import statsmodels.stats.multitest as smt
from mne.stats import bonferroni_correction, fdr_correction

from rpy2.robjects.packages import importr
from rpy2.robjects import FloatVector, r

stats = importr('stats')


""" Handle arguments """
parser = argparse.ArgumentParser(
    description="Subset the data"
)

parser.add_argument("-i", "--input_file",type=argparse.FileType('r'),
    help="Meta data of the entire project",
 )

parser.add_argument("-b", "--balaji_file",type=argparse.FileType('r'),
    help="Meta data of the entire project",
)

parser.add_argument("-o", "--out_file", type=argparse.FileType('w'), 
        help="Output file"
        )

args      = parser.parse_args()
input_df  = pd.read_csv(args.input_file, sep=",")
balaji_df = pd.read_csv(args.balaji_file, sep='\t')

balaji_df = balaji_df.rename(columns = {'Approved.Symbol' : 'gene_name'} )
balaji_df = balaji_df.filter(regex='_fdr$|gene_name')

# Define function to apply FDR correction to a value
def fdr_correct(val):
    p_adjust = r['p.adjust'](FloatVector(np.array(val, dtype=float)), method='BH')
    return np.array(p_adjust)
    
    #reject, corrected_val = fdr_correction([val], method='indep')
    #p_adj = smt.multipletests([val], method='fdr_bh')[1]
    #return corrected_val[0]
    #return(p_adj)


RBP_list= []
for col in input_df.columns[1:]:
     x = col.split('_', 1)[0]  
     if (x!=('Control')) and x not in RBP_list:
        RBP_list.append(x)

for RBP in RBP_list:
    if RBP == "log":
        continue
    input_df[f'{RBP}_pval_adj']    = fdr_correct(input_df[f'{RBP}_log_pval'])
    input_df[f'{RBP}_MW_pval_adj'] = fdr_correct(input_df[f'{RBP}_log_Mann_whitney_pval'])


input_df.to_csv(args.out_file, index=None)

p_adj_new= input_df.filter(regex='pval_adj$|gene_name')
merged_df = pd.merge(p_adj_new, balaji_df, on = "gene_name", how='outer')
print(merged_df.filter(regex='^CTTN|^ADK'))

# for RBP in RBP_list:
         
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
#     if f'{RBP}_vs_CTRL_fdr' in merged_df.columns and f'{RBP}_pval_adj' in merged_df.columns:
#         x = (merged_df[f'{RBP}_vs_CTRL_fdr'])
#         y = (merged_df[f"{RBP}_pval_adj"])
        
#         idx = np.logical_not(np.logical_or(np.isnan(x), np.isnan(y)))
#         x = x[idx]
#         y = y[idx]

#         x = -np.log10(x)
#         y = -np.log10(y)
        
#         if len(x) > 0 and len(y) > 0:
#             a, b = np.polyfit(x, y, 1)
        
#             ax1.scatter(x,y)
#             ax1.plot(x, a*x+b, color='black')

#             ax1.set_xlabel("Old p-val fdr")
#             ax1.set_ylabel("New p_val fdr")
#             ax1.set_title(f"p_val_adj comparison {RBP}")

    # #plot for p-value
    # if f'{RBP}_vs_CTRL_pv' in merged_df.columns and f'{RBP}_log_pval' in merged_df.columns:
    #     pv_x = merged_df[f'{RBP}_vs_CTRL_pv']
    #     pv_y = merged_df[f'{RBP}_log_pval']

    #     idx_pval = np.logical_not(np.logical_or(np.isnan(pv_x), np.isnan(pv_y)))
    #     pv_x = pv_x[idx_pval]
    #     pv_y = pv_y[idx_pval]

    #     pv_x = -np.log10(pv_x)
    #     pv_y = -np.log10(pv_y)

    #     if len(pv_x) > 0 and len(pv_y) > 0:    
    #         pv_a, pv_b = np.polyfit(pv_x, pv_y, 1)
    #         ax2.scatter(pv_x,pv_y)
    #         ax2.plot(pv_x, pv_a*pv_x+pv_b, color='black')
    #         ax2.set_xlabel("-log10(pval) old data")
    #         ax2.set_ylabel("-log10(pval) new data")
    #         ax2.set_title(f"Comparison -log10(pval) {RBP}")

        # plt.savefig(f'p-val-compare/{RBP}_Plot.png')
        # plt.close(fig)
