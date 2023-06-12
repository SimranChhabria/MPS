#!/usr/bin/env python3
#---------------------------------------------------------------------------
#    To get the statistics value from the output of ENCODE.py script 
#---- * Usage : python3 sort_python.py -i {input_file} -o {output_file} -----*



import argparse
import pandas as pd
import numpy as np
import os
import itertools


""" Handle arguments """
parser = argparse.ArgumentParser(
    description="Subset the data"
)

parser.add_argument("-i", "--input_file",type=argparse.FileType('r'),
    help="Meta data of the entire project",
)

parser.add_argument("-d", "--dir_path", type=str,
    help="Directory with all the tsv files"
)

args          = parser.parse_args()
input_file    = args.input_file

df_stats = pd.read_csv(input_file, sep=",")

activated_path = os.path.join(args.dir_path, "activated")
if not os.path.isdir(activated_path):
    os.mkdir(activated_path)

repressed_path = os.path.join(args.dir_path, "repressed")
if not os.path.isdir(repressed_path):
    os.mkdir(repressed_path)

differential_path = os.path.join(args.dir_path, "differential")
if not os.path.isdir(differential_path):
    os.mkdir(differential_path)

# #---------- TPM Counts ------------------#

#Get the t-test of all the genes for that RBP Experiment using TPM values:
RBP_list= []
for col in df_stats.columns[1:]:
     x = col.split('_', 1)[0]  
     if (x!=('Control')) and x not in RBP_list:
        RBP_list.append(x)


#RBP_list = ['CTTN']

print(df_stats)
# Define the log fold change and adjusted p-value thresholds to filter the data
thresholds = {'log2fold': [0], 'adj_pval': [0.01, 0.1, 0.05, 0.001]}

for RBP in RBP_list:
    if RBP == "log":
        continue
    for lfc in thresholds['log2fold']:
        for pval in thresholds['adj_pval']:

            #-------- PART 1: t-test ----------------# 
            
            #-----ACTIVATED TARGETS OF RBP (log2fold > threshold) for t-test p_adj -------------#
            #-----Up regulated genes when the RBP was perturbed (log2fold > threshold) -------#
            # Filter the data based on the current log fold change and adjusted p-value thresholds (Lesser pval, greater significance)
            tactivate_df = df_stats[(df_stats[f'log_{RBP}_log2change'] > lfc) & (df_stats[f'{RBP}_pval_adj'] < pval)]
            
            # Add the gene_name for the filtered rows
            tactivate_df.loc[:, f'{RBP}_gene_name'] = df_stats.loc[tactivate_df.index, 'gene_name']
            
            # Keep only the desired columns in the filtered dataframe
            tactivate_df = tactivate_df.loc[:, [f'{RBP}_gene_name', f'{RBP}_pval_adj', f'log_{RBP}_log2change']]
            
            # Save the filtered data to a CSV file
            file_name = f'{RBP}_log{lfc}_pval{pval}.csv'
            RBP_activated_path = os.path.join(activated_path,"t-test",f"{RBP}")
            if not os.path.isdir(RBP_activated_path):
                os.makedirs(RBP_activated_path)
            tactivate_df.to_csv(os.path.join(RBP_activated_path, file_name), index=None)


            #-----REPRESSED TARGETS OF RBP (log2fold < -(threshold) for t-test p_adj -------------#
            #-----Down regulated genes when the RBP was perturbed log2fold < -(threshold) -------#
            trepressed_df =  df_stats[(df_stats[f'log_{RBP}_log2change'] < -(lfc)) & (df_stats[f'{RBP}_pval_adj'] < pval)]
            # Add the gene_name for the filtered rows
            trepressed_df.loc[:, f'{RBP}_gene_name'] = df_stats.loc[trepressed_df.index, 'gene_name']
            
            # Keep only the desired columns in the filtered dataframe
            trepressed_df = trepressed_df.loc[:, [f'{RBP}_gene_name', f'{RBP}_pval_adj', f'log_{RBP}_log2change']]
            
            # Save the filtered data to a CSV file
            file_name = f'{RBP}_log{lfc}_pval{pval}.csv'
            RBP_repressed_path = os.path.join(repressed_path,"t-test",f"{RBP}")
            if not os.path.isdir(RBP_repressed_path):
                     os.makedirs(RBP_repressed_path)
            trepressed_df.to_csv(os.path.join(RBP_repressed_path, file_name), index=None)


            #-----DIFFERENTIAL EXPRESSED TARGETS OF RBP (log2fold < -(threshold) or log2fold > (threshold) for t-test p_adj -------------#        
            tdifferential_df = df_stats[((df_stats[f'log_{RBP}_log2change'] > lfc) | (df_stats[f'log_{RBP}_log2change'] < -(lfc))) & (df_stats[f'{RBP}_pval_adj'] < pval)]
            
            # Add the gene_name for the filtered rows
            tdifferential_df.loc[:, f'{RBP}_gene_name'] = df_stats.loc[tdifferential_df.index, 'gene_name']
            
            # Keep only the desired columns in the filtered dataframe
            tdifferential_df = tdifferential_df.loc[:, [f'{RBP}_gene_name', f'{RBP}_pval_adj', f'log_{RBP}_log2change']]
            
            # Save the filtered data to a CSV file
            file_name = f'{RBP}_log{lfc}_pval{pval}.csv'
            RBP_differential_path = os.path.join(differential_path,"t-test",f"{RBP}")
            if not os.path.isdir(RBP_differential_path):
                     os.makedirs(RBP_differential_path)
            tdifferential_df.to_csv(os.path.join(RBP_differential_path, file_name), index=None)


           #-------- PART 2: Mann Whitney ----------------# 

           #-----ACTIVATED TARGETS OF RBP (log2fold > threshold) for Mann Whitney p_adj -------------#
            # Filter the data based on the current log fold change and adjusted p-value thresholds
            mactivate_df = df_stats[(df_stats[f'log_{RBP}_log2change'] > lfc) & (df_stats[f'{RBP}_MW_pval_adj'] < pval)]
            
            # Add the gene_name for the filtered rows
            mactivate_df.loc[:, f'{RBP}_gene_name'] = df_stats.loc[mactivate_df.index, 'gene_name']
            
            # Keep only the desired columns in the filtered dataframe
            mactivate_df = mactivate_df.loc[:, [f'{RBP}_gene_name', f'{RBP}_MW_pval_adj', f'log_{RBP}_log2change']]
            
            # Save the filtered data to a CSV file
            file_name = f'{RBP}_log{lfc}_pval{pval}.csv'
            RBP_activated_path = os.path.join(activated_path,"Man_Whitney",f"{RBP}")
            if not os.path.isdir(RBP_activated_path):
                     os.makedirs(RBP_activated_path)
            mactivate_df.to_csv(os.path.join(RBP_activated_path, file_name), index=None)

            #-----REPRESSED TARGETS OF RBP (log2fold < -(threshold) for Man Whitney p_adj -------------#
            mrepressed_df =  df_stats[(df_stats[f'log_{RBP}_log2change'] < -(lfc)) & (df_stats[f'{RBP}_MW_pval_adj'] < pval)]
            # Add the gene_name for the filtered rows
            mrepressed_df.loc[:, f'{RBP}_gene_name'] = df_stats.loc[mrepressed_df.index, 'gene_name']
            
            # Keep only the desired columns in the filtered dataframe
            mrepressed_df = mrepressed_df.loc[:, [f'{RBP}_gene_name', f'{RBP}_MW_pval_adj', f'log_{RBP}_log2change']]
            
            # Save the filtered data to a CSV file
            file_name = f'{RBP}_log{lfc}_pval{pval}.csv'
            RBP_repressed_path = os.path.join(repressed_path,"Mann-Whitney",f"{RBP}")
            if not os.path.isdir(RBP_repressed_path):
                     os.makedirs(RBP_repressed_path)
            mrepressed_df.to_csv(os.path.join(RBP_repressed_path, file_name), index=None)


            #-----DIFFERENTIAL EXPRESSED TARGETS OF RBP (log2fold < -(threshold) or log2fold > (threshold) for Mann-Whitney p_adj -------------#        
            mdifferential_df = df_stats[((df_stats[f'log_{RBP}_log2change'] > lfc) | (df_stats[f'log_{RBP}_log2change'] < -(lfc))) & (df_stats[f'{RBP}_MW_pval_adj'] < pval)]
            
            # Add the gene_name for the filtered rows
            mdifferential_df.loc[:, f'{RBP}_gene_name'] = df_stats.loc[mdifferential_df.index, 'gene_name']
            
            # Keep only the desired columns in the filtered dataframe
            mdifferential_df = mdifferential_df.loc[:, [f'{RBP}_gene_name', f'{RBP}_MW_pval_adj', f'log_{RBP}_log2change']]
            
            # Save the filtered data to a CSV file
            file_name = f'{RBP}_log{lfc}_pval{pval}.csv'
            RBP_differential_path = os.path.join(differential_path,"Mann-Whitney",f"{RBP}")
            if not os.path.isdir(RBP_differential_path):
                     os.makedirs(RBP_differential_path)
            mdifferential_df.to_csv(os.path.join(RBP_differential_path, file_name), index=None)





