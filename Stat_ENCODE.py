#!/usr/bin/env python3
#---------------------------------------------------------------------------
#    To get the statistics value from the output of ENCODE.py script 
#---- * Usage : python3 sort_python.py -i {input_file} -o {output_file} -----*



import argparse
import pandas as pd
import numpy as np
import os
import itertools
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu

""" Handle arguments """
parser = argparse.ArgumentParser(
    description="Subset the data"
)

parser.add_argument("-i", "--input_file",type=argparse.FileType('r'),
    help="Meta data of the entire project",
)

parser.add_argument("-gtf", "--gtf_file", type=argparse.FileType('rb'),
    help=" gtf csv file for gene_id and symbol "
)

parser.add_argument("-o", "--output_file", type=argparse.FileType('w'), 
        help="Output file"
        )

parser.add_argument("-l", "--log2fold_file", type=argparse.FileType('w'), 
        help="Output file"
        )

args          = parser.parse_args()
input_file    = args.input_file
out_file      = args.output_file


TPM_df = pd.read_csv(input_file, sep=",")

# create a dictionary to store the counts of each gene_name
counts = {}
for col in TPM_df.columns[1:]:
    gene_name = col.split('_', 1)[0]    
    # check if the column starts with 'gene_name'
    if gene_name in counts:
        # get the current count for this gene_name
        counts[gene_name] += 1
        count = counts[gene_name]
        # rename the column with a suffix of _count
        new_col = f"{gene_name}_{count}"
        TPM_df.rename(columns={col: new_col}, inplace=True)
        # increment the count for this gene_name     
        #print(f"Incremented count for {gene_name} to {counts[gene_name]}")
    else:
        #print(f"{gene_name} is not in counts")
        # add the gene_name to the counts dictionary
        counts[gene_name] = 1
        count = counts[gene_name]
        new_col = f"{gene_name}_{count}"
        TPM_df.rename(columns={col: new_col}, inplace=True)


# display the new column names


#Get rid of all the only one TPM Experiment no statistical interpretation
new_TPM_df = TPM_df.copy()
for col in TPM_df.columns[1:]:
    gene_name = col.split('_', 1)[0]
    if gene_name + '_2' not in TPM_df.columns:
        new_TPM_df.drop(col, axis=1, inplace=True)


# Read the gtf.csv to dataframe and map the gene ids to gene symbols (source: gtf_genesym.py)
df_gtf = pd.read_csv(args.gtf_file,sep="\t")
df_gtf = df_gtf.rename(columns = {'Approved.Symbol' : 'gene_name'} )

#Map ensemble IDs to their corresponding symbols in the dataframe from the df_gtf
mergedgenesym_df = pd.merge(TPM_df, df_gtf , on = 'gene_name', how="right")
mergedgenesym_df.dropna(subset=['gene_name'],inplace=True)
duplicate_genes = mergedgenesym_df[mergedgenesym_df['gene_name'].duplicated(keep=False)]

mergedgenesym_df = mergedgenesym_df.drop(['Ensembl.ID', 'Entrez.ID','names'], axis=1)
print(mergedgenesym_df)


#mergedgenesym_df.to_csv('CTTN_Raw_Df.csv', index=None)
new_TPM_df = mergedgenesym_df.copy()
logTPM_df  = mergedgenesym_df.copy()

new_TPM_df.to_csv(args.log2fold_file, index=None)

#--- Create copy for log(TPM+1) ---------#
logTPM_df = new_TPM_df.copy()


# #---------- TPM Counts ------------------#

#Get the t-test of all the genes for that RBP Experiment using TPM values:
RBP_list= []
for col in new_TPM_df.columns[1:]:
     x = col.split('_', 1)[0]  
     if (x!=('Control')) and x not in RBP_list:
        RBP_list.append(x)


for RBP in RBP_list:
    col_vals    = [col for col in new_TPM_df.columns if col.startswith(RBP)]
    control_vals = [col for col in new_TPM_df.columns if col.startswith('Control')]
    for index, row in new_TPM_df.iterrows():
        rbp_tpm = row.loc[col_vals].dropna().values.astype(float)
        control_tpm = row.loc[control_vals].dropna().values.astype(float)
        t_stat, p_value = ttest_ind(rbp_tpm,control_tpm,equal_var=False)
        new_TPM_df.loc[index, f'{RBP}_tstat'] = t_stat
        new_TPM_df.loc[index, f'{RBP}_pval']  = p_value

        #Calcute Mann-Whitney Stats:
        if len(rbp_tpm) > 0 and len(control_tpm) > 0:
                 statistic, pvalue = mannwhitneyu(rbp_tpm, control_tpm)
                 new_TPM_df.loc[index,f'{RBP}_Mann_whitney_pval' ] = pvalue
                 new_TPM_df.loc[index,f'{RBP}_Mann_whitney_stat' ] = statistic
        else:
            # handle the case where one or both arrays are empty
            # for example, you could set pvalue and statistic to NaN
                 new_TPM_df.loc[index,f'{RBP}_Mann_whitney_pval' ] = np.nan
                 new_TPM_df.loc[index,f'{RBP}_Mann_whitney_stat' ] = np.nan

        #log2change_tpm = row.loc[col_vals].values.mean()
        #log2change_control = row.loc[control_vals].values.mean()
        log2change_tpm = np.median(row.loc[col_vals].values)
        log2change_control = np.median(row.loc[control_vals].values)
        log_value =  np.log2( (log2change_tpm + 1) / ( log2change_control + 1) )
        new_TPM_df.loc[index,f'{RBP}_log2change' ] = log_value
        
#print(new_TPM_df.head)
#new_TPM_df.to_csv('stat_raw.csv', index=None)

#--------- log(TPM + 1 ) Counts --------------#
logTPM_df = logTPM_df.applymap(lambda x: np.log2(x+1) if isinstance(x, (int, float)) else x)

for RBP in RBP_list:
    col_vals    = [col for col in logTPM_df.columns if col.startswith(RBP)]
    control_vals = [col for col in logTPM_df.columns if col.startswith('Control')]
    for index, row in logTPM_df.iterrows():
        rbp_tpm = row.loc[col_vals].dropna().values.astype(float)
        control_tpm = row.loc[control_vals].dropna().values.astype(float)
        t_stat, p_value = ttest_ind(rbp_tpm,control_tpm,equal_var=False)
        logTPM_df.loc[index, f'{RBP}_log_tstat'] = t_stat
        logTPM_df.loc[index, f'{RBP}_log_pval']  = p_value

        #Calcute Mann-Whitney Stats:
        if len(rbp_tpm) > 0 and len(control_tpm) > 0:
                 statistic, pvalue = mannwhitneyu(rbp_tpm, control_tpm)
                 logTPM_df.loc[index,f'{RBP}_log_Mann_whitney_pval' ] = pvalue
                 logTPM_df.loc[index,f'{RBP}_log_Mann_whitney_stat' ] = statistic
        else:
            # handle the case where one or both arrays are empty
            # for example, you could set pvalue and statistic to NaN
                 logTPM_df.loc[index,f'{RBP}_log_Mann_whitney_pval' ] = np.nan
                 logTPM_df.loc[index,f'{RBP}_log_Mann_whitney_stat' ] = np.nan


        #log2change_tpm = row.loc[col_vals].values.mean()
        #log2change_control = row.loc[control_vals].values.mean()
        
        #According to Balaji's Script: Median
        log2change_tpm = np.median(row.loc[col_vals].values)
        log2change_control = np.median(row.loc[control_vals].values)

        #Since these are log transformed values already, we would subtract the values instead of log_value = np.log2( (log2change_tpm + 1) / ( log2change_control + 1) )
        log_value = log2change_tpm - log2change_control
        logTPM_df.loc[index,f'log_{RBP}_log2change' ] = log_value


#-------Subset only t_test and p_value columns --------#

pval_cols = [col for col in new_TPM_df.columns if col.endswith(('name','pval','stat','log2change'))]
p_val_df = new_TPM_df.loc[:, new_TPM_df.columns.isin(pval_cols)]


log_pval_cols = [col for col in logTPM_df.columns if col.endswith(('name','pval','stat','log2change'))]
log_pval_df = logTPM_df.loc[:, logTPM_df.columns.isin(log_pval_cols)]

#---------- Merge the pandas dataframe --------#

final_df = pd.merge(p_val_df,log_pval_df, on='gene_name', how='outer')
final_df.to_csv(out_file,index=None)




  
