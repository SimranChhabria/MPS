#!/usr/bin/env python3
#---------------------------------------------------------------------------
#    To subset clinical data based on RNA-seq sampleID
#---- * Usage : python3 sort_python.py -i {input_file} -o {output_file} -----*
#
# python ENCODE.py -i /rugpfs/fs0/tavz_lab/scratch/schhabria/Module_Project/shRNA-ENCODE/metadata.tsv  
#    -o /rugpfs/fs0/tavz_lab/scratch/schhabria/Module_Project/shRNA-ENCODE/out.csv 
#    -d /rugpfs/fs0/tavz_lab/scratch/schhabria/Module_Project/shRNA-ENCODE/test-run/ 
#    -a /rugpfs/fs0/tavz_lab/scratch/schhabria/Module_Project/shRNA-ENCODE/out_avg.csv 
#    -gtf /rugpfs/fs0/tavz_lab/scratch/schhabria/Module_Project/geneid-genesym.csv
#--------------------------------------------------------------------------

import argparse
import pandas as pd
import os
import shutil
import mygene
from biomart import BiomartServer
import gzip
from functools import reduce
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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

parser.add_argument("-n", "--newdir_path", type=str,
    help="New directory to store all the renamed the tsv files"
)


parser.add_argument("-gtf", "--gtf_file", type=argparse.FileType('rb'),
    help=" gtf csv file for gene_id and symbol "
)

parser.add_argument("-o", "--output_file", type=argparse.FileType('w'), 
        help="Output file"
        )

parser.add_argument("-a", "--avg_file", type=argparse.FileType('w'), 
        help="Avg Output file"
        )


args          = parser.parse_args()
meta_file     = args.input_file
out_file      = args.output_file

#-- Reading the meta_tsv from ENCORE DATA in the pandas dataframe ---#
meta_df = pd.read_csv(meta_file, sep="\t", header=0)

#-- Filter all the bigWig files -----#
def filter_rows_by_values(df, col, values):
    return df[~df[col].isin(values)]

def empty_col_val_control(df,col) :
    df.loc[df[col].isnull(), col] = "Control"
    return df

meta_df = filter_rows_by_values(meta_df, "File format", ["bigWig"]) 
metacon_df = empty_col_val_control(meta_df, "Experiment target")


if args.dir_path is not None:
    # Get all tsv files in the directory
    files = [f for f in os.listdir(args.dir_path) if f.endswith('.tsv')]
    #print(files)
else:
    print('Please provide a directory path using the -d option')

# Loop through the file list and match each file name with a row in the DataFrame
for file_name in files:
    # Find the row in the DataFrame that matches the file name
    file_namet = file_name.replace('.tsv', '')
    #print(file_namet)
    row = metacon_df[metacon_df['File accession'] == str(file_namet)]
    #print(row)
      
    if not row.empty:
        # get the new name from the corresponding column in the row
        new_name = row.iloc[0]["Experiment target"] + str("_") + row.iloc[0]["Experiment accession"] + ".tsv"
        print(new_name)
       
        src_path = os.path.join(args.dir_path, file_name)
        dst_path = os.path.join(args.newdir_path, file_name)
        print(dst_path)
        shutil.copy(src_path, dst_path)
    # Rename the file using the new name      
        os.rename(dst_path, os.path.join(args.newdir_path, new_name))
        

dfs_list = []
files_new = [f for f in os.listdir(args.newdir_path) if f.endswith('.tsv')]
for new_file in files_new:

     new_name = new_file.replace('.tsv', '')
     new_name = str(new_name) + '_TPM'
     new_name = new_name.replace("-human_", "_")
     df_tmp = pd.read_csv(os.path.join(args.newdir_path, new_file), sep='\t')
     print(df_tmp.head(5))
     df_tmp = df_tmp.loc[:, ['gene_id','TPM']]
     df_tmp = df_tmp.rename(columns={'TPM':new_name})
     dfs_list.append(df_tmp)
     
#Merge the dataframe on gene_id of all the tsv files
merged_df = reduce(lambda left, right: pd.merge(left, right, on='gene_id', how='outer'), dfs_list)
#merged_df.to_csv(out_file, index= None)

numeric_cols = merged_df.columns[merged_df.columns != 'gene_id']
merged_df[numeric_cols] = merged_df[numeric_cols].apply(pd.to_numeric)

# Read the gtf.csv to dataframe and map the gene ids to gene symbols (source: gtf_genesym.py)
df_gtf = pd.read_csv(args.gtf_file,sep=",", header=0)
df_gtf = df_gtf.rename(columns = {'gene_id' : 'gene_id_short'} )


# Split gene_id at the decimal point and take the first element
merged_df['gene_id_short'] = merged_df['gene_id'].str.split('.').str[0]

#Map ensemble IDs to their corresponding symbols in the dataframe from the df_gtf
mergedgenesym_df = pd.merge(merged_df, df_gtf , on = 'gene_id_short', how = 'outer')
mergedgenesym_df.dropna(subset=['gene_name'],inplace=True)
duplicate_genes = mergedgenesym_df[mergedgenesym_df['gene_name'].duplicated(keep=False)]

#print(merged_df[merged_df['gene_id']=='ENSG00000143248.12'])
#print(mergedgenesym_df[mergedgenesym_df['gene_id_short']=='ENSG00000143248'])
#print(duplicate_genes[duplicate_genes['gene_name'] == 'RGS5'])


# For duplicate query terms, take the average of their corresponding values
mergedgenesym_df= mergedgenesym_df.drop(['gene_id', 'gene_id_short'], axis=1)
#mergedgenesym_df.to_csv(out_file,index= None)

#mergedgenesym_df.set_index('gene_name', inplace=True)
df_genenames = mergedgenesym_df.groupby(['gene_name']).mean().reset_index()
print(df_genenames[df_genenames['gene_name']=='NUP35'])
df_genenames.to_csv(out_file, index= None)


#Average of controls and common genes in the merged dataset
df_genenames.set_index('gene_name', inplace=True)
df_avg = df_genenames.groupby(df_genenames.columns.str.split('_').str[0], axis=1).mean()
df_avg.reset_index(inplace=True)

#print(df_avg[df_avg['gene_name']=='NUP35'])
#Get the log2fold change
control = 'Control'
genes = df_avg.columns[df_avg.columns != 'gene_name']
df_l2fc = df_avg[['gene_name']].copy()

for gene in genes:
     if 'Control' in gene:
        continue
     col_name = gene.split('_')[0] + '_log2fc'
     df_l2fc[col_name] = np.log2((df_avg[gene] + 1) /(df_avg[control] + 1) )

df_l2fc.to_csv(args.avg_file,index=None)
print(df_l2fc)

log2fc_df_new = pd.DataFrame(columns=['gene_name', 'log2fold'])
# Loop over the genes and extract the log2fold change value
for gene in genes:
    if 'Control' in gene:
        continue
    col_name = gene + '_log2fc'
    log2fc = df_l2fc.loc[df_l2fc['gene_name'] == gene, col_name]
    if not log2fc.empty:
        log2fc_df_new.loc[len(log2fc_df_new)] = {'gene_name': gene, 'log2fold': log2fc.values[0]}

#log2fc_df_new.to_csv(args.avg_file,index=None)
print(log2fc_df_new)
plt.hist(log2fc_df_new['log2fold'], bins=10,alpha=0.5)
# Add labels and title
plt.xlabel('Log2 Fold Change')
plt.ylabel('Frequency')
plt.title('Histogram of Log2 Fold Change')
plt.savefig('hist5.png')


# ------ Extra code ------------# 

# Initialize a PDF file
# with PdfPages('log2fold_histograms2.pdf') as pdf:
#     # Plot a histogram for each log2fold value
#     for col in df_l2fc.columns[1:]:
#         fig, ax = plt.subplots()
#         ax.hist(df_l2fc[col], bins=10, alpha=0.5)
#         ax.set_title(col)
#         ax.set_xlabel('log2fold')
#         ax.set_ylabel('Frequency')
#         # Add the histogram to the PDF file
#         pdf.savefig(fig)
#         # Close the figure to free up memory
#         plt.close(fig)




# for gene_name in df_l2fc["gene_name"].unique():
#     column_name = str(gene_name) + "_log2fc"
#     if column_name in df_l2fc.columns:
#         plt.hist(df_l2fc[df_l2fc["gene_name"]==gene_name][str(gene_name) + "_log2fc"], alpha=0.5)

# plt.xlabel("log2fc")
# plt.ylabel("Frequency")
# plt.title("Histogram of log2fc for all gene names")
# #plt.legend(loc='upper right', bbox_to_anchor=(0.5, 0.3), fontsize='small')
# plt.savefig('hist3.png')
# plt.show()
# log2fold_values = df_avg['SRFBP1_log2fc']
# plt.hist(log2fold_values, bins=10)
# plt.xlabel('Log2fold values')
# plt.title('Histogram of Log2fold values of SRFBP1_Experiment')
# plt.savefig('hist.png')