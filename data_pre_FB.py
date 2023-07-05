#!/usr/bin/env python3
#---------------------------------------------------------------------------
#    Pre processing TCGA data : expression matrix, z-score, bin expression
#---- * Usage : python data_pre.py -i {clinical_file} -r {exp_file} -s {sample_ID_file} -d {disease_name} -o {out_dir} -----*
#
# python data_pre.py -i /rugpfs/fs0/tavz_lab/scratch/schhabria/Module_Project/brca_tcga/data_clinical_patient.txt  
#    -r /rugpfs/fs0/tavz_lab/scratch/schhabria/Module_Project/brca_tcga/data_mrna_seq_v2_rsem.txt 
#    -s /rugpfs/fs0/tavz_lab/scratch/schhabria/Module_Project/brca_tcga/data_clinical_sample.txt
#    -d "TCGA-BRCA"
#    -o "/rugpfs/fs0/tavz_lab/scratch/schhabria/Module_Project/TCGA_BRCA_processed_data"
#--------------------------------------------------------------------------


import argparse
import pandas as pd
import os
import numpy as np

""" Handle arguments """
parser = argparse.ArgumentParser(
    description="Subset the data"
)

parser.add_argument("-i", "--input_file",type=argparse.FileType('r'),
    help="Patient ID data of the entire project",
)

parser.add_argument("-r", "--rna_file",type=argparse.FileType('r'),
    help="Expression counts of the entire project",
)

parser.add_argument("-z", "--z_file",type=argparse.FileType('r'),
    help="Expression counts of the entire project",
)

parser.add_argument("-d", "--disease_name", type=str,
    help="Directory with all the tsv files"
)

parser.add_argument("-o", "--out_path", type=str,
    help="Directory with all the tsv files"
)

parser.add_argument("-gtf", "--gtf_file", type=argparse.FileType('rb'),
    help=" gtf csv file for gene_id and symbol "
)



args          = parser.parse_args()
data_file     = args.input_file
#sample_file   = args.sample_file
exp_file      = args.rna_file
disease_name  = args.disease_name
out_dir       = args.out_path

print(args.input_file)
print(args.rna_file)
#-- Reading the clincical_data from TCGA BRCA in the pandas dataframe ---#
clin_df = pd.read_csv(data_file, sep=",",low_memory=False)

#data_df = data_df.rename(columns = {'bcr_patient_barcode' : 'PATIENT_ID'})
#print(data_df)


#-- Reading the RNA_seq_rsem continous from TCGA BRCA in the pandas dataframe ---#
exp_df = pd.read_csv(exp_file, sep = "\t", low_memory=False)
# Drop the 0 index containing extra information
exp_df = exp_df.drop(labels=0, axis=0)
exp_df.reset_index(drop=True, inplace=True)
exp_df = exp_df.rename(columns = {'Hybridization REF' : 'Gene_ID'})
print(exp_df)

#--Keeping only the patients mentioned in clinical file---#
sample_list = exp_df.columns
col_list = exp_df.columns.str.split("-").str[:3].str.join("-")
sample_dict = dict(zip(col_list, sample_list))
print(sample_dict)

dis_clin = clin_df.loc[clin_df["bcr_patient_barcode"].isin(col_list)]
dis_clin["SAMPLE_ID"] = dis_clin["bcr_patient_barcode"].map(sample_dict)
print(dis_clin)


#--Keep only the columns for expression data in clinical file --#
patient_list = (dis_clin["SAMPLE_ID"])
columns_to_keep = ["Gene_ID"] + patient_list.tolist()
exp_df = exp_df[columns_to_keep]
print(exp_df)


exp_df["Gene_ID"] = exp_df["Gene_ID"].str.split('|').str[0]
exp_df = (exp_df[exp_df["Gene_ID"] != "?"])
exp_df.reset_index(drop=True, inplace=True)

#-- Convert the dataframe to numberic values
exp_df.iloc[:, 1:] = exp_df.iloc[:, 1:].applymap(lambda x: pd.to_numeric(x, errors='coerce'))
#plt.hist(exp_df["TCGA-3C-AAAU-01A-11R-A41B-07"])
#plt.savefig("Histogram raw counts.png")

#print(exp_df)

#---Read the z-score files -----#
#z_df = pd.read_csv(args.z_file, sep="\t", header=0)

# ------- FUNCTIONS ------------------#

# Function 1 : Create a function which takes the gene_symbol column name 
#            > Define the col name in command line
#            > Remove NA and average of duplicate genes
#            > Keep only unique gene_symbols and remove the extra columns

def dup_genes(exp_data, col_name):
    exp_data.dropna(subset=[col_name],inplace=True)
    duplicate_sym = exp_data[exp_data[col_name].duplicated(keep=False)]
    print ("The number of duplicated genes are", duplicate_sym.shape[0])
    exp_data= exp_data.groupby([col_name]).mean().reset_index()
    #exp_data= exp_data.drop([drop_col], axis=1)
    exp_data = exp_data.drop_duplicates(subset=col_name)
    return(exp_data)

# Function 2 : Create the output dir path
def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print("Directory created:", path)
    else:
        print("Directory already exists:", path)
        
#--------------------------------------**

#-- Average the duplicate genes and remove NAs --#
#exp_fun_df = dup_genes(exp_df, 'Gene_ID')
#print(exp_fun_df)



#-- Get the log transform (base2) of gene expression and z-score calculation --#
log_exp_df = exp_df.applymap(lambda x: np.log2(x+1) if isinstance(x, (int, float)) else x)
#plt.hist(log_exp_df["TCGA-3C-AAAU-01A-11R-A41B-07"])
#plt.savefig("Histogram log transformed counts.png")
#print(log_exp_df)


#--Get the mean and standard deviation of each gene --##
row_mean = log_exp_df.iloc[:, 1:].apply(lambda row: np.mean(row), axis=1)
row_std =  log_exp_df.iloc[:, 1:].apply(lambda row: np.std(row), axis=1)

#--Calculate the z-score for each value in the DataFrame
# -- Note : z_df = is same as data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt

z_df = log_exp_df.copy()  # Make a copy of the original DataFrame
z_df.iloc[:, 1:] = (log_exp_df.iloc[:, 1:].sub(row_mean, axis=0)).div(row_std, axis=0)

# ---- Create discrete bin expression data of sd and mean thresholds --------#
exp_bins = MI_bins_split = 10
# Define the quantile thresholds
quant_th = [0, 0.05, 0.1]

#quant_th = [0.1]
for th in quant_th:
    # Calculate the threshold values
    sd_th = np.quantile(row_mean, th)
    mu_th = np.quantile(row_std, th)
    
    #Get the gene indexes that satisfies the threshold
    g_row = z_df[(z_df.isna().sum(axis=1) == 0) & (row_std > sd_th) & (row_mean > mu_th)].index.tolist()

    # Create a dataFrame of those genes
    th_z_df = z_df.loc[g_row, :]
    print(th_z_df)

    #Get all the gene names satisfying that threshold and number of genes in each bin
    
    genes_vec = list(th_z_df.iloc[:, 0].astype(str))
    bin_value = pd.cut(range(len(genes_vec)), exp_bins, labels=False)
     
    # New dataframe to store the results
    z_bins_df = pd.DataFrame({"Approved.Symbol": genes_vec})
    # Create a dictionary to store the results for each column
    z_bins_dict = {}

    # Iterate over the columns
    for col in th_z_df.columns[1:]:
        sorted_indices = np.argsort(th_z_df[col].values)  # Sort the column values
        sorted_symbols = [genes_vec[i] for i in sorted_indices]
        
        tmp_ = pd.DataFrame({'Approved.Symbol': sorted_symbols})
        tmp_[col] = bin_value.astype('int32')
        z_bins_dict[col] = tmp_.set_index('Approved.Symbol')[col].to_dict()
        
        #z_bins_df = pd.merge(z_bins_df, tmp_, how='outer')

 
    # Create the resulting dataframe from the dictionary
    z_bins_df = pd.DataFrame(z_bins_dict).sort_index().rename_axis('Approved.Symbol').reset_index()
    print(z_bins_df)

    filename_cont = f"{disease_name}_{th}_primary_zscore.csv"
    filename_dis  = f"{disease_name}_{th}_primary_zscore_bins_10.csv"
        
    create_dir(out_dir)
    th_z_df.to_csv(os.path.join(out_dir, filename_cont), index=False)
    z_bins_df.to_csv(os.path.join(out_dir, filename_dis), index= False)


# # #--------Store the clinical data -------------#
# #--Create a list for column of TCGA patient in expression data

# sample_list = th_z_df.columns
# #print(col_list)
# col_list = th_z_df.columns.str.split("-").str[:3].str.join("-")
# sample_dict = dict(zip(col_list, sample_list))


# dis_clin = clin_df.loc[clin_df["bcr_patient_barcode"].isin(col_list)]
# dis_clin["SAMPLE_ID"] = dis_clin["bcr_patient_barcode"].map(sample_dict)


# #--Subset the patient name in from all the clinical data to only disease (BRCA) 
dis_filename = f"{disease_name}_updated_clinical.csv"
dis_clin.to_csv(os.path.join(out_dir,dis_filename), index=False)