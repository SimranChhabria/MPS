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
import shutil
import mygene
import numpy as np
import time
import pyreadr
import random
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.compare import compare_survival
from sksurv.linear_model import CoxPHSurvivalAnalysis
from scipy.stats import chi2
import matplotlib.pyplot as plt



""" Handle arguments """
parser = argparse.ArgumentParser(
    description="Subset the data"
)

parser.add_argument("-g", "--module_gene_file",type=argparse.FileType('r'),
    help="List of genes in the module",
)

parser.add_argument("-a", "--all_gene_file",type=argparse.FileType('r'), default=None,
    help="All the genes names in the experiment: either RDS or create new one",
)

parser.add_argument("-d", "--disease_name", type=str,
    help="Disease name of the analysis"
)

parser.add_argument("-m", "--module_name", type=str,
    help="Module name of the analysis"
)

parser.add_argument("-cr", "--rna_cont_file",type=argparse.FileType('r'),
    help="continous  z-score expression data of rna seq experiment",
)

parser.add_argument("-dr", "--rna_dis_file",type=argparse.FileType('r'),
    help="discrete z-score in bins data of rna seq experiment ",
)

parser.add_argument("-c", "--clinical_file",type=argparse.FileType('r'),
    help="Expression counts of the entire project",
)

parser.add_argument("-s", "--survival_analysis", type=str,
    help="Overall survival or PFS"
)

parser.add_argument("-t", "--survival_time", type=str, default=None,
    help="Overall survival or PFS time"
)

parser.add_argument("-e", "--survival_event", type=str,default=None,
    help="Overall survival or PFS event"
)

parser.add_argument("-o", "--out_path", type=str,
    help="Directory with all the tsv files"
)


parser.add_argument("-gtf", "--gtf_file", type=argparse.FileType('rb'),
    help=" gtf csv file for gene_id and symbol "
)



args                    = parser.parse_args()
module_gene_file        = args.module_gene_file
disease_name            = args.disease_name
module_name             = args.module_name
surv_type               = args.survival_analysis
surv_time               = args.survival_time
surv_event              = args.survival_event

exp_cont_file           = args.rna_cont_file
exp_dis_file            = args.rna_dis_file
out_dir                 = args.out_path
    

#------------- * Read all the inputs *------------#

# List of gene symbols in that module:
gene_list_df = pd.read_csv(module_gene_file, sep=",", header=0)

# Continous expression data and discrete expression data
exp_cont_df = pd.read_csv(exp_cont_file,sep="\t", header=0)
exp_dis_df  = pd.read_csv(exp_dis_file,sep="\t", header=0 )

# Read (rds object) or create all genes list in the experiment:
if args.all_gene_file is not None:
    all_gene_file = args.all_gene_file
    #Read a rds object in python
    result = pyreadr.read_r(all_gene_file.name)
    all_gene_df = result[None]
    print(all_gene_df)
    
else:
    all_gene_df = exp_cont_df.iloc[:, [0]]
    print(all_gene_df)

# Read Clinical file
if args.clinical_file.name.endswith("rds"):
    clin_df = args.clinical_file
    #Read a rds object in python
    result = pyreadr.read_r(clin_df.name)
    clin_df = result[None]
    print(clin_df)
    
else:
    clin_df = pd.read_csv(args.clinical_file,sep="\t", header=0)





#----------------------- FUNCTIONS ------------------------------#


# ---------------------------------
# Function : Calculate MI values
# This function is called in get_MPS_newModule twice. 
# One : To calculate MPS; Two: To shuffle and randomize the values to get the statistically significant values 
# Input : Module Vector (Genes in modules)  and 
# Output: Matrix containing MI values
# ---------------------------------

def calculateMI_v2(x,y):
    rs_cMatrix = np.add(x , y)
    num_genes = np.sum(rs_cMatrix)
    

    log_t1 = np.log((num_genes*x)/(np.sum(x) * rs_cMatrix))
    log_t1[np.isinf(log_t1)] = 0

    log_t2 = np.log((num_genes*y)/(np.sum(y) * rs_cMatrix))
    log_t2[np.isinf(log_t2)] = 0

    val_ret = np.sum((1 / num_genes) * ((x) * log_t1 + (y) * log_t2))
    return val_ret



# ---------------------------------------------
# Function: Module Perturbation Score calculation
# Inputs :  path to input = list of genes in the module
# Output:   Matrix where col = genes; rows= patients with MPS for each module (Activated, Repressed and Differential)
# -----------------------------------------------

def get_MPS_newModule(input_gene_df, all_df, dis, module_name, cont_df, dis_df, pat_df, number_bins=10, pval_thresh = 0.01, MPS_significance = 0, number_rand = 10000, MPS_thresh = 0):

  #------  Step 1: Randomization and Null Distribution -------------------#
  #  Get all the genes in the module and shuffle the genes in the module and genes not in the module. 
  #  Goal : To ensure that high MI is not by chance but is based on the data values
  #  Thus, we shuffle the number of genes in modules with genes not in module, plot null distribution!
  #----* Inputs: 
  #       g_set =  TCGA gene list or experiment gene list ; this will chance depending on the analysed data - It has to be gene symbols! 
  #       inp_g =  genes in one of the module (path_to_input provided by us!)
  #       s_set = sample IDs
  #----* Outputs:
  #     genes_in_mod_vec : List of genes in that module
  #     v_rand_q         : MPS threshold if the MPS_significance is NOT 0 coming from null distribution


       T_1 = time.time()
       #--Get the list of module genes and all genes as characters 
       inp_g       = list(input_gene_df.iloc[:, 0].astype(str))
       g_set       = list(all_df.iloc[:,0].astype(str))
       total_genes = len(g_set)

      
       
       

       #-genes_in_mod = vector of genes present in both TCGA and module
       #-div_f = Total number of genes in the module (sum of 1s in module membership vector)
       genes_in_mod = list(set(inp_g) & set(g_set))
       div_f = int(len(genes_in_mod))
       #-For storing purposes
       genes_in_mod_vec = '|'.join(sorted(set(genes_in_mod)))


       #--Generate fixed random values for null distribution
       random.seed(108)

       #--Define number of times randomization of genes in the module with genes not in the module should happen (10000 times)
       num_rand = number_rand

       #--Get the random total number of genes in each bin
       total_genes_bin = np.random.multinomial(total_genes, [1/number_bins] * number_bins, size=1)
       #--Assign number of genes in the module which is equal to genes in our module
       num_in_path = div_f
       

       #--Count table of col = sample ID, rows = bins and value is number of genes with that bin expression
       path_genes_bin = np.random.multinomial(num_in_path, [1/number_bins] * number_bins, size=num_rand).T

       #--c_vec_r = number of genes in that module and nc_vec_r = number of genes which are not the module for probabilty needed in MI
       c_vec_r = path_genes_bin 
       #print(c_vec_r[:,1])
       
       nc_vec_r = np.tile(total_genes_bin.T, (1, num_rand)) - c_vec_r 
       #print(nc_vec_r[:,1])

       #--Calculate the MI values : for null distribution to prove that high MI is not just a random occasion. Shuffling of the values = num_rand
       #--MI is calculated in this way for computational purposes! v_rand = MI values of null distribution
        
       v_rand = []
       for i in range(num_rand):
             v = calculateMI_v2(c_vec_r[:, i], nc_vec_r[:, i])
             v_rand.append(v)
       
       #--Mean and SD of MI from random distribution
       mu_r = np.mean(v_rand)
       sd_r = np.std(v_rand)

       #--pval_thresh = 99%, 
       #--the v_rand_q = threshold of the MI with mutual information value was within the top 1% of the values in null distribution
       v_rand_q = int(np.quantile(v_rand, (1-pval_thresh)))
       len_rand = len(v_rand)

       
       # --------- STEP 2: MPS Calculation ---------**
       # GOAL: 
       # --* Load the clinical data:
       #     NOTE: We always work on primary tumor data. No normal or control data is used for analysis since the focus is survival analysis
       #             con_g    = tumdat_zscore = Continous expression data used for calculating correlation
       #             dis_g    = tumdata_zscore_bins10 = Discrete expression data (in bins) for MI value
       #             dis_clin = clinical data with survival details
       # --* Get the genes which are in TCGA common genes (g_set) and map it with sample IDs
       # --* Correlation:
       #     Sign of correlation is calculated in continous data to give MI a direction (+ve or -ve) since MI are always +ve
       #            vv = con_g[sam_]               = each sampleID
       #            vv[genes_in_mod,]$bin = 0 or 1 = module membership values
       #            cor(vv)[1,2] = gets the sign of correlation value (-0.044) from cor(vv) -->
       #                                              TCGA-Z5-AAPL-01A-12R-A41B-07         bin
       #  TCGA-Z5-AAPL-01A-12R-A41B-07                   1.00000000                     -0.04458026
       #  bin                                           -0.04458026                      1.00000000
       # NOTE: ---- This is the end of using continuous expression data
       # --* MI : 
       #     c_vec  = vector of length 10 with number of genes in that module in each bin
       #     nc_vec = vector of length 10 with number of genes NOT in that module in each bin
       #     If you add both the vectors, get total number of genes in that bin
       #     Calculate MI and sign of MI
       # --* MPS:
       #     id_x = threshold defined above to filter the values. So ix_0 is the index where MI < threshold (v_rand_q)    
       #     v_z  = calculate MI values z score and filter negative z scores of MI since they are not statistically significant
       #     MPS = v_z*sign(MI values) to get the +ve or -ve direction of MPS since MI is always positive
       # -------------------------------------------------------------------------------------------------------
       
       #-- Intersect common genes in continous data with discrete data and with common genes in all_genes_df
       g_row = np.intersect1d(g_set, np.intersect1d(cont_df.iloc[:, 0], dis_df.iloc[:, 0]))
       
       #-- Sample IDs excluding the 1st column which is gene names 
       s_set = np.intersect1d(cont_df.columns[1:], dis_df.columns[1:])   
       
       # --Ensure the order of the dataframe
       cont_df = cont_df.set_index(cont_df.columns[0])
       dis_df = dis_df.set_index(dis_df.columns[0])

       cont_df = cont_df.loc[g_row, s_set]
       dis_df  = dis_df.loc[g_row, s_set]
       #--total number of genes 
       Nrow    = len(cont_df.iloc[:, 0])
       per_bin = round(Nrow/number_bins)
       
       
       #--Calculate correlation with continous data and MI with discrete data
       out_list = []
       for sample_id in s_set:
               # Calculate sign of correlation on continuous data:
               cvalue_df = cont_df[[sample_id]]
               cvalue_df['bin'] = 0
               cvalue_df.loc[genes_in_mod, 'bin'] = 1
               correlation_matrix = cvalue_df.corr()
               sign_fact = np.sign(correlation_matrix.iloc[0, 1])
               #print(sign_fact)
               

               # Discrete data sample IDs
               dvalue_df = dis_df[[sample_id]]
               #List of number of genes and the names of all the genes in each bin 
               each_mod_genes = [dvalue_df.index[np.where(dvalue_df.values[:, 0] == i)].tolist() for i in range(1, number_bins + 1)]
               df_each_mod = pd.DataFrame({'Bin': range(1, number_bins + 1),
                            'Num_Genes': [len(genes) for genes in each_mod_genes],
                            'Genes': each_mod_genes})
               #Get the number of genes in each module
               c_vec = np.array([len(np.intersect1d(genes_in_mod, genes)) for genes in df_each_mod['Genes']])
               nc_vec = np.array([len(genes) - c for genes, c in zip(each_mod_genes, c_vec)])
               #print(isnan(c_vec))
               mi_ = calculateMI_v2(c_vec, nc_vec)
               mi_sign = sign_fact * mi_
               out_list.append(mi_sign)

       mi_vec = np.array(out_list)
       abs_mi = np.abs(mi_vec)
       #--z-score of MI values
       v_z = (abs_mi - mu_r) / sd_r
       v_z[v_z < 0] = 0
       v_z[np.isnan(v_z)] = 0
       
       #-- id_x = threshold defined above to filter the values. So ix_0 is the index where MI < threshold (v_rand_q)
       if MPS_significance == 1:
            ix_0 = np.where(abs_mi < v_rand_q)[0]
            v_z[ix_0] = 0
       
       #--MPS = MI * sign of correlation giving MPS a direction 
       v_z = v_z*np.sign(mi_vec)

       #--Merge clinical data with MPS values:
       #--TCGA Pancreatic ---#
       #tmp = pd.DataFrame({'SAMPLE_ID': s_set, 'MPS': v_z})
       #--TCGA Breast ----
       tmp = pd.DataFrame({'SAMPLE_ID': s_set, 'MPS': v_z})
       print(pat_df)
       tmp.index = s_set
       
       #tmp = pd.merge(tmp, pat_df, on = "SAMPLE_ID")
       tmp = pd.merge(tmp, pat_df, on = "SAMPLE_ID")

       tmp['module_name'] = module_name
       tmp['module_id'] = 'user_module ' + module_name
       tmp['genes_in_module'] = genes_in_mod_vec
       tmp['num_genes_in_module'] = div_f
       tmp['MPS_groups'] = np.nan
       

       # Get MPS+ and MPS- based on the threshold
       tmp.loc[tmp['MPS'] > MPS_thresh, 'MPS_groups'] = 'MPS+'
       tmp.loc[tmp['MPS'] < (-1 * MPS_thresh), 'MPS_groups'] = 'MPS-'
       
       print(len(g_set))
       print(len(g_row))
       print(cont_df)
       
       return (tmp)

       #tmp.to_csv(out_mod_dir + dis_ + '_MPS.tsv', sep='\t', index=False)




# ---------- PLOTS -------------------##

#------------------------------------
# Function : getSurv 
# Parent Function :plot_Surv
#            Get the survival plots 
#
#--------------------------------------

def getSurv(all_clin,surv_type,surv_time,surv_event,module_name,disease_name,rand_iter = 1000, samp_name = 'MPS+', ctrl_name = 'MPS-'):
           all_clin[surv_time] = all_clin[surv_time].astype(float)
           all_clin[surv_event] = all_clin[surv_event].astype(bool)

        #    fig, ax = plt.subplots()

        #    #Number of MPS+ and MPS-
        #    len_s = np.sum(all_clin['MPS_groups'] == samp_name)
        #    len_c = np.sum(all_clin['MPS_groups'] == ctrl_name)

        #    #Survival analysis
        #    for MPS_type in ("MPS+", "MPS-"):
        #          MPS = all_clin["MPS_groups"] == MPS_type
        #          time, survival_event = kaplan_meier_estimator(
        #             all_clin[surv_event][MPS],
        #             all_clin[surv_time][MPS]
                    
        #         )
        #          plt.step(time, survival_event, where="post",label=f'{MPS_type}')

        #    plt.ylabel("Probability of survival")
        #    plt.xlabel("Time")
        #    plt.legend(loc="best")
        #    plt.savefig(f'{surv_type}_{disease_name}_{module_name}.png')


           structured_array = np.zeros(len(all_clin), dtype=[('event', '?'), ('time', int)])
           structured_array['event'] = all_clin[surv_event].astype(bool)
           structured_array['time']  = all_clin[surv_time].astype(int)

           group_indicator = all_clin.loc[:, 'MPS_groups']
           groups = group_indicator

           chi2, pvalue= compare_survival(structured_array, groups)
           print(pvalue)


           # Convert groups to 0s and 1s
           groups_encoded = np.where(groups == "MPS+", 1, 0)

           # Reshape groups to a 2D array with a single column
           groups_encoded = groups_encoded.reshape(-1, 1)

           # Fit the CoxPHSurvivalAnalysis estimator
           estimator = CoxPHSurvivalAnalysis()
           estimator.fit(groups_encoded, structured_array)
           hazard_ratios = np.exp(estimator.coef_)
           print(hazard_ratios)

           out_df = pd.DataFrame({'pval_surv': pvalue, 'HZR_cox': hazard_ratios})
           


           # Separate the data for MPS+ and MPS-
           mps_plus_data = all_clin[all_clin["MPS_groups"] == "MPS+"]
           mps_minus_data = all_clin[all_clin["MPS_groups"] == "MPS-"]

           # Perform Kaplan-Meier estimation for MPS+
           time_mps_plus, survival_event_mps_plus = kaplan_meier_estimator(mps_plus_data[surv_event].astype(bool),mps_plus_data[surv_time].astype(int))

           # Perform Kaplan-Meier estimation for MPS-
           time_mps_minus, survival_event_mps_minus = kaplan_meier_estimator(mps_minus_data[surv_event],mps_minus_data[surv_time])

           # Plot the survival curves
           plt.step(time_mps_plus, survival_event_mps_plus, where="post", label="MPS+", color="red")
           plt.step(time_mps_minus, survival_event_mps_minus, where="post", label="MPS-", color="blue")

           # Customize the plot as needed (e.g., axes labels, title, legend)
           plt.xlabel("Time")
           plt.ylabel("Survival Probability")
           plt.title(f"Module= {module_name} pval = {pvalue} , HR = {hazard_ratios}")
           plt.legend()
           plt.savefig(f'{surv_type}_{disease_name}_{module_name}.png')

           return (out_df)




res = get_MPS_newModule(gene_list_df, all_gene_df, disease_name, module_name, exp_cont_df, exp_dis_df, clin_df)
print(res['MPS'])
#Clinical data 
all_clin = res.dropna(subset=['MPS_groups'])
all_clin.to_csv("Brain_MPS.csv",index=None)
out_df = getSurv(all_clin,surv_type,surv_time,surv_event,module_name,disease_name)
print(out_df)
print(np.sum(all_clin['MPS_groups']=="MPS+"))
print(np.sum(all_clin['MPS_groups']=="MPS-"))



