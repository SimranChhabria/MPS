# MPS

Command example to run the MPS script <br>

    python MPS_breast_test.py -g exampleModule1.csv \ 
    -d "Breast" -m "Nick-test" \
    -cr <path_to_file>/TCGA-BRCA_0.1_primary_zscore.csv \
    -dr <path_to_file>/TCGA-BRCA_0.1_primary_zscore_bins_10.csv \
    -c <path_to_file>/TCGA-BRCA_clinical.csv \
    -s "PFS" -t "PFS_MONTHS" -e "PFS_STATUS"
