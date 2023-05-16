Short guide as to how figures & tables related to the threshold & LOBICO analyses in Chang et al. 2023 were created

---

Fig. 1E & Fig. 4A: run threshold.py (creates various image file outputs of the form AUC_f?_t0.1_log10.pdf)

Fig. 4B: run threshold.py (creates text file output Samples_t0.1_log10.csv from which Fig. 4B can be created)

Fig. 5B: 
1) Make sure LOBICO & Matlab are installed
2) Run threshold.py (creates AUC_t0.1_log10.csv which contains the binarization threshold)
3) Run LOBICO_crossval.m (runs LOBICO and outputs Validation_*.csv and Testing_*.csv files)
4) run plot_ROC_cross.py (creates pdf images)

Fig. 5C
1) run steps 1)-3) from Fig. 5B
2) run tables_ROC_cross.py (creates ROC_formulae_test_t0.1_log10_limit4.csv); file contains the Pareto optimal formulae (according to the test set) for each drug and training-test data split
3) run gene_importance.py (creates file gene_importance_t0.1_log10_limit4_ew.csv which provides the imporance (weight) for each gene)
4) run plot_gene_importance.py (creates the pdf image)

---

Supplementary Fig. 2: same as for Fig. 4B

Supplementary Fig. 3: same as for Fig. 5B

---

Supplementary Table 3: run threshold.py (creates text file AUC_t0.1_log10.csv from which Suppl. Table 3 can be created)

Supplementary Table 4: run threshold.py and look at text file output Samples_t0.1_log10.csv

Supplementary Table 6: Supplementary Table 6 corresponds to file ROC_formulae_test_t0.1_log10_limit4.csv, see discussion of Fig. 5C above

Supplementary Table 7: Supplementary Table 7 corresponds to file gene_importance_t0.1_log10_limit4_ew.csv, see discussion of Fig. 5C above
