# WT1_TCR

## Overview of the analysis pipeline

&nbsp;
### **PART 1: Train a TCRex model for every WT1 epitope** (src/training)

Step 1: Collect all data from the data/raw folder into one file (collect_data.py)

The raw folder contains all MiXCR files for the two studied WT1 epitopes and two VZV epitopes. The VZV epitopes were not taken into account in this paper, but were used to build VZV-specific TCRex models for another paper [ADD reference]. Since all TCRs for the 4 epitopes were sequences together, this data was used for quality control purposes as explained in the paper. Following 4 folders are present:

- run1: TCR data from the first run for WT1-126, WT1-37 and IE62
- run1_orf18: TCR data from the first run for ORF18
- run2: WT1-37 data filtered on high and low threshold gating. 
- run 3: extra TCR data for WT1-126, WT1-37 aligned with MiXCR

Step 2: Identify all TCRs linked with multiple epitopes (get_shared_data.py)

Step 3: Remove TCRs linked with multiple epitopes from the general data set (remove_shared_data.py)

Step 4: Split the TCR data in epitope-specific files, i.e, one TCR data file for every epitope. It is possible to filter on clone count, however we used all TCR data. 
The final training data is stored in data/training_data

Step 5: Train epitope-specific TCRex models (the two WT1-specific models are available at tcrex.biodatamining.be)

&nbsp;
### **PART 2: Study publicity of training data** (src/public_cdrs)

Step 1: Generate a full dataset with all TCRs from all patients (retrieve_all_tcr_info.py)

Step 2: Retrieve all CDRs which are shared with at least two individuals (get_public_tcrs.py)
results/public_cdrs/public_126.tsv and results/public_cdrs/public_37.tsv give an overview of the public CDR3 beta sequences, their V/J genes and their count for every volunteer. CDR3 beta sequences derived from identical RNA sequences are placed within the square brackets.

Step 3: Assemble info in table (get_public_plot_info.py ) and plot (plot_public_tcrs.R)

&nbsp;
### **PART 3: Clustering of training data** (notebooks/training_data)
Clustering of the training data for WT1-37 and WT1-126, and generating sequence logos was done using notebooks/training_data/clustering_training_data.ipynb

&nbsp;
### **PART 4: Study V/J gene usage of training data** (src/vj_genes)

Step 1: Download background data set 
Download the results from Pogorelyy et al. from the iReceptor gateway at day 0.
More information about the downloaded data is present in data/background/info.tsv

Step 2: Plot the V/J gene usage of the training data and the background data set (vj_frequencies_with_background.py)

Step 3: Calculate background counts (get_background_counts.py)

Step 4: Perform gene enrichment analysis (calculate_gene_enrichment.py) and follow up with multiple testing correction (bh_correction.R)

&nbsp;
### **PART 5: Analyse data from a published AML study using the trained TCRex models** (src/aml_study)

Step 1: The original TCR repertoires can be downloaded from TCRdb (PRJNA510967)

Step 2: Parse the downloaded TCR repertoires (parse_tcrdb_aml.py)

Step 3: Combine TCR data for every volunteer in one file (aml_combine_data.py)

Step 4: Send parsed data through TCRex and make table of TCRex results (tcrex_results.py)
The raw TCRex results are stored in results/aml_study/tcrex_results since these were retrieved with an offline version.
Tables with the identified TCRs and counts are stored in results/aml_study/tcrex_tcrs.tsv and results/aml_study/tcrex_results.tsv.

Step 5: Collect response and tag info into table (aml_info_table.py)
The result is saved in Results/aml_study/info_PRJNA510967.tsv

Step 6: Perform look-up method (lookup_results.py)

Step 7: Allocate WT1-specific TCRs in clusters (aml_clustering.py) and visualize the edges with notebooks/aml_study/viz_aml_clusters.ipynb. Summarize the results with notebooks/aml_study/WT1_specificity_of_clusters.ipynb