# Imports
import datetime
import os
import pandas as pd

import general_functions as gf

# Adjust the working directory
os.chdir('path_to_your_dir')

# Start task
print('Start task on: '+str(datetime.datetime.now()) + '\n')

# Read data
control_data = pd.read_csv('./data/ireceptor_data/ireceptor-public-archive.tsv', sep='\t')

# Parse data
control_data = gf.parse_background(control_data)

# Calculate V gene frequencies
background_v_genes = gf.get_percentages(control_data, 'TRBV_gene', 'background')
background_v_genes.to_csv('./data/parsed_data/ireceptor_data/background_v_genes.tsv', sep='\t')

# Calculate J gene frequencies
background_j_genes = gf.get_percentages(control_data, 'TRBJ_gene', 'background')
background_j_genes.to_csv('./data/parsed_data/ireceptor_data/background_j_genes.tsv', sep='\t')

# End task
print('End task on: '+str(datetime.datetime.now()) + '\n')
