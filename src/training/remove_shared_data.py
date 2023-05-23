# Imports
import datetime
import os
import pandas as pd

import general_functions as gf


# Adjust the working directory
os.chdir('path_to_your_dir')

# Read data sets
all_data = pd.read_csv('./data/parsed_data/all_data.tsv',
                       sep='\t')
shared_data = pd.read_csv('./data/parsed_data/all_shared_data.tsv',
                          sep='\t')

keep_tcrs = pd.read_csv('./data/parsed_data/keep_info.tsv',
                        sep='\t')
keep_tcrs = keep_tcrs.drop(columns=['count'])

# Remove shared sequences from positive data
print('Remove shared sequences \n')
all_data = all_data.set_index(['CDR3_beta'])
shared_data = shared_data.set_index(['CDR3_beta'])
filtered_sequences = gf.remove_shared_sequences(shared_data, all_data,
                                                inplace=False)
filtered_sequences = filtered_sequences.reset_index()

# Add the keep tcr data
final_data = pd.concat([filtered_sequences, keep_tcrs])
print('Number of sequences retained:', final_data.shape[0])

# Save final data set
final_data.to_csv('./data/parsed_data/all_data_without_shared_tcrs.tsv',
                  sep='\t', index=False)
# End task
print('End task on: '+str(datetime.datetime.now()) + '\n')
