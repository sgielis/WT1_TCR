# imports
import datetime
import os

import pandas as pd

# Adjust the working directory
os.chdir('path_to_your_dir')
directory = './data/parsed_data/tcrdb_aml'


def combine_data(data_directory, files):
    # Empty dataframe
    df = pd.DataFrame(columns=['TRBV_gene', 'CDR3_beta',
                               'TRBJ_gene'])
    # Add all files to empty dataframe
    for f in files:
        print('Add file ', f)
        data = pd.read_csv(os.path.join(data_directory, f+'.tsv'),
                           sep='\t')
        df = df.append(data)
        df = df.drop_duplicates()
    return df


# Get all volunteer names
info_df = pd.read_csv(os.path.join(directory, 'info_PRJNA510967.tsv'), sep='\t')
volunteers = set(info_df['volunteer'].tolist())

# Append files for every volunteer
for v in volunteers:
    print('Combining files for volunteer: ', v)
    volunteer_df = info_df[info_df['volunteer'] == v]
    files = volunteer_df['File'].tolist()
    all_data = combine_data(os.path.join(directory, 'PRJNA510967'), files)
    print('Total unique TCRs: ', all_data.shape[0], '\n')
    all_data.to_csv(os.path.join(directory, 'combined_PRJNA510967', v + '.tsv'),
                    sep='\t', index=False)

# End task
print('End task on: '+str(datetime.datetime.now()) + '\n')
