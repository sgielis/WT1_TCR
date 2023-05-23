# Imports
import datetime
import os
import pandas as pd

# Adjust the working directory
os.chdir('path_to_your_dir')

# Read in data sets
# All TCRs (no dropping of duplicates)
all_data = pd.read_csv('./data/parsed_data/all_data.tsv',
                       sep='\t')

# TCRs linked with one epitope
keep_tcrs = pd.read_csv('./data/parsed_data/keep_info.tsv',
                        sep='\t')

# Shared TCRs that must be removed
remove_tcrs = pd.read_csv('./data/parsed_data/remove_tcrs.tsv',
                          sep='\t')


# Extract the volunteer names
all_data['volunteer'] = all_data['file_name'].str.extract(r'[0-9]*[-_]*(DR[0-9]+)[a-z-A-Z-0-9_-]+')
all_data['ID'] = all_data['file_name'].str.replace('_L001_clones.txt', '')


# Remove TCRs linked with wrong epitopes
def action_tcrs(reference, cdr3, epitope):

    if cdr3 in reference['CDR3_beta'].tolist():
        ref_data = reference[reference['CDR3_beta'] == cdr3]
        ref_epitope = ref_data['epitope'].tolist()[0]
        if epitope == ref_epitope:
            action = 'keep'
        else:
            action = 'remove'
    else:
        action = 'keep'

    return action


all_data['action'] = all_data.apply(lambda x: action_tcrs(keep_tcrs, x['CDR3_beta'], x['epitope']), axis=1)
all_data = all_data[all_data['action'] == 'keep']
all_data = all_data.drop(['action'], axis=1)

# Remove shared TCRs
remove_tcrs = remove_tcrs['CDR3_beta'].tolist()


def drop_shared_tcrs(remove_tcrs, cdr3):

    if cdr3 in remove_tcrs:
        action = 'remove'
    else:
        action = 'keep'

    return action


all_data['action'] = all_data.apply(lambda x: drop_shared_tcrs(remove_tcrs, x['CDR3_beta']), axis=1)
all_data = all_data[all_data['action'] == 'keep']

all_data.to_csv('./data/parsed_data/full_filtered_dataset.tsv', index=False)


# End task
print('End task on: '+str(datetime.datetime.now()) + '\n')
