# imports
import datetime
import os
import pandas as pd

# Adjust the working directory
os.chdir('path_to_your_dir')
resultsdir = './results/public_tcrs'

# Read in full dataset
full_data = pd.read_csv('./data/parsed_data/full_filtered_dataset.tsv')

# Keep all CDR3 beta duplicates
data = full_data[full_data.duplicated(subset=['CDR3_beta'], keep=False)].sort_values('CDR3_beta', ascending=True)
data = data.drop(columns=['file_name', 'folder', 'ID', 'action'])
data = data.drop_duplicates()
data = data.rename(columns={'targetSequences': 'RNA'})
data['volunteer_count'] = data['volunteer'] + '(' + data['cloneCount'].astype(str) + ')'


# Concat volunteers for every RNA sequence
def concat_volunteers(data, rna):
    rows = data[data['RNA'] == rna]
    vols = str(rows['volunteer_count'].tolist())
    return vols.replace(',', ':')


data['all_volunteers'] = data['RNA'].apply(lambda x: concat_volunteers(data, x))
data['RNA_info'] = data['RNA'] + data['all_volunteers']


# Get public CDR3
def filter_public_tcrs(dataset):
    # assemble a list of volunteers for each CDR3
    dataset_id = (dataset.groupby(['CDR3_beta'])['volunteer'].apply(', '.join)
                  .reset_index().set_index(['CDR3_beta']))
    dataset_rna = (dataset.groupby(['CDR3_beta'])['RNA'].apply(', '.join)
                   .reset_index().set_index(['CDR3_beta']))
    dataset_info = (dataset.groupby(['CDR3_beta'])['RNA_info'].apply(', '.join)
                    .reset_index().set_index(['CDR3_beta']))
    dataset_V = (dataset.groupby(['CDR3_beta'])['TRBV_gene'].apply(', '.join)
                 .reset_index().set_index(['CDR3_beta']))
    dataset_J = (dataset.groupby(['CDR3_beta'])['TRBJ_gene'].apply(', '.join)
                 .reset_index().set_index(['CDR3_beta']))
    grouped_data = pd.concat([dataset_id, dataset_rna,
                              dataset_info, dataset_V,
                              dataset_J], axis=1)
    # avoid duplication of the same volunteer, e.g vol1,vol1 is replaced by vol1
    grouped_data['volunteer'] = grouped_data['volunteer'].apply(
                              lambda x: x if ',' not in str(x)
                              else ','.join(set(y.strip() for y in x.split(','))))
    grouped_data['RNA'] = grouped_data['RNA'].apply(
                              lambda x: x if ',' not in str(x)
                              else ','.join(set(y.strip() for y in x.split(','))))
    grouped_data['RNA_info'] = grouped_data['RNA_info'].apply(
                              lambda x: x if ',' not in str(x)
                              else ','.join(set(y.strip() for y in x.split(','))))
    grouped_data['TRBV_gene'] = grouped_data['TRBV_gene'].apply(
                              lambda x: x if ',' not in str(x)
                              else ','.join(set(y.strip() for y in x.split(','))))
    grouped_data['TRBJ_gene'] = grouped_data['TRBJ_gene'].apply(
                              lambda x: x if ',' not in str(x)
                              else ','.join(set(y.strip() for y in x.split(','))))
    # retain CDR3s occuring in more than one volunteer
    grouped_data = grouped_data[grouped_data['volunteer'].str.contains(',')]
    grouped_data['nr_volunteers'] = grouped_data.volunteer.str.count(',') + 1
    return grouped_data.sort_values(by='nr_volunteers', ascending=False)


# Collect TCRs for WT-37
collected_37 = data[data['epitope'].isin(['WT1-37'])]

# Public TCRs for WT1-37
public_37 = filter_public_tcrs(collected_37)
public_37 = public_37.reset_index()
public_37['Volunteer_count'] = public_37['RNA_info'].str.replace(r'A|T|G|C', '',
                                                                 regex=True)
public_37 = public_37[['CDR3_beta', 'TRBV_gene', 'TRBJ_gene',
                       'nr_volunteers', 'RNA_info', 'Volunteer_count']]
public_37.to_csv(os.path.join(resultsdir, 'public_37.tsv'),
                 index=False, sep='\t')

size_37 = full_data[full_data['epitope'].isin(['WT1-37'])].drop_duplicates(subset=['CDR3_beta']).shape[0]
print(public_37.shape[0], ' from ', size_37, ' is public.')
print(public_37.shape[0]/size_37, 'percent of WT1-37 is public.')

# Collect TCRs for WT1-126
collected_126 = data[data['epitope'].isin(['WT1-126'])]

# Public TCRs for WT1-126
public_126 = filter_public_tcrs(collected_126)
public_126 = public_126.reset_index()

public_126['Volunteer_count'] = public_126['RNA_info'].str.replace(r'A|T|G|C', '',
                                                                   regex=True)
public_126 = public_126[['CDR3_beta', 'TRBV_gene', 'TRBJ_gene',
                         'nr_volunteers', 'RNA_info', 'Volunteer_count']]
public_126.to_csv(os.path.join(resultsdir, 'public_126.tsv'),
                  index=False, sep='\t')

size_126 = full_data[full_data['epitope'].isin(['WT1-126'])].drop_duplicates(subset=['CDR3_beta']).shape[0]
print(public_126.shape[0], ' from ', size_126, ' is public.')
print(public_126.shape[0]/size_126, 'percent of WT1-126 is public.')

# End task
print('End task on: '+str(datetime.datetime.now()) + '\n')
