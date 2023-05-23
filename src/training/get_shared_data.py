# imports
import datetime
import os
import pandas as pd


# Adjust the working directory
os.chdir('path_to_your_dir')
data = './data/parsed_data/all_data.tsv'

# Identify TCRs with different epitopes
all_tcrs = pd.read_csv(data, sep='\t')
all_tcrs['count'] = all_tcrs['cloneCount'].astype(str)
all_tcrs['sorted_epitope'] = all_tcrs['epitope']

group_V = (all_tcrs.groupby(['CDR3_beta'])['TRBV_gene'].apply(', '.join)
           .reset_index().set_index(['CDR3_beta']))
group_J = (all_tcrs.groupby(['CDR3_beta'])['TRBJ_gene'].apply(', '.join)
           .reset_index().set_index(['CDR3_beta']))
group_sorted_epitope = (all_tcrs.groupby(['CDR3_beta'])['sorted_epitope'].apply(', '.join)
                        .reset_index().set_index(['CDR3_beta']))
group_epitope = (all_tcrs.groupby(['CDR3_beta'])['epitope'].apply(', '.join)
                 .reset_index().set_index(['CDR3_beta']))
group_count = (all_tcrs.groupby(['CDR3_beta'])['count'].apply(', '.join)
               .reset_index().set_index(['CDR3_beta']))
grouped_data = pd.concat([group_epitope,
                          group_sorted_epitope,
                          group_V,
                          group_J, group_count], axis=1)
grouped_data['epitope'] = grouped_data['epitope'].apply(
    lambda x: x if ',' not in str(x) else ','.join(set(y.strip()
                                                   for y in x.split(','))))
shared_data = grouped_data[grouped_data['epitope'].str.contains(',')]
shared_data = shared_data.reset_index()
shared_data.to_csv('./data/parsed_data/all_shared_data.tsv',
                   index=False, sep='\t')


# Calculate the fraction of the two highest clone counts
def get_counts(cdr3, shared_info):

    # Get all entries with shared CDR3
    cdr3_data = shared_info[shared_info['CDR3_beta'] == cdr3]
    # Order the entries descending on clone count
    cdr3_data = cdr3_data.sort_values(by=['cloneCount'], ascending=False)
    # Only retain the entry with the highest clone count for every epitope
    cdr3_data = cdr3_data.drop_duplicates(subset=['epitope'], keep='first')
    # Retrieve all remaining count values for the CDR3
    count_list = cdr3_data[cdr3_data['CDR3_beta'] == cdr3]['cloneCount'].sort_values(ascending=False).tolist()
    # Calculate the fraction between the first and second count
    fraction = count_list[0]/count_list[1]

    # Assign the epitope for which the fraction is at least 100, otherwise return NAN
    if fraction >= 100:
        count = str(count_list[0])
    else:
        count = 'NAN'

    return count


# Get all info on all shared tcrs
shared_tcrs = shared_data['CDR3_beta'].tolist()
shared_info = all_tcrs[all_tcrs['CDR3_beta'].isin(shared_tcrs)]
shared_info['selected_epitope'] = shared_info['CDR3_beta'].apply(lambda x: get_counts(x, shared_info))

# Identify TCRs that may not be removed
selected_epitopes = (shared_info.groupby(['CDR3_beta'])['selected_epitope'].apply(', '.join)
                     .reset_index().set_index(['CDR3_beta']))
selected_epitopes['count'] = selected_epitopes['selected_epitope'].apply(
    lambda x: x if ',' not in str(x) else ','.join(set(y.strip()
                                                   for y in x.split(','))))
selected_epitopes = selected_epitopes.reset_index()

# Save all info of TCRs that need to be retained
keep_tcrs = selected_epitopes[selected_epitopes['count'] != 'NAN']
keep_tcrs.to_csv('./data/parsed_data/keep_tcrs.tsv',
                 index=False, sep='\t')
keep_info = keep_tcrs.merge(all_tcrs, how='inner', on=['CDR3_beta', 'count'])
keep_info = keep_info[['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'cloneCount', 'epitope', 'file_name', 'folder', 'count']]
keep_info.to_csv('./data/parsed_data/keep_info.tsv',
                 index=False, sep='\t')

remove_tcrs = selected_epitopes[selected_epitopes['count'] == 'NAN']
remove_tcrs.to_csv('./data/parsed_data/remove_tcrs.tsv',
                   index=False, sep='\t')

# End task
print('End task on: '+str(datetime.datetime.now()) + '\n')
