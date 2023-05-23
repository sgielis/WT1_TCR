# Imports
import datetime
import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import string

import general_functions as gf

sns.set_style("white")

# Adjust the working directory
os.chdir('path_to_your_dir')

# Read control data
control_data = pd.read_csv('./data/ireceptor_data/ireceptor-public-archive.tsv', sep='\t')

# Parse control data
control_data = gf.parse_background(control_data)
control_data['epitope'] = 'background'

# Read in full dataset
data = pd.read_csv('./data/parsed_data/full_filtered_dataset.tsv')
data = data[['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope']].drop_duplicates().reset_index(drop=True)

# Parse  V/J genes
data['TRBV_gene'] = ('TRBV'+data['TRBV_gene'].apply(
        lambda gene_name: '-'.join([sub_gene.zfill(2) for sub_gene in
                                    gene_name.replace('TRBV', '').split('-')])))
data['TRBJ_gene'] = ('TRBJ' + data['TRBJ_gene'].apply(
        lambda gene_name: '-'.join([sub_gene.zfill(2) for sub_gene in
                                    gene_name.replace('TRBJ', '').split('-')])))

# Concat all data
data = pd.concat([data, control_data], axis=0)

# Get families
data['TRBV_family'] = data['TRBV_gene'].str.split('-').str[0].str.zfill(2)
data['TRBJ_family'] = data['TRBJ_gene'].str.split('-').str[0].str.zfill(2)

# Plot frequencies
f, axes = plt.subplots(2, 2, figsize=(20, 30))

for i, gene in enumerate(['TRBV_gene', 'TRBV_family', 'TRBJ_gene', 'TRBJ_family']):

    gene_freqs = []
    for epitope in ['WT1-37', 'WT1-126', 'background']:
        epitope_counts = data[(data['epitope'] == epitope)][gene].value_counts(normalize=True, sort=False)
        gene_freqs.append(pd.DataFrame({'gene': epitope_counts.index, 'frequency': epitope_counts.values*100,
                                        'epitope': [epitope] * len(epitope_counts)}))
    gene_freqs = pd.concat(gene_freqs).sort_values(by='gene')

    ax = axes[i // 2, i % 2]

    sns.barplot(x='frequency', y='gene', hue='epitope', data=gene_freqs,
                palette='colorblind', ax=ax, hue_order=['WT1-37', 'WT1-126', 'background'])
    ax.set_title('{} frequency'.format(gene.replace('_', '-')), fontsize=20)
    ax.set_ylabel('{}'.format(gene.replace('_', '-')), fontsize=20)
    ax.set_xlabel('frequency (%)', fontsize=20)
    ax.set_yticklabels(gene_freqs['gene'].unique(), rotation='horizontal', size=12)
    ax.legend(loc='upper right', fontsize=15)
    ax.annotate(string.ascii_uppercase[i], xy=(-0.1, 1.1), xycoords='axes fraction', textcoords='offset points',
                fontsize=20, xytext=(0, -15), weight='bold', ha='right', va='top')

plt.savefig('./figures/vj_frequency/vj_frequency_with_background.pdf',
            bbox_inches='tight', dpi=600)
plt.close()

# End task
print('End task on: '+str(datetime.datetime.now()) + '\n')
