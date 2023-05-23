# Imports
import datetime
import os

import math
import numpy as np
import pandas as pd
from scipy import stats as statistics

import general_functions as gf

# Adjust the working directory
os.chdir('path_to_your_dir')


def parse_genes(data):
    data['TRBV_gene'] = ('TRBV' + data['TRBV_gene'].apply(
            lambda gene_name: '-'.join([sub_gene.zfill(2) for sub_gene in
                                        gene_name.replace('TRBV', '').split('-')])))
    data['TRBJ_gene'] = ('TRBJ' + data['TRBJ_gene'].apply(
            lambda gene_name: '-'.join([sub_gene.zfill(2) for sub_gene in
                                        gene_name.replace('TRBJ', '').split('-')])))
    return data


def fill_in_zeros(data):
    for column in list(data):
        data[column] = data[column].apply(lambda x: 0 if math.isnan(x) else x)
    return data


def assemble_data(percentages, background):
    data = pd.concat([percentages, background], axis=1)
    data = fill_in_zeros(data)
    return data


def calculate_p(epitope_count, background_count, epitope_repertoire, background_repertoire):
    if epitope_count >= 2:
        oddsratio, p_value = statistics.fisher_exact([[epitope_count, background_count],
                                                     [epitope_repertoire, background_repertoire]],
                                                     alternative='greater')
    else:
        p_value = np.NaN
    return p_value


def gene_enrichment(data, epitope_repertoire, background_repertoire, epitope):
    counts = epitope + '_counts'
    data['p_value'] = data.apply(lambda row: calculate_p(row[counts], row['background_counts'],
                                                         epitope_repertoire, background_repertoire), axis=1)
    data[epitope + '_repertoire'] = epitope_repertoire
    data['background_repertoire'] = background_repertoire
    return data.sort_values(by=['p_value'])


# Read background dataset
background_v_genes = pd.read_csv('./data/parsed_data/ireceptor_data/background_v_genes.tsv',
                                 sep='\t', index_col=0)
background_j_genes = pd.read_csv('./data/parsed_data/ireceptor_data/background_j_genes.tsv',
                                 sep='\t', index_col=0)
background_repertoire = background_v_genes['background_counts'].sum()

# Read in WT1 data
data = pd.read_csv('./data/parsed_data/full_filtered_dataset.tsv')
data = data[['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope']].drop_duplicates().reset_index(drop=True)
data = parse_genes(data)

# Gene enrichment
for epitope in ['WT1-37', 'WT1-126', 'ORF18', 'IE62']:
    # Get epitope data and size
    epitope_data = data[data['epitope'] == epitope]
    epitope_repertoire = epitope_data.shape[0]
    for gene, background in zip(['TRBV_gene', 'TRBJ_gene'], [background_v_genes, background_j_genes]):
        # Get counts and percentages for epitope and background
        epitope_percentages = gf.get_percentages(epitope_data, gene, epitope)
        all_percentages = assemble_data(epitope_percentages, background)
        # Perform enrichment analysis
        result = (gene_enrichment(all_percentages, epitope_repertoire, background_repertoire, epitope)
                  .rename_axis(gene).reset_index())
        # Save results
        filename = epitope + '_enriched_' + gene + '.tsv'
        result.to_csv(os.path.join('./results/vj_genes/',
                                   filename), na_rep='NA', index=False)

# End task
print('End task on: '+str(datetime.datetime.now()) + '\n')
