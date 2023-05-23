import pandas as pd
import os
import re

from pyteomics import parser


def get_all_folder_names(data_directory):
    # get all zip folders
    all_folders = os.listdir(data_directory)
    r = re.compile(".+zip")
    zip_folders = [x for x in all_folders if r.match(x)]
    # get name of ziped files
    all_folder_names = [x.split('.')[0] for x in zip_folders]
    print('Nr of folders:', len(all_folder_names))
    return all_folder_names


def get_specific_data(all_folder_names, epitope):
    r = re.compile(".+"+epitope+".+")
    specific_folders = [x for x in all_folder_names if r.match(x)]
    return specific_folders


def collect_mixcr_data(data_dir, folder, epitope, remove_list=None):

    print('Collect data for ' + epitope + ' from data set ' + folder)

    files = get_specific_data(os.listdir(os.path.join(data_dir, folder)),
                              epitope)
    if remove_list is not None:
        files = [x for x in files if x not in remove_list]
    print('Number of files: ', len(files))

    # Empty dataframe
    df = pd.DataFrame(columns=['targetSequences', 'cloneId',
                               'cloneCount', 'cloneFraction',
                               'allVHitsWithScore', 'allJHitsWithScore',
                               'minQualCDR3', 'aaSeqCDR3'])
    # Add all files to empty dataframe
    for f in files:
        print('Add data from ', f)
        data = pd.read_csv(os.path.join(data_dir, folder, f), sep='\t',
                           usecols=['targetSequences', 'cloneId',
                                    'cloneCount', 'cloneFraction',
                                    'allVHitsWithScore', 'allJHitsWithScore',
                                    'minQualCDR3', 'aaSeqCDR3'])
        data['epitope'] = epitope
        data['file_name'] = f
        data['folder'] = folder
        df = df.append(data)
    return df


# Remove sequences that are shared between different epitopes
def remove_shared_sequences(sequences1, sequences2, inplace=False):
    if inplace:
        sequences2.drop(sequences1.index, inplace=True, errors='ignore')
    else:
        return sequences2.drop(sequences1.index, errors='ignore')


# Remove sequences with non-amino acids
def _is_amino_acid_sequence(peptide):
    return all(aa in parser.std_amino_acids for aa in peptide)


def _remove_orphon_genes(sequences, V_gene):
    return sequences[~sequences[V_gene].str.contains('or', case=False)]


def _parse_genes(data):
    data['TRBV_gene'] = ('TRBV'+data['TRBV_gene'].apply(
            lambda gene_name: '-'.join([sub_gene.zfill(2) for sub_gene in
                                        gene_name.replace('TRBV', '').split('-')])))
    data['TRBJ_gene'] = ('TRBJ' + data['TRBJ_gene'].apply(
            lambda gene_name: '-'.join([sub_gene.zfill(2) for sub_gene in
                                        gene_name.replace('TRBJ', '').split('-')])))
    return data


# Parse the collected data
def parse_mixcr_data(sequences, clonecount):

    print('Number of starting sequences:', sequences.shape[0], '\n')

    # Retain TRBV chain
    print('Retain TRBV sequences')
    sequences['Gene'] = sequences['allVHitsWithScore'].apply(lambda x: 'beta' if 'TRBV' in x else 'non-beta')
    sequences = sequences[sequences['Gene'] == 'beta']
    print('Number of sequences retained:', sequences.shape[0], '\n')

    # Rename CDR3 beta
    sequences = sequences.rename(columns={'aaSeqCDR3': 'CDR3_beta'})

    print('Remove non-canonical CDRs')
    # Remove non-canonical TCRs: containing non-canonical amino acids.
    sequences = sequences[sequences['CDR3_beta'].apply(_is_amino_acid_sequence)]
    # Remove non-canonical TCRs: not starting with C or ending with F
    start_c = sequences['CDR3_beta'].str.startswith('C', na=False)
    end_f = sequences['CDR3_beta'].str.endswith('F', na=False)
    sequences = sequences[start_c & end_f].reset_index(drop=True)
    print('Number of sequences retained:', sequences.shape[0], '\n')

    # Parse V/J genes
    sequences = _remove_orphon_genes(sequences, 'allVHitsWithScore')
    sequences['TRBV_gene'] = sequences['allVHitsWithScore'].str.replace(r'(TRBV[ABC0123456789\-]+)\*.+', lambda m: m.group(1), regex=True)
    sequences['TRBJ_gene'] = sequences['allJHitsWithScore'].str.replace(r'(TRBJ[0123456789P\-]+)\*.+', lambda m: m.group(1), regex=True)

    # Filter on clone count
    print('Filter on clone count')
    sequences = sequences[sequences['cloneCount'] >= clonecount]
    print('Number of sequences retained:', sequences.shape[0], '\n')

    # Do not drop duplicates before clone count filtering
    sequences = sequences[['targetSequences', 'TRBV_gene', 'CDR3_beta', 'TRBJ_gene',
                           'cloneCount', 'epitope', 'file_name', 'folder']]
    print('Number of sequences retained:', sequences.shape[0], '\n')

    return sequences


def parse_background(control_data):
    # Read background data
    control_data = control_data[['productive', 'v_call', 'j_call', 'junction_aa']]
    control_data = control_data[control_data['productive'] == 'T']
    control_data = control_data.drop_duplicates()

    # Remove non-canonical beta sequences
    start_c = control_data['junction_aa'].str.startswith('C', na=False)
    end_f = control_data['junction_aa'].str.endswith('F', na=False)
    control_data = control_data[start_c & end_f].reset_index(drop=True)
    control_data = control_data[control_data['junction_aa'].apply(_is_amino_acid_sequence)]

    # Remove allele info
    control_data['TRBV_gene'] = control_data['v_call'].str.split('*').str[0]
    control_data['TRBJ_gene'] = control_data['j_call'].str.split('*').str[0]

    # Remove orphon genes and rename columns
    control_data = control_data[['junction_aa', 'TRBV_gene', 'TRBJ_gene']]
    control_data = _remove_orphon_genes(control_data, 'TRBV_gene')
    control_data = control_data.rename(columns={'junction_aa': 'CDR3_beta'})

    # Parse V/J genes
    control_data = _parse_genes(control_data)

    # Drop duplicates
    control_data = control_data.drop_duplicates()
    return control_data


def get_percentages(data, gene, epitope):
    counts = data[gene].value_counts(sort=False).to_frame(name=epitope + '_counts')
    percentages = data[gene].value_counts(normalize=True, sort=False)*100
    percentages = percentages.to_frame(name=epitope + '_percentage')
    return pd.concat([counts, percentages], axis=1).sort_index()
