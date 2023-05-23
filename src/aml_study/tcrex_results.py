# imports
import pandas as pd
import os

# Adjust the working directory
os.chdir('path_to_your_dir')


# Training data
data_37 = './data/parsed_data/tcrex_training_data/WT1-37_count1/20210811_target.tsv'
data_37 = pd.read_csv(data_37, sep='\t')

data_126 = './data/parsed_data/tcrex_training_data/WT1-126_count1/20210811_target.tsv'
data_126 = pd.read_csv(data_126, sep='\t')


def calculate_percentages(directory, all_folders):
    # empty dictionaries
    results = {}
    specific_tcrs = pd.DataFrame()

    for volunteer in all_folders:
        results[volunteer] = {}

        # Get the TCR repertoire of the volunteer
        # Use the 'all_results' folder as TCR repertoire as this contains all TCRs parsed by TCRex
        all_results = pd.read_csv(os.path.join(directory,
                                               volunteer, 'all_results.tsv'),
                                  skiprows=[0, 1, 2, 3, 4, 5], sep='\t')
        
        # Repertoire size
        # Drop results since every TCR is present twice (predictions for 2 epitopes were made)
        size = all_results.drop_duplicates(['TRBV_gene', 'CDR3_beta',
                                            'TRBJ_gene']).shape[0]
        results[volunteer]['TCR_repertoire'] = size

        # Get the WT1-specific TCRs
        filtered_results = pd.read_csv(os.path.join(directory,
                                                    volunteer,
                                                    'filtered_results.tsv'),
                                       skiprows=[0, 1, 2, 3, 4, 5, 6], sep='\t')
        # Collect predicted WT1-37 TCRs
        aml_37 = filtered_results[filtered_results['epitope'] == 'WT1-37_count1']
        aml_37['Volunteer'] = volunteer
        specific_tcrs = pd.concat([specific_tcrs,aml_37], axis=0)
        results[volunteer]['WT1-37 TCRs'] = aml_37.shape[0]
        
        # Collect predicted WT1-126 TCRs
        aml_126 = filtered_results[filtered_results['epitope'] == 'WT1-126_count1']
        aml_126['Volunteer'] = volunteer
        specific_tcrs = pd.concat([specific_tcrs,aml_126], axis=0)
        results[volunteer]['WT1-126 TCRs'] = aml_126.shape[0]

        # Percentage
        results[volunteer]['% WT1-37'] = (aml_37.shape[0]/size)*100
        results[volunteer]['% WT1-126'] = (aml_126.shape[0]/size)*100

        # Overlap tcrex and training data
        overlap_37 = aml_37.merge(data_37, how='inner', on=['CDR3_beta']).shape[0]
        overlap_126 = aml_126.merge(data_126, how='inner', on=['CDR3_beta']).shape[0]

        results[volunteer]['# public WT1-37'] = overlap_37
        results[volunteer]['# public WT1-126'] = overlap_126

    return results, specific_tcrs


# Get list of all files in directory.
directory = './results/tcrex_tcrdb_aml/combined_PRJNA510967/results'
all_folders = os.listdir(directory)
all_folders.remove('.DS_Store')

# Calculate identification percentages
results, specific_tcrs = calculate_percentages(directory, all_folders)
df = pd.DataFrame(results)
df = df.transpose()
df = df[['TCR_repertoire', 'WT1-37 TCRs', 'WT1-126 TCRs', '% WT1-37', '% WT1-126', '# public WT1-37', '# public WT1-126']]
df = df.reset_index()

# Store TCRex results
df.to_csv('./results/tcrex_tcrdb_aml/combined_PRJNA510967/tcrex_results.tsv')
specific_tcrs.to_csv('./results/tcrex_tcrdb_aml/combined_PRJNA510967/tcrex_tcrs.tsv')