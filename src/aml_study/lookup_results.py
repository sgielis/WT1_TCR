# imports
import os
import pandas as pd

# Adjust the working directory
os.chdir('path_to_your_dir')


# training data
data_37 = './data/parsed_data/tcrex_training_data/WT1-37_count1/20210811_target.tsv'
data_37 = pd.read_csv(data_37, sep='\t')

data_126 = './data/parsed_data/tcrex_training_data/WT1-126_count1/20210811_target.tsv'
data_126 = pd.read_csv(data_126, sep='\t')


def look_up(directory, all_folders):
    # empty dictionary/df
    results = {}
    lookups = pd.DataFrame()

    for volunteer in all_folders:
        results[volunteer] = {}

        # Get the TCR repertoire of the volunteer
        # Use the 'all_results' folder as TCR repertoire as this contains all TCRs parsed by TCRex
        all_results = pd.read_csv(os.path.join(directory,
                                               volunteer, 'all_results.tsv'),
                                  skiprows=[0, 1, 2, 3, 4, 5], sep='\t')
        
        # Repertoire size
        # Drop results since every TCR is present twice (predictions for 2 epitopes were made)
        all_results = all_results.drop_duplicates(['TRBV_gene', 'CDR3_beta', 'TRBJ_gene'])
        size = all_results.shape[0]
        
        # Match TCR repertoire with training data to identify WT1-37 TCRs
        tcrs_37 = all_results.merge(data_37, how='inner', on=['CDR3_beta']).drop_duplicates(['CDR3_beta'])
        # Add results to the lookup df
        tcrs_37['Volunteer'] = volunteer
        lookups = pd.concat([lookups,tcrs_37], axis=0)
        # Number of WT1-37 TCRs identified using look-up method
        tcrs_37 = tcrs_37.shape[0]
        
        # Match TCR repertoire with training data to identify WT1-126 TCRs
        tcrs_126 = all_results.merge(data_126, how='inner', on=['CDR3_beta']).drop_duplicates(['CDR3_beta'])
        # Add results to the lookup df
        tcrs_126['Volunteer'] = volunteer
        lookups = pd.concat([lookups,tcrs_126], axis=0)
        # Number of WT1-37 TCRs identified using look-up method
        tcrs_126 = tcrs_126.shape[0]
        
        # Store all results in a dictionary
        results[volunteer]['TCR_repertoire'] = size
        results[volunteer]['WT1-37 TCRs'] = tcrs_37
        results[volunteer]['WT1-126 TCRs'] = tcrs_126
        results[volunteer]['% WT1-37'] = (tcrs_37/size)*100
        results[volunteer]['% WT1-126'] = (tcrs_126/size)*100

    return results,lookups


# Get list of all files in directory.
directory = './results/tcrex_tcrdb_aml/combined_PRJNA510967/results'
all_folders = os.listdir(directory)
all_folders.remove('.DS_Store')

# Perform look-up method
results,lookups = look_up(directory, all_folders)
df = pd.DataFrame(results)
df = df.transpose()
df = df[['TCR_repertoire', 'WT1-37 TCRs', 'WT1-126 TCRs', '% WT1-37', '% WT1-126']]
df = df.reset_index()

# Store look-up results
df.to_csv('./results/tcrex_tcrdb_aml/combined_PRJNA510967/summary_lookup_results.tsv')
lookups.to_csv('./results/tcrex_tcrdb_aml/combined_PRJNA510967/lookup_tcrs.tsv')
