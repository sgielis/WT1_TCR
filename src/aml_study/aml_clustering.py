# Imports
import datetime
import numpy as np
import os
import pandas as pd
from clustcr import Clustering

# Adjust the working directory
os.chdir('path_to_your_dir')
datadir = './results/tcrex_tcrdb_aml/combined_PRJNA510967/results'
results_dir = './results/wt1_visualization/PRJNA510967'
if not os.path.exists(results_dir):
        os.makedirs(results_dir)

# Response info
info = pd.read_csv('./data/parsed_data/tcrdb_aml/info_PRJNA510967.tsv',
                   sep='\t')


# Assemble all data from study
def read_data(datadir, volunteer):
    all_results = pd.read_csv(os.path.join(datadir,
                                           volunteer, 'all_results.tsv'),
                              skiprows=[0, 1, 2, 3, 4, 5], sep='\t')
    all_results = all_results.drop_duplicates(['TRBV_gene',
                                               'CDR3_beta', 'TRBJ_gene'])
    return all_results


def get_response(info, volunteer):
    response = info[info['volunteer'] == volunteer]['Response']
    return response.tolist()[0]


# Get all files of interest
files = os.listdir(datadir)
files.remove('.DS_Store')
# Empty dataframe
df = pd.DataFrame(columns=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene'])
# Add all files to empty dataframe
for volunteer in files:
    # Read in TCRex parsed TCRs
    data = read_data(datadir, volunteer)
    data['Volunteer'] = volunteer
    data['Response'] = get_response(info, volunteer)
    df = df.append(data)
print('Total unique TCRs: ', df.drop_duplicates(subset=['CDR3_beta']).shape[0])
print('Total TCRs: ', df.shape[0])
df.to_csv(os.path.join(results_dir, 'All_tcrs.tsv'), sep='\t')

#  Cluster data
sequences = df['CDR3_beta']
clustering = Clustering(faiss_training_data=sequences,
                        fitting_data_size=len(sequences),
                        max_sequence_size=sequences.str.len().max(),
                        n_cpus=8)


categories = set(df['Response'].tolist())

for cat in categories:
    # cluster per category
    data = df[df['Response'] == cat]['CDR3_beta']
    # Cluster canonical CDR3s
    clustering.batch_precluster(data, name=cat)

for cluster in clustering.batch_cluster(calc_feature_matrix=True):
    results = cluster.clusters_df
    results.to_csv(os.path.join(results_dir, 'clusters.tsv'), sep='\t', index=False)
    clustered_tcrs = set(results['CDR3'].tolist())
    print('Nr of TCR clustered:', len(clustered_tcrs))

matrix = clustering.batch_feature_matrix()
clustering.batch_cleanup()

# Calculate number of clusters combining healthy and cancer patients
overlap = matrix.transpose()
overlap = overlap[overlap['Healthy'] >= 1]
overlap = overlap.replace(0, np.nan)
overlap = overlap.dropna(how='all',subset=['Complete_remission','Relapse'])
print(overlap.shape[0], ' Clusters combining healthy and cancer patients')

# Collect TCRex results
tcrex_dir = os.path.join('./results/tcrex_tcrdb_aml/combined_PRJNA510967/results')
files = os.listdir(tcrex_dir)
files.remove('.DS_Store')
tcrex_results = pd.DataFrame(columns=['TRBV_gene', 'CDR3_beta',
                                      'TRBJ_gene', 'epitope'])
for file in files:
    data = pd.read_csv(os.path.join(tcrex_dir, file, 'filtered_results.tsv'),
                       skiprows=[0, 1, 2, 3, 4, 5, 6], sep='\t')
    data = data[['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope']]
    tcrex_results = tcrex_results.append(data)
tcrex_126 = tcrex_results[tcrex_results['epitope'] == 'WT1-126_count1']
tcrex_126 = tcrex_126['CDR3_beta'].tolist()
tcrex_37 = tcrex_results[tcrex_results['epitope'] == 'WT1-37_count1']
tcrex_37 = tcrex_37['CDR3_beta'].tolist()

# Collect train data
data_37 = pd.read_csv('./data/parsed_data/tcrex_training_data/WT1-37_count1/20210811_target.tsv',
                      sep='\t')
train_37 = data_37['CDR3_beta'].tolist()
data_126 = pd.read_csv('./data/parsed_data/tcrex_training_data/WT1-126_count1/20210811_target.tsv',
                       sep='\t')
train_126 = data_126['CDR3_beta'].tolist()

# Combine TCRex and train CDRs
all_37 = set(train_37 + tcrex_37)
all_126 = set(train_126 + tcrex_126)

# Add WT1-info to clusters_df
results['WT1_37'] = results['CDR3'].apply(lambda x: 1 if x in all_37 else 0)
results['WT1_126'] = results['CDR3'].apply(lambda x: 1 if x in all_126 else 0)


def wt1_specific_clusters(results, epitope):
    clus = results[results[epitope]>0]['cluster'].tolist()
    df = results[results['cluster'].isin(clus)]
    df.to_csv('./results/wt1_visualization/PRJNA510967/'+epitope+'_specific_clusters', sep='\t', index=False)
    return df


# Get WT1-sepcific clusters
wt1_specific_clusters(results, 'WT1_37')
# Get WT1-sepcific clusters
wt1_specific_clusters(results, 'WT1_126')

# Group results per cluster
results = results.groupby(['cluster']).sum()
results['SUM'] = results.sum(axis=1)

# Add pathology info
category_counts = matrix.transpose()
final = pd.concat([category_counts, results], axis=1)

# Clusers with WT1 info
final = final[final['SUM'] > 0]
final.to_csv(os.path.join(results_dir, 'Response_matrix.tsv'), sep='\t')

# Save edges
edges = cluster.export_network(filename=os.path.join(results_dir, 'edges.txt'))

# End task
print('End task on: '+str(datetime.datetime.now()) + '\n')
