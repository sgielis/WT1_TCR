import os

import pandas as pd

# Adjust the working directory
os.chdir('path_to_your_dir')
resultsdir = './results/public_tcrs'

# Calculate publicity for WT1-126
wt1_126 = pd.read_csv(os.path.join(resultsdir, 'public_126.tsv'), sep='\t')
wt1_126 = wt1_126.groupby(['nr_volunteers'])[['CDR3_beta']].count()
wt1_126['epitope'] = 'WT1-126'

# Calculate publicity for WT1-37
wt1_37 = pd.read_csv(os.path.join(resultsdir, 'public_37.tsv'), sep='\t')
wt1_37 = wt1_37.groupby(['nr_volunteers'])[['CDR3_beta']].count()
wt1_37['epitope'] = 'WT1-37'

# Combine info for both epitopes
plot_info = pd.concat([wt1_126, wt1_37], axis=0).reset_index()
df = pd.DataFrame([['5', 0, 'WT1-126'], ['6', 0, 'WT1-126']],
                  columns=['nr_volunteers', 'CDR3_beta', 'epitope'])
plot_info = plot_info.append(df).sort_values(by=['epitope'])
plot_info.to_csv(os.path.join(resultsdir, 'public_counts.tsv'),
                 index=False, sep='\t')
