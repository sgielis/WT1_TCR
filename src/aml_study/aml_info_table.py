# imports
import datetime
import os

import pandas as pd

# Adjust the working directory
os.chdir('path_to_your_dir')
directory = './data/parsed_data/tcrdb_aml'

response = {'SRR8359186': 'Complete_remission',
            'SRR8359187': 'Complete_remission',
            'SRR8359192': 'Healthy',
            'SRR8359193': 'Healthy',
            'SRR8359194': 'Healthy',
            'SRR8359195': 'Complete_remission',
            'SRR8359196': 'Relapse',
            'SRR8359197': 'Relapse',
            'SRR8359198': 'Relapse',
            'SRR8359199': 'Relapse',
            'SRR8359200': 'Relapse',
            'SRR8359201': 'Relapse',
            'SRR8359202': 'Relapse',
            'SRR8359203': 'Relapse'}

tag = {'SRR8359186': 'PT2_CR2_s0_beta',
       'SRR8359187': 'PT3_CR3_s0_beta',
       'SRR8359192': 'HD1_s0_beta',
       'SRR8359193': 'HD2_s0_beta',
       'SRR8359194': 'HD3_s0_beta',
       'SRR8359195': 'PT1_CR1_s0_beta',
       'SRR8359196': 'PT6_REL3_s0_beta',
       'SRR8359197': 'PT5_REL2_s1_beta',
       'SRR8359198': 'PT6_REL3_s2_beta',
       'SRR8359199': 'PT6_REL3_s1_beta',
       'SRR8359200': 'PT4_REL1_s1_beta',
       'SRR8359201': 'PT4_REL1_s0_beta',
       'SRR8359202': 'PT5_REL2_s0_beta',
       'SRR8359203': 'PT4_REL1_s2_beta'}


# Combine response and tag info in pandas df
# Response to df
response_df = pd.DataFrame.from_dict(response, orient='index')
response_df = response_df.rename(columns={0: 'Response'})
# Tag to df
tag_df = pd.DataFrame.from_dict(tag, orient='index')
tag_df = tag_df.rename(columns={0: 'Tag'})
# Combine
info_df = pd.concat([response_df, tag_df], axis=1).reset_index()
# Parse volunteer
info_df['volunteer'] = info_df['Tag'].apply(lambda x: '_'.join(x.split('_')[:-2]))
info_df = info_df.rename(columns={'index': 'File'})
# Save info
info_df.to_csv(os.path.join(directory, 'info_PRJNA510967.tsv'),
               sep='\t', index=False)

# End task
print('End task on: '+str(datetime.datetime.now()) + '\n')
