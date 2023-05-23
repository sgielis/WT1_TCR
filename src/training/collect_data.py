import datetime
import os
import pandas as pd

import general_functions as gf

# Adjust the working directory
os.chdir('path_to_your_dir')
data_dir = './data/raw'


print('Start data collection \n')
# First data set
all_data = pd.DataFrame(columns=['targetSequences', 'cloneId',
                                 'cloneCount', 'cloneFraction',
                                 'allVHitsWithScore', 'allJHitsWithScore',
                                 'minQualCDR3', 'aaSeqCDR3'])

for epitope in ['WT1-37', 'WT1-126', 'IE62']:
    data = gf.collect_mixcr_data(data_dir, 'run1', epitope)
    all_data = all_data.append(data)
    print('')

# ORF18 data
data = gf.collect_mixcr_data(data_dir, 'run1_orf18', 'ORF18')
all_data = all_data.append(data)
print('')

# Second data set
data = gf.collect_mixcr_data(data_dir, 'run2', 'WT137',
                             remove_list=['8_DR18_WT137_IL15a_S4_L001_clones.txt',
                                          '9_DR18_WT137_IL15b_S5_L001_clones.txt'])
all_data = all_data.append(data)
print('')

# Third data set
for epitope in ['WT1-37', 'WT1-126']:
    data = gf.collect_mixcr_data(data_dir, 'run3', epitope)
    all_data = all_data.append(data)
    print('')

print('Parse collected data')
# Use one notation for every epitope
all_data['epitope'] = all_data['epitope'].apply(lambda x: 'WT1-37' if 'WT137' in x else x)

all_data = gf.parse_mixcr_data(all_data, 1)

all_data.to_csv('./data/parsed_data/all_data.tsv',
                sep='\t',  index=False)

# End task
print('End task on: '+str(datetime.datetime.now()) + '\n')
