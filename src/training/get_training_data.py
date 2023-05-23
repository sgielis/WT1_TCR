import datetime
import os
import pandas as pd


# Adjust the working directory
os.chdir('path_to_your_dir')
result_dir = './data/parsed_data/tcrex_training_data/'

data = pd.read_csv('./data/parsed_data/all_data_without_shared_tcrs.tsv',
                   sep='\t')

all_epitopes = set(data['epitope'].tolist())

for epitope in all_epitopes:
    print('Collect data for epitope ', epitope)
    epitope_data = data[data['epitope'] == epitope]
    epitope_data = epitope_data.sort_values(by=['cloneCount'], ascending=False)

    for clone_count in [1, 5, 10, 25, 50, 100]:
        # Filter on clone count
        cl_data = epitope_data[epitope_data['cloneCount'] >= clone_count]
        cl_data = cl_data[['TRBV_gene', 'CDR3_beta', 'TRBJ_gene']]
        cl_data = cl_data.drop_duplicates(subset=['CDR3_beta'], keep='first')

        # Save filtered data
        folder_name = epitope + '_count' + str(clone_count)
        print(folder_name, 'contains', cl_data.shape[0], 'different TCRs')
        epitope_dir = os.path.join(result_dir, folder_name)
        if not os.path.exists(epitope_dir):
            os.makedirs(epitope_dir)
        cl_data.to_csv(os.path.join(epitope_dir, '20210811_target.tsv'),
                       index=False, sep='\t')

# End task
print('End task on: '+str(datetime.datetime.now()) + '\n')
