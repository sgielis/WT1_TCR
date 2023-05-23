import datetime
import os
import pandas as pd


# Adjust the working directory
os.chdir('path_to_your_dir')
data_dir = './data/original_data/tcrdb_aml/'
results_dir = './data/parsed_data/tcrdb_aml'


def parse_tcrdb(data_dir, exp, results_dir):
    fn = exp + '.tsv'
    data = pd.read_csv(os.path.join(data_dir, fn), sep='\t')
    data = data.rename(columns={'AASeq': 'CDR3_beta',
                                'Vregion': 'TRBV_gene',
                                'Jregion': 'TRBJ_gene'})
    for runid in set(data['RunId'].tolist()):
        id_data = data[data['RunId'] == runid]
        id_data = id_data[['TRBV_gene', 'CDR3_beta', 'TRBJ_gene']]
        id_dir = os.path.join(results_dir, exp)
        if not os.path.exists(id_dir):
            os.makedirs(id_dir)
        id_data.to_csv(os.path.join(id_dir, runid + '.tsv'),
                       sep='\t', index=False)

    return


parse_tcrdb(data_dir, 'PRJNA510967', results_dir)

# End task
print('End task on: '+str(datetime.datetime.now()) + '\n')
