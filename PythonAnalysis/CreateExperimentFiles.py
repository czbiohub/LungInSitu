from InSituToolkit.imaging_database import write_experiment
from .helpers import getPositions
import imaging_db.utils.db_utils as db_utils
import os, csv, pickle

# information taken from InSituToolkit notebooks, which happened to use the LungInsitu dataset
# we must pull all image id's that contain "SDG" and make an experiment for them.
# We assume that each file is a single round of imaging

spot_channels = ['Opal_570_low', 'Opal_620', 'Opal_690']
nuc_channel = ['DAPI']
db_credentials = '/Users/andrew.cote/Documents/db_credentials.json'
credentials_str = db_utils.get_connection_str(db_credentials)

'''
We want to save the experiment files in a specific format / directory naming convention, i.e.

in "/Users/andrew.cote/Documents/In-Situ_Transcriptomics/LungInSitu/experiments/" 

The list of experiments reads like

/TH134_E2_B1_assay2/roi1/<experiment_files>
/TH134_E2_B1_assay2/roi2/<experiment_files>



'''

base_path = "/Users/andrew.cote/Documents/In-Situ_Transcriptomics/LungInSitu/experiments/"
list_of_experiments = []

csv_file = open('metadata_lung.csv')
csv_reader = csv.reader(csv_file, delimiter=',')
line_count = 0
for row in csv_reader:
    if line_count > 0 and line_count < 2:

        # extract relevant columns from the csv file
        dataset_id = row[0]
        roi = row[-1]
        full_name = row[2]

        # construct the directory name according to the convention specified:
        idx = full_name.find('assay')
        dir_name = full_name[0:idx + 6]
        save_path = base_path + dir_name
        save_path_exp = save_path + '/roi' + roi + '/'

        if not os.path.exists(save_path):
            os.mkdir(save_path)

        if not(os.path.exists(save_path_exp)):
            os.mkdir(save_path_exp)

            # pull positions from database and write the experiment file
            positions = getPositions(db_credentials, dataset_id)

            write_experiment(db_credentials, save_path_exp, [dataset_id],
                             spot_channels=spot_channels, nuc_channels=nuc_channel,
                             positions=positions
                             )
            print("wrote experiment to " + save_path_exp)
        list_of_experiments.append(save_path_exp)

    line_count += 1

pickle.dump(list_of_experiments, open('list_of_experiments.p'), 'wb')