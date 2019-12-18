from InSituToolkit.imaging_database import write_experiment
import imaging_db.database.db_operations as db_ops
import imaging_db.filestorage.s3_storage as s3_storage
import imaging_db.filestorage.local_storage as local_storage
import imaging_db.utils.db_utils as db_utils
import os

# information taken from InSituToolkit notebooks, which happened to use the LungInsitu dataset
# we must pull all image id's that contain "SDG" and make an experiment for them.
# We assume that each file is a single round of imaging

spot_channels = ['Opal_570_low', 'Opal_620', 'Opal_690']
nuc_channel = ['DAPI']
db_credentials = '/Users/andrew.cote/Documents/db_credentials.json'


# accessing imagingDB and pulling the dataset id's. What we want is a set of unique id's that we then
# write into experiment files

def getIDs(db_credentials, string):
    credentials_str = db_utils.get_connection_str(db_credentials)
    with db_ops.session_scope(credentials_str) as session:
        frames = session.query(db_ops.DataSet)

    set_of_ids = list()
    for f in frames:
        name = f.dataset_serial
        if string in name:
            set_of_ids.append(name)
    return list(set_of_ids)


def getPositions(db_credentials, dataset_identifier):
    credentials_str = db_utils.get_connection_str(db_credentials)

    with db_ops.session_scope(credentials_str) as session:
        frames = session.query(db_ops.Frames) \
            .join(db_ops.FramesGlobal) \
            .join(db_ops.DataSet) \
            .filter(db_ops.DataSet.dataset_serial == dataset_identifier)
    positions = set()
    for f in frames:
        positions.add(f.pos_idx)
    return list(positions)

list_of_datasets = getIDs(db_credentials, 'SDG')

# Construct experiment files for each dataset, finding the positions as well

list_of_datasets = getIDs(db_credentials, 'SDG')

for dataset_identifier in list_of_datasets:
    pos = getPositions(db_credentials, dataset_identifier)

    # Saving experiment file to the repo
    output_dir = "/Users/andrew.cote/Documents/In-Situ_Transcriptomics/LungInSitu/experiments/" + dataset_identifier + "/"

    # first try to make the directory
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Experiment already exists: " + dataset_identifier)
    else:
        # print("Successfully created the directory %s " % output_dir)    # write the experiment
        try:
            write_experiment(db_credentials, output_dir, dataset_identifier,
                         spot_channels=spot_channels, nuc_channels=nuc_channel,
                         positions=pos
                         )
        except :
            None
        else:
            print('Experiment created: ' + dataset_identifier)

