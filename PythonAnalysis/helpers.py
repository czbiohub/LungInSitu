# accessing imagingDB and pulling the dataset id's. What we want is a set of unique id's that we then
# write into experiment files

from InSituToolkit.imaging_database import write_experiment
import imaging_db.database.db_operations as db_ops
import imaging_db.filestorage.s3_storage as s3_storage
import imaging_db.filestorage.local_storage as local_storage
import imaging_db.utils.db_utils as db_utils
import os

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
            .join(db_ops.DataSet)  \
            .filter(db_ops.DataSet.dataset_serial == dataset_identifier) 
    positions = set()
    for f in frames:
        positions.add(f.pos_idx)
    return list(positions)