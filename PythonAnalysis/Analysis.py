from InSituToolkit.imaging_database import write_experiment
import imaging_db.database.db_operations as db_ops
import imaging_db.filestorage.s3_storage as s3_storage
import imaging_db.filestorage.local_storage as local_storage
import imaging_db.utils.db_utils as db_utils
import os, csv, pickle

import numpy as np
from skimage import io
from starfish import Experiment, display, Codebook, ExpressionMatrix, FieldOfView, BinaryMaskCollection, LabelImage
from starfish.image import Filter
from starfish.spots import FindSpots, DecodeSpots, AssignTargets
from starfish.types import Axes, Coordinates, Features, FunctionSource, TraceBuildingStrategies

from InSituToolkit.analysis import save_stack

db_credentials = '/Users/andrew.cote/Documents/db_credentials.json'

# TODO: something clever with loading the experiment names, determining the assay, and dataset ID to get all the positions

# load the list of experiments and iterate over all of them TODO
# list_of_experiments = pickle.load(open('list_of_experiments.p'))
# for path in list_of_experiments[0:1]:


'''
this particular experiment has 12 FOVs, 
RNAscope assay #2:
RNAscope staining 4. NKX2-1: C2 Opal dye 620
RNAscope staining 5. SELENBP1: C3 Opal dye 690
RNAscope staining 6. IGFBP3: C1 Opal dye 570
'''


#TODO: setup a small script that loads up the binary mask as output from ilastik into napari along with some fovs and
# send over to Ashley for her to look over.

'''
Plan of attack:
- create a pipeline in Ilastik that can classify regions as stroma or non-stroma
- Produce a binary mask where non-stroma regions are '1' and stroma regions are '0'
- Logical AND this mask with all three channels on the fovs such that we are only left with non-zero intensities
    in non-stroma regions
- find spots using bandpass filter and blob-detector
- construct codebook according to the type of assay used (will need to read off file-names etc)
- produce the count x cell matrix, where there is only 'one' cell which is non-stroma region
- estimate total area of non-stroma region in pixels (could just sum elements of binary mask matrix)
- produce density of target genes, i.e. number of spots/pixel area, for each channel

Some notes:
- binary mask will be the same across all channels
- Check in with Ashley about whether the stroma and non-stroma classification is valid
- could either use the binary mask to force pixel intensities to be zero, or just use it as a label image for the 
    codebook gene assignment
'''

list_of_datasets = pickle.load(open('list_of_experiments.obj', 'rb'))
dict_of_datasets = pickle.load(open('dict_of_experiments.obj', 'rb'))

for dataset in list_of_datasets:
    path = dict_of_datasets[dataset][0]
    assayNo = dict_of_datasets[dataset][1]

    exp = Experiment.from_json(path + 'experiment.json')

    # get all fovs from the experiment
    fovs = [k for k in exp.keys()]

    # To import images into Ilastik, we want to take the max projection of z-Planes across all color channels,
    # since the stroma tissue flouresces brightly and we want to paint these regions in Ilastik.
    for fov in fovs:
        img_stack = next(exp[fov].get_images(FieldOfView.PRIMARY_IMAGES))

        img_stack_reduced = img_stack.reduce({Axes.ZPLANE, Axes.CH}, func='max')
        if not os.path.exists(path + 'ilastik_images/'):
            os.mkdir(path + 'ilastik_images/')
        save_stack(img_stack_reduced, path + 'ilastik_images/' + fov + '.tif')