from InSituToolkit.imaging_database import write_experiment
import imaging_db.database.db_operations as db_ops
import imaging_db.filestorage.s3_storage as s3_storage
import imaging_db.filestorage.local_storage as local_storage
import imaging_db.utils.db_utils as db_utils
import os, csv, pickle

import numpy as np
from skimage import io
from starfish import Experiment, display, Codebook, ExpressionMatrix, FieldOfView, BinaryMaskCollection, LabelImage, ImageStack
from starfish.image import Filter
from starfish.spots import FindSpots, DecodeSpots, AssignTargets
from starfish.types import Axes, Coordinates, Features, FunctionSource, TraceBuildingStrategies

from InSituToolkit.analysis import save_stack

db_credentials = '/Users/andrew.cote/Documents/db_credentials.json'
root_path = '/Users/andrew.cote/Documents/In-Situ_Transcriptomics/LungInSitu/'

list_of_datasets = pickle.load(open('list_of_experiments.obj', 'rb'))
dict_of_datasets = pickle.load(open('dict_of_experiments.obj', 'rb'))
CODEBOOK = pickle.load(open('codebook.obj', 'rb'))
export_path = root_path + 'ilastik/reduced_fovs/'

"""
Pipeline procedure: 
 - functions should operate on datasets individually
 - Order of operations:
        1. Reduce images for export to Ilastik
        2. Filter images and find spots
        3. Import classifications from Ilastik tif files and generate masks
        4. Assign genes to cells
        5. TODO: calculate area of non-stroma tissue, add up each type of gene, normalize gene density. 
"""

# Objects used commonly in the below functions:

SpotFinder = FindSpots.BlobDetector(min_sigma=1,
                                    max_sigma=10,
                                    num_sigma=10,
                                    threshold=0.01,
                                    measurement_type='mean'
                                   )

glp = Filter.GaussianLowPass(sigma=1)
ghp = Filter.GaussianHighPass(sigma=3)

def reduceImages_for_Ilastik(dataset, dict_of_datasets, export_path) -> None:
    save_path, exp_name, _ = dict_of_datasets[dataset]

    exp = Experiment.from_json(save_path + 'experiment.json')
    exp_name_safe = exp_name.replace('/', '-')

    # get all fovs from the experiment
    fovs = [k for k in exp.keys()]

    # To import images into Ilastik, we want to take the max projection of z-Planes across all color channels,
    # since the stroma tissue flouresces brightly and we want to paint these regions in Ilastik.
    for fov in fovs:
        img_stack = next(exp[fov].get_images(FieldOfView.PRIMARY_IMAGES))
        img_stack_sel = img_stack.sel({Axes.CH: 0})
        img_stack_reduced = img_stack_sel.reduce({Axes.ZPLANE}, func='max')
        save_stack(img_stack_reduced, export_path + exp_name_safe + fov + '.tif')

def produceMaskFromTif(path: str, img_stack: ImageStack) -> BinaryMaskCollection:
    """
    Converts a segmentation TIF file into a binary mask which can be used for assigning spots to cells
    """
    label_image = io.imread(path)

    # we must also get the physical coords from the original to produce the binary mask
    yc = img_stack.xarray.yc.values
    xc = img_stack.xarray.xc.values
    physical_ticks = {Coordinates.Y: yc, Coordinates.X:xc}

    y = img_stack.xarray.y.values
    x = img_stack.xarray.x.values
    pixel_coords = {Axes.Y: y, Axes.X: x}

    label_im = LabelImage.from_label_array_and_ticks(label_image, \
                                                pixel_ticks=pixel_coords, \
                                                physical_ticks=physical_ticks, \
                                                log = img_stack.log)

    mask = BinaryMaskCollection.from_label_image(label_im)
    return mask, label_image


list_of_datasets = pickle.load(open('list_of_experiments.obj', 'rb'))
dict_of_datasets = pickle.load(open('dict_of_experiments.obj', 'rb'))
CODEBOOK = pickle.load(open('codebook.obj', 'rb'))

# could iterate this next line of all datasets
dataset = list_of_datasets[0]

save_path, exp_name, assayNo = dict_of_datasets[dataset]
codebook = Codebook.from_code_array(CODEBOOK[int(assayNo)])
exp = Experiment.from_json(save_path + 'experiment.json')


"""
Pipeline procedure: 
 - functions should operate on datasets individually
 - Order of operations:
        1. Reduce images for export to Ilastik
        2. Filter images and find spots
        3. Import classifications from Ilastik tif files and generate masks
        4. Assign genes to cells
        5. Calculate area of non-stroma tissue, add up each type of gene, normalize gene density. 
"""
gene_counts_across_fovs = []

for fov in exp.fovs():
    # 1 - already completed
    # 2: Filter images and project Zplanes
    img_stack = next(exp[fov].get_images(FieldOfView.PRIMARY_IMAGES))
    img_stack_f1 = ghp.run(img_stack,  in_place=False)
    img_stack_f2 = glp.run(img_stack_f1,  in_place=False)
    img_proj_z = img_stack_f2.reduce({Axes.ZPLANE}, func='max')

    # 2: Find Spots
    spots = SpotFinder.run(img_proj_z)
    decoder = DecodeSpots.SimpleLookupDecoder(codebook=codebook)
    decoded_intensities = decoder.run(spots=spots)

    # 3:
    path = root_path + 'ilastik/classified_fovs/' + exp_name + fov + '.tif'
    mask, label_image = produceMaskFromTif(path, img_stack)

    # 4:
    al = AssignTargets.Label()
    labeled = al.run(mask, decoded_intensities)
    cg = labeled.to_expression_matrix()

    # 5:
    STROMA_ID = 2
    NON_STROMA_ID = 1

    xdiff = np.mean(np.diff(img_stack.xarray.xc.values))
    ydiff = np.mean(np.diff(img_stack.xarray.yc.values))

    area = xdiff * ydiff

    pixel_area = np.sum(label_image == NON_STROMA_ID)
    physical_area = pixel_area * area

    # NOTE: highly likely that 0 == cell_id for non-stroma tissue, yet to be determined however

    gene_counts = {}
    gene_counts['area'] = physical_area

    for target in codebook.indexes['target']:
        gene_counts[target] = cg.loc[0, target].values

    gene_counts_across_fovs.append(gene_counts)
    print(fov + ' completed')