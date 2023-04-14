#! /group/tuominen/anaconda3/bin/python

# coding: utf-8

import numpy as np
import pandas as pd
from neuromaps.images import load_data, load_gifti, annot_to_gifti
from neuromaps.datasets import fetch_annotation
from neuromaps.resampling import resample_images
from neuromaps.nulls import alexander_bloch, burt2020
from neuromaps.parcellate import Parcellater
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
# from statsmodels.stats.multitest import multipletests
from neuromaps.images import relabel_gifti
from neuromaps import parcellate
from neuromaps import transforms 
from neuromaps.stats import compare_images

# get parcellation
path = '/home/lauri/Documents/neuromaps/'

# load in parcellation file
print('loading parcellations')
dk_fsaverage_10k = (path + 'parcellations/atlas-desikankilliany_space-fsaverage_den-10k_hemi-L.label.gii.gz',
                     path + 'parcellations/atlas-desikankilliany_space-fsaverage_den-10k_hemi-R.label.gii.gz')
dk_fsaverage_164k = (path + 'parcellations/atlas-desikankilliany_space-fsaverage_den-164k_hemi-L.aparc-1.annot',
                     path + 'parcellations/atlas-desikankilliany_space-fsaverage_den-164k_hemi-R.aparc-1.annot')
dk_mni = path + 'parcellations/atlas-desikankilliany_space-MNI_res-1mm.nii.gz'

# make sure label IDs are consecutive across hemispheres
dk_fsaverage_10k = relabel_gifti(dk_fsaverage_10k)
dk_fsaverage_164k = annot_to_gifti(dk_fsaverage_164k)  # this does relabel_gift and also converts the annot file to gifti

# make the parcellaters
print('making parcellaters')
parcellater_fs10k = Parcellater(dk_fsaverage_10k, 'fsaverage')
parcellater_fs164k = Parcellater(dk_fsaverage_164k, 'fsaverage')
parcellater_mni = Parcellater(dk_mni, 'MNI152')


# get enigma table and save partial_r as a vector
# download enigma
print('loading enigma maps')
enigmamap = pd.read_csv(path+'data/ENIGMA_S32_partial_correlation_between_cortical_thickness_and_chlorpromazine_equivalents.csv')
enigmamap.drop([68, 69], inplace=True)  # remove the last two rows
enigma_parc = enigmamap['partial_r'].to_numpy()


# create a list of annotations 
print('get a list of annotations')
annotations = list(fetch_annotation(source=['hcps1200',
                                            'margulies2016',
                                            'raichle',
                                            'ding2010', 
                                            'finnema2016', 
                                            'dubois2015',
                                            'dukart2018',
                                            'gallezot2010',
                                            'gallezot2017',
                                            'hillmer2016',
                                            'jaworska2020',
                                            'kaller2017',
                                            'kantonen2020',
                                            'laurikainen2018',
                                            'normandin2015',
                                            'radnakrishnan2018',
                                            'sandiego2015',
                                            'satterthwaite2014',
                                            'sasaki2012',
                                            'savli2012',
                                            'satterthwaite2014',
                                            'smith2017',
                                            'tuominen',
                                            'neurosynth']).keys())

annotations.extend(fetch_annotation(source=['norgaard2021', 'beliveau2017'], space='fsaverage').keys())

# get annotations and parcellate depending on space and density
print('parcellating data')
parcellated = dict([])
for (src, desc, space, den) in annotations:

    annot = fetch_annotation(source=src, desc=desc, space=space, den=den)

    if space == 'MNI152':
        parcellater = parcellater_mni
    elif space == 'fsaverage' and den == '164k':
        parcellater = parcellater_fs164k
    elif space == 'fsLR' and den == '164k':
        space = 'fsaverage'
        annot = transforms.fslr_to_fsaverage(annot, target_density='164k')
        parcellater = parcellater_fs164k
    elif space == 'fsLR' and den != '164k':
        # unfortunately for fsLR-4k we are upsampling to fsaverage-10k to parcellate but it should be fine
        space = 'fsaverage'
        annot = transforms.fslr_to_fsaverage(annot, target_density='10k')
        parcellater = parcellater_fs10k

    parcellated[desc] = parcellater.fit_transform(annot, space=space, ignore_background_data=True)


# loop over all the annotations, calculate rotation and do spin test, add n_perm and remove splicing of annotations when running the whole thing
print('making comparisons')

for a in annotations[0:3]:
    (src, desc, space, den) = a
    if space == 'MNI152':
        
        rotated = burt2020(parcellated[desc], atlas='MNI152', density='1mm',
                       n_perm=100, seed=1234, parcellation=dk_mni)
        
    elif space == 'fsaverage' and den == '164k':
  
        rotated = alexander_bloch(parcellated[desc], atlas='fsaverage', density='164k',
                                n_perm=100, seed=1234, parcellation=dk_fsaverage_164k)
        
    elif space == 'fsLR' and den == '164k':
        
        rotated = alexander_bloch(parcellated[desc], atlas='fsaverage', density='164k',
                                n_perm=100, seed=1234, parcellation=dk_fsaverage_164k)
        
    elif space == 'fsLR' and den != '164k':
    
        rotated = alexander_bloch(parcellated[desc], atlas='fsaverage', density='10k',
                                n_perm=100, seed=1234, parcellation=dk_fsaverage_10k)

    corr, pval = compare_images(parcellated[desc], enigma_parc, nulls=rotated)
    
    print(src, desc)
    print(f'r = {corr:.3f}, p = {pval:.3f}')





