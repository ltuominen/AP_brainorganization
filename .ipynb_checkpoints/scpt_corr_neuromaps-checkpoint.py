"""
Tuominen 2023 antipsychotic exposure
"""

import numpy as np
import pandas as pd
from neuromaps.images import load_data, load_gifti
from neuromaps.datasets import fetch_annotation
from neuromaps.resampling import resample_images
from neuromaps.nulls import alexander_bloch
from neuromaps.parcellate import Parcellater
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

path = "C:/Users/justi/OneDrive - McGill University/MisicLab/collaborations/Tuominen_2023_antipsychotics/"

"""
data
"""

sig_lh = load_gifti(path+'data/lh.sig.mgh.gii')
sig_rh = load_gifti(path+'data/rh.sig.mgh.gii')
mymap = (sig_lh, sig_rh)

"""
list all annotations to compare with
"""

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

"""
get indices for rotating surface data and save e.g.:
rot = alexander_bloch(data=None, atlas='fsaverage', density='164k', n_perm=1000, seed=1234)
np.save(path+'data/alexander-bloch-null_fsaverage_164k.npy', rot)
this way I can reuse the indices instead of generating new rotations for each annotation
"""

rot = alexander_bloch(data=None, atlas='fsaverage', density='164k', n_perm=1000, seed=1234)
np.save(path+'data/alexander-bloch-null_fsaverage_164k.npy', rot)

# indices for rotating surface data
rotated = dict([])
rotated['fsLR-164k'] = np.load(path+'data/alexander-bloch-null_fsLR_164k.npy')
rotated['fsLR-32k'] = np.load(path+'data/alexander-bloch-null_fsLR_32k.npy')
rotated['fsLR-4k'] = np.load(path+'data/alexander-bloch-null_fsLR_4k.npy')
rotated['fsaverage-164k'] = np.load(path+'data/alexander-bloch-null_fsaverage_164k.npy')

nspins = 1000  # number of spins

"""
run correlations
"""

# initialize dictionary to save out later
corrs = {'annotation' : [],
         'pearsonr' : [],
         'pspin' : []}

# initialize null distributions for plotting
nulls = dict([])

for counter, (src, desc, space, den) in enumerate(annotations):

    print(desc+': '+str(counter)+'/'+str(len(annotations)))

    # get the annotation
    target = fetch_annotation(source=src, desc=desc, space=space, den=den)
    newspace = space  # initialize space/densiy of correlation as target space
    newden = den

    # transform annotations to same space
    if not(space=='fsaverage' and den=='164k'):
        mymap_rs, target_rs = resample_images(src=mymap, trg=target,
                                              src_space='fsaverage', trg_space=space,
                                              method='linear', resampling='downsample_only')
        if space == 'MNI152':  # if target map is MNI152 then we're resampling to fsaverage
            newspace = 'fsaverage'
            newden = '164k'
    else:
        mymap_rs = mymap
        target_rs = target

    # get data arrays
    mymap_rs = load_data(mymap_rs)
    target_rs = load_data(target_rs)
    rho = pearsonr(mymap_rs, target_rs)[0]
    corrs['annotation'].append((src, desc))
    corrs['pearsonr'].append(rho)

    # run spin test "manually" to save null distributions for plotting
    n = np.zeros((nspins, ))
    for i in range(nspins):
        n[i] = pearsonr(mymap_rs, target_rs[rotated[newspace+'-'+newden][:, i]])[0]
    nulls[src+'_'+desc] = n
    
    # calculate p-value
    # (we should do a multiple comparisons correction)
    pspin = (1 + sum(abs(n - np.mean(n)) > abs(rho - np.mean(n)))) / (nspins + 1)
    corrs['pspin'].append(pspin)

pd.DataFrame(corrs, index='annotation').to_csv(path+'results/corrs.csv')
pd.DataFrame(nulls).to_csv(path+'results/corrs_null.csv')

"""
plot boxplot
"""

plt.ion()

fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(range(len(annotations)),
           corrs['pearsonr'],
           c=(np.array(corrs['pspin']) < 0.05).astype(int),  # colour points by significance
           s=40)
ax.boxplot(pd.DataFrame(nulls).values, positions=range(len(annotations)))
ax.set_ylabel('pearson r')
ax.set_xticks(range(len(annotations)))
ax.set_xticklabels([item[1] for item in annotations], rotation=90)
fig.tight_layout()
fig.savefig(path+'figures/boxplot_corrs.png')

"""
parcellate data to Desikan Killiany

once I have the appropriate files, I will loop through all the annotations,
and depending on the space of the annotation, parcellate the data to DK
accordingly.
"""

# need to find this file somewhere (I'm sure it exists, I will look for it)
dk_fsaverage_164k = (path + 'data/atlas-desikankilliany-lh.label.gii',
                     path + 'data/atlas-desikankilliany-rh.label.gii')

parcellater = Parcellater(dk_fsaverage_164k, 'fsaverage')
mymap_parc = parcellater.fit_transform(mymap, 'fsaverage', ignore_background_data=True)

