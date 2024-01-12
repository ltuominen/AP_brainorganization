# Molecular, physiological and functional features underlying cortical thinning related to antipsychotic medication use

This is a repository for the code for a project on antipsychotic medication effects on cortical thickness and the underlying features of brain organization.

## Folders

### Data

Raw imaging data or individual level demographic/clinical data is not shared.  

The [data](data/) folder contains:

Group level statistical maps:
+ lh.sig.nii
+ rh.sig.nii
+ lh.gamma.nii
+ rh.gamma.nii

Maps corrected for multiple comparisons:
+ lh.perm.th13.abs.sig.masked.mgh
+ lh.perm.th20.abs.sig.masked.mgh
+ rh.perm.th13.abs.sig.masked.mgh
+ rh.perm.th20.abs.sig.masked.mgh

Parcel-wise partial correlation data from the turku sample and from the [ENIGMA study](https://doi.org/10.1016/j.biopsych.2018.04.023).
+ [turku_partial_r.csv](data/turku_partial_r.csv)
+ [ENIGMA_S32_partial_correlation_between_cortical_thickness_and_chlorpromazine_equivalents.csv](data/ENIGMA_S32_partial_correlation_between_cortical_thickness_and_chlorpromazine_equivalents.csv)

Precalculated neurosynth loadings

### Code

The [code](code/) folder contains the scripts and programs used to analyze the neuroimaging data and create all the figures.
`the effect of lifetime antipsychotic exposure on cortical thickness.ipynb` contains R code to test the effects of lifetime exposure to antipsychotics on cortical thickness in the Turku sample and for sensitivity analyses.
`plot_ap_effects_turku.ipynb` plots group level statistical maps to create [figures 1](figures/Figure1.jpg)  and [supp fig 1](figures/Supp_fig1.jpg)
`neuromaps_analysis.ipynb` does the main comparisons between cortical thickness changes
and brain features in the Turku and ENIGMA samples.
`Plot_neurosynth.ipynb` plots the brain features correlations for figure [2](figures/Figure2.jpg) and figure [4](figures/Figure2.jpg). This also creates [supp fig 2](figures/Supp_fig2.jpg)
`compare_turku_enigma.ipynb` creates fig [3](figures/Figure3.jpg) comparing Turku and ENIGMA samples.
`Neurosynth_analysis.ipynb` makes comparison between antipsychotic related cortical thinning and neurosynth
`Plot_neurosynth.ipynb` plots comparisons between antipsychotic related cortical thinning and neurosynth and creates fig [5](figures/Figure3.jpg)

### Figures

This [folder](figures/) is an output folder for all figures created by the code

### Parcellations

This [folder](pacellations/) contains the desikankilliany parcellation in fsaverage and MNI152 spaces. It also contains precalculated spins used in the analyses

### tables

This [folder](tables/) holds all the tables in the manuscript and supplement
