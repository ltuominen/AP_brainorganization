# Antipsychotics and Brain Organization

This is a repository for the code for a project on antipsychotic medication effects on cortical thickness and the underlying brain organization.

## Folders

### Data

The [data](data/) folder contains group level statistical maps:
+ lh.sig.nii
+ rh.sig.nii
+ lh.gamma.nii
+ rh.gamma.nii

And maps corrected for multiple comparisons:
+ lh.perm.th13.abs.sig.masked.mgh
+ lh.perm.th20.abs.sig.masked.mgh
+ rh.perm.th13.abs.sig.masked.mgh
+ rh.perm.th20.abs.sig.masked.mgh

It also contains parcel-wise partial correlation data from the turku sample and from the [ENIGMA study](https://doi.org/10.1016/j.biopsych.2018.04.023).
+ [turku_partial_r.csv](data/turku_partial_r.csv)
+ [ENIGMA_S32_partial_correlation_between_cortical_thickness_and_chlorpromazine_equivalents.csv](data/ENIGMA_S32_partial_correlation_between_cortical_thickness_and_chlorpromazine_equivalents.csv)

### Code

The [code](code/) folder contains the scripts and programs used to process and analyze the neuroimaging data.
`plot_ap_effects_turku.ipynb` plots group level statistical maps to create [figures 1](figures/Figure1.jpg)  and [supp fig 1](figures/Supp_fig1.jpg)
`compare_turku_enigma.ipynb` creates fig [2](figures/Figure2.jpg).
`neuromaps_analysis_plots.ipynb` does the main comparisons between cortical thickness changes
and brain organization in the Turku and ENIGMA samples. It also create figs [3](figures/Figure3.jpg) and [4](figures/Figure4.jpg), and supp fig [2](figures/Supp_fig1.jpg)

### Figures

This [folder](figures/) is an output folder for all figures created by the code

### Parcellations

This [folder](pacellations/) contains the desikankilliany parcellation in fsaverage and MNI152 spaces. It also contains precalculated spins used in the analyses 
