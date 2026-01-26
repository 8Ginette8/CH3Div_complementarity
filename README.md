
# Complementary (biodiversity) conservation for Switzerland <img src="https://speed2zero.ethz.ch/wp-content/uploads/2023/02/SPEED2ZERO_Logo_trans.png" width="300" align="right">

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15629783.svg)](https://doi.org/10.5281/zenodo.15629783)

This repository contains the whole pipeline used to spatially estimate the importance of complementary biodiversity for conservation prioritization in Switzerland. The Complementarity Indicator (CI) measures the importance of a pixel based on its contribution to taxonomic, functional, and phylogenetic diversity ($\alpha$ component), as well as its distinct composition of species ($\beta$ component), within the context of the [SPEED2ZERO](https://speed2zero.ethz.ch/en/) project.

The CI indicator was generated using the spatial conservation planning software Zonation 5.0 (Moilanen et al. 2022), and biodiversity feature inputs (see SDMapCH data below). Zonation operates through an iterative removal process in which, at each step, the pixel contributing least to overall biodiversity representation is eliminated and the remaining pixel values are recalculated. This produces a hierarchical prioritization of the landscape, with pixels most important for biodiversity conservation receiving the highest ranks.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/28fad961-e6a1-4891-b2be-3631f9ca1e32"
    alt="image"
    style="width:600px; height:auto;"
  />
</p>

> *Spatial distribution of the complementarity biodiversity indicator (CI) across Switzerland. Higher values indicate higher pixel contributions for overall complementarity, lower values imply lesser importance.*

The Core Area Zonation 2 (CAZ2) marginal loss rule algorithm was here applied. This algorithm gives a particular focus on improving the protection of the least represented species, while maintaining a reasonable level of conservation across all species.

Zonation also allows users to assign weights to input features to reflect their relative importance in the prioritization. In our case, we incorporated species-level weights based on each species’ phylogenetic and functional uniqueness to give greater importance to those contributing disproportionately to evolutionary history and ecosystem functioning (Grenié et al. 2017), using the [Open Tree of Life](https://doi.org/10.1111/2041-210X.12593) and the [TraitCH](https://doi.org/10.5281/zenodo.15063844) dataset.

## Requirements

#### Hardware
R (Linux), all scripts can also be run under windows except for the prioritization (Zonation) runs.

#### Input data
Calculating the complementarity indicator (CI) is dependent on species habitat suitability maps derived from the [SDMapCH](https://doi.org/10.1038/s41597-025-06037-x) dataset. **Data is natively at 25 x 25 m resolution**, however, for example feasibility and copyright reason, the available data was uploaded on GitHub at **10 x 10 km resolution**.

For further usages and perfect reproducibility, the SDMapCH native database (v1.3) and associated metadata are openly available on [Dryad](https://doi.org/10.5061/dryad.stqjq2cdk).

Also note that v1.0 of the TraitCH dataset was here used. v2.0 is now available folowing the same link at [Zenodo](https://doi.org/10.5281/zenodo.15063844).

## Scripts description
- `open_r_windows.bat` and `open_r_maclinus.sh`: to open R in the right environment.

From there, you can successively run the scripts:
- `1_CH3Div_OToL_phylogenies.R`: pipeline to generate the phylogenies of our 17 taxonomic groups from the ['Open Tree of Life'](https://tree.opentreeoflife.org/opentree/argus/opentree16.1@ott93302) API. 
- `2_CH3Div_Functional_trees.R`: pipeline to generate the functional trees of our taxonomic groups based on [TraitCH](https://doi.org/10.5281/zenodo.15063844) v1.0.
- `3_CH3Div_Uniqueness_metrics.R`: pipeline to generate the species' indices of functional and phylogenetic uniqueness based on above.
- `4_CH3Div_Zonation5_(linux_only).R`: linux pipeline to run the prioritization (Zonation 5) runs.

## Reference
Adde, A., Rey, PL., Külling, N. et al. SDMapCH: a Comprehensive database of >7,500 modelled species habitat suitability maps for Switzerland. Sci Data 12, 1752 (2025). [10.1038/s41597-025-06037-x](https://doi.org/10.1038/s41597-025-06037-x)

Chauvier, Y., Adde, A., Bergamini, A., Rambold, G., Stofer, S., Graf, N., Gross, A., Blaser, S., Roberts, S. P. M., Potts, S., Casanelles Abella, J., Moretti, M., Nobis, M., Theurillat, J.-P., Hofmann, H., Hartwig, A.-M., Claude, F., Saucy, G., & Altermatt, F. (2025). TraitCH: a multi-taxa functional trait dataset for Switzerland and Europe [Data set]. Zenodo. [10.5281/zenodo.15063844](https://doi.org/10.5281/zenodo.15063844)

Grenié, M., Denelle, P., Tucker, C. M., Munoz, F., & Violle, C. (2017). funrar: An R package to characterize functional rarity. Diversity and Distributions, 23(12), 1365-1371.

Moilanen, A., Lehtinen, P., Kohonen, I., Jalkanen, J., Virtanen, E. A., & Kujala, H. (2022). Novel methods for spatial prioritization with applications in conservation, land use planning and ecological impact avoidance. Methods in Ecology and Evolution, 13(5), 1062-1072.

## Citation
Adde, A., Boussange, V., Chauvier, Y., Dahito, M.-A., Früh, J., Graham, C., Pellissier, L., Zimmermann, N., & Altermatt, F. (2025). Spatial biodiversity indicators and a composite index for conservation prioritization in Switzerland, bioRxiv. [10.1101/2025.06.10.657334](https://doi.org/10.1101/2025.06.10.657334)
