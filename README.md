# custom_completeness
 A small Python tool for recaluclating Kepler survey completeness for custom subsets of stars. This tool was used when calculating survey completeness for occurrence rate calculations in bins of stellar mass in Bergsten et al. (in press). 
 
 Usage
 ---
 The ``Overview.ipynb`` notebook includes several examples on how to generate the survey completeness, including cases for binning by stellar properties, altering disposition score cuts, returning the component elements of completeness, and more.
 
 Acknowledgements
 ---
 Stellar samples can be selected by binning within ranges of the parameters available in the Gaia-Kepler Stellar Properties Catalog from [Berger et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020AJ....159..280B/abstract).

The methodology for calculating survey completeness is the same as in the [``epos``](https://github.com/GijsMulders/epos) package from [Mulders et al. (2018)](https://ui.adsabs.harvard.edu/abs/arXiv:1805.08211). For the detection efficiency, we make use of the Kepler DR25 per-target detection contours from [Burke & Catanzarite (2017)](https://exoplanetarchive.ipac.caltech.edu/docs/KSCI-19111-002.pdf) via [``KeplerPORTs``](https://github.com/nasa/KeplerPORTs). For the vetting efficiency, we follow [Thompson et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018ApJS..235...38T/abstract) using Kepler simulated data products.

 Disclaimers
 ---
 Some properties of some stars in the Gaia-Kepler Stellar Properties Catalog are drawn from distributions in lieu of spectroscopic data -- please see [Berger et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020AJ....159..280B/abstract) for more on the presented stellar properties and their origins.
 
 Recent literature (e.g., [Bryson et al. 2020](https://ui.adsabs.harvard.edu/abs/2020AJ....160..200B/abstract)) has suggested that disposition score cuts (used in evaluating the vetting efficiency) are perhaps unneccesary when treating for reliability in occurrence rate calculations. For more on how reliability may be implemented for the latter, see [Bryson et al. 2020](https://ui.adsabs.harvard.edu/abs/2020AJ....160..200B/abstract) or Bergsten et. al 2022 (in press).
 
 Attribution
 --- 
 For works making use of this tool, please cite Bergsten et al. (in press).
 ```
 Citation coming soon!
 ```
 
 Contact
 ---
 For any questions, requests, or ideas for additional features, please contact Galen Bergsten: gbergsten@arizona.edu
