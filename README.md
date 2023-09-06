# ggACene (global gridded Air Conditioning energy) projections

![framework (2)](https://github.com/giacfalk/ggACene/assets/36954873/3ffea590-b7dd-40ef-a706-b879ab5fa2c5)

This data repository hosts output data for SSPs126, 245, 370 and 585 on the estimated and future projected ownership of residential air conditioning, its energy consumption, and the underlying population (useful to quantify the per-capita average consumption or the headcount of people affected by the cooling gap).

## 

## Running the model
To reproduce the model and generate the dataset from scratch, please refer to the following steps:
 - Download input data from the Zenodo data repository: https://doi.org/10.5281/zenodo.7845126
 - Adjust the path folder in the preamble of the sourcer.R script
 - Run the sourcer.R script to train the ML model and make projections

## Data download
Instead, to directly download the latest version of the (pre-generated) **ggACene** dataset (.ncdf files), access the *output_data* folder at the Zenodo data repository: https://doi.org/10.5281/zenodo.7845126
   
## References
Falchetta, G., De Cian, E., Pavanello, F., & Wing, I. S. (in preparation). Global gridded scenarios of residential cooling energy demand to 2050. *Under review*
