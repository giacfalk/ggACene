# ggACene (global gridded Air Conditioning energy) projections

![framework (2)](https://github.com/giacfalk/ggACene/assets/36954873/3ffea590-b7dd-40ef-a706-b879ab5fa2c5)

This repository contains computer code (in the R programming language) to replicate the ggACene (global gridded Air Conditioning energy) projections dataset, including preparing data, training models, and obtaining gridded projections. For queries: giacomo.falchetta@cmcc.it

### Input data and analysis replication

## Instructions
To reproduce the model and generate the dataset from scratch, please refer to the following steps:
- Download input data "replication_package_input_data.7z" by cloning the following Zenodo data repository [https://doi.org/10.5281/zenodo.12541994](https://doi.org/10.5281/zenodo.12541994)
- Decompress the folder using 7-Zip (https://www.7-zip.org/download.html)
- Open RStudio and adjust the path folder in the sourcer.R script
- Run the sourcer.R script to train the ML model, make projections, and represent result files

### References
Falchetta, G., De Cian, E., Pavanello, F., & Wing, I. S. Inequalities in global residential cooling energy use to 2050. Forthcoming at Nature Communications
