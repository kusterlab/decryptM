![logo](logo.png?raw=true)

-----------------

[![DOI](https://zenodo.org/badge/523340827.svg)](https://zenodo.org/badge/latestdoi/523340827)
[![DOI](https://img.shields.io/badge/DOI-10.1126%2Fscience.ade3925-%23B31B1B)](https://science.org/doi/10.1126/science.ade3925)
-----------------

## What is it?

**decryptM** is a quantitative chemical-proteomics approach to characterize 
pathway engagement of biomolecules (e.g. drugs) in cells. These molecules are 
applied in a dose- or time-dependent manner. The modulation of ten-thousands 
posttranslational modifications (such as phosphorylation, acetylation, 
and ubiquitinylation) are then recorded by mass spectrometry in parallel. Read 
more about the technique in our recent publication.

Here, you can access the **decryptM pipeline**, which is written completely in 
Python. We hope that this pipeline will enable other researchers to implement 
the decryptM approach in their labs more quickly. 

-----------------

We have developed a better way of analyzing dose-reponse curves using p-values for individual curves.

For more information, please visit: [![Static Badge TO CurveCurator](https://img.shields.io/badge/GitHub-CurveCurator-orange)](https://github.com/kusterlab/curve_curator)

-----------------

## Installation:

#### 1. Download the code of this repository and place it next to your data

#### 2. Install the virtual environment manager anaconda to install package dependencies safely (https://www.anaconda.com/).

#### 3. Move to the repository in your shell
```sh
(base)$ cd <path_to_repository>
```

#### 4. Install the decryptM environment
```sh
(base)$ conda env create -f ./environment.yml
```

#### 5. Activate the decryptM environment
```sh
(base)$ conda activate decryptM
(decryptM)$
```

-----------------

## Run pipeline

#### 1. Search your TMT-data with MaxQuant (v.1.6.x.x)
A copy of version 1.6.12.0 can be downloaded via pride.

#### 2. In case of global proteome data, transform the proteinGroups.txt file
```sh
(decryptM)$ python ./decryptM/pg_transform.py <in_path>
```

#### 3. Fill out the parameter toml-file for each dataset
Each dataset comes with a parameter file. This file contains all necessary 
information for each processing batch. Multiple experiments can be processed 
together if they share the same experimental parameter structure. 

For more information about the toml parameters, see the comments in the toml file.  

#### 4. Run the pipeline script
```sh
(decryptM)$ python ./decryptM/pipeline.py -q -p <toml_path>
```
Depending on the number of curves, the execution may take several hours.

-----------------


