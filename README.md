# chrom-tf-kinetic

This repository contains python-based library for setting up and running TF simulations that combines Hi-C predicted structures and Chip-seq predicted binding landscape of TFs


1. To start using the library activate the environment

```
conda env create -f nuc4d
conda activate nuc4d
```

2. The folder src contains the python codes that has the classes and functions you need for building the initial structures needed for running the simulations using HOOMD. The useful files are:

	i. main_with_bindingsites.py -- this file has the classes needed for full architecture simulations 
