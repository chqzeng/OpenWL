## Introduction
this subfolder provides tools needed to collect Sentinel-2 MSI (S2MSI) and Landsat 8/9 (LST8) matchup with in situ water samples, to build models of small water bodies Rrs ---> [Chl]

## Tools
here are a list of tools that provided in this package.

-	GEE S2 satellite data extraction (Level 1 and level 2), with user defined location and time (lat, lon, datetime field in a .csv)
-	GEE LST8 satellite data extraction (level 1 and level 2) ,with user defined location and time (lat, lon, datetime field in a .csv)
-	Climate and meteorological data extraction, with user defined location and time   
-	Machine learning training models (RF, lightGBM, TENSORFLOW) , category classification or numeric regression
-	Py6S for Raileigh Correction of level-1 data.
-	[GLORIA](https://www.nature.com/articles/s41597-023-01973-y) dataset for training example
-	 ...


## How to use

- use the `.ipynb` files and open in google collab to follow the code step by step
- use the `.py` files of the same name as part of your project. (python environment managemnet, conda...)

## Contributors
- [Yulun Wu ](Yulun.Wu@uottawa.ca)
- [Chui Zeng](chqzeng@gmail.com)
