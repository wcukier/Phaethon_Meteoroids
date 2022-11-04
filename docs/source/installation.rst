Installation and Setup
============

The code can be downloaded from the [Github page](https://github.com/wcukier/Phaethon_Meteoroids)

The code requires the following packages and versions

matplotlib==3.6.1
numpy==1.23.3
pandas==1.5.1
rebound==3.21.0
reboundx==3.7.2
scipy==1.9.3
spiceypy==5.1.2
tqdm==4.64.1

These packages can be installed by running the command
`pip install -r requirements.txt`
from the project folder.

SPICE kernels
=============
There are a number of SPICE kernels that are needed for the code to run.
Running the file `downloadKernals.sh` should download all kernals not included in the github.

Output Directories
==================
To create the output directories run the file `createDirs.sh`