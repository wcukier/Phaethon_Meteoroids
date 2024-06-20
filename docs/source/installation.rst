Installation and Setup
======================

Installation
------------
The code can be downloaded from the Github page: https://github.com/wcukier/Phaethon_Meteoroids or downloaded with git

.. code-block:: bash

    git clone https://github.com/wcukier/Phaethon_Meteoroids


Dependencies
------------

The code requires the following packages and versions

 - matplotlib==3.3.4
 - numpy==1.20.1
 - pandas==1.2.4
 - rebound==3.18.0
 - reboundx==3.4.1
 - scipy==1.6.2
 - spiceypy==4.0.2
 - tqdm==4.59.0

These packages can be installed by running the command

.. code-block:: bash

    pip install -r code_reqs.txt

from the project folder.

SPICE kernels
-------------
There are a number of SPICE kernels that are needed for the code to run.
Running the file downloadKernals.sh should download all kernals not included in the github.


.. code-block:: bash

    chmod +x downloadKernals.sh
    ./downloadKernals.sh


If the links in the file fail, please email wolfcukier AT gmail... so I can fix the issue.



Output Directories
------------------
To create the output directories run the file createDirs.sh

.. code-block:: bash

    chmod +x createDirs.sh
    ./createDirs.sh

Data Availability
-----------------
The code output as used in Cukier and Szalay (2023) is available for download at doi.org/10.5281/zenodo.11474474 or https://zenodo.org/uploads/11474474.  Placeing these files in the  :code:`output/cached` directory will allow the data analysis code provided in this repo to read those data.  


