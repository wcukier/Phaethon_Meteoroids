Installation and Setup
======================

Installation
------------
The code can be downloaded from the Github page: https://github.com/wcukier/Phaethon_Meteoroids

Dependencies
------------

The code requires the following packages and versions

- matplotlib==3.6.1
- numpy==1.23.3
- pandas==1.5.1
- rebound==3.21.0
- reboundx==3.7.2
- scipy==1.9.3
- spiceypy==5.1.2
- tqdm==4.64.1

These packages can be installed by running the command

.. code-block:: bash

    pip install -r requirements.txt

from the project folder.

SPICE kernels
-------------
There are a number of SPICE kernels that are needed for the code to run.
Running the file downloadKernals.sh should download all kernals not included in the github.


.. code-block:: bash

    chmod +x downloadKernals.sh
    ./downloadKernals.sh


Note: Some of the links have changed so these dropbox links should have everything
necessary.  Email wolfcukier AT gmail... if the links fail.

small kernels: https://www.dropbox.com/scl/fi/ayk4384rjqr9loa7g6vcu/small_kernels.zip?rlkey=iei7x40aelfo5cuv2cqm0xygx&dl=1
de431_part-1: https://www.dropbox.com/scl/fi/ll2jjh7bnr3vz7m2yk9qd/de431_part-1.bsp?rlkey=gl565mjj0waoas2i0shu3fa32&dl=1
de431_part-2: https://www.dropbox.com/scl/fi/dadjwyn5u70ggxi4fdvip/de431_part-2.bsp?rlkey=ryoepyjg8dyic3s2drl7lvwn1&dl=1


Output Directories
------------------
To create the output directories run the file createDirs.sh

.. code-block:: bash

    chmod +x createDirs.sh
    ./createDirs.sh