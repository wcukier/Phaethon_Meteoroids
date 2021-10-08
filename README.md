
## Setup
Download the neccesary SPICE files from the links in SPICE/meta.tm and add them to your SPICE folder. 

## Execution
to run call the command line argument:
python main.py *run_number* *model_type* *n_particles*

Where *run_number* is an integer between 0 and 99 inclusive 

and *model_type* is 0,1, or 2
0: Standard Creation event -> output in output/novel
1: Violent Creation Event -> output in output/vel
2: Distributed Creation Even -> output in output/distr

*n_particles* is the number of particles to simulate.  100 was used for standard and violent creation events and 10 was used for distributed creation event.

A bash script or wrapper file may be useful to run multiple trials of the script.
