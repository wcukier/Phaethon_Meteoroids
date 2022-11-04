Data Output
===========

The data outputted by the simulation will be put into the appropriate output directory dependent on the model.  The base model will be put in "output/novel", the violent creation model in "output/vel", the cometary creation model in "output/distr".  The young_comet composition models are put in respective directories with "_comet" appended to the directory.

Saved are 5 data products

- betaXX.npy: This is the beta value for each of the particles
- elementsXX.npy: This saves the orbital elements after each year for each particle
- massXX.npy: This is the mass of each particle in grams
- particlesXX.npy: This is the state vector of each particle sampled ~10000 times in a 2 year period at the end of the simulation
- simXX.bin: This is a snapshot of the simulation at the end but before the 2 year period where particles are sampled

Data Processing
---------------
The easiest way to process data is through running the jupyter notebook :code:`scripts/gen_figures.ipynb`