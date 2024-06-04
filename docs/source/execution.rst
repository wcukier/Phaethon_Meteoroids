Execution
=========

The code is run by calling :code:`main.py` with up to 3-4 command line arguments.
 - The first argument is :code:`k`, the batch number
 - The second argument is the model number, which correspond to the model versions as follows
    - 0: The base model, known in the code base as the "novel" model
    - 1: The violent creation model, known in the code base as the "vel" model
    - 2: The cometary creation model, known in the code as the "distr" model
    - Models 3-5 are similar to models 0-2 but with a young_cometary composition as defined by Wilck and Mann (1996) instead of asteroidal
 - The third argument, :code:`n` is the number of particles to simulate per run.  This should evenly divide 10000 for the base and violent creation models and 100000 for the cometary creation model
 - The fourth argument is the age of the stream.  It defaults to 2000

 :code:`k` must be between 0, and 10000 divided by :code:`n` for the base and violent creation models

 :code:`k` must be between 0, and 100000 divided by :code:`n` for the base and cometary creation models

An example execution is below:

.. code-block:: bash

    python main.py 0 1 100 2000

This will run the code using the Violent creation model, using the first 100 particles, with a stream age of 2000 years.

Note: The way the code works, all 10000 particles for the base and violent creation model and all 100000 particles for the cometary creation model must be simulated for proper read_data


If you have access to a HPC cluster using :code:`slurm`, example :code:`.slurm` files have been provided to automate running
many runs.


