# EvoEpi
In this repository we provide the code and the minimum documentation that accompany our publication.
It provides an environment to simulate a simple branching process and SEIR epidemiological model including inhost branching process with bottleneck.

# Installation and run
To be able to reproduce or reuse our code, the following softwares are required:
## Requirement
- Python version 2.7 or higher
- Numpy
- Pandas
- Scipy
- Python-constraint
- multiprocess
- pp,
To install all the requirements automatically via minicondo, simply type the following command:,

      pip -r requirement.txt

## Run our code.

The code is organized in three main parts, all located in the folder `src/`
- The analytical results in `src/analytic`:
It contains two files:
    - `analytic.py`: set of python functions that implement all the analytical results in our paper. it is used in `analytical_seir.py`.
    - main.py: which shows a simple running example how to use the `analytic.py`

to run it, simply type in the director src/analytic:
      python `main.py`

- The branching process in src/bp:


- The epidemiological model: src/seir
