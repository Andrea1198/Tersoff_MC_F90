# Montecarlo simulation of phase transitions with Tersoff potential
## Phase transitions of SiO<sub>2</sub> : 
 * crystal ($\alpha$-quartz)
 * liquid
 * glass

In order to **run** the code, type:
```sh
$ make -(FLAGS) all clean run > output.out
```

<!-- where _FLAGS_ can be DTHREAD=nthreads:
* nthreads is the number of Open-mp parallelization order that the user wants -->

***output.out*** should be substituted with the output filename wanted.

As default steps, the simulation runs:
- 10 steps for initial gdr computation
- 50M + 50M w/ gdr computation steps at **300 K** to stabilize the inital crystal configuration
- 50M + 50M w/ gdr computation steps at **7000 K** to melt the system
- 50M + 50M w/ gdr computation steps at **3000 K** to stabilize the liquid
- 50M + 50M w/ gdr computation steps at **300 K** to freeze it a dn obtain a glass

To change simulation steps go to *files/input/steps.txt* and follow the next directives for the file format (**first column is just for reference**, it is not needed in the file steps.txt)

STEP INDEX | TEMPERATURE (eV) | STEPS | GDR STEPS | GDR UPDATE FREQUENCY|
-|-|-|-|-|
STEP1|0.03|10|10|1
STEP2|0.7|50000000|50000000|10000
STEP3|0.3|50000000|50000000|10000
STEP4|0.03|50000000|50000000|10000

---
---

Montecarlo simulations (from now on: **MC sim.**) is that the code perform a random move and then, with a certain probability, it decides to keep the move or reject it.  
In this case the acc/rej probability is given by the energy difference, which is computed depending on the chosen model.  
In our case we decided to perform computations using Tersoff potential, which is a 3-body potential, much more difficult to handle with respect to pair (couple) potentials, like Lennard-Jones (ref: [Tersoff](https://ui.adsabs.harvard.edu/abs/1988PhRvB..37.6991T/abstract)), which is very effective for Silicon system (SiO<sub>2</sub> in out case)

A better description of equations may be found in the *equations* folder (also there are references for the Python code of Molecular Dynamics approach)[^1].

[^1]: This is my bachelor thesis, so it is written in Italian.
