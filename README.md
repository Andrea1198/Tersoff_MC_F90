# Montecarlo simulation of phase transitions with Tersoff potential
## Phase transitions of SiO<sub>2</sub> : 
 * crystal ($\alpha$-quartz)
 * liquid
 * glass

In order to run the code, type:
```sh
$ make -(FLAGS) all clean run > output.out
```

<!-- where _FLAGS_ can be DTHREAD=nthreads:
* nthreads is the number of Open-mp parallelization order that the user wants -->

output.out should be substituted with the output filename wanted.

As default steps, the simulation runs:
- 10 steps for initial gdr computation
- 50M + 50M w/ gdr computation steps at 300 K to stabilize the inital crystal configuration
- 50M + 50M w/ gdr computation steps at 7000 K to melt the system
- 50M + 50M w/ gdr computation steps at 3000 K to stabilize the liquid
- 50M + 50M w/ gdr computation steps at 300 K to freeze it a dn obtain a glass

To change simulation steps go to files/input/steps.txt and follow the next directives for the file format (first column is just for reference, it is not needed in the file steps.txt)

STEP INDEX | TEMPERATURE (eV) | STEPS | GDR STEPS | GDR UPDATE FREQUENCY|
-|-|-|-|--|
STEP1|0.03|10|10|1
STEP2|0.7|50000000|50000000|10000
STEP3|0.3|50000000|50000000|10000
STEP4|0.03|50000000|50000000|10000