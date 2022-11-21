# Montecarlo simulation of phase transitions with Tersoff potential
## Phase transitions of SiO<sub>2</sub> : 
 * crystal ($\alpha$-quartz)
 * liquid
 * glass

In order to run the code, type:
```sh
$ make -(FLAGS) all clean run > output.out
```

where _FLAGS_ can be DTHREAD=nthreads:
* nthreads is the number of Open-mp parallelization order that the user wants

output.out should be substituted with the output filename wanted.

