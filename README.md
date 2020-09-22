# N-body code

This is a direct sum n-body code for integrating systems such as planetary orbits. 'Direct sum' means that it calculates the gravitational force of each particle on every other particle explictly, so the computational time for a run scales as n^2. The code uses a 4th order Runge-Kutta integration scheme and additonally does adaptive time-stepping so that energy and angular momentum are conserved to some user-specified degree of accuracy.

```
git clone --recursive https://github.com/alexander-mead/nbody
```

The `--recursive` is important in order to simultaneously clone the required library. The code should compile easily with the contained `Makefile`, simply type `>make`, by default it is configured to work for `gfortran` but it should work for other `Fortran` compilers.

The code must then be run specifying: 
(a) An input file 
(b) The simulation time length in years 
(c) Number of dimensions (either 2 or 3, using 2 suppresses all z information)
(c) A specification of a centre-of-mass (CM) boost or not 
(d) An accuracy parameter that determines to what level energy and angular-momentum are conserved
(e) The folder in which to output the results

(a) The input file should be an ascii file with n rows, and columns that are 'mass' 'x' 'y' 'z' 'vx' 'vy' 'vz'. The units of mass are Solar masses, distances are AU and speeds are (2*pi)*AU/year. Each column should be tab separated. The folder 'input' contains some example input files. (e.g., 2_body.dat)

(b) The time length is simply the number of years into the future you wish to integrate the system to (e.g., 10)

(c) Should be set to either 0 or 1. If set to 1 then the particles are all boosted to the centre-of-mass frame both in terms of position and velocity (e.g., 1)

(d) The accuracy parameter governs how energy and angular momentum are conserved, which then determines the time-stepping. I usually set this to 1e-6, which means both quantities are conserved to one part in a million. (e.g., 1e-6)

(e) The folder for the output data to be stored in (e.g., data)

An example run command would therefore be
```
./bin/nbody input/2_body.dat 10 3 1 1e-6 data
```
using the initial condition file `input/2_body.dat` to run a simulation lasting `10` years, in `3` dimensions, boosting to the CM frame, with an accuracy of better than one part in a million for energy and angular momentum conservation. The output data will be foud in `data/particle_1.dat` and `data/particle_2.dat`.

When the code runs it gives a rough % of completion. The output are files `particle_1.dat`, `particle_2.dat`, ... , `particle_n.dat` which are creatred in the specified output folder. The files give time, x, y, z, vx, vy, vz for every particle in the simulation at 1000 linearly-spaced values of 't'.

These data can then be visualised using the gnuplot script `plot.p`, which can do either a static or animated plot. This script can be run by entering gnuplot `>gnuplot` and then type `gnuplot> load 'plot.p'`. The script requires some things to be specified that are listed in each plotting file; things like number of particles, dimension for the plot.
