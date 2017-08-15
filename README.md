#N-body code

This is a direct sum n-body code for integrating systems such as planetary orbits. 'Direct sum' means that it calculates the gravitational force of each particle on every other particle explictly, so the computational time for a run scales as n^2. The code uses a 4th order Runge-Kutta integration scheme and additonally does adaptive time-stepping so that energy and angular momentum are conserved to some user-specified degree of accuracy.

It should compile with any fortran compiler and requires no libraries. I use '>gfortran nbody.f90', which then produces the binary 'a.out'.

It should then be run specifying (a) the input file (b) the simulation time length in years (c) a specification of a centre-of-mass (CM) boost or not (d) an accuracy parameter that determines to what level energy and angular-momentum are conserved to and (e) a folder to output the results to.

(a) The input file should be an ascii file with n rows and columns that are 'mass' 'x' 'y' 'z' 'vx' 'vy' 'vz'. The units of mass are solar masses, distance are AU and speed are (2*pi)*AU/year. Each column should be tab separated. The folder 'input' contains some example input files. (e.g. 2_body.dat)

(b) The time length is simply the number of years into the future you wish to integrate the system to (e.g. 10)

(c) This should be set to either 0 or 1. If set to 1 then the particles are all boosted to the centre-of-mass frame both in terms of position and velocity (e.g. 1)

(d) The accuracy parameter governs how energy and angular momentum are conserved. This then determines the time-stepping. I usually set this to 1e-6, which means both quantities are conserved to one part in a million. (e.g. 1e-6)

(e) The folder for the output data to be stored in (e.g. data)

When the code runs it gives a rough % of completion. The output are files particle_1.dat, particle_2.dat, ... , particle_n.dat which are creatred in the specified output folder. The files give time, x, y, z, vx, vy, vz for every particle in the simulation at 1000 values of 't' that are linearly spaced.

These data can then be visualised using the plotting script 'plot_animate.p', which can do either a static or animated plot. This script can be run by entering gnuplot (type gnuplot) and then doing 'gnuplot> load 'plot_animate.p''. The script requires some things to be speicied which are listed in each plotting file (things like number of particles, dimension for the plot).
