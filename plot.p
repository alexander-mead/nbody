reset

#I use aqua-term on mac and wxt on linux machine. set print=1 to make a postscript file
if(print==-1) set term aqua
if(print==0)  set term wxt
if(print==1)  set term post enh col sol; set output 'orbit.eps'

#L - size of plot in x, y, z in AU
#f - set frequency of rotation for plot. Set to f=0 for plot in inertial frame.
#n - set number of particles in plot (should be the same as the number simulated)
#dim - number of dimensions for plot, can be set to either 2 or 3

set size square

set xlabel 'x/(AU)'
set xrange [-L:L]

set ylabel 'y/(AU)'
set yrange [-L:L]

set zlabel 'z/AU'
set zrange [-L:L]

unset key

filename(n) = sprintf("./data/particle_%d.dat", n)

w=2.*pi*f

if(dim==2) plot for [i=1:n] filename(i) using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)) w l lw 3 noti
if(dim==3) splot for [i=1:n] filename(i) using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)):4 w l lw 3 noti
