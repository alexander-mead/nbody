reset

if(print==-1) set term aqua
if(print==1) set term post enh col sol; set output 'orbit.eps'

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
dim=2

if(dim==2) plot for [i=1:n] filename(i) using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)) w l lw 3 noti
if(dim==3) splot for [i=1:n] filename(i) using 2:3:4 w l lw 3 noti