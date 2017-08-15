reset

#print - use x11 on mac, wxt on linux, set print=1 to create a gif
#if(print==-1) {set term x11 enhanced size 500,500}
if(print==0) {set term aqua enhanced size 700,700}
if(print==-1)  {set term wxt}
if(print==1)  {set term gif animate size 500,500 delay 1 optimize; set output 'movie.gif'}

#L - size of plot in AU
#f - rotation frequency of plot
#dim - number of dimensions, either 2 or 3
#traj - set to 1 to make particle trajectories, 0 otherwise
#n - number of particles in plot (usually this should be the same as the number simulated)

set size square

if(print==0 || print==-1) {set xlabel 'x/(AU)'; set ylabel 'y/(AU)'}
if(print==1) {set format x ''; set format y ''}

set xrange [-L:L]
set yrange [-L:L]
set zrange [-L:L]

unset key

filename(n) = sprintf("./data/particle_%d.dat", n)

w=2.*pi*f

#For other terminals
if(print==0 || print==1) {ptype=7; psize=3}

#For aquaterm
if(print==-1) {ptype=5; psize=2}

nmax=1001 #Number of time outputs to plot
if(animate==0) {nmin=nmax-1}
if(animate==1) {nmin=1}
do for [a=nmin:nmax-1] {

if(animate==0){b=0}
if(animate==1 && traj==0) {b=a}
if(animate==1 && traj==1) {b=0}
if(animate==1 && traj==2) {b=a-20}

#if(dim==2 && traj==0) {
#plot for [i=1:n] filename(i) every ::a::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)) pt ptype ps psize lc i
#}

#if(dim==2 && traj==1) {
if(dim==2){
plot for [i=1:n] filename(i) every ::b::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)) w l lw 2 lc i,\
     for [i=1:n] filename(i) every ::a::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)) pt ptype ps psize lc i
}

#if(dim==2 && traj==2) {
#plot for [i=1:n] filename(i) every ::a-20::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)) w l lw 2 lc i,\
#     for [i=1:n] filename(i) every ::a::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)) pt ptype ps psize lc i
#}

#if(dim==3 && traj==0) {
#splot for [i=1:n] filename(i) every ::a::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)):4 pt ptype ps psize lc i
#}

#if(dim==3 && traj==1) {
if(dim==3){
splot for [i=1:n] filename(i) every ::b::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)):4 w l lw 2 lc i,\
      for [i=1:n] filename(i) every ::a::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)):4 pt ptype ps psize lc i
}

#if(dim==3 && traj==2) {
#splot for [i=1:n] filename(i) every ::a-20::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)):4 w l lw 2 lc i,\
#      for [i=1:n] filename(i) every ::a::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)):4 pt ptype ps psize lc i
#}


if(print==0 || print==-1) {pause .01}

}
