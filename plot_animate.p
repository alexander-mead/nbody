reset

if(print==-1) set term x11
if(print==0)  set term aqua
if(print==1)  set term gif animate size 500,500 delay 1 optimize; set output 'movie.gif'

set size square

if(print==0 || print==-1) set xlabel 'x/(AU)'; set ylabel 'y/(AU)'
if(print==1) set format x ''; set format y ''

set xrange [-L:L]
set yrange [-L:L]
set zrange [-L:L]

unset key

filename(n) = sprintf("./data/particle_%d.dat", n)

w=2.*pi*f

do for [a=0:999] {

if(dim==2 && traj==0) {
plot for [i=1:n] filename(i) every ::a::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)) pt 7 ps 3 lc i
}

if(dim==2 && traj==1) {
plot for [i=1:n] filename(i) every ::a::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)) pt 7 ps 3 lc i,\
for [i=1:n] filename(i) every ::0::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)) w l lw 2 lc i
}

if(dim==3 && traj==0) {
splot for [i=1:n] filename(i) every ::a::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)):4 pt 7 ps 2 lc i
}

if(dim==3 && traj==1) {
splot for [i=1:n] filename(i) every ::a::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)):4 pt 7 ps 2 lc i,\
for [i=1:n] filename(i) every ::0::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)):4 w l lw 2 lc i
}

if(print==0 || print==-1) {pause .01}

}
