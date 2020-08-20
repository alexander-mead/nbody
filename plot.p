reset

#print - use x11 on mac, wxt on linux, set print=1 to create a gif
#if(print==-1) {set term x11 enhanced size 500,500}
if(!exists('print')){print=0}
if(print==0){set term qt enhanced size 700,700}
if(print==1){set term gif animate size 500,500 delay 1 optimize; set output 'movie.gif'}

#L - size of plot in AU
if(!exists('L')){L=1}
print 'Size of plot [au]: L:', L

#f - rotation frequency of plot
if(!exists('f')){f=0}
print 'Rotation frequency [yr^-1]: f:', f

#dim - number of dimensions, either 2 or 3
if(!exists('dim')){dim=2}
print 'Number of dimensions: dim: ', dim

#traj - set to 1 to make particle trajectories, 0 otherwise
if(!exists('traj')){traj=2}
if(traj==0){print 'Trajectories: traj: off'}
if(traj==1){print 'Trajectories: traj: full'}
if(traj==2){print 'Trajectories: traj: limited'}

#n - number of particles in plot (usually this should be the same as the number simulated)
if(!exists('n')){n=2}
print 'Number of particles [n]: ', n

#Animation
if(!exists('animate')){animate=1}
if(animate==0){print 'Static plots'}
if(animate==1){print 'Animating plots'}

#Print some white space
print ''

#Make plot square
set size square

#Make axis labels
if(print==0 || print==-1) {set xlabel 'x / AU'; set ylabel 'y / AU'}
if(print==1) {set format x ''; set format y ''}

#Set ranges for axes
set xrange [-L:L]
set yrange [-L:L]
set zrange [-L:L]

#Remove key
unset key

#Filenames to plot
filename(n) = sprintf("./data/particle_%d.dat", n)

#Convert rotation frequency to angular frequency
w=2.*pi*f

#For other terminals
if(print==0 || print==1) {ptype=7; psize=3}

#For aquaterm
if(print==-1) {ptype=5; psize=2}

#Number of time outputs to plot and loop for animation
nmax=1001
if(animate==0) {nmin=nmax-1}
if(animate==1) {nmin=1}
do for [a=nmin:nmax-1] {

#Do the trajectories
ntraj=20
if(animate==0){b=0}
if(animate==1 && traj==0) {b=a}
if(animate==1 && traj==1) {b=0}
if(animate==1 && traj==2) {b=a-ntraj}

#If doing 2D plot
if(dim==2){
plot for [i=1:n] filename(i) every ::b::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)) w l lw 2 lc i,\
     for [i=1:n] filename(i) every ::a::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)) pt ptype ps psize lc i
}

#If doing 3D plot
if(dim==3){
splot for [i=1:n] filename(i) every ::b::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)):4 w l lw 2 lc i,\
      for [i=1:n] filename(i) every ::a::a using ($2*cos(w*$1)+$3*sin(w*$1)):(-$2*sin(w*$1)+$3*cos(w*$1)):4 pt ptype ps psize lc i
}

if(print==0 || print==-1) {pause .01}

}
