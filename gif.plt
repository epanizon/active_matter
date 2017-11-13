#!/usr/bin/gnuplot
#

reset 
set terminal pngcairo 

system('mkdir -p animation')

do for [i=900:999]{
outfile = sprintf('animation/%4i.png',i)
set output outfile
plot [0:10][0:9] 'conf.xyz' index i u 2:3:4:5 w vectors not, '' index i u 2:3 w p ls 7 lc 1 not
}

