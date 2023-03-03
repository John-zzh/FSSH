#!/usr/bin/env gnuplot

set terminal pdfcairo dashed size 3.37in, 2.08in font "Myriad Pro,10"  fontscale 0.75


set output "toy-FSSH.pdf"

#set key off
set key fixed right top vertical Right noreverse autotitle nobox samplen 1

set xlabel "position x"

set ylabel "Enery"
set ytics nomirror
set mytics 6
set xtics nomirror
set arrow 1 from -10,-0.01 to 10,-0.01 nohead lc -1 dt 2
set arrow 2 from -10,0 to 10,0 nohead lc -1 dt 2

set xrange [-10:10]
set yrange [-0.02:0.02]
set ytics -0.02,0.01,0.02
plot \
  'trj.txt' using 2:4 with linespoints pt 7 ps 0.1 lw 1 lc -1 title "E",\
  # 'trj.txt' using 2:5 with linespoints pt 7 ps 0.2 lw 1 lc -1 title "E1 derivative"
