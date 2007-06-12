# -*- mode: gnuplot -*-

# This is a gnuplot script.

# Invocation:
# gnuplot banner.view.solution.gp

plot "g0.data" using 2:4:($3-$2):($5-$4) with vectors nohead linewidth 2 linetype 5 title "g0.data", \
     "g0.data.cha.dbf" using 2:4:($3-$2):($5-$4) with vectors nohead linewidth 6 linetype 3 title "g0.data.cha.dbf"

set terminal postscript enhanced color solid
set output "banner.solution.ps"
replot

# comment out the following line if you don't have the ps2pdf program to convert ps to pdf
!ps2pdf banner.solution.ps


