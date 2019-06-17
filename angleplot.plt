set xlabel "Time"
set ylabel "Theta"

set xrange[0:60]
set yrange[0:2*pi]

plot 'binary1.dat' u 1:11 w l title "Time vs Theta"

pause -1 "Hit enter to continue"
