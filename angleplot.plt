set xlabel "Theta"
set ylabel "R"

set xrange[0:2*pi]
set yrange[-3:3]

plot 'binary1.dat' u 12:11 w l title "R over Theta Graph"

pause -1 "Hit enter to continue"
