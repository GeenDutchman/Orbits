set xlabel "Theta"
set ylabel "R"

set xrange[0:2*pi]
set xtics (0, 'π' pi, '2π' 2*pi)
set yrange[3.6:4.5]

# plot R vs theta, and shift over some of the later orbits to compare with the first few orbits
plot 'binary1.dat' u ($12):($11) w l, 'binary1.dat' u ($12 - 2*pi):($11) w l

pause -1 "Hit enter to continue"
