set xlabel "Theta"
set ylabel "R"

set xrange[pi:3*pi]
set xtics ('π' pi, '3π' 3*pi)
set yrange[3.6:4.5]

# plot R vs theta, and shift over some of the later orbits to compare with the first few orbits
plot 'binary2.dat' u ($1):($2) w l

pause -1 "Hit enter to continue"
