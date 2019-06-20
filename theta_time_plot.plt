set xlabel "Time"
set ylabel "Theta"

set xrange[0:100]
set yrange[0:2*pi]
set ytics (0, 'π' pi, '2π' 2*pi)

# plot R vs theta, and shift over some of the later orbits to compare with the first few orbits
plot 'binary1.dat' u ($1):($12) w l, 'binary1.dat' u ($1):($12 - 1*(2*pi)) w l

pause -1 "Hit enter to continue"
