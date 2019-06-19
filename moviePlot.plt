set xlabel "X"
set ylabel "Y"
set zlabel "Z"

length=1000
x_min=-3
x_max=3
y_min=-3
y_max=3
z_min=-1
z_max=1

if (ARGC > 0) {
    length=ARG1
}
print length

if (ARGC > 1) {
    x_min=ARG2 - 1
    x_max=ARG3 + 1
    y_min=ARG4 - 1
    y_max=ARG5 + 1
    z_min=ARG6 - 1
    z_max=ARG7 + 1
}

set xrange[x_min:x_max]
set yrange[y_min:y_max]
set zrange[z_min:z_max]
set size square

do for [i=0:length] {
splot 'binary1.dat' u 2:3:4 every ::i::i ps 1 pt 7 title "Neutron Star", 'binary1.dat' u 5:6:7 every ::i::i ps 3 pt 7 title "Black Hole 1", 'binary1.dat' u 8:9:10 every ::i::i ps 3 pt 7 title "Black Hole 2"
}

pause -1 "Hit enter to continue"
