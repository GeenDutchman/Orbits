set xlabel "X"
set ylabel "Y"
set zlabel "Z"
set xrange[-3:3]
set yrange[-3:3]
set zrange[-0.5:0.5]
# set size square

length=ARG1

print length

do for [i=0:length] {
splot 'binary1.dat' u 2:3:4 every ::i::i ps 1 pt 7 title "Neutron Star", 'binary1.dat' u 5:6:7 every ::i::i ps 3 pt 7 title "Black Hole 1", 'binary1.dat' u 8:9:10 every ::i::i ps 3 pt 7 title "Black Hole 2"
}

pause -1 "Hit enter to continue"