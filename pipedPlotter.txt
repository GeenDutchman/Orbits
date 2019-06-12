set size square
set xrange [-2.0:2.0]
set yrange [-2.0:2.0]
set zrange [-2.0:2.0]
set xlabel 'X'
set ylabel 'Y'
set zlabel 'Z'

print "hello"
set terminal png

do for [i=0:10] {
set output "".i."test.png"
splot '-' u 1:2:3 every ::i::i pt 7 ps 3 title sprintf("time: %f", i*0.001)
}
#set output "test.png"
#splot '-' u 1:2:3 pt 7 ps 3

pause -1 "Hit enter to continue"
