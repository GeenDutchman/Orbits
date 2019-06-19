set xlabel "X"
set ylabel "Y"
set zlabel "Z"

x_min=-3
x_max=3
y_min=-3
y_max=3
z_min=-1
z_max=1

if (ARGC > 0) {
    x_min=ARG1 - 1
    x_max=ARG2 + 1
    y_min=ARG3 - 1
    y_max=ARG4 + 1
    z_min=ARG5 - 1
    z_max=ARG6 + 1
}

CUBED_VIEW=1
if (CUBED_VIEW == 1) {
    # define function
    max(x, y) = (x > y ? x : y)

    # find the biggest single dimension
    x_dim=max(abs(x_min), x_max)
    y_dim=max(abs(y_min), y_max)
    z_dim=max(abs(z_min), z_max)
    xyz_dim=max(max(x_dim, y_dim), z_dim)

    # assign x    
    x_min=-1 * xyz_dim
    x_max=xyz_dim
    # assign y    
    y_min=-1 * xyz_dim
    y_max=xyz_dim
    # assign z    
    z_min=-1 * xyz_dim
    z_max=xyz_dim
}

set xrange[x_min:x_max]
set yrange[y_min:y_max]
set zrange[z_min:z_max]
set size square

splot 'binary1.dat' u 2:3:4 ps 1 pt 7 title 'Neutron Star', 'binary1.dat' u 5:6:7 ps 3 pt 7 title 'Black Hole 1', 'binary1.dat' u 8:9:10 ps 3 pt 7 title 'Black Hole 2'

pause -1 "Hit enter to continue"
