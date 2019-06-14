#!/bin/bash


function pre_checks {
    binaryNewton_exists=$( ls | grep -c 'binaryNewton.py' )
    if [ $binaryNewton_exists == 0 ]; then
        (>&2 echo 'This script cannot be run in this directory.') 
        (>&2 echo 'Run in a directory with "binaryNewton.py"')
        exit 1
    fi

    python_exists=$( which python3 )
    if [ $? != 0  ]; then
        (>&2 echo 'This script requires "python3" to be installed.')
        (>&2 echo 'Please install "python3" and then rerun.')
        exit 1
    fi
}

function display_animation {
    if [ -e 'binary1.dat' ]; then
        echo 'Displaying animation'
        BIN_LEN=$(wc -l < binary1.dat)
        NUM_COMMENTS=$(grep -c "#" binary1.dat)
        # echo $BIN_LEN
        # echo $NUM_COMMENTS
        # echo $((BIN_LEN - NUM_COMMENTS))
        gnuplot -c 'binPlot.txt' $((BIN_LEN - NUM_COMMENTS))
    else
        (>&2 echo 'No data file')
    fi

    # how to make a movie
    ## ffmpeg -framerate 100 -i ./images/%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p all_orbit.mp4

}


function main {
    write_success=$( python3 binaryNewton.py $@ > binary1.dat )
    if [ $? != 0 ]; then
        (>&2 echo 'An error occured during python execution!')
        (>&2 echo 'Logging the error')
        date >> errors.log.txt
        cat binary1.dat >> errors.log.txt
        (>&2 echo 'Removing the errant file')
        rm ./binary1.dat
    else
        echo 'The data was produced and written successfully.'
        display_animation
    fi

}

pre_checks
main "$@"
