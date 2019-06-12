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
    if [ -e 'binary1.txt' ]; then
        echo 'Displaying animation'
        gnuplot 'binPlot.txt'
    else
        (>&2 echo 'No data file')
    fi
}


function main {
    write_success=$( python3 binaryNewton.py $@ > binary1.txt )
    if [ $? != 0 ]; then
        (>&2 echo 'An error occured during python execution!')
        (>&2 echo 'Logging the error')
        date >> errors.log.txt
        cat binary1.txt >> errors.log.txt
        (>&2 echo 'Removing the errant file')
        rm ./binary1.txt
    else
        echo 'The data was produced and written successfully.'
        display_animation
    fi

}

pre_checks
main "$@"
