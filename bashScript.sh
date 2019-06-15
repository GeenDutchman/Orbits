#!/bin/bash

function file_exists() {
    echo -n "Checking for \"$1\" in this directory..."
    exists=$( ls | grep -c $1 )
    if [ $exists == 0 ]; then
        echo "not found!!"
        return 1
    fi
    echo "found."
}

function pre_checks {
    # binaryNewton_exists=$( ls | grep -c 'binaryNewton.py' )
    # if [ $binaryNewton_exists == 0 ]; then
    #     (>&2 echo 'This script cannot be run in this directory.') 
    #     (>&2 echo 'Run in a directory with "binaryNewton.py"')
    #     exit 1
    # fi

    declare -a files=("binaryNewton.py" "RK.py")
    for file_name in ${files[@]}
    do
        file_exists $file_name
        if [ $? != 0 ]; then
            (>&2 echo 'Please run in the directory with this file')
            exit 1
        fi
    done

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
        gnuplot -c 'moviePlot.plt' $((BIN_LEN - NUM_COMMENTS))
    else
        (>&2 echo 'No data file')
    fi

    # how to make a movie
    ## ffmpeg -framerate 100 -i ./images/%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p all_orbit.mp4

}

function display_plate {
    if [ -e 'binary1.dat' ]; then
        echo "Displaying plate"
        gnuplot -c 'platePlot.plt'
    else
        (>&2 echo 'No data file')
    fi
}


function main {
    write_success=$( python3 binaryNewton.py $@ > binary1.dat )
    if [ $? != 0 ]; then
        (>&2 echo 'An error occured during python execution!')
        err_msg=$(grep '!!' binary1.dat)
        (>&2 echo $err_msg)
        (>&2 echo 'Logging the error')
        date >> errors.log.txt
        echo $err_msg >> errors.log.txt
        (>&2 echo 'Removing the errant file')
        rm ./binary1.dat
    else
        echo 'The data was produced and written successfully.'
        display_animation
        display_plate
    fi

}

pre_checks
main "$@"
