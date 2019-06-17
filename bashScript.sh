#!/bin/bash

# only print to stdout and stderr if not equal to "0"
NO_TEE="0"

err_file="errors.log.txt"
#err_flag=0
log_file="hist.log.txt"
#log_flag=0

function tee_print() {
    if  { { [ "$1" == "-nt" ] || [ "$1" == "--no-tee" ]; } && shift; } || [ "$NO_TEE" != "0" ]; then
        # echo "no tee"
        if [ "$1" == "-p" ] || [ "$1" == "--prob" ]; then
            shift 
            echo $@ >> /dev/stderr
        else
            echo $@
        fi
    else # just echo normalish
        # echo "normal echo"
        if [ "$1" == "-p" ] || [ "$1" == "--prob" ]; then     
            shift
            date_string = $(date)
            echo -n $date_string >> $err_file
            echo $@ | tee -a $log_file $err_file
        else
            echo $@ | tee -a $log_file
        fi
    fi

}

function file_exists() {
    tee_print -nt -n "Checking for \"$1\" in this directory..."
    exists=$( ls | grep -c $1 )
    if [ $exists == 0 ]; then
        tee_print -nt "not found!!"
        return 1
    fi
    tee_print -nt "found."
}

function pre_checks {
    declare -a files=("binaryNewton.py" "RK.py")
    for file_name in ${files[@]}
    do
        file_exists $file_name
        if [ $? != 0 ]; then
            tee_print -nt -p "Please run in the directory with '$file_name'"
            exit 1
        fi
    done

    declare -a programs=("python3" "gnuplot")
    for prog_name in ${programs[$@]}
    do
        python_exists=$( which $prog_name )
        if [ $? != 0  ]; then
            tee_print -p "This script requires '$prog_name' to be installed."
            tee_print -p "Please install '$prog_name' and then rerun."
            exit 1
        fi
    done
}

function display_animation {
    if [ -e 'binary1.dat' ]; then
        tee_print -nt 'Displaying animation'
        BIN_LEN=$(wc -l < binary1.dat)
        NUM_COMMENTS=$(grep -c "#" binary1.dat)
        # echo $BIN_LEN
        # echo $NUM_COMMENTS
        # echo $((BIN_LEN - NUM_COMMENTS))
        gnuplot -c 'moviePlot.plt' $((BIN_LEN - NUM_COMMENTS))
    else
        tee_print -p 'No data file for the animation!'
    fi

    # how to make a movie
    ## ffmpeg -framerate 100 -i ./images/%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p all_orbit.mp4

}

function display_plate {
    if [ -e 'binary1.dat' ]; then
        tee_print -nt "Displaying plate"
        gnuplot -c 'platePlot.plt'
    else
        tee_print -p 'No data file for the plate!'
    fi
}

function display_angplate {
    if [ -e 'binary1.dat' ]; then
        tee_print -nt "Displaying r vs theta"
        gnuplot -c 'angleplot.plt'
    else
        tee_print -p 'No data file for the plate!'
    fi
}


function main {
    date >> $log_file
    write_success=$( python3 binaryNewton.py $@ > binary1.dat )
    if [ $? != 0 ]; then
        tee_print -p 'An error occured during python execution!'
        err_msg=$(grep '!!' binary1.dat)
        tee_print -nt 'Logging the error'
        tee_print -p $err_msg        
        # date >> errors.log.txt
        # echo $err_msg >> errors.log.txt
        tee_print -nt 'Removing the errant file'
        rm ./binary1.dat
    else
        tee_print "$@"
        tee_print 'The data was produced and written successfully.'
        #display_animation
        display_plate
	    display_angplate
    fi

}

pre_checks
main "$@"
