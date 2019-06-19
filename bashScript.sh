#!/bin/bash

# only print to stdout and stderr if not equal to "0"
NO_TEE="0"
ADJUST_MIN_MAX="1"

err_file="errors.log.txt"
#err_flag=0
log_file="hist.log.txt"
#log_flag=0
data_file="binary1.dat"

function tee_print() {
    if  { { [ "$1" == "-nt" ] || [ "$1" == "--no-tee" ]; } && shift; } || [ "$NO_TEE" != "0" ]; then
        # echo "no tee"
        if [ "$1" == "-p" ] || [ "$1" == "--prob" ]; then
            shift 
            echo -e "$@" >> /dev/stderr
        else
            echo -e "$@"
        fi
    else # just echo normalish
        # echo "normal echo"
        if [ "$1" == "-p" ] || [ "$1" == "--prob" ]; then     
            shift
            date_string=$(date)
            echo -e -n $date_string ": " >> $err_file
            echo -e "$@" | tee -a $log_file $err_file
        else
            echo -e "$@" | tee -a $log_file
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
    if [ -e "$data_file" ]; then
        tee_print -nt 'Displaying animation'
        # find length of data and ignore comments
        BIN_LEN=$( wc -l < $data_file )
        NUM_COMMENTS=$(grep -c "#" $data_file)
        # get the min and max of x,y,z of star
        IFS=' ' read -a MIN_MAXs <<< $( tail -n 1 "$data_file")
        MIN_MAXs=("${MIN_MAXs[@]:1}")
        # if we don't want to use it, then empty MIN_MAXs
        if [ ADJUST_MIN_MAX == "0" ]; then
            MIN_MAXs=""
        fi
        gnuplot -c 'moviePlot.plt' $((BIN_LEN - NUM_COMMENTS)) ${MIN_MAXs[@]}
    else
        tee_print -p 'No data file for the animation!'
    fi

    # how to make a movie
    ## ffmpeg -framerate 100 -i ./images/%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p all_orbit.mp4

}

function display_plate {
    if [ -e "$data_file" ]; then
        tee_print -nt "Displaying plate"
        # get the min and max of x,y,z of star
        IFS=' ' read -a MIN_MAXs <<< $( tail -n 1 "$data_file")
        MIN_MAXs=("${MIN_MAXs[@]:1}")
        # if we don't want to use it, then empty MIN_MAXs
        if [ ADJUST_MIN_MAX == "0" ]; then
            MIN_MAXs=""
        fi
        gnuplot -c 'platePlot.plt' ${MIN_MAXs[@]}
    else
        tee_print -p 'No data file for the plate!'
    fi
}

function display_angplate {
    if [ -e "$data_file" ]; then
        tee_print -nt "Displaying r vs theta"
        gnuplot -c 'r_theta_plot.plt'
    else
        tee_print -p 'No data file for the plate!'
    fi
}

function display_time_angplate {
    if [ -e "$data_file" ]; then
        tee_print -nt "Displaying theta vs time"
        gnuplot -c 'theta_time_plot.plt'
    else
        tee_print -p 'No data file for the plate!'
    fi
}



function main {

    if [ $# == 0 ]; then
        date >> $log_file
        tee_print "Running with the default parameters. Example:"
        tee_print "    --star -x 2 -y 0 -vy 1 --tstep 1.0e-2 --tmax 60 --mratio 1 --sep 2"
        write_success=$( python3 binaryNewton.py > "$data_file")
    else
        # declare -a ARG_ARRAY
        # this is for processing the exceptions
        for arg in "$@";
        do
            case $arg in
                -h|--help)
                    tee_print -nt $( python3 binaryNewton.py --help )
                    return 0
                    ;;
                -d|--default)
                    tee_print -nt $( python3 binaryNewton.py --default )
                    return 0
                    ;;
            esac
        done
        date >> $log_file
        tee_print "$@"
        write_success=$( python3 binaryNewton.py "$@" > "$data_file" )
    fi
    # echo $write_success
    # if it ran incorrectly
    if [ $? != 0 ]; then
        tee_print -p 'An error occured during python execution!'
        err_msg=$(grep '!!' "$data_file")
        tee_print -nt 'Logging the error'
        tee_print -p $err_msg        
        # date >> errors.log.txt
        # echo $err_msg >> errors.log.txt
        tee_print -nt 'Removing the errant file'
        rm ./"$data_file"
    else # it ran correctly
        tee_print 'The data was produced and written successfully.'
        display_animation
        display_plate
	    display_angplate
	    display_time_angplate
    fi

}

pre_checks
main "$@"
