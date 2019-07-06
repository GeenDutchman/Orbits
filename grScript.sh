#!/bin/bash

# only print to stdout and stderr if not equal to "0"
NO_TEE="0"
ADJUST_MIN_MAX="1"

err_file="errors.log.txt"
err_flag=0
log_file="hist.log.txt"
#log_flag=0
data_file="binary1.dat"

function tee_print() {
    if  { { [ "$1" == "-nt" ] || [ "$1" == "--no-tee" ]; } && shift; } || [ "$NO_TEE" != "0" ]; then
        # echo "no tee"
        if [ "$1" == "-p" ] || [ "$1" == "--prob" ]; then
            shift 
            printf "$@" >> /dev/stderr
        else
            printf "$@"
        fi
    else # just echo normalish
        # echo "normal echo"
        if [ "$1" == "-p" ] || [ "$1" == "--prob" ]; then     
            shift
            if [ $err_flag == 0 ]; then
                date_string=$(date)
                printf "$date_string\n" | tee -a $err_file
                err_flag=1
            fi
            printf "$@" | tee -a $log_file $err_file
        else
            printf "$@" | tee -a $log_file
        fi
    fi

}

function file_exists() {
    tee_print -nt "Checking for \"$1\" in this directory..."
    exists=$( ls | grep -c $1 )
    if [ $exists == 0 ]; then
        tee_print -nt "not found!!\n"
        return 1
    fi
    tee_print -nt "found.\n"
}

function pre_checks {
    declare -a files=("relativisticBinary.py" "RK.py" "opt.py")
    for file_name in ${files[@]}
    do
        file_exists $file_name
        if [ $? != 0 ]; then
            tee_print -nt -p "Please run in the directory with '$file_name'\n"
            exit 1
        fi
    done

    declare -a programs=("python3" "gnuplot")
    for prog_name in ${programs[$@]}
    do
        python_exists=$( which $prog_name )
        if [ $? != 0  ]; then
            tee_print -p "This script requires '$prog_name' to be installed.\n"
            tee_print -p "Please install '$prog_name' and then rerun.\n"
            exit 1
        fi
    done
}

function display_animation {
    if [ -e "$data_file" ]; then
        tee_print -nt "Displaying animation\n"
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
        tee_print -p "No data file for the animation!\n"
    fi

    # how to make a movie
    ## ffmpeg -framerate 100 -i ./images/%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p all_orbit.mp4

}

function display_plate {
    if [ -e "$data_file" ]; then
        tee_print -nt "Displaying plate\n"
        # get the min and max of x,y,z of star
        IFS=' ' read -a MIN_MAXs <<< $( tail -n 1 "$data_file")
        MIN_MAXs=("${MIN_MAXs[@]:1}")
        # if we don't want to use it, then empty MIN_MAXs
        if [ ADJUST_MIN_MAX == "0" ]; then
            MIN_MAXs=""
        fi
        gnuplot -c 'platePlot.plt' ${MIN_MAXs[@]}
    else
        tee_print -p "No data file for the plate!\n"
    fi
}

function display_angplate {
    if [ -e "$data_file" ]; then
        tee_print -nt "Displaying r vs theta\n"
        gnuplot -c 'r_phi_plot.plt'
    else
        tee_print -p "No data file for the plate!\n"
    fi
}

function main {

    if [ $# == 0 ]; then
        date >> $log_file
        tee_print "Running with the default parameters. Example:\n"
        tee_print "\t--star -x 5000 -y 0 -vy 0.014 --tstep 1.0e-2 --tmax 8.5e6 --mratio 1 --sep -x 100\n"
        write_success=$( python3 relativisticBinary.py > "$data_file")
    else
        # if not tee, then don't log the time
        if [ "$NO_TEE" == "0" ]; then
            date >> $log_file
        fi
        # this is for processing the exceptions
        for arg in "$@";
        do
            case $arg in
                -h|--help)
                    python3 relativisticBinary.py --help
                    return 0
                    ;;
                -d|--default)
                    python3 relativisticBinary.py --default
                    return 0
                    ;;
            esac
            tee_print "%s " "$arg"
        done
        tee_print "\n"
        write_success=$( python3 relativisticBinary.py "$@" > "$data_file" )
    fi
    # echo $write_success
    # if it ran incorrectly
    if [ $? != 0 ]; then
        tee_print -p "An error occured during python execution!\n"
        err_flag=1
        err_msg=$(grep '!!' "$data_file")
        tee_print -nt "Logging the error\n"
        if [ -n "$err_msg" ]; then
            tee_print -p "$err_msg\n"
        else
            tee_print -p "No error message found!!\n"
        fi 
        tee_print -nt "Removing the errant file\n"
        rm ./"$data_file"
    else # it ran correctly
        tee_print "The data was produced and written successfully.\n"
        min_max_info=$( tail -n 2 "$data_file")
        tee_print "$min_max_info\n"
        # display_animation
        # display_plate
        # display_angplate
        tee_print "Analyzing data for precession\n"
        precession_analysis=$( python3 opt.py )
        if [ $? != 0 ]; then
            tee_print -p "$precession_analysis\n"
            return 1
        else
            tee_print "$precession_analysis\n"
        fi
    fi

}

pre_checks
main "$@"
