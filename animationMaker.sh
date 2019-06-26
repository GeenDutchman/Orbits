#!/bin/bash

PRINT_IT="1"
ADJUST_MIN_MAX="1"
data_file="binary1.dat"

function tee_print() {
    if [ $PRINT_IT == "1" ]; then
        printf "$@"
    fi
}

function file_exists() {
    tee_print "Checking for \"$1\" in this directory..."
    exists=$( ls | grep -c $1 )
    if [ $exists == 0 ]; then
        tee_print "not found!!\n"
        return 1
    fi
    tee_print "found.\n"
}

function pre_checks {
    declare -a files=("animationPlot.plt" "images" "$data_file")
    for file_name in ${files[@]}
    do
        file_exists $file_name
        if [ $? != 0 ]; then
            tee_print "Please run in the directory with '$file_name'\n"
            exit 1
        fi
    done

    declare -a programs=("python3" "gnuplot" "ffmpeg")
    for prog_name in ${programs[$@]}
    do
        python_exists=$( which $prog_name )
        if [ $? != 0  ]; then
            tee_print "This script requires '$prog_name' to be installed.\n"
            tee_print "Please install '$prog_name' and then rerun.\n"
            exit 1
        fi
    done
}

function make_animation {
    if [ -e "$data_file" ]; then
        tee_print "Making animation\n"
        # find length of data and ignore comments
        BIN_LEN=$( wc -l < $data_file )
        NUM_COMMENTS=$(grep -c "#" $data_file)
        # get the min and max of x,y,z of star
        IFS=" " read -a MIN_MAXs <<< $( tail -n 1 "$data_file")
        MIN_MAXs=("${MIN_MAXs[@]:1}")
        # if we don't want to use it, then empty MIN_MAXs
        if [ ADJUST_MIN_MAX == "0" ]; then
            MIN_MAXs=""
        fi
	# delete old images
	rm images/*
	# run gnuplot
	tee_print "Running gnuplot to make images\n"
        gnuplot -c 'animationPlot.plt' $((BIN_LEN - NUM_COMMENTS)) ${MIN_MAXs[@]}
    else
        tee_print "No data file for the animation!\n"
	exit 1
    fi

    # how to make a movie
    ## ffmpeg -framerate 100 -i ./images/%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p all_orbit.mp4

}

function main {
    make_animation
    # make the movie
    tee_print "Running ffmpeg to put images together as a movie.\n"
    ffmpeg -framerate 100 -i ./images/%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p orbit.mp4
    tee_print "Done with making movie\n"
}

pre_checks
main
