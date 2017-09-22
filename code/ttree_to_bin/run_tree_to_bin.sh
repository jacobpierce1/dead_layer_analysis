#!/bin/bash

# define this to disable overwrites 
PROTECT=1

# takes one arg: the name of the .root file to be converted. 
arg=""


# check number of args supplied 
if [[ $# > 2 ]]
then
    echo $#
    echo "ERROR in run: too many arguments"
    exit -1
elif [[ $# < 1 ]]
then
    echo "ERROR: didnt supply .root file to convert."
    exit -1
else

    # check that a .root was supplied 
    if [[ "${1##*.}" = "root" ]]
    then
	prefix=$(basename $1)
	prefix=${prefix%.*}
	echo 'INFO: processing ' $prefix
    else
	echo ${1##*.}
	echo "ERROR: file supplied doesn't have .root suffix."
	exit -1
    fi

    # disable protection if requested.
    if [[ $# == 2 && $2 == "-0" ]]
    then
	echo "INFO: disabling protection for this file."
	PROTECT=0
    fi
fi





# in principle u could supply different output path here and the cpp program would work as expected. 
OUTPUT_DIR=../extracted_ttree_data/$prefix/


# give warning if this dir exists
if [[ -d $OUTPUT_DIR ]]
then
    if [[ $PROTECT == 1 ]]
    then
	echo "ERROR: directory already exists and PROTECT is enabled, exiting." 
	exit -1
    else
	echo "INFO: directory already exists and PROTECT is disabled. contents will be overwritten."
    fi
fi


mkdir -p $OUTPUT_DIR




# run it
root -q "ttree_to_bin.cpp(\"$1\", \"${OUTPUT_DIR}\")"

