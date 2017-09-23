#!/bin/bash

# this script exports PYTHONPATH so that it finds our libraries. 


main()
{
    libs=$(myreadlink "python_libs/" )
    echo $libs
    PYTHONPATH=$PYTHONPATH":"$libs
    export PYTHONPATH
}



# extract absolute directory path 
# https://stackoverflow.com/questions/17577093/how-do-i-get-the-absolute-directory-of-a-file-in-bash
function myreadlink()
{
    (
	cd $(dirname $1)         # or  cd ${1%/*}
	echo $PWD/$(basename $1) # or  echo $PWD/${1##*/}
    )
}


main
