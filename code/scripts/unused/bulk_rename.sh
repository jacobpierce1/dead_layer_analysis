# rename a lot of files in the same dir that have the same prefix.

#!/bin/bash

if [ $# -ne 3 ]
then
    echo 'ERROR, wrong num args. \nUSAGE: bulk_rename DIR OLDPREFIX NEWPREFIX'
    exit -1
fi


dir=$1
oldprefix=$2
newprefix=$3

oldprefix_length=$(( ${#oldprefix} + 1 ))

cd $dir

for name in $oldprefix*
do
    base=$(echo "$name" | cut -c"$oldprefix_length"- )
    newname=$newprefix$base
    mv $name $newname
done

