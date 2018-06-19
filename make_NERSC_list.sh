#!/bin/bash
# Script for making filelists work on NERSC
# Script to fill ~/hbbISR/filelists with .txt files representing the dids to be run over.
# Each .txt is filled with the exact path(s) to the local / localGroupDisk copy of the sample files inside the container.

if [ ${#} -lt 1 ]; then
    echo "usage: ${0} files to be used"
    exit -1
fi

for dir in "$@"
do
  cd $dir
  for file in *
  do
    path="$(pwd)"
    touch ~/hbbISR/filelists/$dir\.txt
    echo $path/$file >> ~/hbbISR/filelists/$dir\.txt
  done
  cd ..
done
