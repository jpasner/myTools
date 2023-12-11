#!/bin/bash
# Script which takes in a did list and uses rucio to download each one
while IFS='' read -r line || [[ -n "$line" ]]; do
    rucio get $line
done < "$1"
