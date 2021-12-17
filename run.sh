#!/bin/bash

#for entry in "input/zhang-instances/"*
for entry in "input/zhang-instances/z200-400-13660.gcc" "input/zhang-instances/z300-800-3196.gcc"
do
    echo "$entry"
    output=$(./ldda $entry | tail -n 3)
    echo "$output"
done




