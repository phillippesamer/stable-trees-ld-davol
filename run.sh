#!/bin/bash

#for entry in "input/zhang-instances/"*
#for entry in "input/zhang-instances/z200-400-13660.gcc" "input/zhang-instances/z300-800-3196.gcc"
for entry in "input/zhang-instances/z100-300-897.gcc" "input/zhang-instances/z100-500-1247.gcc" "input/zhang-instances/z200-400-13660.gcc" "input/zhang-instances/z200-600-1797.gcc" "input/zhang-instances/z200-800-3196.gcc"  "input/zhang-instances/z300-600-31000.gcc" "input/zhang-instances/z300-800-3196.gcc"
do
    echo "$entry"
    #output=$(./ldda $entry | tail -n 1 >> "$entry-new.out")
    output=$(./ldda $entry | tail -n 1 >> "lp-new.out")
    echo "$output"
done




