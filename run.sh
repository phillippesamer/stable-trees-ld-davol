#!/bin/bash

##for entry in "input/zhang-instances/z200-400-13660.gcc" "input/zhang-instances/z300-800-3196.gcc"
##for entry in "input/zhang-instances/z100-300-897.gcc" "input/zhang-instances/z100-500-1247.gcc" "input/zhang-instances/z200-400-13660.gcc" "input/zhang-instances/z200-600-1797.gcc" "input/zhang-instances/z200-800-3196.gcc"  "input/zhang-instances/z300-600-31000.gcc" "input/zhang-instances/z300-800-3196.gcc"
#for entry in "input/zhang-instances/"*
#do
#    echo "$entry"
#    #output=$(./ldda $entry | tail -n 1 >> "lp-new.out")
#    output=$(./ldda $entry >> "$entry""_xp1.out")
#    echo "$output"
#done

zhang_all=( 
"input/zhang-instances/z50-200-199.gcc" 
"input/zhang-instances/z50-200-398.gcc" 
"input/zhang-instances/z50-200-597.gcc" 
"input/zhang-instances/z50-200-995.gcc" 
"input/zhang-instances/z50-200-type2-3903.gcc" 
"input/zhang-instances/z50-200-type2-4877.gcc" 
"input/zhang-instances/z50-200-type2-5864.gcc" 
"input/zhang-instances/z100-300-448.gcc" 
"input/zhang-instances/z100-300-897.gcc" 
"input/zhang-instances/z100-300-1344.gcc" 
"input/zhang-instances/z100-300-type2-8609.gcc" 
"input/zhang-instances/z100-300-type2-10686.gcc" 
"input/zhang-instances/z100-300-type2-12761.gcc" 
"input/zhang-instances/z100-500-1247.gcc" 
"input/zhang-instances/z100-500-2495.gcc" 
"input/zhang-instances/z100-500-3741.gcc" 
"input/zhang-instances/z100-500-6237.gcc" 
"input/zhang-instances/z100-500-12474.gcc" 
"input/zhang-instances/z100-500-type2-24740.gcc" 
"input/zhang-instances/z100-500-type2-30886.gcc" 
"input/zhang-instances/z100-500-type2-36827.gcc" 
"input/zhang-instances/z200-400-13660.gcc" 
"input/zhang-instances/z200-400-17089.gcc" 
"input/zhang-instances/z200-400-20469.gcc" 
"input/zhang-instances/z200-600-1797.gcc" 
"input/zhang-instances/z200-600-3594.gcc" 
"input/zhang-instances/z200-600-5391.gcc" 
"input/zhang-instances/z200-600-34504.gcc" 
"input/zhang-instances/z200-600-42860.gcc" 
"input/zhang-instances/z200-600-50984.gcc" 
"input/zhang-instances/z200-800-3196.gcc" 
"input/zhang-instances/z200-800-6392.gcc" 
"input/zhang-instances/z200-800-9588.gcc" 
"input/zhang-instances/z200-800-15980.gcc" 
"input/zhang-instances/z200-800-62625.gcc" 
"input/zhang-instances/z200-800-78387.gcc" 
"input/zhang-instances/z200-800-93978.gcc" 
"input/zhang-instances/z300-600-31000.gcc" 
"input/zhang-instances/z300-600-38216.gcc" 
"input/zhang-instances/z300-600-45310.gcc" 
"input/zhang-instances/z300-800-3196.gcc" 
"input/zhang-instances/z300-800-59600.gcc" 
"input/zhang-instances/z300-800-74500.gcc" 
"input/zhang-instances/z300-800-89300.gcc" 
"input/zhang-instances/z300-1000-4995.gcc" 
"input/zhang-instances/z300-1000-9990.gcc" 
"input/zhang-instances/z300-1000-14985.gcc" 
"input/zhang-instances/z300-1000-96590.gcc" 
"input/zhang-instances/z300-1000-120500.gcc" 
"input/zhang-instances/z300-1000-144090.gcc" 
)

zhang_type1=( 
"input/zhang-instances/z50-200-199.gcc" 
"input/zhang-instances/z50-200-398.gcc" 
"input/zhang-instances/z50-200-597.gcc" 
"input/zhang-instances/z50-200-995.gcc" 
"input/zhang-instances/z100-300-448.gcc" 
"input/zhang-instances/z100-300-897.gcc" 
"input/zhang-instances/z100-300-1344.gcc" 
"input/zhang-instances/z100-500-1247.gcc" 
"input/zhang-instances/z100-500-2495.gcc" 
"input/zhang-instances/z100-500-3741.gcc" 
"input/zhang-instances/z100-500-6237.gcc" 
"input/zhang-instances/z100-500-12474.gcc" 
"input/zhang-instances/z200-600-1797.gcc" 
"input/zhang-instances/z200-600-3594.gcc" 
"input/zhang-instances/z200-600-5391.gcc" 
"input/zhang-instances/z200-800-3196.gcc" 
"input/zhang-instances/z200-800-6392.gcc" 
"input/zhang-instances/z200-800-9588.gcc" 
"input/zhang-instances/z200-800-15980.gcc" 
"input/zhang-instances/z300-800-3196.gcc" 
"input/zhang-instances/z300-1000-4995.gcc" 
"input/zhang-instances/z300-1000-9990.gcc" 
"input/zhang-instances/z300-1000-14985.gcc" 
)

for entry in "${zhang_type1[@]}";
do
    echo "$entry"
    #output=$(./ldda $entry | tail -n 1 >> "lp-new.out")
    output=$(./ldda $entry >> "$entry""_xp1.out")
    echo "$output"
done


