#!/bin/bash

for entry in "input/zhang-instances/z100-300-"*
do
    echo "$entry"

    n=2
    while [ $n -le 3 ]
    do
        echo -n $n" "
        output=$(./kstab $entry $n | tail -n 1)
        echo $output
        n=$(( n+1 ))
    done
done
