#!/bin/bash

c=0
while read seq; do
    echo ">"$c"_shuffled"
    echo $seq
    let c=c+1
done < $1
