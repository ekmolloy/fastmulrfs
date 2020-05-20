#!/bin/bash

tools="../python-tools"
fastmulrfs="../external/FastRFS/build/FastRFS"

python $tools/preprocess_multrees_v3.py \
    -i g_trees-mult.trees

if [ ! -e $fastmulrfs ]; then
    echo "Need to get external dependencies!"
else
    $fastmulrfs \
        -i g_trees-mult-for-fastrfs.trees \
        -o fastmulrfs.tree &> fastmulrfs.log
fi
