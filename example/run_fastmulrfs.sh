#!/bin/bash

preprocess="../python-tools/preprocess_multrees_v3.py"
fastmulrfs="../external/FastRFS/build/FastRFS"

python $preprocess \
    -i g_trees-mult.trees \
    --verbose

if [ ! -e $fastmulrfs ]; then
    echo "Need to get external dependencies!"
else
    $fastmulrfs \
        -i g_trees-mult-for-fastrfs.trees \
        -o fastmulrfs.tree &> fastmulrfs.log
fi
