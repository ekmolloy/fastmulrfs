#!/bin/bash

tools="../python-tools"
fastmulrfs="../external/FastRFS/build/FastRFS"

python $tools/map_species_to_gene_simphy.py \
    -i g_trees.trees

python $tools/preprocess_multrees_for_fastrfs.py \
    -i g_trees-s2g.trees \
    -a g_trees-s2g-map.txt

if [ ! -e $fastmulrfs ]; then
    echo "Need to get external dependencies!"
else
    $fastmulrfs \
        -i g_trees-s2g-for-fastrfs.trees \
        -o fastmulrfs.tree &> fastmulrfs.log
fi
