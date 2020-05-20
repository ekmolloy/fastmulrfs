#!/bin/bash

mulrfscorer="../external/MulRF1.2/executables/MulRFScorerMac"
checkv1="../python-tools/check_mulrf_scores_v1.py"
checkv2="../python-tools/check_mulrf_scores_v2.py"
checkv3="../python-tools/check_mulrf_scores_v3.py"


true_rfs=( 1052 12548)

for i in 1 2; do
    j=$[i-1]
    true_rf=${true_rfs[$j]}

    data=$(term=ansi python $checkv1 -s s_tree_${i}.trees \
                                     -g g_trees_${i}-s2g.trees \
                                     -a g_trees_${i}-s2g-map.txt \
                                     -x $mulrfscorer)
    esti_rf=$(echo $data | awk '{print $1}')
    if [ $true_rf == $esti_rf ]; then
        echo "Version 1 passed test $i."
    else
        echo "Version 1 failed test $i, because"
        echo "    $data"
    fi

    data=$(term=ansi python $checkv2 -s s_tree_${i}.trees \
                                     -g g_trees_${i}-s2g.trees \
                                     -a g_trees_${i}-s2g-map.txt \
                                     -x $mulrfscorer)
    esti_rf=$(echo $data | awk '{print $1}')
    if [ $true_rf == $esti_rf ]; then
        echo "Version 2 passed test $i."
    else
        echo "Version 2 failed test $i, because"
        echo "    $data"
    fi

    data=$(term=ansi python $checkv3 -s s_tree_${i}.trees \
                                     -g g_trees_${i}-mult.trees \
                                     -x $mulrfscorer)
    esti_rf=$(echo $data | awk '{print $1}')
    if [ $true_rf == $esti_rf ]; then
        echo "Version 3 passed test $i."
    else
        echo "Version 3 failed test $i, because"
        echo "    $data"
    fi
done
