FastMulRFS Example
==================

**Step 1: Preprocess gene family trees. This transforms each MUL-tree into singly-labeled (and potentially unresolved) tree, but we have shown that this does not impact the optimal solution to the MulRF supertree optimization problem.**
```
python ../python-tools/preprocess_multrees_v3.py \
    -i g_trees-mult.trees
```

**Step 2: Run FastRFS (see install instructions [here](../external/README.md)), a heuristic for estimating a RF supertree.**
```
../external/FastRFS/build/FastRFS \
        -i g_trees-mult-for-fastrfs.trees \
        -o fastmulrfs.tree &> fastmulrfs.log
```

**Also see [this](run_fastmulrfs.sh) bash script.**
