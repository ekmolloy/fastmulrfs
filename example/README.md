FastMulRFS Example
==================

**Step 1:** Preprocess gene family trees. NOTE: This transforms each MUL-tree into singly-labeled (and potentially unresolved) tree, but this does not impact the optimal solution to the MulRF supertree optimization problem, as shown [here](https://doi.org/10.1101/835553).
```
python ../python-tools/preprocess_multrees_v3.py \
    -i g_trees-mult.trees \
    --verbose
```

**Step 2:** Run FastRFS, a heuristic for estimating a RF supertree. NOTE: First download and build FastRFS, as described [here](../external/README.md))
```
../external/FastRFS/build/FastRFS \
        -i g_trees-mult-for-fastrfs.trees \
        -o fastmulrfs.tree &> fastmulrfs.log
```

**Also see [this](run_fastmulrfs.sh) bash script.**
