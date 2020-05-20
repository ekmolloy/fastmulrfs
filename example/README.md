FastMulRFS Example
==================
```
git clone https://github.com/ekmolloy/fastmulrfs.git
cd example
```

**Step 1:** Preprocess gene family trees. This transforms each MUL-tree into singly-labeled (and potentially unresolved) tree, but this does not impact the optimal solution to the MulRF supertree optimization problem, as shown [here](https://doi.org/10.1101/835553). NOTE: If you haven't already, first install [TreeSwift](https://github.com/niemasd/TreeSwift).
```
python ../python-tools/preprocess_multrees_v3.py \
    -i g_trees-mult.trees \
    --verbose
```

**Step 2:** Run FastRFS, a heuristic for estimating a RF supertree. NOTE: If you haven't already, first install [FastRFS](https://github.com/pranjalv123/fastrfs), as described [here](../external/README.md).
```
../external/FastRFS/build/FastRFS \
        -i g_trees-mult-for-fastrfs.trees \
        -o fastmulrfs.tree &> fastmulrfs.log
```

**Also see [this](run_fastmulrfs.sh) bash script.**
