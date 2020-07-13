FastMulRFS Example
==================
**Step 0:** Clone this repository.
```
git clone https://github.com/ekmolloy/fastmulrfs.git
cd example
```

**Step 1:** Preprocess gene family trees. This transforms each MUL-tree into singly-labeled (and potentially unresolved) tree. Importantly, the optimal solution to the MulRF supertree optimization problem for the MUL-trees is the optimal solution to the RF supertree optimization problem for the preprocessed MUL-trees ([Theorem 2](https://doi.org/10.1093/bioinformatics/btaa444)). NOTE: If you haven't already, first install [TreeSwift](https://github.com/niemasd/TreeSwift).
```
python ../python-tools/preprocess_multrees_v3.py \
    -i g_trees-mult.trees \
    -o g_trees-mult-for-fastrfs.trees \
    --verbose
```

**Step 2:** Run FastRFS, a heuristic for estimating a RF supertree, on the preprocessed gene family trees. NOTE: If you haven't already, first install [FastRFS](https://github.com/pranjalv123/fastrfs), as described [here](../external/README.md).
```
../external/FastRFS/build/FastRFS \
        -i g_trees-mult-for-fastrfs.trees \
        -o fastmulrfs.tree &> fastmulrfs.log
```

**Also see [this](run_fastmulrfs.sh) bash script.**
