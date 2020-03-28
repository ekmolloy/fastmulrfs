FastMulRFS Example
==================

**Step 1: Transform MUL-trees into singly-labeled trees.**
```
python ../python-tools/preprocess_multrees_for_fastrfs.py \
    -i g_trees-s2g.trees \
    -a g_trees-s2g-map.txt
```

**Step 2: Run FastRFS (see install instructions [here](../external/README.md)).**
```
../external/FastRFS/build/FastRFS \
        -i g_trees-s2g-for-fastrfs.trees \
        -o fastmulrfs.tree &> fastmulrfs.log
```

**Also see [this](run_fastmulrfs.sh) bash script.**
