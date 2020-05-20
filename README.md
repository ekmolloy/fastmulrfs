FastMulRFS
==========
[FastMulRFS](https://doi.org/10.1101/835553) is a pipeline for estimating species trees from multi-copy gene trees, based on the Robinson-Foulds supertree problem for MUL-trees ([MulRF](https://doi.org/10.1186/1748-7188-8-28)). See this [example](example/README.md).

PYTHON DEPENDENCIES
-------------------
+ Python 2.7 or later (including Python 3)
+ [DendroPy](https://www.dendropy.org) if using version 1
+ [TreeSwift](https://github.com/niemasd/TreeSwift) if using version 2 or 3 (recommended)

OTHER DEPENDENCIES (see install instructions [here](external/README.md))
------------------
+ [FastRFS](https://github.com/ekmolloy/fastrfs)
+ [ASTRAL](https://github.com/smirarab/astral)

CITATION
--------
```
@article{molloy2020fastmulrfs,
  author = {Molloy, Erin K. and Warnow, Tandy},
  title = {FastMulRFS: Fast and accurate species tree estimation under generic gene duplication and loss models},
  year = {2020},
  journal = {bioRxiv},
  elocation-id = {835553},
  publisher = {Cold Spring Harbor Laboratory},
  doi = {10.1101/835553}
}
```

LICENSE
-------
See [here](LICENSE.txt) for details.
