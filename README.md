FastMulRFS
==========
[FastMulRFS](https://doi.org/10.1093/bioinformatics/btaa444) is a pipeline for estimating species trees from multi-copy gene trees, based on the Robinson-Foulds supertree problem for MUL-trees ([MulRF](https://doi.org/10.1186/1748-7188-8-28)). See this [example](example/README.md).

PYTHON DEPENDENCIES
-------------------
+ Python 2.7 or later (including Python 3)
+ [DendroPy](https://www.dendropy.org) if using version 1
+ [TreeSwift](https://github.com/niemasd/TreeSwift) if using version 2 or 3 (recommended)

OTHER DEPENDENCIES (see install instructions [here](external/README.md))
------------------
+ [FastRFS](https://github.com/ekmolloy/fastrfs)
+ [ASTRAL](https://github.com/smirarab/astral)

CITATIONS
--------
Citation for FastMulRFS:
```
@article{molloy2020fastmulrfs,
  author = {Molloy, Erin K. and Warnow, Tandy},
  title = {{FastMulRFS: Fast and accurate species tree estimation under generic gene duplication and loss models}},
  year = {2020},
  journal = {Bioinformatics},
  volume = {36},
  number = {Supplement_1},
  pages = {i57--i65},
  doi = {10.1093/bioinformatics/btaa444}
}
```

Citations for external dependencies:
```
@article{vachaspati2017fastrfs,
  author = {Pranjal Vachaspati and Tandy Warnow},
  title = {{FastRFS: fast and accurate Robinson-Foulds Supertrees using constrained exact optimization}},
  year = {2017},
  journal = {Bioinformatics},
  volume = {33},
  number = {5},
  pages = {631--639},
  doi = {10.1093/bioinformatics/btw600}
}
```
```
@article{zhang2018astral3,
  author = {Chao Zhang and Maryam Rabiee and Erfan Sayyari and Siavash Mirarab},
  title = {{ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees}},
  year = {2018},
  journal = {BMC Bioinformatics},
  volume = {19},
  number = {Suppl 6},
  pages = {153},
  doi = {10.1186/s12859-018-2129-y}
}
```

LICENSE
-------
See [here](LICENSE.txt) for details.
