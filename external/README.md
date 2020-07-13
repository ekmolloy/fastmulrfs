How to get external dependencies on Mac or Linux systems
========================================================

Step 0: Get FastMulRFS
----------------------
```
git clone https://github.com/ekmolloy/fastmulrfs.git
cd fastmulrfs/external
```

Step 1: Get PhyloKit
--------------------
```
git clone https://github.com/pranjalv123/phylokit.git
cd phylokit/
mkdir build
cd build/
cmake -DCMAKE_INSTALL_PREFIX=../../install ../src
make install
cd ../../
```

Step 2: Get PhyloNaut
---------------------
**Step 2a: Download phylonaut.**
```
git clone https://github.com/ekmolloy/phylonaut.git
```

**Step 2b: Build phylonaut.**
```
cd phylonaut
mkdir build
cd build/
cmake -DCMAKE_INSTALL_PREFIX=../../install ../src
make install
cd ../../
```

Step 3: Get FastRFS
---------------------
**Step 3a: Download FastRFS.**
```
git clone https://github.com/ekmolloy/FastRFS.git
```

**Step 3b: Build FastRFS.**
```
cd FastRFS/
mkdir build
cd build/
cmake ../src
make
```

Step 3: Get ASTRAL-III
----------------------
```
git clone https://github.com/smirarab/ASTRAL.git
cd ASTRAL
unzip Astral*zip
mv Astral ..
cd ../../../../example
```

***Go to [example](../example/README.md).***
