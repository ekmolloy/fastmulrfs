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
git checkout 90763d0ebc32e0a4eceea5d190e796e7110e610d
mkdir build
cd build/
cmake -DCMAKE_INSTALL_PREFIX=../../install ../src
make install
cd ../../
```

Step 2: Get PhyloNaut
---------------------
```
git clone https://github.com/ekmolloy/phylonaut.git
cd phylonaut
mkdir build
cd build/
cmake -DCMAKE_INSTALL_PREFIX=../../install ../src
make install
cd ../../
```

Step 3: Get FastRFS
---------------------
```
git clone https://github.com/ekmolloy/FastRFS.git
cd FastRFS/
mkdir build
cd build/
cmake ../src
make
```

Step 4: Get ASTRAL-III
----------------------
```
git clone https://github.com/smirarab/ASTRAL.git
cd ASTRAL
unzip Astral*zip
mv Astral ..
cd ../../../../example
```

***Go to [example](../example/README.md).***
