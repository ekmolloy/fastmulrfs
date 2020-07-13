How to get external dependencies
====================================

```
git clone https://github.com/ekmolloy/fastmulrfs.git
cd fastmulrfs/external
```

Step 1: Get PhyloKit
----------------------
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
git clone https://github.com/pranjalv123/phylonaut.git
```

**Step 2b: Open file `phylonaut/src/wASTRAL.cpp` and comment out the following lines using `//` sign.**
```
wASTRAL.cpp:56: undefined reference to `Timer::start(
wASTRAL.cpp:74: undefined reference to `Timer::stop(
wASTRAL.cpp:78: undefined reference to `Timer::start(
wASTRAL.cpp:81: undefined reference to `Timer::stop(
wASTRAL.cpp:85: undefined reference to `Timer::start(
wASTRAL.cpp:89: undefined reference to `Timer::stop(
wASTRAL.cpp:94: undefined reference to `Timer::start(
wASTRAL.cpp:98: undefined reference to `Timer::stop(
```

**Step 2c: Open file `phylonaut/src/CMakeLists.txt` and add the following text below line 10**
```
include_directories("../../install/include")
link_directories("../../install/lib")
```

**Step 2d: Build phylonaut.**
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

**Step 3d: Build FastRFS.**
```
cd FastRFS/
mkdir build
cd build/
cmake ../src
make
```

Step 3: Get ASTRAL
------------------
```
git clone https://github.com/smirarab/ASTRAL.git
cd ASTRAL
unzip Astral*zip
mv Astral ..
cd ../../../../example
```

***Go to [example](../example/README.md).***
