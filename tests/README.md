How to run tests
============
```
git clone https://github.com/ekmolloy/fastmulrfs.git
```

**Step 1:** Install [MulRF](http://genome.cs.iastate.edu/CBL/MulRF).
``` 
cd external
wget http://genome.cs.iastate.edu/CBL/MulRF/Downloads/NewVersion/MulRF1.2.zip
unzip MulRF1.2.zip
rm MulRF1.2.zip 
chmod a+x  MulRF1.2/executables/MulRF*
cd ../tests
```

**Step 2:** Run tests.
```
./run_tests.sh
```

NOTE: These tests pass when the last equation in [Lemma 13](https://doi.org/10.1093/bioinformatics/btaa444) holds.
