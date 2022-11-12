# LD-DAVOL
## Lagrangean Decomposition dual ascent + volume algorithm

Lagrangean Decomposition based algorithm to compute dual bounds for minimum weight stable spanning trees in a graph - aka minimum spanning tree under conflict constraints (MSTCC) problem. Main references about this reformulation and the resulting algorithm design:

- _Towards stronger Lagrangean bounds for stable spanning trees_. 10th International Network Optimization Conference (INOC), Aachen, Germany, 2022. Published in a volume of OpenProceedings (open access). DOI: [10.48786/inoc.2022.06 ](http://dx.doi.org/10.48786/inoc.2022.06)
- _Polyhedral results and stronger Lagrangean bounds for stable spanning trees_. Optimization Letters, 2022, in press (open access). DOI: [10.1007/s11590-022-01949-8](https://doi.org/10.1007/s11590-022-01949-8)

### To get COIN-OR::LEMON:

Download lemon-1.3.1 source code at https://lemon.cs.elte.hu/trac/lemon

Unpack the file (e.g. on /opt) and open a terminal on that directory

```
mkdir build
cd build
cmake ../
make
make check
[optional] sudo make install
```
If LEMON is not installed (skipping the last step above), we need to add -I flags to the makefile indicating where to find the corresponding headers.


### To get COIN-OR::Vol:

```
git clone https://github.com/coin-or/Vol 
cd Vol
./configure -C
make
make install
```

Optionally, we may test COIN-OR::Vol by compiling and running the uncapacitated facility location problem example:

```
cd Vol/examples/VolUfl
make
./ufl
```


### To compile and run LD-davol:

```
git clone https://github.com/phillippesamer/stable-trees-ld-davol
cd stable-trees-ld-davol
```

One must edit the Makefile initial definitions of variables VOL_PATH and GRB_PATH to reflect the root folder where COIN-OR::Vol and Gurobi are installed, respectively. Finally, we compile and run LD-davol:

```
make
./main [input file path] [optimum (optional)]
```
