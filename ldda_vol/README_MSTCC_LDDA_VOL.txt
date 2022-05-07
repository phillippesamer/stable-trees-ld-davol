-------------------------------------------------------------------------------
-- To get COIN-OR::LEMON:

download lemon-1.3.1 source code at https://lemon.cs.elte.hu/trac/lemon

unpack the file (e.g. on /opt) and open a terminal on that directory

mkdir build

cd build

cmake ../

make

make check

[optional] sudo make install


-------------------------------------------------------------------------------
-- To get COIN-OR::Vol:

git clone https://github.com/coin-or/Vol 

cd Vol

./configure -C

make

make install


-------------------------------------------------------------------------------
-- To compile and run the uncapacitated facility location problem example:

cd Vol/examples/VolUfl

make

./ufl


-------------------------------------------------------------------------------
-- To compile and run LD-davola:

git clone https://github.com/phillippesamer/stable-trees-ld-davola

cd stable-trees-ld-davola

edit the Makefile initial definitions of variables VOL_PATH and GRB_PATH to 
reflect the root folder where COIN-OR::Vol and Gurobi are installed,
respectively

make

./main [input file path] [optimum (optional)]
