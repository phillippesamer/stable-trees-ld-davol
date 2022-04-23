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
-- To compile and run MSTCC LDDA_Vol:

git clone https://github.com/phillippesamer/mstcc-ldda

cd mstcc-ldda/ldda_vol/

edit the Makefile initial definitions of variables VOL_PATH and GRB_PATH to 
reflect the root folder where COIN-OR::Vol and Gurobi are installed,
respectively

make

./main [input file path] [optimum (optional)]
