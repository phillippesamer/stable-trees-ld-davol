
-- To get COIN-OR::Vol

git clone https://github.com/coin-or/Vol 
cd Vol
./configure -C
make
make install


-- To compile and run the uncapacitated facility location problem example:

cd Vol/examples/VolUfl
make
./ufl


-- To setup MSTCC LDDA_Vol

assuming that building and running the UFL example above worked fine,
make a copy of the VolUfl folder on the same directory (Vol/examples/)

once inside the new folder, download all MSTCC LDDA_Vol files
edit the Makefile preamble, replacing the corresponding lines as follows
(NB! please use the path to your current Gurobi installation where relevant):

DRIVER = main

# CHANGEME: This should be the name of your executable
EXE = $(DRIVER)

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS =  graph.o io.o kstab_model.o mstcc_model.o ldda.o ldda_vol.o $(DRIVER).o

# CHANGEME: Additional libraries
ADDLIBS = -L/opt/gurobi950/linux64/lib/ -lgurobi_g++5.2 -lgurobi95 -lm -lemon

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS = -I/opt/gurobi950/linux64/include/ 

# CHANGEME: Directory to the sources for the (example) problem definition
# files
SRCDIR = .


-- To make and run LDDA_Vol

make
./main [input file]
