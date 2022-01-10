CC 		= g++ -Wall -O3 -m64
CC_DEBUG	= g++ -Wall -O3 -m64 -Wextra -g -DDEBUG -lefence
FILES_CC	= graph.cpp io.cpp kstab_model.cpp mstcc_model.cpp ldda.cpp main.cpp
BINARY 		= ldda

UNAME_S := $(shell uname -n)
ifeq ($(UNAME_S),ii3102747)
	GRB_INCLUDE = -I/scratch/gurobi912/linux64/include/
    GRB_LINK    = -L/scratch/gurobi912/linux64/lib/ -lgurobi_g++5.2 -lgurobi91 -lm -lemon
else
	GRB_INCLUDE = -I/opt/gurobi950/linux64/include/
    GRB_LINK    = -L/opt/gurobi950/linux64/lib/ -lgurobi_g++5.2 -lgurobi95 -lm -lemon
endif

all:   clean compile
debug: clean compile_debug

clean:
	find . -name '*.o' -exec rm -f '{}' ';'
	rm -f $(BINARY);

compile:
	$(CC)  -o $(BINARY)  $(FILES_CC)  $(GRB_INCLUDE)  $(GRB_LINK)

compile_debug:
	$(CC_DEBUG)  -o $(BINARY)  $(FILES_CC)  $(GRB_INCLUDE)  $(GRB_LINK)
