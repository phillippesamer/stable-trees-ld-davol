CC		    = g++ -Wall -O3 -m64 -lemon
CC_DEBUG    = g++ -Wall -O0 -m64 -lemon -g
FILES_CC	= io.cpp model.cpp main.cpp
BINARY		= ldla

UNAME_S := $(shell uname -n)
ifeq ($(UNAME_S),ii3102747)
	GRB_INCLUDE = -I/scratch/gurobi912/linux64/include/
    GRB_LINK    = -L/scratch/gurobi912/linux64/lib/ -lgurobi_g++5.2 -lgurobi91 -lm
else
	GRB_INCLUDE = -I/opt/gurobi912/linux64/include/
    GRB_LINK    = -L/opt/gurobi912/linux64/lib/ -lgurobi_g++5.2 -lgurobi91 -lm
endif

all:	clean compile

clean:
	find . -name '*.o' -exec rm -f '{}' ';'
	rm -f $(BINARY);

compile:
	$(CC)  -o $(BINARY)  $(FILES_CC)  $(GRB_INCLUDE)  $(GRB_LINK)

debug:
	$(CC_DEBUG)  -o $(BINARY)  $(FILES_CC)  $(GRB_INCLUDE)  $(GRB_LINK)
