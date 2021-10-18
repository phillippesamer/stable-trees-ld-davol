CC		= g++ -Wall -O3 -m64

GRB_INCLUDE_NYBOKS      = -I/opt/gurobi811/linux64/include/
GRB_LINK_NYBOKS   = -L/opt/gurobi811/linux64/lib/ -lgurobi_g++5.2 -lgurobi81 -lm

FILES_CC	= io.cpp model.cpp main.cpp

BINARY		= kstab

all:	clean compile

clean:
	find . -name '*.o' -exec rm -f '{}' ';'
	rm -f $(BINARY);

compile:
	$(CC)  -o $(BINARY)  $(FILES_CC)  $(GRB_INCLUDE_NYBOKS)  $(GRB_LINK_NYBOKS)
