
UNAME_S := $(shell uname -n)
ifeq ($(UNAME_S),ii3102747)
    VOL_PATH    = /scratch/volume
    GRB_PATH    = /scratch/gurobi912/linux64
    GRB_INCLUDE = -I$(GRB_PATH)/include/
    GRB_LINK    = -L$(GRB_PATH)/lib/ -lgurobi_g++5.2 -lgurobi91 -lm -lemon
else
    ifeq ($(UNAME_S),pingvin)
        VOL_PATH    = /opt/Vol
        GRB_PATH    = /opt/gurobi951/linux64/
        GRB_INCLUDE = -I$(GRB_PATH)/include/
        GRB_LINK    = -L$(GRB_PATH)/lib/ -lgurobi_g++5.2 -lgurobi95 -lm -lemon
    else
        VOL_PATH    = /opt/volume
        GRB_PATH    = /opt/gurobi950/linux64/
        GRB_INCLUDE = -I$(GRB_PATH)/include/
        GRB_LINK    = -L$(GRB_PATH)/lib/ -lgurobi_g++5.2 -lgurobi95 -lm -lemon
    endif
endif

# In this case we define DRIVER to be the name of the executable and filename
# without extension
DRIVER = main

# CHANGEME: This should be the name of your executable
EXE = $(DRIVER)

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS =  graph.o io.o kstab_model.o mstcc_model.o cut_generator.o ldda.o ldda_vol.o $(DRIVER).o

# CHANGEME: Additional libraries
ADDLIBS = $(GRB_LINK)

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS = $(GRB_INCLUDE)

# CHANGEME: Directory to the sources for the (example) problem definition
# files
SRCDIR = .


##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile the      #
#  COIN-OR Vol package.                                                  #
##########################################################################

COIN_HAS_PKGCONFIG = TRUE
COIN_CXX_IS_CL = #TRUE

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -O3 -m64 -pipe -DNDEBUG -Wparentheses -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wno-unknown-pragmas -Wno-long-long   -DVOL_BUILD

# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,$(VOL_PATH)/lib

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
ifeq ($(COIN_HAS_PKGCONFIG), TRUE)
  INCL = `PKG_CONFIG_PATH=$(VOL_PATH)/lib64/pkgconfig:$(VOL_PATH)/lib/pkgconfig:$(VOL_PATH)/share/pkgconfig: pkg-config --cflags vol`
else
  INCL = 
endif
INCL += $(ADDINCFLAGS)

# Linker flags
ifeq ($(COIN_HAS_PKGCONFIG), TRUE)
  LIBS = `PKG_CONFIG_PATH=$(VOL_PATH)/lib64/pkgconfig:$(VOL_PATH)/lib/pkgconfig:$(VOL_PATH)/share/pkgconfig: pkg-config --libs vol`
else
  ifeq ($(COIN_CXX_IS_CL), TRUE)
    LIBS = -link -libpath:`$(CYGPATH_W) $(VOL_PATH)/lib` libVol.lib 
  else
    LIBS = -L$(VOL_PATH)/lib -lVol 
  endif
endif

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla $(LIBS) $(ADDLIBS)

clean:
	rm -rf $(EXE) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<


.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`

.c.o:
	$(CC) $(CFLAGS) $(INCL) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<


.c.obj:
	$(CC) $(CFLAGS) $(INCL) -c -o $@ `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`
