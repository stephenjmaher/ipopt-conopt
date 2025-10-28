# Copyright (C) 2005, 2010 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change the four CHANGEME entries below.   #
##########################################################################

# CHANGEME: This should be the name of your executable
EXE = hs071_cpp

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS = hs071_main.o \
	hs071_nlp.o \
	$(CONOPT_ROOT)/src/IpoptToConoptCallbacks.o

# CHANGEME: Additional libraries
ADDLIBS =

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile Ipopt.   #
##########################################################################

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -O2 -DNDEBUG 

# additional C++ Compiler options for linking
CXXLINKFLAGS = -Wl,--rpath -Wl,/home/maher/Documents/GAMS/Development/conopt/ipopt-conopt/external/conopt-linux-x86_64/lib -Wl,--rpath -Wl,/home/maher/Documents/GAMS/Development/conopt/ipopt-conopt/external/ipopt/install/lib

# CONOPT paths
CONOPT_ROOT = /home/maher/Documents/GAMS/Development/conopt/ipopt-conopt
CONOPT_LIB = $(CONOPT_ROOT)/external/conopt-linux-x86_64/lib
CONOPT_INC = $(CONOPT_ROOT)/external/conopt-linux-x86_64/include
IPOPT_LIB = $(CONOPT_ROOT)/external/ipopt/install/lib
IPOPT_INC = $(CONOPT_ROOT)/external/ipopt/install/include/coin-or
SHIM_INC = $(CONOPT_ROOT)/src

# Include directories
INCL = -I$(CONOPT_INC) -I$(SHIM_INC) -I$(SHIM_INC)/Ipopt -I$(IPOPT_INC) $(ADDINCFLAGS)

# Linker flags
LIBS = -L$(CONOPT_LIB) -L$(IPOPT_LIB) -lconopt -lipopt -lstdc++ -lm -ldl

all: $(EXE)

.SUFFIXES: .cpp .o

$(EXE): $(OBJS)
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $(OBJS) $(ADDLIBS) $(LIBS)

clean:
	rm -rf $(EXE) $(OBJS) ipopt.out conopt.out

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<

# Special rule for IpoptToConoptCallbacks.cpp
$(CONOPT_ROOT)/src/IpoptToConoptCallbacks.o: $(CONOPT_ROOT)/src/IpoptToConoptCallbacks.cpp
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<
