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

# CONOPT paths and build options (must be specified: make CONOPT_ROOT=/path/to/ipopt-conopt CONOPT_LIB=/path/to/conopt/lib IPOPT_LIB=/path/to/ipopt/lib [OPT=dbg|opt])
ifndef CONOPT_ROOT
$(error CONOPT_ROOT is not set. Please specify: make CONOPT_ROOT=/path/to/ipopt-conopt CONOPT_LIB=/path/to/conopt/lib IPOPT_LIB=/path/to/ipopt/lib [OPT=dbg|opt])
endif
ifndef CONOPT_LIB
$(error CONOPT_LIB is not set. Please specify: make CONOPT_ROOT=/path/to/ipopt-conopt CONOPT_LIB=/path/to/conopt/lib IPOPT_LIB=/path/to/ipopt/lib [OPT=dbg|opt])
endif
ifndef IPOPT_LIB
$(error IPOPT_LIB is not set. Please specify: make CONOPT_ROOT=/path/to/ipopt-conopt CONOPT_LIB=/path/to/conopt/lib IPOPT_LIB=/path/to/ipopt/lib [OPT=dbg|opt])
endif

# OPT defaults to 'opt' if not specified
OPT ?= opt
ifneq ($(OPT),dbg)
ifneq ($(OPT),opt)
$(error OPT must be either 'dbg' or 'opt'. Please specify: make OPT=dbg or make OPT=opt)
endif
endif

CONOPT_EXTERNAL_LIB = $(CONOPT_ROOT)/external/conopt-linux-x86_64/lib
CONOPT_INC = $(CONOPT_ROOT)/external/conopt-linux-x86_64/include
IPOPT_INC = $(CONOPT_ROOT)/external/ipopt/install/include/coin-or
SHIM_INC = $(CONOPT_ROOT)/src

# C++ Compiler options (set based on OPT)
ifeq ($(OPT),dbg)
CXXFLAGS = -g -O0 -DDEBUG
else
CXXFLAGS = -O2 -DNDEBUG
endif

# additional C++ Compiler options for linking
CXXLINKFLAGS = -Wl,--rpath -Wl,$(CONOPT_LIB) -Wl,--rpath -Wl,$(IPOPT_LIB)

# Include directories
INCL = -I$(CONOPT_INC) -I$(SHIM_INC) -I$(SHIM_INC)/Ipopt -I$(IPOPT_INC) $(ADDINCFLAGS)

# Linker flags
LIBS = -L$(CONOPT_EXTERNAL_LIB) -L$(IPOPT_LIB) -lconopt -lipopt -lstdc++ -lm -ldl

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
