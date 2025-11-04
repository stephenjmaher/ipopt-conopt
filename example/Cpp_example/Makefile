# Copyright (C) 2003, 2010 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change the four CHANGEME entries below.   #
##########################################################################

# CHANGEME: This should be the name of your executable
EXE = cpp_example

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS =  cpp_example.o \
	MyNLP.o \
	$(IP2CO_DIR)/IpoptToConoptCallbacks.o

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

# CONOPT paths and build options
# All paths must be explicitly specified:
#   make CONOPT_DIR=/path/to/conopt IP2CO_DIR=/path/to/ipopt-conopt/src IPOPT_DIR=/path/to/ipopt [OPT=dbg|opt]
# Only check for these variables if not running clean target (or if clean is not the only target)
ifeq ($(filter clean,$(MAKECMDGOALS)),)
ifndef CONOPT_DIR
$(error CONOPT_DIR is not set. Please specify: make CONOPT_DIR=/path/to/conopt ...)
endif
ifndef IP2CO_DIR
$(error IP2CO_DIR is not set. Please specify: make IP2CO_DIR=/path/to/ipopt-conopt/src ...)
endif
ifndef IPOPT_DIR
$(error IPOPT_DIR is not set. Please specify: make IPOPT_DIR=/path/to/ipopt ...)
endif
endif

# OPT defaults to 'opt' if not specified
OPT ?= opt
ifneq ($(OPT),dbg)
ifneq ($(OPT),opt)
$(error OPT must be either 'dbg' or 'opt'. Please specify: make OPT=dbg or make OPT=opt)
endif
endif

# C++ Compiler options (set based on OPT)
ifeq ($(OPT),dbg)
CXXFLAGS = -g -O0 -DDEBUG
else
CXXFLAGS = -O2 -DNDEBUG
endif

# additional C++ Compiler options for linking
CXXLINKFLAGS = -Wl,--rpath -Wl,$(CONOPT_DIR)/lib -Wl,--rpath -Wl,$(IPOPT_DIR)/lib

# Include directories
INCL = -I$(CONOPT_DIR)/include -I$(IP2CO_DIR) -I$(IP2CO_DIR)/Ipopt -I$(IPOPT_DIR)/include/coin-or $(ADDINCFLAGS)

# Linker flags
LIBS = -L$(CONOPT_DIR)/lib -L$(IPOPT_DIR)/lib -lconopt -lipopt -lstdc++ -lm -ldl

all: $(EXE)

.SUFFIXES: .cpp .o

$(EXE): $(OBJS)
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $(OBJS) $(ADDLIBS) $(LIBS)

clean:
	rm -rf $(EXE) $(OBJS) ipopt.out conopt.out

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<

# Special rule for IpoptToConoptCallbacks.cpp
$(IP2CO_DIR)/IpoptToConoptCallbacks.o: $(IP2CO_DIR)/IpoptToConoptCallbacks.cpp
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<
