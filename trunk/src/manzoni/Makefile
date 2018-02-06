##
# Manzoni
#
# Author: Cedric Hernaslteens <cedric.hernalsteens@cern.ch>
# 
# European Organization for Nuclear Research
#
# Copyright (c) 2010+ CERN. All rights reserved.
#
##

# Manzoni
# Global Makefile

#
# General variables
EXEC=manzoni
EXEC_DEBUG=manzonid
VERSION=1.0

#
# Platform dependent options
PLATFORM := $(shell uname -m)
OS := $(shell uname)

# Paths and libraries
## MATHEMATICA
ifeq ($(OS), Darwin)
MATH_KERNEL = /Applications/Mathematica.app/Contents/MacOS/MathKernel
MATH_PATH = /Applications/Mathematica.app/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64/CompilerAdditions
endif
ifeq ($(OS), Linux)
MATH_KERNEL = /afs/cern.ch/project/parc/math90/bin/MathKernel
MATH_PATH = /afs/cern.ch/project/parc/math90/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions
endif
## XERCES
XERCES_PATH = xerces-c-3.1.1
## GSL
GSL_PATH = gsl-1.9
## PUTTING IT TOGETHER
LIBRARIES_INC=-isystem$(MATH_PATH)
LIBRARIES_LIB=-L$(MATH_PATH)
ifeq ($(OS), Linux)
LIBS_PATH = /afs/cern.ch/user/c/chernals/public/libs
LIBRARIES_INC= -isystem$(MATH_PATH) -I$(LIBS_PATH)/$(XERCES_PATH)/include -I$(LIBS_PATH)/$(GSL_PATH)/include
LIBRARIES_LIB= -L$(MATH_PATH) -L$(LIBS_PATH)/$(XERCES_PATH)/lib -L$(LIBS_PATH)/$(GSL_PATH)/lib
endif
LIB_PATH= -L/usr/local/lib $(LIBRARIES_LIB)
INC_PATH= -I. -I./include -I/include -I /usr/include -I/usr/local/include $(LIBRARIES_INC)

#
# Configuration variables
BUFFER_SIZE = 512
DEFS= -DMATH_KERNEL_PATH="$(MATH_KERNEL)" -DVERSION=$(VERSION) -DBUFFER_SIZE=$(BUFFER_SIZE) -DINTERMEDIATE_TURNS_FLAG -DINVARIANT_FLAG -DDAMPER_FLAG
DEFS_DEBUG= -DMATH_KERNEL_PATH="$(MATH_KERNEL)" -DVERSION=$(VERSION) -DBUFFER_SIZE=$(BUFFER_SIZE) -DINVARIANT_FLAG -DDEBUG_FLAG -DDAMPER_FLAG -DITERATOR_TO_FILE_FLAG -DINTERMEDIATE_TURNS_FLAG #-DXML_PARSER_DEBUG #-DBZ_DEBUG

#
#
# Compiler
ifeq ($(OS), Darwin)
CC=g++ -stdlib=libstdc++
LINKER=g++ -stdlib=libstdc++ -lstdc++
else
CC=g++
LINKER=g++
endif
MATHTM=${MATH_PATH}/mprep
OPTIONS= -O3 -funroll-loops -ffast-math -std=gnu++98
WARNINGS= -m64 -Wall -Wextra -Wconversion -Wshadow -Wold-style-cast -Wformat-y2k -Wformat-nonliteral -Wformat-security -Winit-self -Wuninitialized -Wswitch-default -Wswitch-enum -Wunused -Wunused-parameter -Wunknown-pragmas -Wundef -Wtrigraphs -Wctor-dtor-privacy -Woverloaded-virtual
ifeq ($(OS), Darwin)
OPT= $(OPTIONS) $(WARNINGS)
else
OPT= $(OPTIONS) -m64 -Wall -Wextra -Wtrigraphs -Wunused -Wswitch-enum -Wno-pragmas -Wno-unknown-pragmas -Wno-long-long
endif

COMP_ROOT_FLAGS = `root-config --cflags`
ifeq ($(OS), Darwin)
# To overcome the problem with the standard library flags
COMP_ROOT_FLAGS = -pthread -m64 -I/usr/local/Cellar/root/5.34.14/include/root
endif

COMP_FLAGS= $(COMP_ROOT_FLAGS) $(OPT) -c $(INC_PATH) $(DEFS)
COMP_FLAGS_DEBUG= $(COMP_ROOT_FLAGS) $(OPT) -c $(INC_PATH) $(DEFS_DEBUG)
COMP= $(CC) $(COMP_FLAGS)
COMP_NO_WARNING= $(CC) $(OPTIONS) -c $(INC_PATH) $(DEFS)
COMP_DEBUG= $(CC) $(COMP_FLAGS_DEBUG)
COMP_DEBUG_NO_WARNING= $(CC) $(OPTIONS) -c $(INC_PATH) $(DEFS_DEBUG)
LINK_ROOT_FLAGS= `root-config --glibs`
ifeq ($(OS), Darwin)
LINK_ROOT_FLAGS=-L/usr/local/Cellar/root/5.34.14/lib/root -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -lm -ldl
endif
ifeq ($(OS), Darwin)
LINK_MATH_FLAGS= -framework Foundation -lMLi4
else
LINK_MATH_FLAGS= -lML64i3
endif
LINK_FLAGS= $(LINK_ROOT_FLAGS) $(LINK_MATH_FLAGS) -lgsl -lgslcblas -lxerces-c
LINK= $(LINKER) $(LIB_PATH) $(LINK_FLAGS)

# Dependencies
EXEC_DEPENDENCIES = main.o \
                    Manzoni.o \
                    Iterator.o \
									  Integrator.o \
                    ParametersEvolution.o \
                    flows/Flow.o \
                    flows/HenonFlow.o \
										flows/SymplecticFlow.o \
                    dataProcessing/procDistrHisto.o \
                    dataProcessing/procDistrProfile.o \
                    dataProcessing/procGraph.o \
                    dataProcessing/procIslands.o \
                    dataProcessing/procTracking.o \
										dataProcessing/procMMA.o \
                    utils/signals.o \
		    						utils/Logger.o \
                    utils/Singleton.o \
		    						utils/Arguments.o \
										utils/MathLink.o \
										utils/XMLInputParser.o \
										distributions/ParticlesDistribution.o \
										distributions/TrackingDistribution.o \
										distributions/ParticlesSet.o \
                    distributions/AADistribution.o \
							      integrators/BaseIntegrator.o \
							      integrators/Pendulum.o \
							      integrators/InterpolatingHenon4.o
		            
EXEC_DEPENDENCIES_DEBUG = main.od \
                    Manzoni.od \
                    Iterator.od \
									  Integrator.od \
                    ParametersEvolution.od \
                    flows/Flow.od \
                    flows/HenonFlow.od \
										flows/SymplecticFlow.od \
                    dataProcessing/procDistrHisto.od \
                    dataProcessing/procDistrProfile.od \
                    dataProcessing/procGraph.od \
                    dataProcessing/procIslands.od \
                    dataProcessing/procTracking.od \
										dataProcessing/procMMA.od \
                    utils/signals.od \
		    						utils/Logger.od \
                    utils/Singleton.od \
		    						utils/Arguments.od \
										utils/MathLink.od \
										utils/XMLInputParser.od \
										distributions/ParticlesDistribution.od \
										distributions/TrackingDistribution.od \
										distributions/ParticlesSet.od \
                    distributions/AADistribution.od \
							      integrators/BaseIntegrator.od \
							      integrators/Pendulum.od \
							      integrators/InterpolatingHenon4.od

# Rules
#
.PHONY: clean cleand cleanall count install
# # Exec
$(EXEC): $(EXEC_DEPENDENCIES)
	$(LINK) -o $(EXEC) $(EXEC_DEPENDENCIES)

$(EXEC_DEBUG): $(EXEC_DEPENDENCIES_DEBUG)
	$(LINK) -o $(EXEC_DEBUG) $(EXEC_DEPENDENCIES_DEBUG)

# # Objects
%.o: %.cc %.h
	$(COMP) -o $@ $<
	
%.od: %.cc %.h
	$(COMP_DEBUG) -o $@ $<
	
# # MathLink templates
utils/MathLink.cc: utils/MathLink.tm
	$(MATHTM) $? -o $@
utils/MathLink.o: utils/MathLink.cc utils/MathLink.h
	$(COMP_NO_WARNING) -o $@ $<
utils/MathLink.od: utils/MathLink.cc utils/MathLink.h
	$(COMP_DEBUG_NO_WARNING) -o $@ $<
	
# # Clean
clean:
	rm -Rf *.o* */*.o* */*/*.o*
	
cleand:
	rm -Rf *.od* */*.od* */*/*.od*

cleanall:
	rm -Rf *.o* */*.o* */*/*/o* *.od* */*.od* */*/*/od* *.pdf *.dat logs/* $(EXEC) $(EXEC_DEBUG)

# # Count lines
count:
	python ./script/count

# # Install on lxplus
ifeq ($(OS), Linux)
install:
		cp manzoni ~/public/bin/manzoni && cp manzonid ~/public/bin/manzoni/manzoni-debug
endif
