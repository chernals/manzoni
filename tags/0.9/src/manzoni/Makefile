# Manzoni
# Global Makefile

#
# General variables
EXEC=manzoni
EXEC_DEBUG=manzonid
VERSION=0.9

#
# Platform dependent options
PLATFORM := $(shell uname -m)
OS := $(shell uname)

#
# Configuration variables
BUFFER_SIZE = 512
DEFS= -DVERSION=$(VERSION) -DBUFFER_SIZE=$(BUFFER_SIZE) -DINTERMEDIATE_TURNS_FLAG -DINVARIANT_FLAG #-DDAMPER_FLAG
DEFS_DEBUG= -DVERSION=$(VERSION) -DBUFFER_SIZE=$(BUFFER_SIZE) -DINVARIANT_FLAG -DDEBUG_FLAG -DDAMPER_FLAG -DITERATOR_TO_FILE_FLAG -DINTERMEDIATE_TURNS_FLAG #-DXML_PARSER_DEBUG #-DBZ_DEBUG

#
# Paths and libraries
BLITZ_PATH = blitz-mercurial
XERCES_PATH = xerces-c-3.1.1
GSL_PATH = gsl-1.9
ifeq ($(OS), Darwin)
LIBS_PATH = /usr/local/lib
LIBRARIES_INCLUDE=
LIBRARIES_LIB=
LIB_PATH= -L/usr/local/lib
endif
ifeq ($(OS), Linux)
LIBS_PATH = /afs/cern.ch/user/c/chernals/public/libs
LIBRARIES_INCLUDE= -I$(LIBS_PATH)/$(BLITZ_PATH) -I$(LIBS_PATH)/$(XERCES_PATH)/include -I$(LIBS_PATH)/$(GSL_PATH)/include
LIBRARIES_LIB=-L$(LIBS_PATH)/$(BLITZ_PATH)/lib -L$(LIBS_PATH)/$(XERCES_PATH)/lib -L$(LIBS_PATH)/$(GSL_PATH)/lib
LIB_PATH = $(LIBRARIES_LIB)
endif
INC_PATH= -I. -I./include -I/include -I /usr/include -I/usr/local/include $(LIBRARIES_INCLUDE)

#
# Compiler
CC=g++
LINKER=g++
OPTIONS= -O3 -funroll-loops -ffast-math -std=gnu++98 -Wno-long-long
WARNINGS= -Wall -Wextra -Wconversion -Wshadow -Wold-style-cast -Wformat-y2k -Wformat-nonliteral -Wformat-security -Winit-self -Wuninitialized -Wswitch-default -Wswitch-enum -Wunused -Wunused-parameter -Wunknown-pragmas -Wundef -Wtrigraphs -Wctor-dtor-privacy -Woverloaded-virtual
ifeq ($(PLATFORM), x86_64)
OPT= $(OPTIONS) -m64 -Wall -Wextra -Wtrigraphs -Wunused -Wswitch-enum -Wno-pragmas -Wno-unknown-pragmas
else
OPT= $(OPTIONS) $(WARNINGS) #-pedantic
endif
COMP_ROOT_FLAGS = `root-config --cflags`
COMP_FLAGS= $(COMP_ROOT_FLAGS) $(OPT) -c $(INC_PATH) $(DEFS)
COMP_FLAGS_DEBUG= $(COMP_ROOT_FLAGS) $(OPT) -c $(INC_PATH) $(DEFS_DEBUG)
COMP= $(CC) $(COMP_FLAGS)
COMP_DEBUG= $(CC) $(COMP_FLAGS_DEBUG)
LINK_ROOT_FLAGS= `root-config --glibs`
LINK_FLAGS= $(LINK_ROOT_FLAGS) -lblitz -lgsl -lgslcblas -lxerces-c
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
                    utils/signals.o \
		    						utils/Logger.o \
                    utils/Singleton.o \
		    						utils/Arguments.o \
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
                    utils/signals.od \
		    						utils/Logger.od \
                    utils/Singleton.od \
		    						utils/Arguments.od \
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
# # Objects
%.o: %.cc %.h
	echo "Compiling in performance mode"
	$(COMP) -o $@ $<

%.od: %.cc %.h
	echo "Compiling in debug mode"
	$(COMP_DEBUG) -o $@ $<

# # Exec
$(EXEC): $(EXEC_DEPENDENCIES)
	$(LINK) -o $(EXEC) $(EXEC_DEPENDENCIES)

$(EXEC_DEBUG): $(EXEC_DEPENDENCIES_DEBUG)
	$(LINK) -o $(EXEC_DEBUG) $(EXEC_DEPENDENCIES_DEBUG)

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
install:
	cp manzoni ~/public/bin && cp manzonid ~/public/bin/manzoni-debug
