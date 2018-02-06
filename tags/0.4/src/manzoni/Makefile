# Manzoni
# Global Makefile

#
# General variables
EXEC=manzoni
EXEC_DEBUG=manzonid
VERSION=0.4

#
# Platform dependent options
PLATFORM := $(shell uname -m)

#
# Configuration variables
BUFFER_SIZE = 512
DEFS= -DVERSION=$(VERSION) -DBUFFER_SIZE=$(BUFFER_SIZE) -DDAMPER_FLAG #-DINTERMEDIATE_TURNS_FLAG
DEFS_DEBUG= -DVERSION=$(VERSION) -DBUFFER_SIZE=$(BUFFER_SIZE) -DDEBUG_FLAG -DDAMPER_FLAG -DITERATOR_TO_FILE_FLAG -DINTERMEDIATE_TURNS_FLAG #-DXML_PARSER_DEBUG #-DBZ_DEBUG

#
# Paths and libraries
LIBS_PATH = /afs/cern.ch/user/c/chernals/public/libs
BLITZ_PATH = blitz-0.9
XERCES_PATH = xerces-c-3.1.1
GSL_PATH = gsl-1.9
LIBRARIES_INCLUDE= -I$(LIBS_PATH)/$(BLITZ_PATH)/include -I$(LIBS_PATH)/$(XERCES_PATH)/include -I$(LIBS_PATH)/$(GSL_PATH)/include
LIBRARIES_LIB= -L$(LIBS_PATH)/$(BLITZ_PATH)/lib -L$(LIBS_PATH)/$(XERCES_PATH)/lib -L$(LIBS_PATH)/$(GSL_PATH)/lib
INC_PATH= -I. -I./include -I/include -I /usr/include -I/usr/local/include $(LIBRARIES_INCLUDE)
LIB_PATH= -L/lib64 -L/lib -L/usr/lib64 -L/usr/lib -L/usr/local/lib64 -L/usr/local/lib64 $(LIBRARIES_LIB)

#
# Compiler
CC=g++
OPTIONS= -std=gnu++98 -O3 -funroll-loops -ffast-math
WARNINGS= -Wall -Wextra -Wconversion -Wshadow -Wold-style-cast -Wformat-y2k -Wformat-nonliteral -Wformat-security -Winit-self -Wuninitialized -Wswitch-default -Wswitch-enum -Wunused -Wunused-parameter -Wunknown-pragmas -Wundef -Wtrigraphs -Wctor-dtor-privacy -Woverloaded-virtual
ifeq ($(PLATFORM), x86_64)
OPT= $(OPTIONS) -m64 -Wall -Wextra -Wtrigraphs -Wunknown-pragmas -Wunused -Wswitch-enum
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
LINK= $(CC) $(LIB_PATH) $(LINK_FLAGS)

# Dependencies
EXEC_DEPENDENCIES = main.o \
                    Manzoni.o \
                    Flow.o \
                    AlessandroFlow.o \
                    ParticlesSet.o \
                    Iterator.o \
                    XMLInputParser.o \
                    ParametersEvolution.o \
                    dataProcessing/procDistrHisto.o \
                    dataProcessing/procDistrProfile.o \
                    dataProcessing/procGraph.o \
                    dataProcessing/procIslands.o \
                    dataProcessing/procTracking.o \
                    utils/signals.o \
		            utils/Logger.o \
                    utils/Singleton.o \
		            utils/Arguments.o
		            
EXEC_DEPENDENCIES_DEBUG = main.od \
                    Manzoni.od \
                    Flow.od \
                    AlessandroFlow.od \
                    ParticlesSet.od \
                    Iterator.od \
                    XMLInputParser.od \
                    ParametersEvolution.od \
                    dataProcessing/procDistrHisto.od \
                    dataProcessing/procDistrProfile.od \
                    dataProcessing/procGraph.od \
                    dataProcessing/procIslands.od \
                    dataProcessing/procTracking.od \
                    utils/signals.od \
		            utils/Logger.od \
                    utils/Singleton.od \
		            utils/Arguments.od
# Rules
#
# # Objects
%.o: %.cc %.h
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

cleanall:
	rm -Rf *.o* */*.o* */*/*/o* *.pdf *.dat $(EXEC) $(EXEC_DEBUG)

# # Count lines
count:
	python ./script/count
