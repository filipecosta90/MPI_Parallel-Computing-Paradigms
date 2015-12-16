################################################################################
# Makefile for PCP 
#
# by (CPD-Minho)
################################################################################

SHELL = /bin/sh

MPI = mpi_reduce
MPI_GATHER = mpi_gather
MPICXX = mpic++
CXX = mpic++

BIN = bin
BIN_MPI = eth_reduce_by_node
BIN_MPI_GATHER = eth_gather

CXXFLAGS   = -O3 -Wall -Wextra -std=c++11 -fopenmp  
LIBS = -L/share/apps/mpiP-3.4.1/ -lmpiP -lbfd -liberty -lunwind -lm

SRC_DIR = src
BIN_DIR = bin
BUILD_DIR = build
SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(patsubst src/*.cpp,build/*.o,$(SRC))
DEPS = $(patsubst build/*.o,build/*.d,$(OBJ))

vpath %.cpp $(SRC_DIR)


################################################################################
# Rules
################################################################################

.DEFAULT_GOAL = all

$(BUILD_DIR)/%.d: %.cpp
	$(CXX) -M $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LIBS)

$(BUILD_DIR)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LIBS)

$(BIN_DIR)/$(BIN_MPI): $(BUILD_DIR)/$(MPI).o $(BUILD_DIR)/$(MPI).d 
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(BUILD_DIR)/$(MPI).o $(LIBS)

$(BIN_DIR)/$(BIN_MPI_GATHER): $(BUILD_DIR)/$(MPI_GATHER).o $(BUILD_DIR)/$(MPI_GATHER).d 
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(BUILD_DIR)/$(MPI_GATHER).o $(LIBS)

checkdirs:
	@mkdir -p build 
	@mkdir -p src
	@mkdir -p timing
	@mkdir -p bin

all: checkdirs $(BIN_DIR)/$(BIN_MPI) $(BIN_DIR)/$(BIN_MPI_GATHER)

clean:
	rm -f $(BUILD_DIR)/* $(BIN_DIR)/* 	
