################################################################################
# Makefile for PCP -- team laptop 
#
# by (CPD-Minho)
################################################################################

SHELL = /bin/sh

MPI_OPENMP = mpi_openmp
MPICXX = mpic++
CXX = mpic++

BIN = bin
BIN_MPI_OPENMP = eth_mpi_openmp

CXXFLAGS   = -O3 -Wall -Wextra -std=c++11 -fopenmp  
LIBS = 

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

$(BIN_DIR)/$(BIN_MPI_OPENMP): $(BUILD_DIR)/$(MPI_OPENMP).o $(BUILD_DIR)/$(MPI_OPENMP).d 
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(BUILD_DIR)/$(MPI_OPENMP).o $(LIBS)

checkdirs:
	@mkdir -p build 
	@mkdir -p src
	@mkdir -p timing
	@mkdir -p bin

all: checkdirs $(BIN_DIR)/$(BIN_MPI_OPENMP)

clean:
	rm -f $(BUILD_DIR)/* $(BIN_DIR)/* 	
