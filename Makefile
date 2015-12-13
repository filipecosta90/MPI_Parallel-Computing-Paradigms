################################################################################
# Makefile for PCP 
#
# by (CPD-Minho)
################################################################################

SHELL = /bin/sh

MPI = mpi_scatter
PAR = par
CXX = mpicc

BIN = bin
BIN_MPI = pcp_tp2_mpi
BIN_PAR = pcp_tp1_par

CXXFLAGS   = -O3 -Wall -Wextra -std=c++11 -fopenmp -g
LIBS =  #-L/share/apps/mpiP-gcc4.9-openmpi-1.6.3/ -lmpiP -lbfd -liberty -lunwind -lm
#S-lbfd
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
	$(CXX) -M $(CXXFLAGS) $(INCLUDES) $< -o $@

$(BUILD_DIR)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

$(BIN_DIR)/$(BIN_PAR): $(BUILD_DIR)/$(PAR).o $(BUILD_DIR)/$(PAR).d 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(BUILD_DIR)/$(PAR).o $(LIBS)

$(BIN_DIR)/$(BIN_MPI): $(BUILD_DIR)/$(MPI).o $(BUILD_DIR)/$(MPI).d 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(BUILD_DIR)/$(MPI).o $(LIBS)

checkdirs:
	@mkdir -p build 
	@mkdir -p src
	@mkdir -p timing
	@mkdir -p bin

all: checkdirs $(BIN_DIR)/$(BIN_PAR) $(BIN_DIR)/$(BIN_MPI) 

clean:
	rm -f $(BUILD_DIR)/* $(BIN_DIR)/* 	
