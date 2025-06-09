SRC_DIR=src
INPUTS=./files
# SRC=$(SRC_DIR)/BoundingBox.cpp $(SRC_DIR)/ex04_funcs.cpp $(SRC_DIR)/main.cpp
SRC=$(SRC_DIR)/BoundingBox.cpp $(SRC_DIR)/main.cpp
BIN=lbm
ARGS=$(INPUTS)/params_Re40.dat
# ARGS=$(INPUTS)/params_Re100_wing.dat
# ARGS=$(INPUTS)/params_Re500.dat

CC=g++
ERRFLAGS=-Wall -pedantic
DEBUG=-g -Og
STD=c++17
CFLAGS=$(ERRFLAGS) -std=$(STD) -O3


default: $(SRC)
	$(CC) $(CFLAGS) -o $(BIN) $(SRC)

debug: $(SRC)
	$(CC) $(CFLAGS) $(DEBUG) -o $(BIN) $(SRC)

run:
	./$(BIN) $(ARGS)

clean:
	rm -rf $(BIN)
	rm -rf $(BIN).exe
	rm -rf *.vtk

vtkclean:
	rm -rf *.vtk

memcheck:
	valgrind --leak-check=full --track-origins=yes -s ./$(BIN) $(ARGS)
