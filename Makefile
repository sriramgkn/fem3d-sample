#Makefile for compiling FEM3D project

TARGET = fem3d

#Compilation flags can be chosen differently for optimization or debug mode
CFLAGS = -Wall -O3 -std=c++17 -I/usr/include/python3.8 -lpython3.8
#-g -pg

#here we include the locations of necessary libraries
INCLUDE = ../../libs/eigen
INCLUDE2 = -lcomplex_bessel -lgfortran

all: fem3d.cpp
	g++ -o $(TARGET) $(TARGET).cpp -I $(INCLUDE) $(CFLAGS) $(INCLUDE2)

clean:
	rm $(TARGET)
