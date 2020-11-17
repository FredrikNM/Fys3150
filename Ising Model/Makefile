# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= mpic++ -fopenmp -g -O3 -std=c++11 
# Here we define the library functions we nee
LIB = #-larmadillo -llapack -lblas -lmpi

# -------------------------------------------------------------------------------- #
# Here we define the name of the executable
PROG = runme
${PROG} :	mainprogram.o 	matrix_allocate.o
			${CPPflags} mainprogram.o 	matrix_allocate.o ${LIB} -o ${PROG}

mainprogram.o	:	mainprogram.cpp 
			${CPPflags}	-c	mainprogram.cpp

matrix_allocate.o :		 matrix_allocate.cpp
					${CPPflags}	-c	matrix_allocate.cpp
