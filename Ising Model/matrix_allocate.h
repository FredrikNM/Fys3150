#ifndef MATRIX_ALLOCATE_H
#define MATRIX_ALLOCATE_H

# include <iostream>
# include <random>
# include <omp.h>

// Setup for regular matrix
void matrix(int**& Mptr, int n_spins);
void matrix_relase(int** Mptr, int n_spins);


// Setup for our try of a faster boundary matrix
void boundary_matrix(void***& Mpointer, int n_spins); // A matrix where the edges are pointer. Described in detailed at the initialize function
void boundary_matrix_relase(void*** Mptr, int n_spins);

// Setup for locks with pointers at the outer indexes
void boundary_locks(void***& count_lock, int n_spins);
void boundary_locks_relase(void*** count_lock, int n_spins);



#endif // MATRIX_ALLOCATE_H