
# include "matrix_allocate.h"


// Uniform random number from 0 -> 1

std::random_device rdev2;      // Random device
std::mt19937_64 randng2(rdev2());  // Choosing generator
std::uniform_real_distribution<double> u_cont_d2;  // Setting up a uniform continuous distribution


// Setup for spin matrix with double pointers
void matrix(int**& Mptr, int n_spins)
{
  Mptr = (int **) malloc(n_spins*sizeof(int*));

  for (int x = 0; x < n_spins; x++){
    Mptr[x] = (int *) malloc(n_spins*sizeof(int));
  }
}


// Free the memory used by matrix
void matrix_relase(int** Mptr, int n_spins)
{
  for (int x = 0; x < n_spins; x++)
  {
    delete [] Mptr[x];
  }
  delete [] Mptr;
}



// Setup for spin matrix with boundary as pointer to the n-2 x n-2
void boundary_matrix(void***& Mpointer, int n_spins)
{
  Mpointer = (void ***) malloc((n_spins+2) * sizeof(void *));
  for (int i=0; i<n_spins+2; i++)
  {
    Mpointer[i] = (void **) malloc((n_spins+2) * sizeof(void *));
    for (int j=0; j<n_spins+2; j++)
    {
      if (i == 0 || i == (n_spins+1) || j == 0 || j == (n_spins+1))
      {
        Mpointer[i][j] = (void *) malloc(sizeof(void *));
      }
      else
      {
        Mpointer[i][j] = (int *) malloc(sizeof(int *));
      }
    }
  }
}


// Free the memory used by matrix
void boundary_matrix_relase(void*** Mptr, int n)
{
  for (int x = 0; x < n; x++)
  {
    delete [] Mptr[x];
  }
  delete [] Mptr;
}





void boundary_locks(void***& count_lock, int n_spins)
{
  // void*** count_lock;

  count_lock = (void ***) malloc((n_spins+2) * sizeof(void *));
  for (int i=0; i<n_spins+2; i++)
  {
    count_lock[i] = (void **) malloc((n_spins+2) * sizeof(void *));
    for (int j=0; j<n_spins+2; j++)
    {
      if (i == 0 || i == (n_spins+1) || j == 0 || j == (n_spins+1))
      {
        count_lock[i][j] = (void *) malloc(sizeof(omp_lock_t *));
      }
      else
      {
        count_lock[i][j] = (omp_lock_t *) malloc(sizeof(omp_lock_t ));
      }
    }
  }
}


// Free the memory used by matrix
void boundary_locks_relase(void*** count_lock, int n_spins)
{

  for (int i = 0; i < n_spins+2; i++)
  {
  	for (int j = 0; j < n_spins+2; j++)
  	{
  	  omp_destroy_lock(&*((omp_lock_t *)(count_lock[i][j])));
  	}
  }
}




