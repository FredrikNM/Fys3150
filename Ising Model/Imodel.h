#ifndef IMODEL_H
#define IMODEL_H

# include <iostream>
# include <iomanip>
# include <string>
# include <random>
# include <fstream>
// # include <cmath>
// # include <ctime>


// inline int get_periodic_index(int index, int n_spins);
int get_periodic_index(int index, int n_spins);
void initialize(int **spins, int n_spins, double& E, double& M, int order);
void initialize_boundary_matrix(void*** Mpointer, int n_spins, double& E, double& M, int order);
void initialize_boundary_locks(void*** count_lock, int n_spins);


void MMC(int n_spins, int MC_samples, double Temp, double expectation[5], int order);
void MMC_boundary_matrix(int n_spins, int MC_samples, double Temp, double expectation[5], int order);
void MMC_boundary_matrix_with_locks(int n_spins, int MC_samples, double Temp, double expectation[5], int order, int cores);
void MMC_boundary_matrix_parallel_spin(int n_spins, int MC_samples, double Temp, double expectation[5], int order, int cores);

void get_analytical_solutions(double Temp, int N);
void print_results(double expectation[5], int n_spins, int MC_samples, double temp);
void write_to_file(double expectation[5], int MC_samples, int n_spins, std::string output_filename);
void print_boundary_lattice(void*** spins, int n_spins);
void print_lattice(int** spins, int n_spins);


#endif // IMODEL_H