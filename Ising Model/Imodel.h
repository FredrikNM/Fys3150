#ifndef IMODEL_H
#define IMODEL_H

# include <iostream>
# include <iomanip>
# include <string>
# include <random>
# include <fstream>
# include <sys/stat.h> 
# include <omp.h>
# include "matrix_allocate.h"


std::random_device rdev;
std::mt19937_64 randng(rdev());
std::uniform_real_distribution<double> u_cont_d(0.0, 1.0);



inline int get_periodic_index(int index, int n_spins);
void initialize(int **spins, int n_spins, double& E, double& M, int order);
void initialize_boundary_matrix(void*** Mpointer, int n_spins, double& E, double& M, int order);
void initialize_boundary_locks(void*** count_lock, int n_spins);



void MMC(int n_spins, int MC_samples, double Temp, double expectation[5], int order,
		 long long& changes_accepted);


void MMC_boundary_matrix(int n_spins, int MC_samples, double Temp, double* expectation, int order, 
					     long long& changes_accepted, int method);


void MMC_boundary_matrix_with_locks(int n_spins, int MC_samples, double Temp, double expectation[5],
									int order, int cores, int print_every_x_MC_cycle, int probability_dist,
									std::ofstream* value_outfile, std::ofstream* prob_outfile);


void MMC_boundary_matrix_parallel_spin(int n_spins, int MC_samples, double Temp, double expectation[5],
									   int order, int cores, int print_every_x_MC_cycle, int probability_dist,
									   std::ofstream* value_outfile, std::ofstream* prob_outfile);



void get_analytical_solutions(double Temp);

void print_results(double expectation[5], int n_spins, int MC_samples, double temp);


void print_results(double expectation[5], int MC_samples, int n_spins, double n_x_n_inv, double Temp,
                   double inv_Temp, double inv_d_Temp);

void write_to_file(double expectation[5], int MC_samples, int n_spins, double n_x_n_inv, double Temp,
                   double inv_Temp, double inv_d_Temp, int changes_accepted, std::ofstream* outfile);


void print_boundary_lattice(void*** spins, int n_spins);
void print_lattice(int** spins, int n_spins);

// TESTS
void periodic_test();
void check_boundary_matrix();
void rng_initialization();

#endif // IMODEL_H