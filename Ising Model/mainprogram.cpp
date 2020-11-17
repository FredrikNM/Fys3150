

# include <omp.h>
# include "matrix_allocate.h"
# include "Imodel.h"

// omp_lock_t lock[n_x_n];

int main(int argc, char *argv[])
{


  //////////// TESTS //////////////

  // std::cout << " Matrix Initialize \n";
  // for(int x = 0; x < n_spins; x++) 
  // {
  //   for(int y = 0; y < n_spins; y++)
  //   {
  //     std::cout << " " << spins[x][y] <<  " ";
  //   }
  // std::cout << std::endl;
  // }
  // std::cout << get_periodic_index(-1,n_spins) << std::endl;

  /////////////////////////////////



  double T_0, T_n, dN_T;
  // int MC_samples, n_temps, configuration, time_pointer_matrix, order, cores, method_to_use, n_spin_size;
  int MC_samples, n_temps, configuration, time_pointer_matrix, order, cores, n_methods_size, n_spin_size;
  int * n_spins;
  int * n_methods;
  std::string output_file;

  // Read commandline arguments if it is provided
  if (argc > 1) 
  {
    if (argc > 1 && argc <= 7)
    {
      std::cout << "Bad Usage: " << argv[0] << 
        " read in output file, N in NxN spins, MC cycles, initial"
        " and final temperature, tempurate step, order or not (0 for "
        " not order, 1 for ordered or 2 for both) and number of cores like this : "
        " ./runme filename.txt 2 100000 2 2.4 0.05 1 4. You have a maximum"
        " of " << omp_get_max_threads() << " Cores available" << std::endl;
      exit(1);
    }
    // else
    // {
    //   output_file = argv[1]; n_spins = atoi(argv[2]); MC_samples = atoi(argv[3]);
    //   T_0 = atof(argv[4]); T_n = atof(argv[5]); dN_T = atof(argv[6]); 
    //   order = atoi(argv[8]); cores = atoi(argv[8]);
    //   int n_temps = (T_n - T_0)/ dN_T; if (n_temps == 0) n_temps = 1; // Setting number of temperatures to run at least to 1
    // }
  }

  // If no commandline arguments are given, ask for them
  else
  {
    std::cout << "Press 1 and Enter if you want to set your own configurations, "
    "or press 0 and Enter if you want the system to run with 2x2 spins, "
    "1 million MC cycles , temprature = 1, all spins start value = 1 "; std::cin >> configuration;

    if (configuration == 1) 
    {

      std::cout << "Type min temp then hit Enter: "; std::cin >> T_0;
      std::cout << "Type max temp then hit Enter: "; std::cin >> T_n;
      std::cout << "Please type stepsize between the temperatures, and hit Enter: "; std::cin >> dN_T;
      std::cout << "How many different lattice size do you want to simulate : "; std::cin >> n_spin_size;
      n_spins = new int [n_spin_size];
      for (int i = 0; i < n_spin_size; i ++){std::cout << "Please enter 'N' in NxN spins for "
        "number "+std::to_string(i+1)+" out of "+std::to_string(n_spin_size)+" : "; std::cin >> n_spins[i];}
      std::cout << "Should the spins be ordered ? 0 and Enter for No, 1 and Enter for Yes, 2 and Enter for both : "; std::cin >> order;
      int n_temps = (T_n - T_0)/ dN_T; if (n_temps == 0) n_temps = 1; // Setting number of temperatures to run at least to 1
      std::cout << "Thats gonna be " << n_temps << " different temperatures. "
      "Enter the number of MC cycles: "; std::cin >> MC_samples;
      std::cout << "How many cores you want to divide the labour on when it is possible to parallelize ?"
      " You have " << omp_get_max_threads() << " available : "; std::cin >> cores;
      std::cout << " \n\n First enter how many methods you want out of the ones listed below, and then enter method number"
      "\n \n 0 -> no parallelization. \n \n 1 -> correct way of parallelization where each core lock the "
      "spin and neighbouring spins while checking to flip it or not."
      "\n \n 2 -> parallelization where each core dont care about if the neighbouring "
      "spin is beeing worked on by another core. Can cause bias in the results."
      "\n \n 3 -> splitting the numbers of MC cycles on cores available, where each "
      "start their own system and do a total of ((MC cycles)/(cores used)) MC cycle each "
      "then averaging over all runs. (Similar to model averaging in statistics). \n \n "; std::cin >> n_methods_size;
      n_methods = new int [n_methods_size];
      for (int i = 0; i < n_methods_size; i ++){std::cout << "Enter method number "+std::to_string(i+1)+" out"
      " of "+std::to_string(n_methods_size)+" : " ; std::cin >> n_methods[i];}


    }
    else
    {
      // Temperature variables (min, max, steps, n_temps = how many of the steps to take)
      T_0 = 1; T_n = 2, dN_T = 1; n_temps = 1;
      // N in the NxN lattice, (ordered spins = 1 means all spins starts value = 1, order = 0 gives random spins 1 or -1 for each spin)
      n_spins = new int [1]; n_spins[0] = 2; n_spin_size = 1; order = 0; 
      // MC_samples, and cores to use
      MC_samples = 1000000; cores = 1; //cores = omp_get_max_threads();
      // This is the different methods of parallelization or not. 
      // method_to_use = 0 -> no parallelization.
      // method_to_use = 1 -> correct way of parallelization where each core lock the 
      //                      spin and neighbouring spins while checking to flip it or not.
      // method_to_use = 2 -> parallelization where each core dont care about if the neighbouring 
      //                      spin is beeing worked on by another core. Can cause bias in the results.
      // method_to_use = 3 -> splitting the numbers of MC cycles on cores available, where each 
      //                      start their own system and do a total of ((MC cycles)/(cores used)) MC cycle each
      //                      then averaging over all runs. (Similar to model averaging in statistics).
      n_methods = new int [1]; n_methods[0] = 0; n_methods_size = 1;
    }
  }


  // Set number of threads
  #define NUM_THREADS cores

  // Just a precaution in case n_temps is set to large
  if (n_temps > (T_n - T_0)/dN_T ){n_temps = int (T_n - T_0)/dN_T ;}
  if (n_temps == 0) n_temps = 1; // Setting number of temperatures to run at least to 1
  // Setting up tempratures to run
  double T[n_temps];
  if (n_temps >= 1){ for (int i = 0; i < n_temps; i++) T[i] = T_0 + dN_T*i; }
  else{ T[1] = {T_0}; }




  // double expectation[5] = {0,0,0,0,0};
  // double t0 = omp_get_wtime();
  // // MMC_boundary_matrix_with_locks(n_spins[0], MC_samples, T[0], expectation, order, cores);
  // MMC_boundary_matrix(n_spins[0], MC_samples, T[0], expectation, order);
  // double time = omp_get_wtime()-t0;
  // printf("(%f)\n", time);


  // int total_spins = n_spins[0]*n_spins[0];

  // double norm = 1.0/((double) (MC_samples));
  // double AllSpins = 1.0/((double) total_spins);
  // double x1 = expectation[0]*norm;
  // double x2 = expectation[1]*norm;
  // double x3 = expectation[2]*norm;
  // double x4 = expectation[3]*norm;
  // double x5 = expectation[4]*norm;

  // double HeatCapacity = (x2- x1*x1)*AllSpins/T[0]/T[0];
  // double chi = (x4 - x5*x5)*AllSpins/T[0];

  // std::cout << std::setprecision(8) << "    E : " << x1 * AllSpins;
  // std::cout << std::setprecision(8) << "    E2 : " << x2 * AllSpins;
  // std::cout << std::setprecision(8) << "    M : " << x5 * AllSpins;
  // std::cout << std::setprecision(8) << "    M^2 : " << x4 * AllSpins;
  // std::cout << std::setprecision(8) << "    C_v : " << HeatCapacity;
  // std::cout << std::setprecision(8) << "    X : " << chi << std::endl;
  // get_analytical_solutions(T[0], n_spins[0]);




  for (int t = 0; t < n_temps; t++)
  {
    for (int n = 0; n < n_spin_size; n++)
    {
      for(int m = 0; m < n_methods_size; m++)
      {
        // Avoiding problems occuring with parallization with locks when same core is
        // trying to lock the same index(i,j) twice, which occure in a 2x2 lattice
        int temp_method = n_methods[m];
        if (n_methods[m] == 1 && n_spins[n] < 3){temp_method = 0;}

        if (temp_method == 0)
        {
        
          double t0 = omp_get_wtime();  // Start timing of method
          double expectation[5] = {0,0,0,0,0};  // Setting up array to save E, E^2, M, M^2, and |M|
          MMC_boundary_matrix(n_spins[n], MC_samples, T[t], expectation, order);
          double time = omp_get_wtime()-t0; // Save time used

          // Print results to commandline
          std::cout << "\n Run of method "+std::to_string(temp_method)+" with lattice size"
          " "+std::to_string(n_spins[n])+"x"+std::to_string(n_spins[n])+" and"
          " "+std::to_string(MC_samples)+" MC samples with temperature "+std::to_string(T[t])+" "
          "took "+std::to_string(time)+" seconds \n" << std::endl;
          print_results(expectation, n_spins[n], MC_samples, T[t]);

        } // end method_to_use == 0
        else if (temp_method == 1)
        {

          double t0 = omp_get_wtime();  // Start timing of method
          double expectation[5] = {0,0,0,0,0};  // Setting up array to save E, E^2, M, M^2, and |M|
          MMC_boundary_matrix_with_locks(n_spins[n], MC_samples, T[t], expectation, order, cores);
          double time = omp_get_wtime()-t0; // Save time used

          // Print results to commandline
          std::cout << "\n Run of method "+std::to_string(temp_method)+" with lattice size"
          " "+std::to_string(n_spins[n])+"x"+std::to_string(n_spins[n])+" and"
          " "+std::to_string(MC_samples)+" MC samples with temperature "+std::to_string(T[t])+" "
          "took "+std::to_string(time)+" seconds \n" << std::endl;
          print_results(expectation, n_spins[n], MC_samples, T[t]);

        } // end method_to_use == 1
        else if (temp_method == 2)
        {

          double t0 = omp_get_wtime();  // Start timing of method
          double expectation[5] = {0,0,0,0,0};  // Setting up array to save E, E^2, M, M^2, and |M|
          MMC_boundary_matrix_parallel_spin(n_spins[n], MC_samples, T[t], expectation, order, cores);
          double time = omp_get_wtime()-t0; // Save time used

          // Print results to commandline
          std::cout << "\n Run of method "+std::to_string(temp_method)+" with lattice size"
          " "+std::to_string(n_spins[n])+"x"+std::to_string(n_spins[n])+" and"
          " "+std::to_string(MC_samples)+" MC samples with temperature "+std::to_string(T[t])+" "
          "took "+std::to_string(time)+" seconds \n" << std::endl;
          print_results(expectation, n_spins[n], MC_samples, T[t]);

        } // end method_to_use == 2
        else if (temp_method == 3)
        {
          double t0 = omp_get_wtime();  // Start timing of method
          double expectation[5] = {0,0,0,0,0};  // Setting up array to save E, E^2, M, M^2, and |M|
          #pragma omp parallel for reduction(+:expectation)
          for (int c = 0; c < NUM_THREADS; c++)
          {
            MMC_boundary_matrix(n_spins[n], MC_samples/NUM_THREADS, T[t], expectation, order);
          }
          double time = omp_get_wtime()-t0; // Save time used

          // Print results to commandline
          std::cout << "\n Run of method "+std::to_string(temp_method)+" with lattice size"
          " "+std::to_string(n_spins[n])+"x"+std::to_string(n_spins[n])+" and"
          " "+std::to_string(MC_samples)+" MC samples with temperature "+std::to_string(T[t])+" "
          "took "+std::to_string(time)+" seconds \n" << std::endl;
          print_results(expectation, n_spins[n], MC_samples, T[t]);

        } // end method_to_use == 3


      } // end n_methods
    } // end n_spins_size
  } // end n_temps


  return 0;
}
















// Metropolis Monte Carlo
void MMC_boundary_matrix_with_locks(int n_spins, int MC_samples, double Temp, double expectation[5], int order, int cores)
{
  std::random_device rd1;
  std::mt19937_64 rng1(rd1());
  std::uniform_real_distribution<double> uniform1(0.0, 1.0);

  double E = 0, M = 0;
  int changes_accepted = 0;
  void*** spins;

  boundary_matrix(spins, n_spins);
  initialize_boundary_matrix(spins, n_spins, E, M, order);
  print_boundary_lattice(spins, n_spins);


  int all_dE[5] = {-8, -4, 0, 4, 8};
  double exp_dE[5];

  for (int k=0; k<5 ; k++) exp_dE[k] = exp(-all_dE[k]/Temp);

  int rand_i, rand_j;
  int n_x_n = n_spins*n_spins/cores;
  int total_spins = n_spins*n_spins;
  double n_x_n_inv = 1./n_x_n;
  int dE = E;
  int m, n;




  void*** count_lock;
  boundary_locks(count_lock, n_spins);
  initialize_boundary_locks(count_lock, n_spins);



  // std::ofstream outfile;
  // outfile.open("energy_count.txt");
  // outfile << std::setw(15) << "MC";
  // outfile << std::setw(15) << "E";
  // outfile << std::setw(15) << "E^2";
  // outfile << std::setw(15) << "|M|";
  // outfile << std::setw(15) << "M^2";
  // outfile << std::setw(15) << "C_v";
  // outfile << std::setw(15) << "X" << std::endl;



  #define NUM_THREADS cores
  #pragma omp parallel num_threads(NUM_THREADS) shared(m) private(n)
  {
    for (m = 1; m<=MC_samples; m++)
    {
      n = 0;
      #pragma omp critical

      while (n < n_x_n)
      {

        int rand_i = (int) (uniform1(rng1)*(double)n_spins) + 1;
        int rand_j = (int) (uniform1(rng1)*(double)n_spins) + 1;

        if (omp_test_lock(&*((omp_lock_t *)(count_lock[rand_i][rand_j]))) && omp_test_lock(&*((omp_lock_t *)(count_lock[rand_i+1][rand_j])))
          && omp_test_lock(&*((omp_lock_t *)(count_lock[rand_i-1][rand_j]))) && omp_test_lock(&*((omp_lock_t *)(count_lock[rand_i][rand_j-1]))))
        {

        dE = 2* *((int*)(spins[rand_i][rand_j])) *
                (*((int*)(spins[rand_i-1][rand_j])) + *((int*)(spins[rand_i+1][rand_j])) + 
                 *((int*)(spins[rand_i][rand_j-1])) + *((int*)(spins[rand_i][rand_j+1])));

        if (dE <= 0)
        {
          *((int*)(spins[rand_i][rand_j])) *= -1; 
          M += (double) (2* *((int*)(spins[rand_i][rand_j])));
          E += (double) dE;
          changes_accepted += 1;
        }
        else
        {
          double current_exp_dE = 0.0;

          for (int i = 0; i < 5; i++)
          {
            if (dE == all_dE[i])
            {
                current_exp_dE = exp_dE[i];
                break;
            }
          }

          if (uniform1(rng1) <= current_exp_dE)
          {
              *((int*)(spins[rand_i][rand_j])) *= -1; 
              M += (double) (2* *((int*)(spins[rand_i][rand_j])));
              E += (double) dE;
              changes_accepted += 1;
          }
        }

        omp_unset_lock(&*((omp_lock_t *)(count_lock[rand_i][rand_j])));
        omp_unset_lock(&*((omp_lock_t *)(count_lock[rand_i+1][rand_j])));
        omp_unset_lock(&*((omp_lock_t *)(count_lock[rand_i-1][rand_j])));
        omp_unset_lock(&*((omp_lock_t *)(count_lock[rand_i][rand_j+1])));
        omp_unset_lock(&*((omp_lock_t *)(count_lock[rand_i][rand_j-1])));

        n = n + 1;

        }
      }

      expectation[0] += E; expectation[1] += E*E;
      expectation[2] += M; expectation[3] += M*M; 
      expectation[4] += fabs(M);

    }
  }  // Slutt prallel 
  boundary_matrix_relase(spins, n_spins);
  boundary_locks_relase(count_lock, n_spins);
  // std::cout << "\n acc states : " << changes_accepted << " divided on MC_samples and n_spins : " << (double) changes_accepted / ((double) MC_samples*n_x_n ) << std::endl;
  // outfile.close();

}















// Metropolis Monte Carlo
void MMC_boundary_matrix_parallel_spin(int n_spins, int MC_samples, double Temp, double expectation[5], int order, int cores)
{
  std::random_device rd1;
  std::mt19937_64 rng1(rd1());
  std::uniform_real_distribution<double> uniform1(0.0, 1.0);

  double E = 0, M = 0;
  int changes_accepted = 0;
  void*** spins;

  boundary_matrix(spins, n_spins);
  initialize_boundary_matrix(spins, n_spins, E, M, order);
  // print_boundary_lattice(spins, n_spins);

  int all_dE[5] = {-8, -4, 0, 4, 8};
  double exp_dE[5];

  for (int k=0; k<5 ; k++) exp_dE[k] = exp(-all_dE[k]/Temp);

  int rand_i, rand_j;
  int n_x_n = n_spins*n_spins/cores;
  int total_spins = n_spins*n_spins;
  double n_x_n_inv = 1./n_x_n;
  int dE = E;
  int m, n;


  #define NUM_THREADS cores
  #pragma omp parallel num_threads(NUM_THREADS) shared(m) private(n)
  {
    for (m = 1; m<=MC_samples; m++)
    {
      n = 0;
      #pragma omp critical

      while (n < n_x_n)
      {

        int rand_i = (int) (uniform1(rng1)*(double)n_spins) + 1;
        int rand_j = (int) (uniform1(rng1)*(double)n_spins) + 1;

        dE = 2* *((int*)(spins[rand_i][rand_j])) *
                (*((int*)(spins[rand_i-1][rand_j])) + *((int*)(spins[rand_i+1][rand_j])) + 
                 *((int*)(spins[rand_i][rand_j-1])) + *((int*)(spins[rand_i][rand_j+1])));

        if (dE <= 0)
        {
          *((int*)(spins[rand_i][rand_j])) *= -1; 
          M += (double) (2* *((int*)(spins[rand_i][rand_j])));
          E += (double) dE;
          changes_accepted += 1;
        }
        else
        {
          double current_exp_dE = 0.0;

          for (int i = 0; i < 5; i++)
          {
            if (dE == all_dE[i])
            {
                current_exp_dE = exp_dE[i];
                break;
            }
          }

          if (uniform1(rng1) <= current_exp_dE)
          {
              *((int*)(spins[rand_i][rand_j])) *= -1; 
              M += (double) (2* *((int*)(spins[rand_i][rand_j])));
              E += (double) dE;
              changes_accepted += 1;
          }
        }

        n = n + 1;
      }

      expectation[0] += E; expectation[1] += E*E;
      expectation[2] += M; expectation[3] += M*M; 
      expectation[4] += fabs(M);

    }
  }  // Slutt prallel 
  boundary_matrix_relase(spins, n_spins);
  // std::cout << "\n acc states : " << changes_accepted << " divided on MC_samples and n_spins : " << (double) changes_accepted / ((double) MC_samples*n_x_n ) << std::endl;
  // outfile.close();

}

















// Metropolis Monte Carlo
void MMC_boundary_matrix(int n_spins, int MC_samples, double Temp, double expectation[5], int order)
{
  std::random_device rd1;
  std::mt19937_64 rng1(rd1());
  std::uniform_real_distribution<double> uniform1(0.0, 1.0);

  double E = 0, M = 0;
  int changes_accepted = 0;
  void*** spins;

  boundary_matrix(spins, n_spins);
  initialize_boundary_matrix(spins, n_spins, E, M, order);
  print_boundary_lattice(spins, n_spins);


  int all_dE[5] = {-8, -4, 0, 4, 8};
  double exp_dE[5];

  for (int k=0; k<5 ; k++) exp_dE[k] = exp(-all_dE[k]/Temp);

  int rand_i, rand_j;
  int n_x_n = n_spins*n_spins;
  double n_x_n_inv = 1./n_x_n;
  int dE = E;

  for (int m = 1; m<=MC_samples; m++)
  {
    for (int n = 0; n<n_x_n; n++)
    {

      int rand_i = (int) (uniform1(rng1)*(double)n_spins) + 1;
      int rand_j = (int) (uniform1(rng1)*(double)n_spins) + 1;

      dE = 2* *((int*)(spins[rand_i][rand_j])) *
              (*((int*)(spins[rand_i-1][rand_j])) + *((int*)(spins[rand_i+1][rand_j])) + 
               *((int*)(spins[rand_i][rand_j-1])) + *((int*)(spins[rand_i][rand_j+1])));

      if (dE <= 0)
      {
        *((int*)(spins[rand_i][rand_j])) *= -1; 
        M += (double) (2* *((int*)(spins[rand_i][rand_j])));
        E += (double) dE;
        changes_accepted += 1;
      }
      else
      {
        double current_exp_dE = 0.0;

        for (int i = 0; i < 5; i++)
        {
          if (dE == all_dE[i])
          {
              current_exp_dE = exp_dE[i];
              break;
          }
        }

        if (uniform1(rng1) <= current_exp_dE)
        {
            *((int*)(spins[rand_i][rand_j])) *= -1; 
            M += (double) (2* *((int*)(spins[rand_i][rand_j])));
            E += (double) dE;
            changes_accepted += 1;
        }
      }
    }
    expectation[0] += E; expectation[1] += E*E;
    expectation[2] += M; expectation[3] += M*M; 
    expectation[4] += fabs(M);

  }
  boundary_matrix_relase(spins, n_spins);

}
















// Metropolis Monte Carlo
void MMC(int n_spins, int MC_samples, double Temp, double expectation[5], int order)
{


  std::random_device rd2;
  std::mt19937_64 rng2(rd2());
  std::uniform_real_distribution<double> uniform(0.0, 1.0);


  double E, M;
  long long changes_accepted = 0;
  int** spins;
  E = 0;
  M = 0;

  matrix(spins, n_spins);
  initialize(spins, n_spins, E, M, order);


  int all_dE[5] = {-8, -4, 0, 4, 8};
  double exp_dE[5];

  for (int k=0; k<5 ; k++) exp_dE[k] = exp(-all_dE[k]/Temp);

  int dE, rand_i, rand_j, n_x_n;
  n_x_n = n_spins*n_spins;
  dE = E;

  for (int m = 0; m<MC_samples; m++)
  {
    for (int j = 0; j<n_x_n; j++)
    {

      int rand_i = (int) (uniform(rng2)*(double)n_spins);  // Finding random indexes
      int rand_j = (int) (uniform(rng2)*(double)n_spins);

      dE = 2*spins[rand_i][rand_j] *
            (spins[get_periodic_index(rand_i-1,n_spins)][rand_j] + spins[get_periodic_index(rand_i+1,n_spins)][rand_j] + 
             spins[rand_i][get_periodic_index(rand_j-1,n_spins)] + spins[rand_i][get_periodic_index(rand_j+1,n_spins)]); 

      if (dE <= 0)
      {

        spins[rand_i][rand_j] *= -1; 
        M += (double) (2 * spins[rand_i][rand_j]);
        E += (double) dE;
        changes_accepted += 1;

      }

      else
      {

        double current_exp_dE = 0.0;
        for (int i = 0; i < 5; i++)
        {
          if (dE == all_dE[i])
          {
            current_exp_dE = exp_dE[i];
            break;
          }
        }

        if (uniform(rng2) <= current_exp_dE)
        {
          spins[rand_i][rand_j] *= -1; 
          M += (double) (2 * spins[rand_i][rand_j]);
          E += (double) dE;
          changes_accepted += 1;
        }

      }
    }
    expectation[0] += E; expectation[1] += E*E;
    expectation[2] += M; expectation[3] += M*M; 
    expectation[4] += fabs(M);
  }

  matrix_relase(spins, n_spins);
}




// If-else function for periodic boundary conditions
inline int get_periodic_index(int index, int n_spins)
{
  return  (index < 0) ? n_spins - 1: (index >= n_spins) ? 0 : index;
}




// Initialize the energy, spin matrix and magnetization
void initialize(int **spins, int n_spins, double& E, double& M, int order)
{
  std::random_device rd3;
  std::mt19937_64 randng(rd3());
  std::uniform_real_distribution<double> u_cont_d(0.0, 1.0);

  for(int x=0; x < n_spins; x++) 
  {
    for (int y=0; y < n_spins; y++)
    {

      // If temperature below 1.5 we set all the spins to 1, if not set them randomly to 1 or -1
      spins[x][y] = (order == 1) ? 1 : (u_cont_d(randng) > 0.5) ? 1 : -1;
      // Magnetization
      M += (double) spins[x][y];

    }
  }

  // Energy
  for(int x=0; x < n_spins; x++)
  {
    for (int y=0; y < n_spins; y++)
    {
      E -= spins[x][y]*
      (spins[get_periodic_index(x-1,n_spins)][y] + spins[x][get_periodic_index(y-1,n_spins)]); 
    }
  }

}






void initialize_boundary_matrix(void*** Mpointer, int n_spins, double& E, double& M, int order)
{

  std::random_device rd4;
  std::mt19937_64 randng(rd4());
  std::uniform_real_distribution<double> u_cont_d(0.0, 1.0);


  // For loops for setting up the matrix
  for (int i=0; i<n_spins+2; i++)
  {
    for (int j=0; j<n_spins+2; j++)
    {

      // Setting the values in the square inside
      if (i>0 && i<n_spins+1 && j>0 && j<n_spins+1 )
      {
      // If temperature below 1.5 we set all the spins to 1, if not set them randomly to 1 or -1
        *((int*)(Mpointer[i][j])) = (order == 1) ? 1 : (u_cont_d(randng) > 0.5) ? 1 : -1;
        M +=  (double) *((int *)(Mpointer[i][j])) ;
      }
      // Creating the pointers
      else if (i == 0 && j != 0 && j != n_spins+1)    // != n_spins+1 is used to avoid setting the corners
      {
        Mpointer[i][j] = &*((int *)(Mpointer[n_spins][j]));   // Syntax is horrible to watch but this might be faster
      }
      else if (i == n_spins+1 && j != 0 && j != n_spins+1)
      {
        Mpointer[i][j] = &*((int *)(Mpointer[1][j])); // & is memory address of *, pointer which is void, so since it is void
      }                         // we use (), insde the () we have the void function, which is 
      else if (j == 0 && i != 0 && i != n_spins+1)    // (int*), int pointer, to matrix element [i][j], Mpointer[1][j].
      {
        Mpointer[i][j] = &*((int *)(Mpointer[i][n_spins]));
      }
      else if (j == n_spins+1 && i != 0 && i != n_spins+1)
      {
        Mpointer[i][j] = &*((int *)(Mpointer[i][1]));
      }

    }
  }  // End of loops for setting up the matrix

  // Energy
  for(int x=1; x <=n_spins; x++)
  {
    for (int y=1; y <=n_spins; y++)
    {
      E -= *((int*)(Mpointer[x][y]))*
      (*((int*)(Mpointer[x-1][y])) + *((int*)(Mpointer[x][y-1])));
    }
  }

}







void initialize_boundary_locks(void*** count_lock, int n_spins)
{
  for(int i=0; i<n_spins+2; i++)
    {
    for(int j=0; j<n_spins+2; j++)
    {
      if (i>0 && i<n_spins+1 && j>0 && j<n_spins+1 )
      {
        omp_init_lock(&*(omp_lock_t*) (count_lock[i][j]));
      }
      //////// Creating the pointers ////////
      else if (i == 0 && j != 0 && j != n_spins+1)
      {
        count_lock[i][j] = &*((omp_lock_t *)(count_lock[n_spins][j]));
      }
      else if (i == n_spins+1 && j != 0 && j != n_spins+1)
      {
        count_lock[i][j] = &*((omp_lock_t *)(count_lock[1][j])); // & is memory address of *, pointer which is void, so since it is void
      }                         // we use (), insde the () we have the void function, which is 
      else if (j == 0 && i != 0 && i != n_spins+1)    // (omp_lock_t*), omp_lock_t pointer, to matrix element [i][j], count_lock[1][j].
      {
        count_lock[i][j] = &*((omp_lock_t *)(count_lock[i][n_spins]));
      }
      else if (j == n_spins+1 && i != 0 && i != n_spins+1)
      {
        count_lock[i][j] = &*((omp_lock_t *)(count_lock[i][1]));
      }
    }
  }
}







void print_results(double expectation[5], int n_spins, int MC_samples, double temp)
{
  int total_spins = n_spins*n_spins;

  double norm = 1.0/((double) (MC_samples));
  double AllSpins = 1.0/((double) total_spins);
  double x1 = expectation[0]*norm;
  double x2 = expectation[1]*norm;
  double x3 = expectation[2]*norm;
  double x4 = expectation[3]*norm;
  double x5 = expectation[4]*norm;

  double HeatCapacity = (x2- x1*x1)*AllSpins/temp/temp;
  double chi = (x4 - x5*x5)*AllSpins/temp;

  std::cout << std::setprecision(8) << "    E : " << x1 * AllSpins;
  std::cout << std::setprecision(8) << "    E2 : " << x2 * AllSpins;
  std::cout << std::setprecision(8) << "    M : " << x5 * AllSpins;
  std::cout << std::setprecision(8) << "    M^2 : " << x4 * AllSpins;
  std::cout << std::setprecision(8) << "    C_v : " << HeatCapacity;
  std::cout << std::setprecision(8) << "    X : " << chi << std::endl;
  get_analytical_solutions(temp, n_spins);

}

void print_boundary_lattice(void*** spins, int n_spins)
{

  ////////////////////////////////////////
  /////////// PRINT OUT LATTICE ////////// 
  std::cout << std::endl;
  for(int x = 0; x < n_spins+2; x++) 
  {
    for(int y = 0; y < n_spins+2; y++)
    {
      std::cout << " " << *((int*)(spins[x][y])) <<  " ";
    }
  std::cout << std::endl;
  }
  std::cout << std::endl;
  ////////////////////////////////////////
  ////////////////////////////////////////



}


void print_lattice(int** spins, int n_spins)
{

  ////////////////////////////////////////
  /////////// PRINT OUT LATTICE ////////// 
  std::cout << std::endl;
  for(int x = 0; x < n_spins; x++) 
  {
    for(int y = 0; y < n_spins; y++)
    {
      std::cout << " " << spins[x][y] <<  " ";
    }
  std::cout << std::endl;
  }
  std::cout << std::endl;
  ////////////////////////////////////////
  ////////////////////////////////////////

}



void get_analytical_solutions(double Temp, int N)
{
    // double one_over_total_spins = Temp/(N*N);
    // double k_B = 1.0;
    // double beta = 1.0/(k_B*Temp);

    // double Z = 4.0*cosh(8*beta) + 12.0;
    // double E = -32.0*sinh(8*beta)/Z;
    // double E2 = 256.0*sinh(8*beta)/Z;
    // double M = (8.0*exp(8.0*beta) + 16.0)/Z;
    // double M2 = 32.0*(exp(8.0*beta) + 1)/Z;
    // double C_v = (E2 - E*E)/(k_B*Temp*Temp);
    // double X = (M2 - M*M)/(k_B*Temp);

    // std::cout << "<E>   = " << E*one_over_total_spins << std::endl;
    // std::cout << "<E^2> = " << E2*one_over_total_spins << std::endl;
    // std::cout << "<M>   = " << M*one_over_total_spins << std::endl;
    // std::cout << "<M^2> = " << M2*one_over_total_spins << std::endl;
    // std::cout << "C_v   = " << C_v*one_over_total_spins << std::endl;
    // std::cout << "X     = " << X*one_over_total_spins << std::endl << std::endl;


    double one_over_total_spins = 1.0/(2.0*2.0);
    double k_B = 1.0;
    double beta = 1.0/(k_B*Temp);

    double Z = 4.0*cosh(8*beta) + 12.0;
    double E = -32.0*sinh(8*beta)/Z;
    double E2 = 256.0*sinh(8*beta)/Z;
    double M = (8.0*exp(8.0*beta) + 16.0)/Z;
    double M2 = 32.0*(exp(8.0*beta) + 1)/Z;
    double C_v = (E2 - E*E)/(k_B*Temp*Temp);
    double X = (M2 - M*M)/(k_B*Temp);


    std::cout << std::setprecision(8) << "    E   = " << E*one_over_total_spins;
    std::cout << std::setprecision(8) << "    E^2 = " << E2*one_over_total_spins;
    std::cout << std::setprecision(8) << "    M  = " << M*one_over_total_spins;
    std::cout << std::setprecision(8) << "    M^2 = " << M2*one_over_total_spins;
    std::cout << std::setprecision(8) << "    C_v   = " << C_v*one_over_total_spins;
    std::cout << std::setprecision(8) << "    X     = " << X*one_over_total_spins << std::endl;

}







void write_to_file(double expectation[5], int MC_samples, int n_spins, std::string output_filename)
{

  std::ofstream outfile;
  outfile.open("energy_count.txt");
  outfile << std::setw(15) << "MC";
  outfile << std::setw(15) << "E";
  outfile << std::setw(15) << "E^2";
  outfile << std::setw(15) << "|M|";
  outfile << std::setw(15) << "M^2";
  outfile << std::setw(15) << "C_v";
  outfile << std::setw(15) << "X" << std::endl;



  // double norm = 1.0/((double) (m));
  // double x1 = expectation[0]*norm;
  // double x2 = expectation[1]*norm;
  // double x3 = expectation[2]*norm;
  // double x4 = expectation[3]*norm;
  // double x5 = expectation[4]*norm;

  // double HeatCapacity = (x2- x1*x1)*n_x_n_inv/Temp/Temp;
  // double chi = (x4 - x5*x5)*n_x_n_inv/Temp;


  // outfile << std::setw(15) << std::setprecision(8) << m;
  // outfile << std::setw(15) << std::setprecision(8) << (double) x1* n_x_n_inv;
  // outfile << std::setw(15) << std::setprecision(8) << (double) x2* n_x_n_inv;
  // outfile << std::setw(15) << std::setprecision(8) << (double) x5* n_x_n_inv;
  // outfile << std::setw(15) << std::setprecision(8) << (double) x4* n_x_n_inv;
  // outfile << std::setw(15) << std::setprecision(8) << (double) HeatCapacity;
  // outfile << std::setw(15) << std::setprecision(8) << (double) chi << std::endl;



}



