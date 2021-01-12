

/////// Creating a code that was efficient and easy to use with lots of alternatives turned out     /////////////
/////// to make the code terrible long. Something might have been avioded, but most of it are things      ///////
/////// that had to be done. It is what it is.                                                              /////


/////// To run the program put all files in a folder and write     make     on the commandline in your shell ////
/////// make use g++ version for compiling and if you do not have it open make file change line 4 to something //
/////// like      CPPflags= c++ -fopenmp -O3 -std=c++11                                                      ////
/////// After compiling write ./runme in the commandline and you will get alternatives where one is you  ////////
/////// setting your own configurations, and the other one is the configurations set between line 71 and ////////
/////// 93 in this file. If you want to configure through commandline write      ./runme 1               ////////
/////// and you will get an error where what to be put in is written                                     ////////




# include "Imodel.h"




int main(int argc, char *argv[])
{


  // Code for making directories to save the different results in
  // One for probability distribution of the energy, and one for stats
  errno = 0;
  int dir_result_prob = mkdir("probdist/", 0755);
  if(dir_result_prob != 0 && errno != EEXIST){
    printf("Unable to create directory\n"); 
    exit(1);   
  }
  else{
    printf("Folder prodist exist or have been created\n"); 
  }

  errno = 0;
  int dir_result_values = mkdir("value/", 0755);
  if(dir_result_values != 0 && errno != EEXIST){
    printf("Unable to create directory\n"); 
    exit(1);   
  }
  else{
    printf("Folder value exist or have been created\n"); 
  }


  double T_0, T_n, dN_T;  // Minumum temp to initialize system with, T_0, maximum, T_n, and steps starting from T_0, dN_t
  int MC_samples, n_temps, configuration, order, threads, method, n_spin_size, probability_dist, print_every_x_MC_cycle;
  std::ofstream* prob_outfile;
  std::ofstream* value_outfile;
  int * n_spins;
  double epsilon = 0.000001;




  // EXAMPLE SET UP THAT RUNS IF NO INPUT IS GIVEN

  // Temperature variables (min, max, steps, n_temps = how many of the steps to take)
  T_0 = 2; T_n = 3, dN_T = 0.2; n_temps = 1;

  // N in the NxN lattice, (ordered spins = 1 means all spins starts value = 1, order = 0 gives random spins 1 or -1 for each spin)
  n_spins = new int [1]; n_spins[0] = 20; n_spin_size = 1; order = 0; 

  // Values to print out. If probability_dist is set to 1 Energy values after 600000 are saved (will slow down the speed)
  // print_x_MC_cycle you can set to a value at which you want Energy, Magnetization, ect. something like 10000 or 100000
  // with 0 meaning save values only after all MC cylces are done
  probability_dist = 0; print_every_x_MC_cycle = 0;

  // MC_samples, and threads to use
  MC_samples = 1000000; threads = omp_get_max_threads();

  // This is the different methods of parallelization or not. 
  // method_to_use = 0 -> no parallelization.
  // method_to_use = 1 -> correct way of parallelization where each thread lock the 
  //                      spin and neighbouring spins while checking to flip it or not.
  // method_to_use = 2 -> parallelization where each thread dont care about if the neighbouring 
  //                      spin is beeing worked on by another thread. Can cause bias in the results.
  // method_to_use = 3 -> splitting the numbers of MC cycles on threads available, where each 
  //                      start their own system and do a total of ((MC cycles)/(threads used)) MC cycle each
  //                      then averaging over all runs. (Similar to model averaging in statistics).
  method = 3;





  // Read commandline arguments if it is provided
  if (argc > 1) 
  {
    if (argc > 1 && argc <= 5)
    {

      std::cout << "\n\nBad Usage: " << argv[0] << 
      " read N as in NxN spins, MC cycles to the as the power of 10, initial"
      " and final temperature, tempurate step : \n\n"
      "           ./runme 20 6 2 2.4 0.05 \n\n If you want you can after the 0.05 input 0 for unorder"
      " or 1 for order, then a 0 or 1 for saving Energy levels after 600000 MC_cycles(will slow down computations)."
      " After that if you want the different values saved for every 10^x MC cycles enter your x (4 or 5 is probably good)"
      " while a 0 is just saving a summary at after all MC cycles is done for each temperature."
      " Threads to run is next, where your max is "
      "" << omp_get_max_threads() << ", and one of the following methods "
      " \n "
      "\n \n 0 -> no parallelization. \n \n 1 -> correct way of parallelization where each thread lock the "
      "spin and neighbouring spins while checking to flip it or not."
      "\n \n 2 -> parallelization where each thread dont care about if the neighbouring "
      "spin is beeing worked on by another thread. Can cause bias in the results."
      "\n \n 3 -> splitting the numbers of MC cycles on threads available, where each "
      "start their own system and do a total of ((MC cycles)/(threads used)) MC cycle each "
      "then averaging over all runs. (Similar to model averaging in statistics). \n \n"
      "           ./runme 20 6 2 2.4 0.05 0 0 4 4 3 " << std::endl;
      exit(1);

    }
    else
    { 

      n_spins[0] = atoi(argv[1]); MC_samples = pow(10, atoi(argv[2]));
      T_0 = atof(argv[3]); T_n = atof(argv[4]); dN_T = atof(argv[5]); 
      if (argc >= 7) order = atoi(argv[6]);
      if (argc >= 8) probability_dist = atoi(argv[7]);
      if (argc >= 9) print_every_x_MC_cycle = (atoi(argv[8]) == 0) ? 0: pow(10, atoi(argv[8]));
      if (argc >= 10) threads = atoi(argv[9]);
      if (argc == 11) method = atoi(argv[10]);
      n_temps = (double) ((T_n+epsilon+dN_T) - T_0)/ dN_T;
      if (n_temps == 0) n_temps = 1; // Setting number of temperatures to run at least to 1      
      
      //////// TESTS //////////
      periodic_test();
      check_boundary_matrix();
      rng_initialization();
      /////   END TESTS   /////

    }
  }

  // If no commandline arguments are given, ask for them
  else
  {
    if (method != 0)
    {
    std::cout << "\n Press 1 and Enter if you want to set your own configurations. \n "
    "Press 0 and Enter if you want the system to run with "+std::to_string(n_spins[0])+"x"+std::to_string(n_spins[0])+" "
    "spins, "+std::to_string(MC_samples)+" MC cycles , temprature =  "+std::to_string(T_0)+", on "+std::to_string(threads)+""
    " threads\n";
    std::cin >> configuration;
    } 
    else
    {
    std::cout << "\n Press 1 and Enter if you want to set your own configurations. \n "
    "Press 0 and Enter if you want the system to run with "+std::to_string(n_spins[0])+"x"+std::to_string(n_spins[0])+" "
    "spins, "+std::to_string(MC_samples)+" MC cycles , temprature =  "+std::to_string(T_0);
    std::cin >> configuration;      
    }

    if (configuration == 1) 
    {

      std::cout << "Type min temp then hit Enter: "; std::cin >> T_0;
      std::cout << "Type max temp then hit Enter: "; std::cin >> T_n;
      std::cout << "Please type stepsize between the temperatures, and hit Enter: "; std::cin >> dN_T;
      std::cout << "How many different lattice size do you want to simulate : "; std::cin >> n_spin_size;
      n_spins = new int [n_spin_size];
      for (int i = 0; i < n_spin_size; i ++){std::cout << "Please enter 'N' in NxN spins for "
        "number "+std::to_string(i+1)+" out of "+std::to_string(n_spin_size)+" : "; std::cin >> n_spins[i];}
      std::cout << "Should the spins be ordered ? 0 and Enter for No, 1 and Enter for Yes: "; std::cin >> order;
      n_temps = ((T_n+epsilon+dN_T) - T_0)/ dN_T; if (n_temps == 0) n_temps = 1; // Setting number of temperatures to run at least to 1
      std::cout << "Thats gonna be " << n_temps << " different temperatures. "
      "Enter the number of MC cycles: "; std::cin >> MC_samples;
      std::cout << "How many threads you want to divide the labour on when it is possible to parallelize ?"
      " You have " << omp_get_max_threads() << " available : "; std::cin >> threads;
      std::cout << " \n\n Finaly enter the method of choice"
      "\n \n 0 -> no parallelization. \n \n 1 -> correct way of parallelization where each thread lock the "
      "spin and neighbouring spins while checking to flip it or not."
      "\n \n 2 -> parallelization where each thread dont care about if the neighbouring "
      "spin is beeing worked on by another thread. Can cause bias in the results."
      "\n \n 3 -> splitting the numbers of MC cycles on threads available, where each "
      "start their own system and do a total of ((MC cycles)/(threads used)) MC cycle each "
      "then averaging over all runs. (Similar to model averaging in statistics). \n \n "; std::cin >> method;
    }
    //////// TESTS //////////
    periodic_test();
    check_boundary_matrix();
    rng_initialization();
    /////   END TESTS   /////
  }




  // Just some precautions in case some values are set wrongly we set them to some values that will run

  // Checking temps
  if (n_temps > ((T_n+epsilon+dN_T) - T_0)/dN_T ){n_temps = int (T_n - T_0)/dN_T ;}
  if (n_temps == 0) n_temps = 1; // Setting number of temperatures to run at least to 1
  double T[n_temps];
  if (n_temps >= 1){ for (int i = 0; i < n_temps; i++) T[i] = T_0 + dN_T*i; }
  else{ T[1] = {T_0}; }
  // Cheking threads
  if (threads > omp_get_max_threads()) threads == omp_get_max_threads();
  // Set number of threads
  #define NUM_THREADS threads
  // Cheking spins
  for (int k = 0; k < n_spin_size; k++)
  {
    if (n_spins[k] < 1)n_spins[k] = 2;
  }
  // Cheking ordered or not
  if (order != 0 && order != 1) order = 0;
  // Cheking method
  if (method != 0 && method != 1 && method != 2 && method != 3) method = 3;
  // Checking probability distribution
  if (probability_dist != 0 && probability_dist != 1) probability_dist = 0;
  // Checking probability distribution
  if (print_every_x_MC_cycle < 1.1 || print_every_x_MC_cycle > MC_samples ) print_every_x_MC_cycle = 0;






  // Loop over different Lattice sizes
  for (int n = 0; n < n_spin_size; n++)
  {


  // Create a summary file if stats is only saved after completed MMC run
  if (print_every_x_MC_cycle == 0 && n_temps > 1)
  {

    prob_outfile = new std::ofstream;
    value_outfile = new std::ofstream;
    prob_outfile->open("probdist/Probdist_for_"+std::to_string(n_spins[n])+"_spins_and_"
                      "for_different_temps_"+std::to_string(MC_samples)+"_MC_samples.txt");
    value_outfile->open("value/Values_for_"+std::to_string(n_spins[n])+"_spins_and_"
                      "for_different_temps_"+std::to_string(MC_samples)+"_MC_samples.txt");
  }

    // Loop over different temperatures
    for (int t = 0; t < n_temps; t++)
    {


      // Create new file for each temperature if stats every x MC cycle is to be stored
      if (print_every_x_MC_cycle != 0 || n_temps == 1)
      {
        prob_outfile = new std::ofstream;
        value_outfile = new std::ofstream;
        prob_outfile->open("probdist/Probdist_for_"+std::to_string(n_spins[n])+"_spins"
                          "_and_temp_"+std::to_string(T[t])+"_"+std::to_string(MC_samples)+"_MC_samples.txt");
        value_outfile->open("value/Values_for_"+std::to_string(n_spins[n])+"_spins_and_temp"
                            "_"+std::to_string(T[t])+"_"+std::to_string(MC_samples)+"_MC_samples.txt"); 
      }


      // Avoiding some issues. Not all methods can do everything.
      int temp_method = method;
      // Since the model averaging can not print out average values for all the chains while they are running
      // we change the methods to method 1 where it is possible to save values at ever x MC cycle.
      if ((method == 3 && probability_dist == 1) || (method == 3 && print_every_x_MC_cycle != 0)){temp_method = 1;}
      // Avoiding problems occuring with parallization in method 1 with locks when same thread is
      // trying to lock the same index(i,j) twice, which occure in a 2x2 lattice
      if (method == 1 && n_spins[n] < 3){temp_method = 0;}



      // START OF METHOD ONE -> NO PARALLIZATION
      if (temp_method == 0)
      {
      

        double t0 = omp_get_wtime();  // Start timing of method
        double expectation[5] = {0,0,0,0,0};  // Setting up array to save E, E^2, M, M^2, and |M|
        long long changes_accepted = 0;


        // Make header to file
        if (t == 0){write_to_file(expectation, 0, n_spins[n], 0, 0, 0, 0, 0, value_outfile);}
        else if(print_every_x_MC_cycle != 0){write_to_file(expectation, 0, n_spins[n], 0, 0, 0, 0, 0, value_outfile);}


        MMC_boundary_matrix(n_spins[n], MC_samples, T[t], expectation, order, 
                            changes_accepted, temp_method);
        // MMC(n_spins[n], MC_samples, T[t], expectation, order, 
        //                     changes_accepted);


        write_to_file(expectation, MC_samples, n_spins[n], (double) 1./(n_spins[n]*n_spins[n]), T[t],
                      (double) 1./(T[t]), (double) 1./T[t]/T[t], changes_accepted, value_outfile);

        double time = omp_get_wtime()-t0; // Time used


        // Print results to commandline
        print_results(expectation, MC_samples, n_spins[n], (double) 1/(n_spins[n]*n_spins[n]), T[t],
                      (double) 1/(T[t]), (double) 1/T[t]/T[t]);


        std::cout << "\n Run of method "+std::to_string(temp_method)+" with lattice size"
        " "+std::to_string(n_spins[n])+"x"+std::to_string(n_spins[n])+" and"
        " "+std::to_string(MC_samples)+" MC samples with temperature "+std::to_string(T[t])+" "
        "took "+std::to_string(time)+" seconds \n" << std::endl;


      } // end method_to_use == 0




      // START OF METHOD TWO -> THE RIGHT WAY TO USE PARALLIZATION
      else if (temp_method == 1)
      {


        double t0 = omp_get_wtime();  // Start timing of method
        double expectation[5] = {0,0,0,0,0};  // Setting up array to save E, E^2, M, M^2, and |M|


        // Make header to file
        if (t == 0){write_to_file(expectation, 0, n_spins[n], 0, 0, 0, 0, 0, value_outfile);}
        else if(print_every_x_MC_cycle != 0){write_to_file(expectation, 0, n_spins[n], 0, 0, 0, 0, 0, value_outfile);}


        MMC_boundary_matrix_with_locks(n_spins[n], MC_samples, T[t], expectation, 
                                        order, threads, print_every_x_MC_cycle, probability_dist,
                                        value_outfile, prob_outfile);


        double time = omp_get_wtime()-t0; // Time used

        // Print results to commandline
        std::cout << "\n Run of method "+std::to_string(temp_method)+" with lattice size"
        " "+std::to_string(n_spins[n])+"x"+std::to_string(n_spins[n])+" and"
        " "+std::to_string(MC_samples)+" MC samples with temperature "+std::to_string(T[t])+" "
        "took "+std::to_string(time)+" seconds \n" << std::endl;

      } // end method_to_use == 1




      // // START OF METHOD THREE -> THE BIASED WAY TO USE PARALLIZATION
      else if (temp_method == 2)
      {


        double t0 = omp_get_wtime();  // Start timing of method
        double expectation[5] = {0,0,0,0,0};  // Setting up array to save E, E^2, M, M^2, and |M|


        // Make header to file
        if (t == 0){write_to_file(expectation, 0, n_spins[n], 0, 0, 0, 0, 0, value_outfile);}
        else if(print_every_x_MC_cycle != 0){write_to_file(expectation, 0, n_spins[n], 0, 0, 0, 0, 0, value_outfile);}


        MMC_boundary_matrix_parallel_spin(n_spins[n], MC_samples, T[t], expectation,
                                          order, threads, print_every_x_MC_cycle, probability_dist,
                                          value_outfile, prob_outfile);

        double time = omp_get_wtime()-t0; // Time used

        // Print results to commandline
        std::cout << "\n Run of method "+std::to_string(temp_method)+" with lattice size"
        " "+std::to_string(n_spins[n])+"x"+std::to_string(n_spins[n])+" and"
        " "+std::to_string(MC_samples)+" MC samples with temperature "+std::to_string(T[t])+" "
        "took "+std::to_string(time)+" seconds \n" << std::endl;

      } // end method_to_use == 2




      // // START OF METHOD FOUR -> THE MODEL AVERAGING WAY
      else if (temp_method == 3)
      {

  
        double t0 = omp_get_wtime();  // Start timing of method
        double expectation[5] = {0,0,0,0,0};  // Setting up array to save E, E^2, M, M^2, and |M|
        long long changes_accepted = 0;

  
        // Make header to file
        if (t == 0){write_to_file(expectation, 0, n_spins[n], 0, 0, 0, 0, 0, value_outfile);}
        else if(print_every_x_MC_cycle != 0){write_to_file(expectation, 0, n_spins[n], 0, 0, 0, 0, 0, value_outfile);}
  

        // OMP keeps track of the values we need (expectation) for the loops which each thread is running
        #pragma omp parallel for reduction(+:expectation)
        for (int c = 0; c < NUM_THREADS; c++)
        {
          MMC_boundary_matrix(n_spins[n], MC_samples/NUM_THREADS, T[t], expectation, order, 
                              changes_accepted, temp_method);

        }

        double time = omp_get_wtime()-t0; // Time used
        write_to_file(expectation, MC_samples, n_spins[n], (double) 1./(n_spins[n]*n_spins[n]), T[t],
                      (double) 1./(T[t]), (double) 1./T[t]/T[t], changes_accepted, value_outfile);

        // Print results to commandline
        print_results(expectation, MC_samples, n_spins[n], (double) 1/(n_spins[n]*n_spins[n]), T[t],
                      (double) 1/(T[t]), (double) 1/T[t]/T[t]);

        std::cout << "\n Run of method "+std::to_string(temp_method)+" with lattice size"
        " "+std::to_string(n_spins[n])+"x"+std::to_string(n_spins[n])+" and"
        " "+std::to_string(MC_samples)+" MC samples with temperature "+std::to_string(T[t])+" "
        "took "+std::to_string(time)+" seconds \n" << std::endl;

      } // end method_to_use == 3



    } // end n_spins_size
  // Close files
  prob_outfile->close();
  value_outfile->close();
  } // end n_temps

  return 0;
}














// Metropolis Monte Carlo
void MMC_boundary_matrix_with_locks(int n_spins, int MC_samples, double Temp, double expectation[5], 
                                    int order, int threads, int print_every_x_MC_cycle, int probability_dist,
                                    std::ofstream* value_outfile, std::ofstream* prob_outfile)
{

  // Initialize variables for Energy, Magnetization, Flips that is accepted by the MMC, and a tripple void 
  // pointer that is beeing used as a bigger matrix where the outer values are pointer to the inner lattize 
  // on the opposite side to be used for implementation of boundary condition in the Ising model. 

  double E = 0, M = 0;
  int changes_accepted = 0;
  void*** spins;


  // Create the bigger matrix with the Lattice "inside", and initialize values either ordered or unordered

  boundary_matrix(spins, n_spins);
  initialize_boundary_matrix(spins, n_spins, E, M, order);
  if (n_spins <= 10){print_boundary_lattice(spins, n_spins);}  // Print out the matrix if it is not to big


  // Keep all possible values of Energy changes in a array to reduce computations

  int all_dE[5] = {-8, -4, 0, 4, 8};
  double exp_dE[5];
  for (int k=0; k<5 ; k++) exp_dE[k] = exp(-all_dE[k]/Temp);


  // Initializing variables

  int rand_i, rand_j; // Random integers to be choosen as index in the Lattice when checking to flip a spin or not
  int dE = E; // Changes in energy
  int m, n; // Counters

  // Some more pre-calculations

  int n_x_n = n_spins*n_spins/threads; // Number of spins to check for each tread per MC cycle
  double n_x_n_inv = 1./(n_spins*n_spins);  // Number of total spins inverse, used in calculation of Energy, Magnetization, ect...
  double inv_Temp = 1./Temp;
  double inv_d_Temp = 1./Temp/Temp;

  // Initializing matrix of locks with boundary edges

  void*** count_lock;
  boundary_locks(count_lock, n_spins);
  initialize_boundary_locks(count_lock, n_spins);


  // Making locks for writing to file, so only one thread do at a time 

  omp_lock_t value_lock, prob_lock;
  omp_init_lock(&value_lock);
  omp_init_lock(&prob_lock);



  // First one is probably not needed but we define it anyway
  #define NUM_THREADS threads
  // Initializing OMP's parallelization, set number of threads, say that the variable they share is
  // m in how many MC cycles that is done, and setting each flip to be private (that is why we divide NxN on threads)
  #pragma omp parallel num_threads(threads) shared(m) private(n)
  {
    for (m = 1; m<MC_samples; m++)
    {
      n = 0; // Counter for spins checked by each thread

      // Say critical so that whenever a new MC cycle is started every thread has to be ready
      // and done with other operations
      #pragma omp critical


      // Loop over lattice size divided by active threads for each thread

      while (n < n_x_n)
      {

        // Finding a random spin to try and flip

        int rand_i = (int) (u_cont_d(randng)*(double)n_spins) + 1;
        int rand_j = (int) (u_cont_d(randng)*(double)n_spins) + 1;


        // Checking that no other thread is working on that spin found in Lattice(rand_i, rand_j)
        // and neighbours used for calculating to flip or not. If one of them is locked omp_test_lock
        // ensures that the thread is picking a new rand_i, rand_j by not doing what is in the if test

        if (omp_test_lock(&*((omp_lock_t *)(count_lock[rand_i][rand_j]))) && omp_test_lock(&*((omp_lock_t *)(count_lock[rand_i+1][rand_j])))
          && omp_test_lock(&*((omp_lock_t *)(count_lock[rand_i-1][rand_j]))) && omp_test_lock(&*((omp_lock_t *)(count_lock[rand_i][rand_j-1]))))
        {


        // Calculating Energy changed by flipping the spin

        dE = 2* *((int*)(spins[rand_i][rand_j])) *
                (*((int*)(spins[rand_i-1][rand_j])) + *((int*)(spins[rand_i+1][rand_j])) + 
                 *((int*)(spins[rand_i][rand_j-1])) + *((int*)(spins[rand_i][rand_j+1])));


        // Accept if the energy is getting lower

        if (dE <= 0)
        {
          *((int*)(spins[rand_i][rand_j])) *= -1;   // Update spin
          M += (double) (2* *((int*)(spins[rand_i][rand_j])));  // Update Magnetization
          E += (double) dE; // Update Energy
          changes_accepted += 1;
        }


        // If energy is not getting lower draw a random number uniformly between 0 and 1 and accept by a given probability

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

          if (u_cont_d(randng) <= current_exp_dE)
          {
              *((int*)(spins[rand_i][rand_j])) *= -1;    // Update spin
              M += (double) (2* *((int*)(spins[rand_i][rand_j])));  // Update Magnetization
              E += (double) dE; // Update Energy
              changes_accepted += 1;
          }
        }

        // Unlock the locks so the spins are open for other threads
        omp_unset_lock(&*((omp_lock_t *)(count_lock[rand_i][rand_j])));
        omp_unset_lock(&*((omp_lock_t *)(count_lock[rand_i+1][rand_j])));
        omp_unset_lock(&*((omp_lock_t *)(count_lock[rand_i-1][rand_j])));
        omp_unset_lock(&*((omp_lock_t *)(count_lock[rand_i][rand_j+1])));
        omp_unset_lock(&*((omp_lock_t *)(count_lock[rand_i][rand_j-1])));

        n = n + 1;  // Add that a spin is checked

        } // end of lock check
      } // end of a spin if locks are open



      // Update values after a sweep through the lattice

      expectation[0] += E; expectation[1] += E*E;
      expectation[2] += M; expectation[3] += M*M; 
      expectation[4] += fabs(M);


      // If values is to be printed every x cycle in the MC check 

      if (print_every_x_MC_cycle != 0 && m % print_every_x_MC_cycle == 0)
      {

        omp_set_lock(&value_lock);
        write_to_file(expectation, m, n_spins, n_x_n_inv, Temp, inv_Temp, inv_d_Temp, changes_accepted, value_outfile);
        omp_unset_lock(&value_lock);

      }


      // Check if prbability density is to be printed out 

      if (m > 600000 && probability_dist == 1)
      {
        omp_set_lock(&prob_lock);
        *prob_outfile << std::setw(15) << std::setprecision(8) << (double) E * n_x_n_inv << std::endl;
        omp_unset_lock(&prob_lock);

      }
    } //  End of a MC Cycle
  }  // End of MMC run 


  // If values where to be saved only once, we write to file now

  if (print_every_x_MC_cycle == 0)
  {
    write_to_file(expectation, m, n_spins, n_x_n_inv, Temp, inv_Temp, inv_d_Temp, changes_accepted, value_outfile);
  }


  // Print out some results to commandline

  print_results(expectation, MC_samples, n_spins, n_x_n_inv, Temp,
                inv_Temp, inv_d_Temp);


  // Release Memory

  boundary_matrix_relase(spins, n_spins);
  boundary_locks_relase(count_lock, n_spins);

}








// Metropolis Monte Carlo
void MMC_boundary_matrix_parallel_spin(int n_spins, int MC_samples, double Temp, double expectation[5],
                                       int order, int threads, int print_every_x_MC_cycle, int probability_dist,
                                       std::ofstream* value_outfile, std::ofstream* prob_outfile)
{

  // Initialize variables for Energy, Magnetization, Flips that is accepted by the MMC, and a tripple void 
  // pointer that is beeing used as a bigger matrix where the outer values are pointer to the inner lattize 
  // on the opposite side to be used for implementation of boundary condition in the Ising model. 

  double E = 0, M = 0;
  int changes_accepted = 0;
  void*** spins;


  // Create the bigger matrix with the Lattice "inside", and initialize values either ordered or unordered

  boundary_matrix(spins, n_spins);
  initialize_boundary_matrix(spins, n_spins, E, M, order);
  if (n_spins <= 10){print_boundary_lattice(spins, n_spins);}   // Print out the matrix if it is not to big


  // Keep all possible values of Energy changes in a array to reduce computations

  int all_dE[5] = {-8, -4, 0, 4, 8};
  double exp_dE[5];
  for (int k=0; k<5 ; k++) exp_dE[k] = exp(-all_dE[k]/Temp);


  // Initializing variables

  int rand_i, rand_j; // Random integers to be choosen as index in the Lattice when checking to flip a spin or not
  int dE = E; // Changes in energy
  int m, n; // Counters


  // Some more pre-calculations

  int n_x_n = n_spins*n_spins/threads; // Number of spins to check for each tread per MC cycle
  double n_x_n_inv = 1./(n_spins*n_spins);  // Number of total spins inverse, used in calculation of Energy, Magnetization, ect...
  double inv_Temp = 1./Temp;
  double inv_d_Temp = 1./Temp/Temp;


  // Making locks for writing to file, so only one thread do at a time 

  omp_lock_t value_lock, prob_lock;
  omp_init_lock(&value_lock);
  omp_init_lock(&prob_lock);


  // First one is probably not needed but we define it anyway
  #define NUM_THREADS threads
  // Initializing OMP's parallelization, set number of threads, say that the variable they share is
  // m in how many MC cycles that is done, and setting each flip to be private (that is why we divide NxN on threads)
  #pragma omp parallel num_threads(NUM_THREADS) shared(m) private(n)
  {
    for (m = 1; m<MC_samples; m++)
    {

      // Say critical so that whenever a new MC cycle is started every thread has to be ready
      // and done with other operations
      #pragma omp critical


      // Loop over lattice size divided by active threads for each thread

      for(n = 0; n < n_x_n; n++)
      {


        // Calculating Energy changed by flipping the spin

        int rand_i = (int) (u_cont_d(randng)*(double)n_spins) + 1;
        int rand_j = (int) (u_cont_d(randng)*(double)n_spins) + 1;


        // Calculating Energy changed by flipping the spin

        dE = 2* *((int*)(spins[rand_i][rand_j])) *
                (*((int*)(spins[rand_i-1][rand_j])) + *((int*)(spins[rand_i+1][rand_j])) + 
                 *((int*)(spins[rand_i][rand_j-1])) + *((int*)(spins[rand_i][rand_j+1])));


        // Accept if the energy is getting lower

        if (dE <= 0)
        {
          *((int*)(spins[rand_i][rand_j])) *= -1;     // Update spin
          M += (double) (2* *((int*)(spins[rand_i][rand_j])));  // Update Magnetization
          E += (double) dE; // Update Energy
          changes_accepted += 1;
        }


        // If energy is not getting lower draw a random number uniformly between 0 and 1 and accept by a given probability

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

          if (u_cont_d(randng) <= current_exp_dE)
          {
              *((int*)(spins[rand_i][rand_j])) *= -1;     // Update spin
              M += (double) (2* *((int*)(spins[rand_i][rand_j])));  // Update Magnetization
              E += (double) dE; // Update Energy
              changes_accepted += 1;
          }
        }

      }

      // Update values after a sweep through the lattice

      expectation[0] += E; expectation[1] += E*E;
      expectation[2] += M; expectation[3] += M*M; 
      expectation[4] += fabs(M);


      // If values is to be printed every x cycle in the MC check 

      if (print_every_x_MC_cycle != 0 && m % print_every_x_MC_cycle == 0)
      {

        omp_set_lock(&value_lock);
        write_to_file(expectation, m, n_spins, n_x_n_inv, Temp, inv_Temp, inv_d_Temp, changes_accepted, value_outfile);
        omp_unset_lock(&value_lock);

      }


      // Check if prbability density is to be printed out 

      if (m > 600000 && probability_dist == 1)
      {
        omp_set_lock(&prob_lock);
        *prob_outfile << std::setw(15) << std::setprecision(8) << (double) E * n_x_n_inv << std::endl;
        omp_unset_lock(&prob_lock);

      }
    } // End of a MC cycle
  }  // Slutt prallel 



  // If values where to be saved only once, we write to file now

  if (print_every_x_MC_cycle == 0)
  {
    write_to_file(expectation, m, n_spins, n_x_n_inv, Temp, inv_Temp, inv_d_Temp, changes_accepted, value_outfile);
  }


  // Print out some results to commandline

  print_results(expectation, MC_samples, n_spins, n_x_n_inv, Temp,
                inv_Temp, inv_d_Temp);


  // Release Memory

  boundary_matrix_relase(spins, n_spins);

}









// Metropolis Monte Carlo
void MMC_boundary_matrix(int n_spins, int MC_samples, double Temp, double* expectation,
                         int order, long long& changes_accepted, int method)
{

  // Initialize own random number generator, since slowed down the speed if different
  // threads is using the same for some reason
  std::random_device rdev1;
  std::mt19937_64 randng1(rdev1());
  std::uniform_real_distribution<double> u_cont_d1(0.0, 1.0);

  // Initialize variables for Energy, Magnetization, and a tripple void 
  // pointer that is beeing used as a bigger matrix where the outer values are pointer to the inner lattize 
  // on the opposite side to be used for implementation of boundary condition in the Ising model. 

  double E = 0, M = 0;
  void*** spins;


  // Create the bigger matrix with the Lattice "inside", and initialize values either ordered or unordered

  boundary_matrix(spins, n_spins);
  initialize_boundary_matrix(spins, n_spins, E, M, order);
  if (n_spins <= 10 && method == 0){print_boundary_lattice(spins, n_spins);}  // Print out the matrix if it is not to big
 

  // Keep all possible values of Energy changes in a array to reduce computations

  int all_dE[5] = {-8, -4, 0, 4, 8};
  double exp_dE[5];
  for (int k=0; k<5 ; k++) exp_dE[k] = exp(-all_dE[k]/Temp);


  // Initializing variables

  int rand_i, rand_j; // Random integers to be choosen as index in the Lattice when checking to flip a spin or not
  int dE = E; // Changes in energy
  int m; // Counters
  int n_x_n = n_spins*n_spins; // Pre-calculations number of spins in the lattice



  // Start MC cycles

  for (m = 1; m<MC_samples; m++)
  {


    // Loop over lattice size divided by active threads for each thread

    for (int n = 0; n<n_x_n; n++)
    {


      // Finding a random spin to try and flip

      int rand_i = (int) (u_cont_d1(randng1)*(double)n_spins) + 1;
      int rand_j = (int) (u_cont_d1(randng1)*(double)n_spins) + 1;


      // Calculating Energy changed by flipping the spin

      dE = 2* *((int*)(spins[rand_i][rand_j])) *
              (*((int*)(spins[rand_i-1][rand_j])) + *((int*)(spins[rand_i+1][rand_j])) + 
               *((int*)(spins[rand_i][rand_j-1])) + *((int*)(spins[rand_i][rand_j+1])));


      // Accept if the energy is getting lower

      if (dE <= 0)
      {
        *((int*)(spins[rand_i][rand_j])) *= -1;      // Update spin
        M += (double) (2* *((int*)(spins[rand_i][rand_j])));  // Update Magnetization
        E += (double) dE; // Update Energy
        changes_accepted += 1;
      }


      // If energy is not getting lower draw a random number uniformly between 0 and 1 and accept by a given probability

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

        if (u_cont_d1(randng1) <= current_exp_dE)
        {
            *((int*)(spins[rand_i][rand_j])) *= -1;      // Update spin
            M += (double) (2* *((int*)(spins[rand_i][rand_j])));  // Update Magnetization
            E += (double) dE; // Update Energy
            changes_accepted += 1;
        }

      } // End of a spin
    } // End of a MC cycle

    // Update values after a sweep through the lattice

    expectation[0] += E; expectation[1] += E*E;
    expectation[2] += M; expectation[3] += M*M; 
    expectation[4] += fabs(M);
  } // End of MMC


  // Relase memory

  boundary_matrix_relase(spins, n_spins);

}









// Metropolis Monte Carlo
void MMC(int n_spins, int MC_samples, double Temp, double expectation[5], int order, long long& changes_accepted)
{


  // Initialize variables for Energy, Magnetization, and a double pointer to int to create a matrix

  double E, M;
  int** spins;
  E = 0;
  M = 0;


  // Initialize the matrix with ordered or not ordered values, as well as calculate Energy and Magnetization

  matrix(spins, n_spins);
  initialize(spins, n_spins, E, M, order);
  if (n_spins <= 10){print_lattice(spins, n_spins);}  // Print out the matrix if it is not to big


  // Keep all possible values of Energy changes in a array to reduce computations

  int all_dE[5] = {-8, -4, 0, 4, 8};
  double exp_dE[5];
  for (int k=0; k<5 ; k++) exp_dE[k] = exp(-all_dE[k]/Temp);


  // Initializing variables

  int rand_i, rand_j; // Random integers to be choosen as index in the Lattice when checking to flip a spin or not
  int dE = E; // Changes in energy
  int m; // Counters
  int n_x_n = n_spins*n_spins; // Pre-calculations number of spins in the lattice


  // Start MC cycles

  for (m = 1; m<MC_samples; m++)
  {


    // Loop over lattice size divided by active threads for each thread

    for (int n = 0; n<n_x_n; n++)
    {


      // Finding a random spin to try and flip

      int rand_i = (int) (u_cont_d(randng)*(double)n_spins);
      int rand_j = (int) (u_cont_d(randng)*(double)n_spins);


      // Calculating Energy changed by flipping the spin

      dE = 2*spins[rand_i][rand_j] *
            (spins[get_periodic_index(rand_i-1,n_spins)][rand_j] + spins[get_periodic_index(rand_i+1,n_spins)][rand_j] + 
             spins[rand_i][get_periodic_index(rand_j-1,n_spins)] + spins[rand_i][get_periodic_index(rand_j+1,n_spins)]); 


      // Accept if the energy is getting lower

      if (dE <= 0)
      {

        spins[rand_i][rand_j] *= -1;       // Update spin
        M += (double) (2 * spins[rand_i][rand_j]);  // Update Magnetization
        E += (double) dE; // Update Energy
        changes_accepted += 1;

      }


      // If energy is not getting lower draw a random number uniformly between 0 and 1 and accept by a given probability

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

        if (u_cont_d(randng) <= current_exp_dE)
        {
          spins[rand_i][rand_j] *= -1;       // Update spin
          M += (double) (2 * spins[rand_i][rand_j]);  // Update Magnetization
          E += (double) dE; // Update Energy
          changes_accepted += 1;
        }

      } // End of a spin
    } // End of a MC cycle


    // Update values after a sweep through the lattice

    expectation[0] += E; expectation[1] += E*E;
    expectation[2] += M; expectation[3] += M*M; 
    expectation[4] += fabs(M);

  }

  // Relase memory

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
      E -= (double) spins[x][y]*
      (spins[get_periodic_index(x-1,n_spins)][y] + spins[x][get_periodic_index(y-1,n_spins)]); 
    }
  }

}







// This is our try to reduce time with boundary conditions. It is best described by this :

// 0 p p p 0 // Where the p         //  0          p_to_2_0    p_to_2_1     p_to_2_2    0           //
// p 1 1 1 p // are pointers to     //  p_to_0_2   0_0         0_1          0_2         p_to_0_0    //
// p 1 1 1 p // the (n-2)x(n-2)     //  p_to_1_2   1_0         1_1          1_2         p_to_1_0    //
// p 1 1 1 p // squareinside        //  p_to_2_2   2_0         2_1          2_2         p_to_2_0    //
// 0 p p p 0 // like this :         //  0          p_to_0_0    p_to_0_1     p_to_0_2    0           //

void initialize_boundary_matrix(void*** Mpointer, int n_spins, double& E, double& M, int order)
{

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
      E -= (double) *((int*)(Mpointer[x][y]))*
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

      //////// Creating the locks ////////
      if (i>0 && i<n_spins+1 && j>0 && j<n_spins+1 )
      {
        omp_init_lock(&*(omp_lock_t*) (count_lock[i][j]));
      } // locks made


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
      } // outer pointers to locks made


    }
  }
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



void get_analytical_solutions(double Temp)
{

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


    std::cout << "                ANALYTICAL SOLUTION         " << std::endl;
    std::cout << std::setprecision(8) << "    E   = " << E*one_over_total_spins;
    std::cout << std::setprecision(8) << "    E^2 = " << E2*one_over_total_spins;
    std::cout << std::setprecision(8) << "    M  = " << M*one_over_total_spins;
    std::cout << std::setprecision(8) << "    M^2 = " << M2*one_over_total_spins;
    std::cout << std::setprecision(8) << "    C_v   = " << C_v*one_over_total_spins;
    std::cout << std::setprecision(8) << "    X     = " << X*one_over_total_spins << std::endl;

}





void print_results(double expectation[5], int MC_samples, int n_spins, double n_x_n_inv, double Temp,
                   double inv_Temp, double inv_d_Temp)
{


  double MC_inv = 1.0/((double) (MC_samples));

  double system_mean_E = expectation[0] * MC_inv;  // Energy divided on MC cycles 
  double system_mean_E_sq = expectation[1] * MC_inv; // Energy squared divided on MC cycles
  double system_mean_M = expectation[2] * MC_inv;  //  Magnetization divided on MC cycles 
  double system_mean_M_sq = expectation[3] * MC_inv; //  Magnetization squared divided on MC cycles 
  double system_mean_abs_M = expectation[4] * MC_inv;  // Absolute magnetization divided on MC cycles 
  double var_E = ( system_mean_E_sq - system_mean_E*system_mean_E ) * n_x_n_inv; // Variation of the energy per spin
  double var_abs_M = (system_mean_M_sq - system_mean_abs_M*system_mean_abs_M )* n_x_n_inv; // Variation of the Magnetization per spin
  double HeatCapacity = var_E * inv_d_Temp;
  double chi = var_abs_M * inv_Temp;

  std::cout << "                CALCULATED VALUES         " << std::endl;
  std::cout << std::setprecision(8) << "    E : " << system_mean_E * n_x_n_inv;
  std::cout << std::setprecision(8) << "    E2 : " << system_mean_E_sq * n_x_n_inv;
  std::cout << std::setprecision(8) << "    M : " << system_mean_abs_M * n_x_n_inv;
  std::cout << std::setprecision(8) << "    M^2 : " << system_mean_M_sq * n_x_n_inv;
  std::cout << std::setprecision(8) << "    C_v : " << HeatCapacity;
  std::cout << std::setprecision(8) << "    X : " << chi << std::endl;



  if (n_spins == 2){get_analytical_solutions(Temp);}

}





void write_to_file(double expectation[5], int MC_samples, int n_spins, double n_x_n_inv, double Temp,
                   double inv_Temp, double inv_d_Temp, int changes_accepted, std::ofstream* outfile)
{

  // Setup of file, header line
  if (MC_samples == 0)
  {
    *outfile << std::setw(10) << "Num spins";
    *outfile << std::setw(18) << "Tempreature";
    *outfile << std::setw(18) << "MC cycles";
    *outfile << std::setw(18) << "Mean energy";  // mean E
    *outfile << std::setw(18) << "Mean abs mag";  // mean |M|
    *outfile << std::setw(18) << "Mean mag";  // Mean M
    *outfile << std::setw(18) << "Energy var ";  // (E^2 - (mean E)^2) 
    *outfile << std::setw(18) << "Mean abs var";  // (|M|^2 - (mean |M|)^2)
    *outfile << std::setw(18) << "Heat cap";  // C_v
    *outfile << std::setw(18) << "Mag susc";  // X
    *outfile << std::setw(18) << "Acc states" << std::endl;
  }

  else
  {
    double MC_inv = 1.0/((double) (MC_samples));

    double system_mean_E = expectation[0] * MC_inv;  // Energy divided on MC cycles 
    double system_mean_E_sq = expectation[1] * MC_inv; // Energy squared divided on MC cycles
    double system_mean_M = expectation[2] * MC_inv;  //  Magnetization divided on MC cycles 
    double system_mean_M_sq = expectation[3] * MC_inv; //  Magnetization squared divided on MC cycles 
    double system_mean_abs_M = expectation[4] * MC_inv;  // Absolute magnetization divided on MC cycles 

    double var_E = (system_mean_E_sq - system_mean_E*system_mean_E) * n_x_n_inv; // Variation of the energy per spin
    double var_abs_M = (system_mean_M_sq - system_mean_abs_M*system_mean_abs_M) * n_x_n_inv; // Variation of the Magnetization per spin
    double HeatCapacity = var_E * inv_d_Temp;
    double chi = var_abs_M * inv_Temp;

    *outfile << std::setw(10) << std::setprecision(8) << n_spins;  // "Num spins"
    *outfile << std::setw(18) << std::setprecision(8) << Temp; // "Tempreature"
    *outfile << std::setw(18) << std::setprecision(8) << MC_samples; // "MC cycles"
    *outfile << std::setw(18) << std::setprecision(8) << (double) system_mean_E * n_x_n_inv; // "Mean energy"
    *outfile << std::setw(18) << std::setprecision(8) << (double) system_mean_abs_M * n_x_n_inv; // "Mean abs mag"
    *outfile << std::setw(18) << std::setprecision(8) << (double) system_mean_M * n_x_n_inv; // "Mean mag
    *outfile << std::setw(18) << std::setprecision(8) << (double) var_E; // "Energy var "
    *outfile << std::setw(18) << std::setprecision(8) << (double) var_abs_M; // "Mean abs var"
    *outfile << std::setw(18) << std::setprecision(8) << (double) HeatCapacity;  // "Heat cap"
    *outfile << std::setw(18) << std::setprecision(8) << (double) chi; // "Mag susc"
    *outfile << std::setw(18) << std::setprecision(8) << (double) changes_accepted << std::endl; // "Acc states"
  }

}









void periodic_test()
{

  bool check1 = get_periodic_index(-1, 20) == 19;
  bool check2 = get_periodic_index(20, 20) == 0;

  if (!check1 || !check2){
      if (!check1) std::cout << "Lower periodic boundary condition FAILED, CHECK get_periodic_index function" << std::endl;
      if (!check2) std::cout << "Upper periodic boundary condition FAILED, CHECK get_periodic_index function" << std::endl;
  }

  else std::cout << "\n Periodic boundary conditions test has passed \n " << std::endl;

}



void check_boundary_matrix()
{

  // Checking a 10 times 10 matrix if the boundary conditions is working
  double E = 0, M = 0;
  void*** mat_check;
  int n_spins = 5;
  boundary_matrix(mat_check, n_spins);
  initialize_boundary_matrix(mat_check, n_spins, E, M, 0);

  int error = 0;

  std::cout << std::endl;
  for(int x = 1; x < n_spins+1; x++) 
  {
    for(int y = 1; y < n_spins+1; y++)
    {
      if (x == 1)
      {
        if ( *((int*)(mat_check[n_spins+1][y])) !=  *((int*)(mat_check[x][y]))) error += 1;
      }
      if (x == n_spins)
      {
        if ( *((int*)(mat_check[0][y])) !=  *((int*)(mat_check[x][y]))) error += 1;
      }

      if (y == 1)
      {
        if ( *((int*)(mat_check[x][n_spins+1])) !=  *((int*)(mat_check[x][y]))) error += 1;
      }
      if (y == n_spins)
      {
        if ( *((int*)(mat_check[x][0])) !=  *((int*)(mat_check[x][y]))) error += 1;
      }
    }
  }

  print_boundary_lattice(mat_check, n_spins);


  if (error == 0) std::cout << "\n Periodic boundary conditions test 2 has passed \n " << std::endl;
  else std::cout << "\n Something went wrong when initializing boundary matrix\n " << std::endl;


}




void rng_initialization()
{

  int checks = 1000;
  int not_different_numbers = 0;
  int wrong_numbers = 0;
  double value_1 = u_cont_d(randng);

  for (int k = 0; k < checks; k++)
  {
    double value_2 = u_cont_d(randng);
    if (value_2 == value_1) not_different_numbers += 1;
    if (value_2 > 1 && value_2 < 0) wrong_numbers += 1;
  }

  if (wrong_numbers != 0 || not_different_numbers > 5)
  {
    std::cout << " RNG SEEMS TO HAVE SOME ERRORS " << std::endl;
  }
  else std::cout << "\n RNG test passed \n" << std::endl;
}




