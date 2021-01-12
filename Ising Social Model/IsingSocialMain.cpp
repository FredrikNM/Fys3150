# include <iomanip>
# include <iostream>
# include <fstream>
# include <random>
# include <omp.h>
# include <algorithm>
# include <vector>
# include <string>
# include <sys/stat.h> 



// g++ -fopenmp IsingSocialMain.cpp


// using namespace std;
std::random_device rdev;
std::mt19937_64 randng(rdev());
std::uniform_real_distribution<double> u_cont_d(0.0, 1.0);


void vector(int*& Vptr, int n_agents);
void initialize(int n_agents, int order, double consentration, int* agents_system);
void MC(std::string filename, int mc_cycles, int n_agents, int* agents_system, double& M, int print_mc, double p_rules, int print_all_agents, int boundary);


int main(int argc, char *argv[])
{


    int order = 0;  // All agents starting with same opinion, or not. 0 means each agent is initialized with a 50% chance of beliving 1 or -1
    int systems_to_run;     // If several systems is runned on different threads, this tells each thread how many to run. Else it is set to 1
    int n_agents;   // Number of agents to initialize the system with
    int mc_cycles;  // How many MC cycles to run
    int print_mc = 0;   // Set to 10 means write to file every 10 MC cycle, set to 1000 means write to file every 1000 mc cycle. 0 means only last result.
    int print_all_agents = 1;   // 1 will write out the system, and set to 0 will write out the magnetization
    int boundary = 0;   // 0 to not apply boundary conditions, 1 to apply them
    double p_rules = 1;     // Set to 1 to follow the rules. NB : if set to 0.000002 then every 0.00002 agents will choose opinion 1 with a 50% chance and -1 with a 50% chance.
    double consentration = 0;   // Set the initial concentration of agents beliving/meaning the same
    std::string filename_start;

    n_agents = atoi(argv[1]);
    if (argc >= 3){order = atoi(argv[2]);}   // If input bigger than 0, you will print out results at each step you input
    if (argc >= 4){print_mc = atoi(argv[3]);}     // If input bigger than 0, you will print out results at each step you input
    if (argc >= 5){print_all_agents = atoi(argv[4]);}     // Input 1 to print out all spins, or input 0 to print out magnetization
    if (argc >= 6){consentration = atof(argv[5]);}    //  Input a concentration. 1 means all 1, and 0.9 means 90% 1's
    if (argc >= 7){p_rules = atof(argv[6]);}    // Probability of not following the neighbour rules, and an agent just randomly choose what to be
    if (argc >= 8){boundary = atoi(argv[7]);}   // Apply boundary conditions if set to 1

    mc_cycles = 100*n_agents;   // Found to give stable results



    // Code for making directory to save the results in
    errno = 0;
    int dir_result_prob = mkdir("Results/", 0755);
    if(dir_result_prob != 0 && errno != EEXIST){
    printf("Unable to create directory\n"); 
    exit(1);   
    }
    else{
    printf("Folder prodist exist or have been created\n"); 
    }

    // Setting file name according to Magnetization written out or not
    int fileID = (int) (u_cont_d(randng)*(double)1000); // Setting a random fileID to avoid overwritting
    if (print_mc == 0)
    {
        filename_start = "MResult";
        systems_to_run = (int) 1000 / omp_get_max_threads();    // Dividing 1000 systems on the number of threads available
        omp_set_num_threads(omp_get_max_threads());
    }
    else
    {
        filename_start = "Result";
        systems_to_run = 1;
        omp_set_num_threads(1);
    }
    // Setting rest of filename to agents, MC cyles and fileID
    std::string filename = "Results/"+filename_start + std::to_string(n_agents) + "Agents" + std::to_string(mc_cycles) + "McCycles" + std::to_string(fileID) + ".txt";
    std::ifstream ifile;
    ifile.open(filename);



    std::ofstream outfile;
    outfile.open(filename, std::ios_base::app); // append instead of overwrite in file


    // // Makes locks so only one thread can write to file at a time
    omp_lock_t value_lock;
    omp_init_lock(&value_lock);


    #pragma omp parallel
    {   // Start parallelization
    for (int n = 0; n < systems_to_run; n++)
    {   // Run all systems

        int *agents_system;
        double M = 0;
        vector(agents_system, n_agents);
        initialize(n_agents, order, consentration, agents_system);
        MC(filename, mc_cycles, n_agents, agents_system, M, print_mc, p_rules, print_all_agents, boundary);


        if (print_mc == 0)
        {
            if (print_all_agents == 0)
            {
                for(int i=1; i < n_agents+1; i++)    // Start with 1 since we initialize a system with index -1 = 0
                {
                    M += (int) agents_system[i];    // Calculationg magnetization
                }
                omp_set_lock(&value_lock);
                outfile << M << "\n";   // Write final magnetization to file
                omp_unset_lock(&value_lock);
            }
            if (print_all_agents == 1)
            {
                omp_set_lock(&value_lock);
                for(int i=1; i < n_agents+1; i++)    // Start with 1 since we initialize a system with index -1 = 0
                {
                    outfile << std::setw(3) << agents_system[i];
                }
                outfile << "\n"; 
                omp_unset_lock(&value_lock);
            }
        }
        delete [] agents_system;    // Release Memory
        std::cout << "Thread nr : " << omp_get_thread_num()+1 << " running number : " << n+1 << " out of " <<  systems_to_run << std::endl;    // Print out which thread is running

    }   // Finished running systems
    }   // End of parallelization
    outfile.close();

}







// function to initialize agent_vector
void vector(int*& Vptr, int n_agents)
{
    Vptr = (int *) malloc((n_agents+1)*sizeof(int));
}





// Initialize the agent_system
void initialize(int n_agents, int order, double consentration, int* agents_system)
{


    if (consentration == 0)
        {   // Initialize randomly from a uniform distribution
        for(int i = 0; i < (n_agents+1); i++) 
        {
            // order = 0 gives random ordering, and 1 gives all equal to -1
            agents_system[i] = (order == 1) ? 1 : (u_cont_d(randng) > 0.5) ? 1 : -1;
        }   // End of initialization
    }
    else
    {   // Initialize with concentration
        if (order == 0) // Randomized order
        {
            std::vector<int> v;
            for (int i = 0; i < n_agents; i++)
            {
                v.push_back(i);
            }
            std::shuffle(v.begin(), v.end(), randng);   // Make a shuffled index
            agents_system[0] = 0;
            for (int i = 0; i < n_agents; i++)
            {
                agents_system[v[i]+1] = (i+1 < consentration*n_agents) ? -1 : 1;
            }
        }   // End of randomizing order
        else
        {   // Initialize with a cluster
            agents_system[0] = 0;
            for (int i = 0; i < n_agents; i++)
            {
                agents_system[i+1] = (i+1 < consentration*n_agents) ? -1 : 1;
            }
        }   // End of cluster
    }   // End of concentration
}   // End Initialize








// Run MC-cycles
void MC(std::string filename, int mc_cycles, int n_agents, int* agents_system, double& M, int print_mc, double p_rules, int print_all_agents, int boundary)
{

    std::ofstream outfile;
    if (print_mc != 0){outfile.open(filename, std::ios_base::app);} // append instead of overwrite

    if (p_rules == 1)
    {   // Start MC cylces
        for (int j = 0; j<mc_cycles; j++)
        {

            if (boundary == 0)
            {   // No boundary conditions
                for (int i = 0; i<n_agents; i++)
                {

                    // Pick a random index
                    int rand_i = (int) (u_cont_d(randng)*(double)(n_agents-1)) + 1 ;

                    // Updating rules
                    if (agents_system[rand_i]*agents_system[rand_i+1]==1)
                    {
                        agents_system[rand_i-1] = agents_system[rand_i];
                        agents_system[rand_i+2] = agents_system[rand_i];
                    }
                    if (agents_system[rand_i]*agents_system[rand_i+1]==-1)
                    {
                        agents_system[rand_i-1] = agents_system[rand_i+1];
                        agents_system[rand_i+2] = agents_system[rand_i];
                    }            
                }
            }   // No boundary conditions loop finished
            else
            {   // With boundary conditions
                for (int i = 0; i<n_agents; i++)
                {   // Applying rules n times


                    // Pick a random index
                    int rand_i = (int) (u_cont_d(randng)*(double)(n_agents-1)) + 1 ;

                    // Updating rules for boundary conditions
                    if(rand_i == 1)
                    {
                        if (agents_system[rand_i]*agents_system[rand_i+1]==1)
                        {
                            agents_system[n_agents] = agents_system[rand_i];
                            agents_system[rand_i+2] = agents_system[rand_i];
                        }
                        if (agents_system[rand_i]*agents_system[rand_i+1]==-1)
                        {
                            agents_system[n_agents] = agents_system[rand_i+1];
                            agents_system[rand_i+2] = agents_system[rand_i];
                        }   
                    }
                    else if(rand_i == n_agents-1)
                    {
                        if (agents_system[rand_i]*agents_system[rand_i+1]==1)
                        {
                            agents_system[rand_i-1] = agents_system[rand_i];
                            agents_system[1] = agents_system[rand_i];
                        }
                        if (agents_system[rand_i]*agents_system[rand_i+1]==-1)
                        {
                            agents_system[rand_i-1] = agents_system[rand_i+1];
                            agents_system[1] = agents_system[rand_i];
                        }   
                    }
                    else if(rand_i == n_agents)
                    {
                        if (agents_system[rand_i]*agents_system[1]==1)
                        {
                            agents_system[rand_i-1] = agents_system[rand_i];
                            agents_system[2] = agents_system[rand_i];
                        }
                        if (agents_system[rand_i]*agents_system[1]==-1)
                        {
                            agents_system[rand_i-1] = agents_system[1];
                            agents_system[2] = agents_system[rand_i];
                        }   
                    }
                    else if (rand_i != 1 && rand_i != n_agents && agents_system[rand_i]*agents_system[rand_i+1]==1)
                    {
                        agents_system[rand_i-1] = agents_system[rand_i];
                        agents_system[rand_i+2] = agents_system[rand_i];
                    }
                    else if (rand_i != 1 && rand_i != n_agents && agents_system[rand_i]*agents_system[rand_i+1]==-1)
                    {
                        agents_system[rand_i-1] = agents_system[rand_i+1];
                        agents_system[rand_i+2] = agents_system[rand_i];
                    }
                }   // Finised applying rules n times
            }   // With boundary conditions finished




        // Save magnetization or spins every step(= print_mc) cylce
        if (print_mc != 0 && j % print_mc == 0)
        {

            if (print_all_agents == 0)
            {
                M = 0;
                for(int i=1; i < n_agents+1; i++)    // Start with 1 since we initialize a system with index -1 = 0
                {
                    M += (double) agents_system[i];
                }
                outfile << M <<"\n";
            }
            else
            {
                for(int i=1; i < n_agents+1; i++)    // Start with 1 since we initialize a system with index -1 = 0
                {
                    outfile << std::setw(3) << agents_system[i];
                }
                outfile << "\n";
            }
        }   // End of writing to file
        } // End of a cycle
    }   // End for loop of MC cycle



    else
    {   // Start MC cycles with a probability = p_rules to of not following rules
        for (int j = 0; j<mc_cycles; j++)
        {


            if (boundary == 0)
            {   // No boundary conditions
                for (int i = 0; i<n_agents; i++)
                {

                    // Pick a random index
                    int rand_i = (int) (u_cont_d(randng)*(double)(n_agents-1)) + 1 ;

                    // Updating rules
                    if (agents_system[rand_i]*agents_system[rand_i+1]==1)
                    {
                        agents_system[rand_i-1] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                        agents_system[rand_i+2] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                    }
                    if (agents_system[rand_i]*agents_system[rand_i+1]==-1)
                    {
                        agents_system[rand_i-1] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i+1] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                        agents_system[rand_i+2] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                    }            
                }
            }   // No boundary conditions loop finished

            else
            {   // With boundary conditions
                for (int i = 0; i<n_agents; i++)
                {   // Applying rules n times


                    // Pick a random index
                    int rand_i = (int) (u_cont_d(randng)*(double)(n_agents-1)) + 1 ;

                    // Updating rules for boundary conditions
                    if(rand_i == 1)
                    {
                        if (agents_system[rand_i]*agents_system[rand_i+1]==1)
                        {
                            agents_system[n_agents] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                            agents_system[rand_i+2] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                        }
                        if (agents_system[rand_i]*agents_system[rand_i+1]==-1)
                        {
                            agents_system[n_agents] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i+1] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                            agents_system[rand_i+2] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                        }   
                    }
                    else if(rand_i == n_agents-1)
                    {
                        if (agents_system[rand_i]*agents_system[rand_i+1]==1)
                        {
                            agents_system[rand_i-1] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                            agents_system[1] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                        }
                        if (agents_system[rand_i]*agents_system[rand_i+1]==-1)
                        {
                            agents_system[rand_i-1] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i+1] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                            agents_system[1] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                        }   
                    }
                    else if(rand_i == n_agents)
                    {
                        if (agents_system[rand_i]*agents_system[1]==1)
                        {
                            agents_system[rand_i-1] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                            agents_system[2] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                        }
                        if (agents_system[rand_i]*agents_system[1]==-1)
                        {
                            agents_system[rand_i-1] = (u_cont_d(randng) > p_rules) ? agents_system[1] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                            agents_system[2] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                        }   
                    }
                    else if (rand_i != 1 && rand_i != n_agents && agents_system[rand_i]*agents_system[rand_i+1]==1)
                    {
                        agents_system[rand_i-1] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                        agents_system[rand_i+2] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                    }
                    else if (rand_i != 1 && rand_i != n_agents && agents_system[rand_i]*agents_system[rand_i+1]==-1)
                    {
                        agents_system[rand_i-1] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i+1] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                        agents_system[rand_i+2] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                    }
                }   // Finised applying rules n times
            }   // With boundary conditions finished




        // Save magnetization or spins every step(= print_mc) cylce
        if (print_mc != 0 && j % print_mc == 0)
        {
            if (print_all_agents == 0)
            {
                M = 0;
                for(int i=1; i < n_agents+1; i++)   // Start with 1 since we initialize a system with index -1 = 0
                {
                    M += (double) agents_system[i];
                }
                outfile << M <<"\n";
            }
            else
            {
                for(int i=1; i < n_agents+1; i++)    // Start with 1 since we initialize a system with index -1 = 0
                {
                    outfile << std::setw(3) << agents_system[i];
                }
                outfile << "\n";
            }
        }   // End of writing to file

        }   // End of a cycle
    }   // End of MC cycles with probability of breaking the rules
}   // End MC cycles




