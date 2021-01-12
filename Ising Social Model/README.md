
Our last project of the semester. Using the Ising model to see social dynamics. Pretty cool, and we got a good feedback.


Since it is a lot of different commandline inputs that can be handled here it is:




To compile write :
  
```
g++ -fopenmp IsingSocialMain.cpp
```
  
a.out is created  

To run you will need to input one argument, how big of a system (#n spins)  

./a.out 10 

is 10 spins.
   
   <br>

Then there are several other applications you can input :
   <br>

###### Second argument is if you want the system ordered or not :  

./a.out 10 1   -> is a systems with 10 spins, started with all 1's  
./a.out 10 0   -> is a systems with 10 spins, each started -1 or 1 drawn from the unifrom distribution.  
   
   <br>
   
###### Third argument and fourth argument is to decide if you want to save the system to file at each step inputed and how    
3 = steps, 4 = magnetization or spins (0 for magnetizations and 1 for spins)
      
if steps is set to 0 we make a 1000 systems and divide it on numbers of threads available on your computer
  
./a.out 10 0 10 1  -> 10 spins, unordered, where each spin is saved to file every 10'th MC cycle  
./a.out 10 0 0 0  -> 10 spins, unordered, and now we run 1000 systems and save the final magnetizations from all the systems to one file  
./a.out 10 0 0 1  -> 10 spins, unordered, and now we run 1000 systems and save the final spins from all the systems to one file  

  <br>
  
###### Fifth argument is concentration level  

cB * N unordered  :  
./a.out 10 0 10 0 0.9  -> 10 = spins, 0 = unordered, 10 = save every 10'th MC cycle, 0 = Magnetization, 0.9 = 90% spins started with -1  

cB * N with the concentration cB starting in a cluster for index 0 to index cB * N :  
./a.out 10 1 10 0 0.9  -> 10 = spins, 1 = ordered, 10 = save every 10'th MC cycle, 0 = Magnetization, 0.9 = 90% spins started with -1  
  
  <br>
  
###### Sixth argument is chance of not following the rules  
  
./a.out 10 0 10 0 0.9 0.000002  -> 
10 = spins, 0 = unordered, 10 = save every 10'th MC cycle,   
0 = Magnetization, 0.9 = 90% spins started with 1, 0.000002 = 0.00002% chance of not following the rule
  
or with all spins starting with -1  :  
./a.out 10 0 10 0 1 0.000002  -> 
10 = spins, 0 = unordered, 10 = save every 10'th MC cycle,   
0 = Magnetization, 1 = all spins -1, 0.000002 = 0.00002% chance of not following the rule
  
  <br>
  
###### Seventh argument is whether to apply periodic or free boundary conditions 
   
./a.out 10 0 10 0 0.9 1 1  -> 
10 = spins, 0 = unordered, 10 = save every 10'th MC cycle,   
0 = Magnetization written out, 0.9 = 90% spins started with 1, 1 = following the rules, 1 = periodic boundary conditions
  
./a.out 10 0 10 0 0.9 1 0  -> 
10 = spins, 0 = unordered, 10 = save every 10'th MC cycle,   
0 = Magnetization written out, 0.9 = 90% spins started with 1, 1 = following the rules, 0 = free boundary conditions



