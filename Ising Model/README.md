  
Download all files and write make in the cmd/console.
    
# Run

Then you write :
```
./runme
```

and the program gives you alternatives. 

# Another possibility is to give commandline inputs

Write 
```
./runme 1
```
which will give you an error and tell you what you can input.

# Parallelization

It did go a lot work in to parallelize the code, with 3 different methods of doing it. In the work to implement it time got lost and the code ended up with some memory leak. On the other hand, if you would like to see how you can use OMP to parallelize your code, the different "void MMC" methods is well explained and contains most of the variables inside so it should be easy to follow as well. The methods are also well explained in the PDF

# The report

We got some critic for some of the math, not talking specifically about the Metropolis Monte Carlo method, that there is other/better methods for the PBC, "Wrong statement about the Boltzmann dist. The Boltzman dist. gives the probability for a system state, e.g. P(state i) -- it does not give the probability for an energy P(E), which is what your histograms show.", also that some of the results needed more discussion.



