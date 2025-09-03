# Senior-Thesis-stuff
Files for Transdimensional Bayesian inversion of 1D geophysical resistivity surveys with a Markov Chain Monte Carlo approach

Step 1: Synthetic Data Generation
  - Main script: main_inversion.m
    - Generates and saves synthetic data for a virtual Schlumberger array field experiment based on user-set parameters
  - Subscripts:
    - calculateRho1D: The forward model for resistivity measurements via Schlumberger array. Outputs apparent resisitivity "measurements" given subsurface structure properties
    - makeLambda: Makes a 'lambda matrix' based on electrode spacings which is required for calculateRho1D
    - subStructGen: A library of virtual subsurface structures from which to generate synthetic data.

Step 2: Inversion
  - Main script: main_inversion.m
    - Loads a data file (as generated in step 1) and performs an MCMC inversion on it to produce a solution ensemble 
  - Subscripts:
    - chooseOption: Used within mcmcAlgorithm. Controls the choice of what way the proposed solution is edited in a given step.
    - genericSln: A class that contains all the information about a solution and all the necessary methods to edit itself
    - mcmcAlgorithm: The overarching file that does the inversion process. The 'inversion' master script is really just a place for a user to change options and parameter bounds, and to save the results after mcmcAlgorithm is done, but mcmcAlgorithm is the bulk of the process.

Step 3: Analysis
     - Main script: main_analysis.m

Step 4: Making figures
     - The figures for Mills, Rudolph, and Lekic can be reproduced using the scripts in Step4FigurePlottingScripts

Sequencing Ensemble Solutions
     - A Python language jupyter notebook is in the Sequencing/ subdirectory and can be used to sequence a subset of solutions from an ensemble.

