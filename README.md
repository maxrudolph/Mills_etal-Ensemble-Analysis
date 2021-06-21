# Senior-Thesis-stuff
Files for Transdimensional Bayesian inversion of 1D geophysical resistivity surveys with a Monte Carlo Markov Chain approach

Step 1: Synthetic Data Generation
  - Master Script: createSyntheticData
    - Generates and saves synthetic data for a virtual Schlumberger array field experiment based on user-set parameters
  - Subscripts:
    - calculateRho1D: The forward model for resistivity measurements via Schlumberger array. Outputs apparent resisitivity "measurements" given subsurface structure properties
    - makeLambda: Makes a 'lambda matrix' based on electrode spacings which is required for calculateRho1D
    - subStructGen: A library of virtual subsurface structures from which to generate synthetic data.

Step 2: Inversion
  - Master script: inversion
    - Loads a data file (as generated in step 1) and performs an MCMC inversion on it to produce a solution ensemble 
  - Subscripts:
    - chooseOption: Used within mcmcAlgorithm. Controls the choice of what way the proposed solution is edited in a given step.
    - genericSln: A class that contains all the information about a solution and all the necessary methods to edit itself
    - mcmcAlgorithm: The overarching file that does the inversion process. The 'inversion' master script is really just a place for a user to change options and parameter bounds, and to save the results after mcmcAlgorithm is done, but mcmcAlgorithm is the bulk of the process.



calculatedModel: a class made for ease of ensemble analysis and plotting; used in ensembleAnalysis3.

doSaving: Just a save statement for when we were generating a lot of large ensembles;
matlab won't allow you to put a save statement in a parfor loop but if you hide that save statement in a function, it works!

ensembleAnalysis: outdated, use ensembleAnalysis3 now



genericMedium is a class which creates a self-contained object used for keeping track of, updating, and checking
all of the information contained in a proposed earth model (layers with resistivities, that model's associated 
misfit, etc.

...to be continued

MillsSeniorThesisMain - the primary script for now
