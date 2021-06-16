# Senior-Thesis-stuff
Files for Transdimensional Bayesian inversion of 1D geophysical resistivity surveys with a Monte Carlo Markov Chain approach

Step 1: Synthetic Data Generation
  - Master Script: Create Synthetic Data
    - Generates and saves synthetic data for a virtual Schlumberger array field experiment based on user-set parameters
  - Subscripts:
    - calculateRho1D: The forward model for resistivity measurements via Schlumberger array. Outputs apparent resisitivity "measurements" given subsurface structure properties
    - makeLambda: Makes a 'lambda matrix' based on electrode spacings which is required for calculateRho1D
    - subStructGen: A library of virtual subsurface structures from which to generate synthetic data.

Step 2: Inversion



calculatedModel: a class made for ease of ensemble analysis and plotting; used in ensembleAnalysis3.

doSaving: Just a save statement for when we were generating a lot of large ensembles;
matlab won't allow you to put a save statement in a parfor loop but if you hide that save statement in a function, it works!

ensembleAnalysis: outdated, use ensembleAnalysis3 now



genericMedium is a class which creates a self-contained object used for keeping track of, updating, and checking
all of the information contained in a proposed earth model (layers with resistivities, that model's associated 
misfit, etc.

...to be continued

MillsSeniorThesisMain - the primary script for now
