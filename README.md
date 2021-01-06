# Senior-Thesis-stuff
Files for Transdimensional Bayesian inversion of 1D geophysical resistivity surveys with a Monte Carlo Markov Chain approach

calculatedModel: a class made for ease of ensemble analysis and plotting; used in ensembleAnalysis3.

calculateRho1D: the script for the forward model in this inversion. Outputs apparent resistivity "measurements" given depths and 
resistivities for layers in a subsurface environment

createSyntheticData: Generates synthetic data for a virtual Schlumberger array field experiment given 
parameters like min and max electrode spacings, number of measurements, etc. Adds a specified amount of random Gaussian noise

...to be continued

MillsSeniorThesisMain - the primary script. Follows this structure:




Step 1: Load previous ensemble or create new one. 
A new ensemble can be created by running an inversion on either synthetic data or actual data loaded in (Note: not set up for this yet)

Step 2: Analyze ensemble with ensembleAnalysis function.