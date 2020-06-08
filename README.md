## Introduction to Q-Factor sideband subtraction at GlueX
[Original Paper - Multivariate side-band subtraction using probabilistic event weights](https://arxiv.org/pdf/0809.2548.pdf)


Calculating Q-factors in the reaction <img src="https://render.githubusercontent.com/render/math?math=\gamma p\rightarrow\pi^0\eta p \rightarrow 4\gamma p"> at GlueX.
Q-factors is an event-by-event multivariate sideband subtraction technique. The only requirement is the knowledge of the signal and background distribution of some discriminating variable. 
1. First the nearest neighbors, under some set of phase space variables, is found for a given event
2. Distribution of the discriminating variable is filled with the nearest neighbors
3. Fit the above distribution with the known/assumed signal and bkg distribution
4. Calculate Q-factor as the signal fraction
5. Do this for all entries


## Requirements
- Minuit is not threadsafe so this code requires TMinuit2. I believe the ROOT version needed for this is atleast 6.19.
- All combinatons from all events that passed your selections should be kept in a root file. Each combination is an entry. At GlueX, this is the flat tree format from the DSelector.  

## Code
There are a few programs working together:
1. getInitParams.C does the initial fit to the full discriminating variable's distribution. Save the fit parameters to a file which is then read in later
2. main.h/C is the Q-factor program. QFactorAnalysis class is defined here. QFactorAnalysis class has multiple methods which load the data and fitted parameters from getInitParams and spawns the threads to do the analysis
3. makeDiagnosticHists loads all the q-factor results which were added all together using root's hadd function. Various plots are made with the q-factor weighting. 
4. run.py drives the main program and makeDiagnosticHists. Most of the important variables are configued here and modifies main directly using search and replace. Might be slightly dangerous? main is also compiled at this stage, directories are cleaned to reduce confusion. 
5. convertROOTtoPNG.C is needed since outputting histograms into image formats, in a multithreaded way, might cause errors. Maybe due to some blocking issues. Anyways, we can save them all as root files and then convert them to pngs or whatever
6. helperFunc.h contains a bunch of auxilliary functions that help run the code. The most important ones to set are the signal/bkg distributions and the parameter limits and degrees of freedom.
7. getUniquenessWeights.C is being implemented to extract the uniqueness tracking weights, where each combo is weighted by the number of combos that passed event selections in a given event. This is done post DSelector for ease of use. The current algorithm in the Q facotr analysis simply checks to see if the combo has been seen before, if it hasnt then it is included. How we define combo can be difficult to assume, so an alternative method can be insightful. i.e. should we be looking at uniqueness of the spectroscopic photons in this case?

Things you need to configure:
1. Most importantly the settings at the top of the run.py file needs to configured properly.
2. The distribution for the signal and bkg must be given. "main" program uses fitFunc, background, and signal function definitions defined in helperFunc.h. A bernstein polynomial bkg is taken with a Gaussian signal. The parameter limits must be set also be set in helperFunc.h. Bernstein polynoimal is useful since it can be made non-negative for any domain.


## Outputs
- fitResults folder is created to hold information about the initialization parameters. These parameters are obtained from a fit to the reference distribution using the entire dataset
- histograms folder contains the saved Q-factor histograms for various combo entries. 
- logs folder contains information about the processing times and chiSqs for each combo. From here we can try to improve performance.


## Functionalities to consider
There are a bit of technical issues which this package tries to tackle. 
- How to include Accidental subtraction and sideband subtraction (if we wanted to)? Should we apply it before the calculation of the q-factors or use the accidental weights during the calculation. I would argue for the latter since there might be some double counting if we dont. 
- Standardization of the phase space variables is needed. Range and standard deviation standardization is implemented
- Convergence of a fit requires a decent initialization. Currently, the full distribution of the discriminating variable is fitted and the relevant parameters are scaled down to match the number of neighbors
- How does the initialization of the fit function change the results? Original paper suggests doing 3 fits per iteration where there is 100% bkg, 50/50% bkg/sig, and 100% signal. Implemented.
- For each combo in an entire dataset, is there a criteria for what it's potential neighbors can be? For each event, there can be multiple combinations that pass our selection criteria. How should we deal with these inter-event combinations? Currently there is a implementation that looks for spectroscopically unique combos only. Inter-event combinations cannot use the same 4 photons. It is not immediately obvious what determines the uniqueness of a combo and how to pick the appropriate one among combos of the same identity. Instead of implemnting it this way, we can try to weigh each combo by 1/N where N is the number of combos that pass event selections for a given event.
- What type of fit functions can one use? Function defintions for Gaussian, double-Gaussian, Breit Wigner has been defined for the signal distribution. The bkg distribution is currently a bernstein polynomial. The current logic goes like this. The efficiency of the detector should be relatively flat if the nearest neighbors all cluster closely in some specific region of the detector. Under this scenario the signal distribution should take more of a Gaussian shape. That would really depend on how large the sample size is and what the phase space variables are.
- The detector that the particle was detected in also matters. There is resolution differences between a pair of photons detected in the BCAL, in the FCAL, and where one photon ends up in the BCAL and the other ends up in the FCAL. Currently, the data is separated in the DSelector stage and run independently. If we can do this, it is greatly beneficial since Q-factor analysis is at least N^2 in complexity.
- What is best way to find the nearest neighbors. Currently, a priority queue is set up which has O(logn) to insert and to remove the smallest element. If doing only one fit per entry then finding neighbors takes 90%+ of the time. Neighbor finding (distance calculating + sorting) is inherently vectorizable.
- How do we know which phase space variables to select? Currently, the program is able to take in a set of variables and just scan through all possible permutations and subsets. Left to the analyzer to decide from there. Good idea to turn down the number of nentries to consider and do a qualitative scan if things take too long...
- Does k random neighbors make sense as a check? Looking for the nearest neighbors makes sense. Looking at the furthest neighbors might not make much sense as a check.  The performance of the Q-factors would take contibutions from the fit function and from the neighbors it selects. We might be able to decouple the performance and look directly at the performance of the fit function if we use k random neighbors. A flag is set up to do this.
- There is a lot of statistics at GlueX, especially in the channel I am interested in. Currently, Q-factors calcuation is multithreaded (require Minuit2 to be thread safe). Aggregating the results in done by another program, makeDiagnosticHists to allow this parallelization. Future works might be in expanding this to HPC job submission and or include GPU compatibiilty. Accelerating the neighbor finding (most time consuming) would probably have the most benefits. Particularly useful if bootstrapping to find errors since there would be proportionally less data transfers between host and device. Almost have a working version and then profiling the acceleration would be the next step. 
