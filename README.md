## Introduction to Q-Factor sideband subtraction at GlueX
[Original Paper - Multivariate side-band subtraction using probabilistic event weights](https://arxiv.org/pdf/0809.2548.pdf)


Calculating Q-factors in the reaction <img src="https://render.githubusercontent.com/render/math?math=\gamma p\rightarrow\pi^0\eta p \rightarrow 4\gamma p"> at GlueX.
Q-factors is an event-by-event multivariate sideband subtraction technique. The only requirement is the knowledge of the signal and background distribution of some discriminating variable. 
1. First the nearest neighbors, under some set of phase space variables, is found for a given event
2. Distribution of the discriminating variable is filled with the nearest neighbors
3. Fit the above distribution with the known/assumed signal and bkg distribution
4. Calculate Q-factor as the signal fraction

There are a bit of caveats which this package tries to tackle. 
- How to include Accidental subtraction and sideband subtraction (if we wanted to)? Should we apply it before the calculation of the q-factors or use the accidental weights during the calculation
- Standardization of the phase space variables is needed. Range and standard deviation standardization is implemented
- Convergence of a fit requires a decent initialization. Currently, the full distribution of the discriminating variable is fitted and the relevant parameters are scaled down to match the number of neighbors
- How does the initialization of the fit function change the results? Original paper suggests doing 3 fits per iteration where there is 100% bkg, 50/50% bkg/sig, and 100% signal. Implemented.
- For each combo in an entire dataset, is there a criteria for what it's potential neighbors can be? For each event, there can be multiple combinations that pass our selection criteria. How should we deal with these inter-event combinations? Currently there is a implementation that looks for spectroscopically unique combos only. Inter-event combinations cannot use the same 4 photons.
- What type of fit functions can one use? Function defintions for Gaussian, double-Gaussian, Breit Wigner has been defined for the signal distribution. The bkg distribution is currently a polynomial. The current logic goes like this. The efficiency of the detector should be relatively flat if the nearest neighbors all cluster closely in some specific region of the detector. Under this scenario the signal distribution should take more of a Gaussian shape. That would really depend on how large the sample size is and what the phase space variables are.
- The detector that the particle was detected in also matters. There is resolution differences between a pair of photons detected in the BCAL, in the FCAL, and where one photon ends up in the BCAL and the other ends up in the FCAL. Currently, the data is separated in the DSelector stage and run independently. If we can do this, it is greatly beneficial since Q-factor analysis is at least N^2 in complexity.
- There is a lot of statistics at GlueX, especially in the channel I am interested in. Currently, Q-factors calcuation is multithreaded (require Minuit2 to be thread safe). Aggregating the results in done by another program, makeDiagnosticHists to allow this parallelization. Future works might be in expanding this to HPC job submission and or include GPU compatibiilty. I have begun looking into accelerating the neighbor finding (most time consuming) using CUDA, have a working version but need to profile the performance to see what the CPU to GPU ratio shoud be. 
- What is best way to find the nearest neighbors. Currently, a priority queue is set up which has O(logn) to insert and to remove the smallest element. Offloading this to a GPU could be beneficial if the sequentail part for the CPU is relatively long.
- How do we know which phase space variables to select? Currently, the program is able to take in a set of variables and just scan through all possible permutations and subsets. Left to the analyzer to decide from there.

## Requirements
Minuit is not threadsafe so this code requires TMinuit2. I believe the ROOT version needed for this is atleast 6.19.
All combinatons from all events that passed your selections should be kept in a root file. Each combination is an entry. At GlueX, this is the flat tree format from the DSelector.  

## Code
There are a few programs working together:
1. getInitParams.C does the initial fit to the full discriminating variable's distribution. Save the fit parameters to a file which is then read in later
2. main.h/C is the Q-factor program. QFactor class is in main.h and there are bunch of helper functions in helperFunc.C. The inputs are not entirely decoupled so some checks should be done to make sure there is consistency of setting in between all the files. QFactorAnalysis class has multiple methods which load the data and fitted parameters from getInitParams and spawns the threads to do the analysis
3. makeDiagnosticHists aggregates all the results from the main program. Various plots are made with the q-factor weighting
4. run.py drives the main program and makeDiagnosticHists. Most of the important variables are configued here and modifies main directly using search and replace, which can be slightly dangerous. run.py also modifies makeDiagnosticHists in a similar way to keep things consistent. main is also compiled at this stage, directories are cleaned to reduce confusion. 
5. convertROOTtoPNG.C is needed since outputting histograms into image formats, in a multithreaded way, causes errors. Maybe due to some blocking issues. Anyways, we can save them all as root files and then convert them to pngs or whatever

Things you need to configure:
1. Most importantly the settings at the top of the run.py file needs to configured properly.
2. The distribution for the signal and bkg must be given. "main" program uses fitFunc, background, and signal function definitions defined in helperFunc.h. Typically a polynomial bkg is taken with a Gaussian signal.  
3. The fit parameters and the fit initializations need to be set in getInitParams and main.h
