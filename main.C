#include "main.h"

void QFactorAnalysis::loadTree(string rootFileLoc, string rootTreeName){
	cout << "Loading root file and tree" << endl;
	TFile* dataFile=new TFile((rootFileLoc).c_str());
	dataFile->GetObject((rootTreeName).c_str(),dataTree);
	// Get the total number of entries and potentially overwrite it if we want to have a shorter run
	total_nentries = (Long64_t)dataTree->GetEntries();
	if (!override_nentries){
		nentries=dataTree->GetEntries();
	}
	cout << "Chosen Total Entries: " << nentries << endl;
}

void QFactorAnalysis::loadFitParameters(string fitLocation){
	cout << "Loading the fit parameters" << endl;
	double eventRatioSigToBkg = 1;
	// -----------------------------------------------------
	// -----------------------------------------------------
	//                         LOAD IN THE FITTED PARAMETERS TO THE FULL DISTRIBUTION
	// -----------------------------------------------------
	// -----------------------------------------------------
        // Import all the fitted values from getInitParams
	string varName;
	double varVal;
	ifstream inFile;
	inFile.open(fitLocation.c_str());
        std::vector<string> initParNames;
        cout << "RootFile(treeName)(fileTag): " << rootFileLoc << "(" << rootTreeName << ")(" << fileTag << ")" << endl;
	while (inFile >> varName >> varVal){
            if (varName.at(0) != '#'){ //not all output from getInitParams is useful. They start with a # 
                initPars.push_back(varVal);
                initParNames.push_back(varName);
                cout << varName << ": " << varVal << endl;
            }
            if (varName == "#eventRatioSigToBkg"){
                eventRatioSigToBkg = varVal;
            }
	}
	cout << "eventRatioSigToBkg: " << eventRatioSigToBkg << endl;
	inFile.close();

	// ----------------------------------------
	// Set the fit constants for eta or pi0 full fits
	// ----------------------------------------
        // We need to have a global scale factor since in the fit to the full distribution of the discriminating variable there are just more events
        // If we use a normalized Gaussian then the amplitude is simply scaled by the kDim/total
	double scaleFactor = (double)kDim/total_nentries;
	cout << "scaleFactor: " << scaleFactor << endl;
	
        // scaling signal distribution
        for ( auto iVar : sigVarsNeedScaling) {
            initPars[iVar] = scaleFactor*initPars[iVar];
        }
        // scaling bkg distribution
        for ( auto iVar : bkgVarsNeedScaling) {
            initPars[iVar] = scaleFactor*initPars[iVar];
        }
        cout << "Scaled initialization parameters, scaled down to kDim" << endl;
        for ( int iVar=0; iVar<numDOFbkg+numDOFsig; ++iVar ){
            cout << initParNames[iVar] << ": " << initPars[iVar] << endl; 
        }

        // We will 3 iterations. Not sure if I am doing this correctly
        // The goal would be to use 100 bkg, 50/50, 100% signal. We can scale the amplitudes by a certain factor related to eventRatioSigToBkg
	cout << "------------- SETTING UP | 100% | 50%/50% | 100% | parameter initialization -----------------" << endl;
        redistributeFactorSig = { 0, (1+1/eventRatioSigToBkg)/2, (1+1/eventRatioSigToBkg) };
        redistributeFactorBkg = { (1+eventRatioSigToBkg), (1+eventRatioSigToBkg)/2, 0 };

        // we can set up the parameter limits now
        parLimits.initPars = initPars;
        parLimits.eventRatioSigToBkg = eventRatioSigToBkg;
        parLimits.setupParLimits();
        parLimits.printParLimits();
}


void QFactorAnalysis::loadData(){
	cout << "Loading the data into arrays" << endl;
	// -----------------------------------------------------
	// -----------------------------------------------------
	//                         LOAD IN THE DATA
	// -----------------------------------------------------
	// -----------------------------------------------------
	// Create variables to hold the data as we read in the data from the tree
	double discrimVar;
        double phaseSpaceVar[dim];
        double sideBandVar;
	double AccWeight;
	ULong64_t spectroscopicComboID;

        // vars we will use to fill but not use directly
	ULong64_t eventNumber;
	double uniqueComboID;
	bool isUniqueEtaB;
	bool isUniquePi0B;
	bool isUniquePi0EtaB;
        double utWeight;

	// Set branch addresses so we can read in the data
        dataTree->SetBranchAddress(s_discrimVar.c_str(), &discrimVar);
	parseVarString parse(varString);
	parse.parseString();
        if ( parse.varStringSet.size() != dim ) { cout << "Uh-oh something went wrong. varString size not the same as dim" << endl; }
	for (int iVar=0; iVar<dim; ++iVar){
            cout << "Setting " << parse.varStringSet[iVar] << " to phaseSpaceVar index " << iVar << endl; 
            dataTree->SetBranchAddress(parse.varStringSet[iVar].c_str(),&phaseSpaceVar[iVar]);
        }
	dataTree->SetBranchAddress(s_discrimVar.c_str(), &discrimVar);
	dataTree->SetBranchAddress(s_sideBandVar.c_str(), &sideBandVar);
	dataTree->SetBranchAddress(s_accWeight.c_str(),&AccWeight);
	dataTree->SetBranchAddress("spectroscopicComboID",&spectroscopicComboID);
        
        // vars we will use to fill but not use directly
	dataTree->SetBranchAddress("event",&eventNumber);
	dataTree->SetBranchAddress("uniqueComboID",&uniqueComboID);
	dataTree->SetBranchAddress("isNotRepeated_eta",&isUniqueEtaB);
	dataTree->SetBranchAddress("isNotRepeated_pi0",&isUniquePi0B);
	dataTree->SetBranchAddress("isNotRepeated_pi0eta",&isUniquePi0EtaB);
        dataTree->SetBranchAddress("uniqunessTrackingWeights",&utWeight);

        double sbWeight;
	
	// We will use a ientry to keep track of which entries we will get from the tree. We will simply use ientry when filling the arrays.  
	for (Long64_t ientry=0; ientry<nentries; ientry++)
	{
		dataTree->GetEntry(ientry);
		discrimVars.push_back(discrimVar);
	        for (int iVar=0; iVar<dim; ++iVar){
                    phaseSpaceVars[iVar].push_back(phaseSpaceVar[iVar]);
                }
                sideBandVars.push_back(sideBandVar);
	        AccWeights.push_back(AccWeight);
                utWeights.push_back(utWeight);
                getSBWeight(sideBandVar,&sbWeight,weightingScheme);
		sbWeights.push_back(sbWeight);
	        spectroscopicComboIDs.push_back(spectroscopicComboID);
	}

	if ( verbose_outputDistCalc ) {
	    cout << "Before standarization the first nentries of phaseSpaceVar[0]" << endl;
	    for ( int ientry=0 ; ientry < nentries; ientry++){
	        cout << phaseSpaceVars[0][ientry] << endl;
	        //cout << cosTheta_X_cms[ientry] << endl;//"," << phi_X_cms[ientry] << endl;
	    }
	}

	// -----------------------------------------------------
	// -----------------------------------------------------
	//                         STANDARDIZE THE DATA
        //                         Range vs StdDev standardization
        //                         Not much difference so far but the option is left here
	// -----------------------------------------------------
	// -----------------------------------------------------
	//
	standardizeArray standarizationClass;
	for (int iVar=0; iVar<dim; ++iVar){
            phaseSpaceVars[iVar].push_back(phaseSpaceVar[iVar]);

            if(standardizationType=="range"){
                standarizationClass.rangeStandardization(phaseSpaceVars[iVar],nentries);
            }
            else if (standardizationType=="std"){
                standarizationClass.stdevStandardization(phaseSpaceVars[iVar],nentries);
            }
        }

	if ( verbose_outputDistCalc ) {
	    cout << "After standarization the first nentries of phaseSpaceVar[0]" << endl;
	    for ( int ientry=0 ; ientry < nentries; ientry++){
	        cout << phaseSpaceVars[0][ientry] << endl;
	        //cout << cosTheta_X_cms[ientry] << endl;//"," << phi_X_cms[ientry] << endl;
	    }
	}

	
	// phasePoint1 will consider all events from lowest to largest since these will be our attached q values. phasePoint2 on the other hand will only look at a subset of the events where the
	// elements must be spectroscopically distinct, i.e. the 4 photons in consideration are different. Not sure if this is the value I should consider or is it better to do a pair of maps
	// where I am tracking the two photon pairs that make up the eta and pi0.         
        if(nRndRepSubset!=0 && nRndRepSubset<nentries){
            std::set<Int_t> rndRepSubset; //probably should make the subset have unique elements. So sampling without replacement
            srand (time(NULL)); // set a random seed
            while(rndRepSubset.size() < nRndRepSubset){
                rndRepSubset.insert(rand()%nentries);
            }
            phasePoint2PotentialNeighbor.assign(rndRepSubset.begin(),rndRepSubset.end());
        }
        else {
	    set<Int_t> setUsedSpectroscopicIDs;
	    for (Int_t ientry=0; ientry<nentries; ientry++){ 
                if (s_uniquenessTracking=="default"){
	            if ( setUsedSpectroscopicIDs.find( spectroscopicComboIDs[ientry] ) == setUsedSpectroscopicIDs.end() ) {
	                setUsedSpectroscopicIDs.insert( spectroscopicComboIDs[ientry] );
	                phasePoint2PotentialNeighbor.push_back(ientry);
	            }
                }
                else if (s_uniquenessTracking=="weighted"){
                    // if we use the weighted setting we will allow any neighbor to be possible
	            phasePoint2PotentialNeighbor.push_back(ientry);
                }
                else{ 
                    cout << "uniqueness tracking setting not set up properly!" << endl;
                    exit(0);
                }
	    }
        }
        cout << "\n\n--------------------------------------" << endl;
	cout << phasePoint2PotentialNeighbor.size() << "/" << nentries << " are used as potential neighbors" << endl;
        cout << "Some indicies that will be used as potential neighbros: " << endl;
        //for (auto iPhasePoint2=0; iPhasePoint2 < phasePoint2PotentialNeighbor.size(); ++iPhasePoint2){
        for (auto iPhasePoint2=0; iPhasePoint2 < 20; ++iPhasePoint2){
            cout << phasePoint2PotentialNeighbor[iPhasePoint2] << " ";
        }
        cout << endl;

}

// This will be a thread that is totally independent and will be spawned in the main process after we have loaded all the data
void QFactorAnalysis::runQFactorThreaded(){
	// make sure the global variables are read in correctly
	cout << "kDim: " << kDim << endl;
	cout << "numberEventsToSavePerProcess: " << numberEventsToSavePerProcess << endl;
	cout << "nProcess: " << nProcess << endl;
	cout << "seedShift: " << seedShift << endl;
	cout << "nentries: " << nentries << endl;
	cout << "override_nentries: " << override_nentries << endl;
	cout << "verbose: " << verbose << endl;

	ROOT::EnableThreadSafety();
	// [=] refers to a capture list which is used by this lambda expression. The lambda gets a copy of all the local variables that it uses when it is created. If we
	// just use [] we will get an error since the lambda will have no idea what these variables are	
	auto f = [=](int iProcess){
		int iFit = 0;
		int iThread = iProcess;
		// ----------------------------------------
		// Open up a root file to save the q-factors and other diagnostics to it
		// ---------------------------------------
                // Initializing some variables we can track during the q-value extraction
        	double chiSq;
        	double chiSq_pi0;
        	ULong64_t flatEntryNumber;
        	double bestChiSq;
		double worstChiSq;
        	double best_qvalue;
                double qvalueBS_std=0;

                // Saving the results along with some diagnostics
		TBranch* b_sbWeight;
		TBranch* b_flatEntryNumber;
        	TFile *resultsFile = new TFile(("logs/"+fileTag+"/results"+to_string(iThread)+".root").c_str(),"RECREATE");
        	TTree* resultsTree = new TTree("resultsTree","results");
        	resultsTree->Branch("flatEntryNumber",&flatEntryNumber,"flatEntryNumber/l");
        	resultsTree->Branch("qvalue",&best_qvalue,"qvalue/D");
                resultsTree->Branch("qvalueBS_std",&qvalueBS_std,"qvalueBS_std/D");
        	resultsTree->Branch("bestChiSq",&bestChiSq,"bestChiSq/D");
        	resultsTree->Branch("worstChiSq",&worstChiSq,"worstChiSq/D");
        	//resultsTree->Branch("combostd",&comboStd,"combostd/D");
        	//resultsTree->Branch("combostd2",&comboStd2,"combostd2/D");
        	cout << "Set up branch addresses" << endl;

		// Define some needed variables like canvases, histograms, and legends
        	auto legend_init = new TLegend(0.1,0.7,0.4,0.9);
        	auto legend_fit = new TLegend(0.1,0.7,0.4,0.9);
        	auto legend_qVal = new TLegend(0.1,0.7,0.4,0.9);

    		TCanvas *allCanvases = new TCanvas(("anyHists"+to_string(iThread)).c_str(),"",1440,900);
    		TCanvas *allCanvases_badFit = new TCanvas(("anyHists_badFit"+to_string(iThread)).c_str(),"",1440,900);
		cout << "Creating canvas " << iThread << endl;
        	TLine* discrimVarLine;
		TLine* qSigLine;
		TLine* qBkgLine;
                TLine* qValLine;
		TH1F* discriminatorHist;
		//TH1F* discriminatorHist2;

		// defining some variables we will use in the main loop to get the distances and then the q-values
		map<double, int> mapDistToEntry;
		set<double> distances;
		double phasePoint1[dim];
		double phasePoint2[dim];
		double distance;
        	double qvalue;

        	distSort_kNN distKNN(kDim);
        	pair<double,int> newPair;

		// opening a file to write my log data to
    		ofstream logFile;
    		logFile.open(("logs/"+fileTag+"/processLog"+to_string(iThread)+".txt").c_str());
		
		// Determine what events each thread should run
		int batchEntries = (int)nentries/nProcess; // batchEntries the size of the batch
		int lowest_nentry = iThread*batchEntries;
		int largest_nentry;
        	if (iThread!=(nProcess-1)) {
        	    largest_nentry  = (iThread+1)*batchEntries;
        	}
        	else {
        	    largest_nentry = nentries; 
        	}
		cout << "nentries we will use for this process: " << lowest_nentry << ", " << largest_nentry << endl;

		// randomly select some events to write histograms for 
		set<int> selectRandomIdxToSave;
		int randomEvent;
		srand(iThread+seedShift);
		for (int i=0; i<numberEventsToSavePerProcess; i++){
			// basically randomly sample a uniform number between (lowest_nentry, largest_nentry)
			randomEvent = rand() % (int)batchEntries; // batchEntries is the size of the batch
			randomEvent += lowest_nentry; // shift by the lowest entry of the batch
			selectRandomIdxToSave.insert( randomEvent );
		}
        	cout << "randomly selected some events to save" << endl;


		std::vector<double> binRangeEta;
        	std::vector<double> fitRangeEta;

		binRangeEta={200,0.25,0.85};
		fitRangeEta={0.3,0.8};

        	// Getting the scaling of the flat bkg is a bit tricky. Have to count how much bins we have in our fit range. kDim divided by the bins in fit range is the height of the flat function.
        	double numBinsInFit = (fitRangeEta[1]-fitRangeEta[0])/((binRangeEta[2]-binRangeEta[1])/binRangeEta[0]);
        	double binSize=((binRangeEta[2]-binRangeEta[1])/binRangeEta[0]);

        	if (verbose) {cout << "Eta hist range: " << binRangeEta[0] << ", " << binRangeEta[1] << ", " << binRangeEta[2] << endl;}
        	if (verbose) {cout << "Eta fit range: " << fitRangeEta[0] << ", " << fitRangeEta[1] << endl;}
        	
                int numFits;
                if(redistributeBkgSigFits){ 
                    numFits=3;
                }
                else{
                    numFits=1;
                }
		TF1* fit = new TF1(("fit"+to_string(iThread)).c_str(),fitFunc,fitRangeEta[0],fitRangeEta[1],numDOFbkg+numDOFsig);
		TF1* bkgFit = new TF1(("bkgFit"+to_string(iThread)).c_str(),background,fitRangeEta[0],fitRangeEta[1],numDOFbkg);
		TF1* sigFit = new TF1(("sigFit"+to_string(iThread)).c_str(),signal,fitRangeEta[0],fitRangeEta[1],numDOFsig);
		discriminatorHist = new TH1F(("discriminatorHist"+to_string(iThread)).c_str(),"",binRangeEta[0],binRangeEta[1],binRangeEta[2]);
                TH1F* dHist_qvaluesBS = new TH1F(("qvaluesBS"+to_string(iThread)).c_str(),"Bootstrapped Q-Factors",100,0,1);
                std::vector<TF1*> initFits;
                TF1* initFit;
		for ( int iFit = 0; iFit < numFits; ++iFit){
                    initFit = new TF1(("initFit"+to_string(iFit)+to_string(iThread)).c_str(),fitFunc,fitRangeEta[0],fitRangeEta[1],numDOFbkg+numDOFsig);
                    initFits.push_back(initFit);
                }
		// the main loop where we loop through all events in a double for loop to calculate dij. Iterating through all j we can find the k nearest neighbors to event i.
		// Then we can plot the k nearest neighbors in the discriminating distribution which happens to be a double gaussian and flat bkg. Calculate Q-Value from event
		// i's discriminating variable's value. This value will be plugged into the signal PDF and the total PDF. The ratio of these two values are taken which is the Q-Value.
		//logFile << std::fixed << std::setprecision(6);
		int randomEntry;
		int saveN_badEvents=1;
		int savedN_badEvents=0;
        	for (int ientry=lowest_nentry; ientry<largest_nentry; ientry++){ 
                        TFile *qHistsFile;
        		if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end()) {
        	            qHistsFile = new TFile(("histograms/"+fileTag+"/qValueHists_"+to_string(ientry)+".root").c_str(),"RECREATE");
                        }
                        dHist_qvaluesBS->Reset();
		        Double_t par[numDOFbkg+numDOFsig]; 
		        Double_t parBest[numDOFbkg+numDOFsig]; // saving the best parameters to plot them later

                        // Outputting the progress of each thread
			double unshiftedEntry = (double)(ientry-lowest_nentry);
                        if(batchEntries<30){ cout << "Due to how we output the progress of the threads, we must have > 20 nentries PER nProcess\nEXITING\nEXITING" << endl; exit(0);}
			int percentOfBatchEntries = (int)(batchEntries/20);
			if ( (ientry-lowest_nentry) % percentOfBatchEntries == 0) { 
				cout << "(Process " << iProcess << ") Percent done: " << (int)round( unshiftedEntry/(largest_nentry-lowest_nentry)*100 ) << "%" << endl;
		       	}

			flatEntryNumber=ientry;
			auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
			auto duration_beginEvent = std::chrono::high_resolution_clock::now();
			if(verbose) { logFile << "Starting event " << ientry << "/" << largest_nentry << " ---- Global Time: " << duration2 << "ms" << endl; }
			
			// clear data from previous events
			mapDistToEntry.clear();
			distances.clear();
			
			for ( int iVar=0; iVar<dim; ++iVar ){
				phasePoint1[iVar] = phaseSpaceVars[iVar][ientry];
			}

                        vector<double> qvalues; qvalues.reserve(nBS);
                        int nPotentialNeighbors=(int)phasePoint2PotentialNeighbor.size();
                        vector<int> phasePoint2PotentailNeighbor_BS;
                        phasePoint2PotentailNeighbor_BS.reserve(nPotentialNeighbors);
                        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); 
                        // RNG using generators and distributions are really cheap. Is threadsafe also whereas rand() is not
                        // https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers
                        // Mersenne twistor requires large state storage but has the highest quality + period. 
                        // Lagged Fibonoci Generator liek ranlux24 has a much smaller period but should be much faster
                        // https://books.google.com/books?id=zi08DQAAQBAJ&pg=PA237&lpg=PA237&dq=ranlux24+%22period%22&source=bl&ots=N4HXolKGDr&sig=ACfU3U1JfZcQczw-xoAv-lyHTOCOf4Xpjw&hl=en&sa=X&ved=2ahUKEwiosZP8-_npAhW2RTABHatDDe8Q6AEwAnoECAsQAQ#v=onepage&q=ranlux24%20%22period%22&f=false
                        // mt19937 has a period around 10^6000 whereas ranlux24 has a period around 10^171. ranlux24 produces 24 bit integers. 24 bits ~ 17M which is much larger than the number of potential neighbors
                        std::uniform_int_distribution<int> distribution(0,nPotentialNeighbors-1); //if we ask for a range of 0 to 10 it will include 10. phasePoint2PotentialNeighbor has a zero index. Need to subtract by 1
                        //static thread_local std::mt19937 generator(seed); 
                        static thread_local std::ranlux24_base generator(seed); 

                        // ---------------------
                        // NOW FIND NEIGHBORS AND EXTRACT Q-FACTORS.
                        // ---------------------
                        for (int iBS=0; iBS<nBS+1; ++iBS){ 
                            // resetting some variables
                            discriminatorHist->Reset();
			    bestChiSq=DBL_MAX;
			    worstChiSq=DBL_MIN;
			    duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
			    if(verbose){logFile << "\tBegin bootstrapping potential neighbors: " << duration2 << "ms" << endl; }
                            phasePoint2PotentailNeighbor_BS.clear();
                            if (iBS==nBS){
                                phasePoint2PotentailNeighbor_BS = phasePoint2PotentialNeighbor;
                            }
                            else{ // resampling neighbors with replacement
                                for(int iNeighbor=0; iNeighbor<nPotentialNeighbors; ++iNeighbor){
                                    phasePoint2PotentailNeighbor_BS.push_back(phasePoint2PotentialNeighbor[distribution(generator)]);
                                }
                            }
			    duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
			    if(verbose){logFile << "\tBegin finding neighbors: " << duration2 << "ms" << endl; }
                            // What if we wanted to look at k random neighbors?
                            if (doKRandomNeighbors){
			        for (int jentry=0; jentry<kDim;++jentry) {  
			              randomEntry = rand() % nentries;
			              distKNN.insertPair(make_pair(1, randomEntry) ); //just using 1 as a distance. Doesnt matter anyways
			        }
                            }
                            else {
			        for (int jentry : phasePoint2PotentailNeighbor_BS) {  
			        	if ( verbose_outputDistCalc ) { cout << "event i,j = " << ientry << "," << jentry << endl;} 
			        	for ( int iVar=0; iVar<dim; ++iVar ){
			        	   	phasePoint2[iVar] = phaseSpaceVars[iVar][jentry];
			        	}
			        	if (spectroscopicComboIDs[jentry] != spectroscopicComboIDs[ientry]){
			        	        distance = calc_distance(dim,phasePoint1,phasePoint2,verbose_outputDistCalc);
			        	        distKNN.insertPair(make_pair(distance,jentry));
			        	}
			        }
                            }
			    duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
			    if(verbose){logFile << "\tFound neighbors: " << duration2 << "ms" << endl; }
			    if (distKNN.kNN.size() != kDim){ cout << "size of distKNN is not equal to kDim! size,kDim="<< distKNN.kNN.size() << "," << kDim 
			        << "\n    -- if size is 1 less than kDim it is probably because kDim=nentries and event i cannot be a neighbor to itself" << 
			        "\n    -- if size != kDim it could also mean that the number of spectroscopically unique neighbors reduces the number of poential neighbors below kDim" << endl;}
                            if(verbose_outputDistCalc){ cout << "These are our neighbors" << endl; }
			    while ( distKNN.kNN.empty() == false ){
			            newPair = distKNN.kNN.top();
			            distKNN.kNN.pop();
			    	    // ATLEST FOR NOW I WILL NOT WORRY ABOUT TRACKING THE UNIQUE ETA COMBOS WHEN FILLING HERE BECAUSE THAT MIGHT NOT BE GOOD IN THE FIRST PLACE
			    	    // WHEN WE FILL THE HISTOGRAMS IN THE END WE CAN DO THIS BUT IF WE USE IT NOW IT WHEN FILLING THE NEIGHBORS HISTOGRAM IT MIGHT BE BAD?
			    	    // IF ANYTHING WE SHOULD JUST SKIP THESE NON UNIQUE COMBINATIONS TO SAVE TIME (IF WE WERE ONLY PLOTTING M(ETA) BUT SINCE WE PLOT ALL THE DISTRIBUTIONS LIKE
			    	    // M(PI0ETA) WE CANT DO THIS
                                    double weight;
                                    if (weightingScheme==""){ weight=1; }
                                    if (weightingScheme=="as"){ weight=AccWeights[newPair.second]; }
                                    if (weightingScheme=="as*bs"){ weight=AccWeights[newPair.second]*sbWeights[newPair.second]; }
                                    if (s_uniquenessTracking=="weighted"){
                                        weight=weight*utWeights[newPair.second];
                                    }
			            discriminatorHist->Fill(discrimVars[newPair.second],weight);

                                    if(verbose_outputDistCalc){
			                cout << "(" << newPair.first << ", " << newPair.second << ", " << discrimVars[newPair.second] << ")" << endl; 
                                    }
			    }
			    
			    duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
			    if(verbose){logFile <<	"\tFilled neighbors: " << duration2 << "ms" << endl;}
			    
			    // /////////////////////////////////////////
			    // Calcuclate q-value
			    // /////////////////////////////////////////
			    //// We use a normalized gaussian and a flat function. 
			    
                            // saving the initializations to see if they sort of make sense
			    for ( int iFit = 0; iFit < numFits; ++iFit){
                                std::vector<double> initParsVaryPercentSig = initPars;
                                if(redistributeBkgSigFits){
                                    for ( auto sigVar : sigVarsNeedScaling ){ 
                                        initParsVaryPercentSig[sigVar] = initParsVaryPercentSig[sigVar]*redistributeFactorSig[iFit];
                                    }
                                    for ( auto bkgVar : bkgVarsNeedScaling ){ 
                                        initParsVaryPercentSig[bkgVar] = initParsVaryPercentSig[bkgVar]*redistributeFactorBkg[iFit];
                                    }
                                }
                                // saving the initialization to plot it later, just checking if the initializations look reasonable
			    	initFits[iFit]->SetParameters(&initParsVaryPercentSig[0]);
                                fit->SetParameters(&initParsVaryPercentSig[0]);
                                for (int iPar=0; iPar<parLimits.numDOF; ++iPar){
                                    fit->SetParLimits(iPar,parLimits.lowerParLimits[iPar], parLimits.upperParLimits[iPar]);
                                }
			        duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
			        if(verbose){logFile <<	"\tinitialized init fits: " << duration2 << "ms" << endl;}
			    	//{
			    	//	// have to use a mutex here or else we get some weird error
			    	//	//R__LOCKGUARD(gGlobalMutex);
			    	//	discriminatorHist->Fit(("fit"+to_string(iThread)).c_str(),"RQBNL"); // B will enforce the bounds, N will be no draw
			    	//}
			    	discriminatorHist->Fit(fit,"RQNL"); // B will enforce the bounds, N will be no draw
			    	//discriminatorHist->Fit(("fit"+to_string(iThread)).c_str(),"RQBNL"); // B will enforce the bounds, N will be no draw
			        duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
			        if(verbose){logFile <<	"\tFitted fits: " << duration2 << "ms" << endl;}
			    	fit->GetParameters(par);
			    	bkgFit->SetParameters(par);
			    	sigFit->SetParameters(&par[numDOFbkg]);
			    	qvalue=sigFit->Eval(discrimVars[ientry])/fit->Eval(discrimVars[ientry]);
			    	
			    	// now that the q-value is found we can get the chiSq and save the parameters with the best chiSq
			    	chiSq = fit->GetChisquare()/(fit->GetNDF());
			    	//if (verbose) { logFile << "\tcurrent ChiSq, best ChiSq: " << chiSq << ", " << bestChiSq << endl; }
			    	if (chiSq < bestChiSq){
			    		best_qvalue = qvalue;
			    		bestChiSq=chiSq;
			    		for (int i=0; i < sizeof(par)/sizeof(Double_t); ++i){
			    			parBest[i]=par[i];
			    		}
			    	} 
			    	if (chiSq > worstChiSq){
			    		worstChiSq = chiSq;
			    	}
			        duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
			        if(verbose){logFile <<	"\tFitted hist with some initialization: " << duration2 << "ms" << endl;}
			    } // close iFit loop
                            if(iBS<nBS){  
                                qvalues.push_back(best_qvalue);
                                dHist_qvaluesBS->Fill(best_qvalue);
                            }
                            if (iBS==nBS || saveBShistsAlso){ //show all histograms or only the true q-factor
        		        // /////////////////////////////////////////
        		        // Drawing histogram 
        		        // /////////////////////////////////////////
        		        if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end()) {
        		                // Here we draw the histograms that were randomly selected
        		                allCanvases->Clear();
        		                legend_init->Clear();
        		                legend_fit->Clear();
                                        legend_qVal->Clear();
                                        TPad *pad1 = new TPad("pad1", "",0.0,0,0.5,1.0);
                                        TPad *pad2 = new TPad("pad2", "",0.5,0.5,1.0,1);
                                        TPad *pad3 = new TPad("pad3", "",0.5,0,1.0,0.5);
                                        pad1->Draw();
                                        pad2->Draw();
                                        pad3->Draw();
        
                                        pad1->cd();
        		        	discriminatorHist->SetTitle(("BEST VALS:  QValue="+to_string(best_qvalue)+"  ChiSq="+to_string(bestChiSq)).c_str() );
        		        	discrimVarLine = new TLine(discrimVars[ientry],0,discrimVars[ientry],kDim);
        		        	discrimVarLine->SetLineColor(kOrange);
        		        	
                	          	fit->SetParameters(parBest);
        		          	bkgFit->SetParameters(parBest);
        		          	sigFit->SetParameters(&parBest[numDOFbkg]);
                	          	fit->SetLineColor(kRed);
          		          	bkgFit->SetFillColor(kMagenta);
                	          	bkgFit->SetLineColor(kMagenta);
          		          	bkgFit->SetFillStyle(3004);
          		          	sigFit->SetFillColor(kBlue);
                	          	sigFit->SetLineColor(kBlue);
          		          	sigFit->SetFillStyle(3005);
        
        		        	qSigLine = new TLine(0,sigFit->Eval(discrimVars[ientry]),binRangeEta[2],sigFit->Eval(discrimVars[ientry]));
        		        	qBkgLine = new TLine(0,bkgFit->Eval(discrimVars[ientry]),binRangeEta[2],bkgFit->Eval(discrimVars[ientry]));
        		        	
                	          	qSigLine->SetLineColor(kBlue);
        		        	qSigLine->SetLineStyle(9);
                	          	qBkgLine->SetLineColor(kMagenta);
        		        	qBkgLine->SetLineStyle(9);
        
        
                                        legend_fit->AddEntry(fit,"Tot");
                                        legend_fit->AddEntry(bkgFit,"Bkg");
                                        legend_fit->AddEntry(sigFit,"Sig");
                                        legend_fit->AddEntry(discrimVarLine,"DiscrimVarVal");
                                        legend_fit->AddEntry(qSigLine,"SigVal @ DiscrimVarVal");
                                        legend_fit->AddEntry(qBkgLine,"BkgVal @ DiscrimVarVal");
        
                	          	discriminatorHist->Draw();
                	          	discrimVarLine->Draw("same");
        
                	          	fit->Draw("SAME");
                	          	bkgFit->Draw("SAME FC");
                	          	sigFit->Draw("SAME FC");
                	        	qBkgLine->Draw("SAME");
                	        	qSigLine->Draw("SAME");
                                        legend_fit->Draw();
                	          	drawText(parBest,numDOFbkg+numDOFsig,"par",sigFit->Eval(discrimVars[ientry]),bkgFit->Eval(discrimVars[ientry]),fit->Eval(discrimVars[ientry]));
                	        	// INTERESTING, IF I WERE TO SAVE THE CANVAS AS A ROOT FILE I GET AN ERROR IF I PUT THE SAME HISTOGRAM ON TWO DIFFERENT PADS. SEEMS LIKE THE CANVAS
                	        	// SAVES A TList OF HISTS+TF1'S AND IF THERE ARE MULTIPLE CALLS TO THE SAME HISTOGRAM IT MIGHT DELETE THE HISTOGRAM AFTER SEEING IT FOR THE FIRST TIME AND
                	        	// THEN IT WOULD NOT BE ABLE TO FIND THE HISTOGRAM AGAIN THE SECOND TIME AROUND. WE HAVE TO CLONE THE HISTOGRAM FIRST AND THEN SAVE THE ROOT FILE SO THE CANVAS
                	        	// ARE TWO DIFFERENT ELEMENTS.
                                        pad2->cd();
                	        	TH1F* clonedHist = (TH1F*) discriminatorHist->Clone();
                                        clonedHist->SetTitle("Initializations");
                	          	clonedHist->Draw();
        
                                        if (redistributeBkgSigFits){
                	        	    initFits[0]->SetLineColor(kBlue-4);
                	        	    initFits[0]->Draw("SAME");
                	        	    initFits[1]->SetLineColor(kRed-3);
                	        	    initFits[1]->Draw("SAME");
                	        	    initFits[2]->SetLineColor(kGreen+2);
                	        	    initFits[2]->Draw("SAME");
                                            legend_init->AddEntry(initFits[0],"100% bkg");
                                            legend_init->AddEntry(initFits[1],"50/50 bkg/sig");
                                            legend_init->AddEntry(initFits[2],"100% sig");
                                        }
                	        	TF1* initFit4 = new TF1(("initFit4"+to_string(iThread)).c_str(),fitFunc,fitRangeEta[0],fitRangeEta[1],numDOFbkg+numDOFsig);
                	        	initFit4->SetParameters(initPars[0],initPars[1],initPars[2],initPars[3],initPars[4]);
                	        	initFit4->SetLineColor(kOrange+1);
                	        	initFit4->Draw("SAME");
                                        legend_init->AddEntry(initFit4,"Scaled Full Fit");
                                        legend_init->Draw();
        
                                        pad3->cd();
                                        if (iBS==nBS){
                                            qvalueBS_std=calculateStd(nBS,&qvalues[0]);
                                            dHist_qvaluesBS->SetTitle(("STD: "+to_string(qvalueBS_std)).c_str());
                                        }
                                        else{
                                            dHist_qvaluesBS->SetTitle("Bootstrapped Q-Factors");

                                        }
                                        dHist_qvaluesBS->Draw();
        		                qValLine = new TLine(best_qvalue,0,best_qvalue,nBS);
        		                qValLine->SetLineColor(kOrange);
                                        legend_qVal->AddEntry(qValLine,"True Q Value");
                                        qValLine->Draw("SAME");
                                        legend_qVal->AddEntry(dHist_qvaluesBS,"Bootstrapped Q Values");
                                        legend_qVal->Draw();
        		        	// need to save as a root file first then convert to pngs or whatever. Seems like saveas doesnt like threaded since the processes might make only one png converter
        		        	// or whatever and maybe if multiple threads calls it then a blocking effect can happen
                                        qHistsFile->cd();
                                        if(iBS==nBS){
        		                    //allCanvases->SaveAs(("histograms/"+fileTag+"/Mass-event"+std::to_string(ientry)+".root").c_str());
                                            allCanvases->Write(("Mass-event"+std::to_string(ientry)).c_str());
                                            qHistsFile->Close();
                                        }
                                        else {
        		        	        //allCanvases->SaveAs(("histograms/"+fileTag+"/Mass-event"+std::to_string(ientry)+"BS"+to_string(iBS)+".root").c_str());
                                                allCanvases->Write(("Mass-event"+std::to_string(ientry)+"BS"+to_string(iBS)).c_str());
                                        }
        		        	    
        		        	if(verbose){
        		        	     duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
        		        	     logFile << "\tSaved this histogram since it was randomly selected: " << duration2 <<  "ms" << endl;
        		        	}
			        } 
                            } // closes condition to save histograms
                        } // finishes nBS loop
                        if(numFits>1){
			    if(verbose){logFile << "\tDelta b/w best and worst chiSq = " << to_string(bestChiSq-worstChiSq) << ": " << duration2 << "ms" << endl; }
                        }
                        else{
			    if(verbose){logFile << "\tBest chiSq = " << to_string(bestChiSq) << ": " << duration2 << "ms" << endl; }
                        }
			double sbWeight = sbWeights[ientry];
			resultsTree->Fill();
			
		}
		delete fit;
		delete bkgFit;
		delete sigFit;
        	resultsFile->cd();
        	resultsTree->Write();
        	cout << "nentries: " << nentries << endl;
        	
		// Finish the log files by including an elapsed time and finally closing the file
		if (verbose){
			auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
			logFile << "Total time: " << duration2 << " ms" << endl;
			logFile << "Time Per Event: " << duration2/nentries  << " ms" << endl;
		}
		logFile.close();
	};


        // Now that we have the lambda function we can start to spawn threads
	cout << "Launching " << nProcess << " threads 1 second apart!" << endl;
	vector<thread> threads;
	for ( int iThread=0; iThread<nProcess; ++iThread){
		cout << "(Thread " << iThread << ") is starting" << endl;
		threads.emplace_back( [f, iThread] { f(iThread); } );
                sleep(1);
		//threads[iThread] = std::thread(QFactorAnalysis::staticEntryPoint, this, iThread);
	}
	for (auto&& t : threads) t.join(); // join waits for completion
	cout << "Threads have completed running!" << endl;
}

int main( int argc, char* argv[] ){
	ROOT::EnableThreadSafety();
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); // For thread safety we need this
	TH1::AddDirectory(kFALSE);

	// This suppresses all the "info" messages from roofit. Like saying that it will use numeric integrator for the generic pdf we define
	// https://root-forum.cern.ch/t/suppressing-info-messages/14642/6
	//RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);

	// -----------------------------------------------------
	// -----------------------------------------------------
	//                          PARSE THE COMMAND LINE ARGUMENTS
	// -----------------------------------------------------
	// -----------------------------------------------------
        int kDim= atoi(argv[1]);
        std::string varString=argv[2];
        std::string standardizationType=argv[3];
        std::string fitLocation=argv[4];
        bool redistributeBkgSigFits;
        if(atoi(argv[5])==1){
           redistributeBkgSigFits=true; 
        }
        else{
            redistributeBkgSigFits=false;
        }
        bool doKRandomNeighbors;
        if(atoi(argv[6])==1){
           doKRandomNeighbors=true; 
        }
        else{
            doKRandomNeighbors=false;
        }
        int numberEventsToSavePerProcess=atoi(argv[7]);
	int nProcess=atoi(argv[8]);
        int seedShift= atoi(argv[9]);
        Long64_t nentries= atoi(argv[10]);
        int nRndRepSubset=atoi(argv[11]);
        int nBS=atoi(argv[12]);
        bool saveBShistsAlso;
        if ( atoi(argv[13]) ==1 ){ saveBShistsAlso=true; }
        else { saveBShistsAlso=false; }
        bool override_nentries;
        cout << "Num Vars="<< dim << endl;
        if ( atoi(argv[14]) ==1 ){ override_nentries=true; }
        else{ override_nentries=false; }
        bool verbose;
        if ( atoi(argv[15])==1 ){ verbose=true;}
        else{ verbose=false;}
        cout << "----------------------------" << endl;
        cout << "kDim: " << kDim << endl;
        cout << "nProcess: " << nProcess << endl;
        cout << "numberEventsToSavePerProcess: " << numberEventsToSavePerProcess << endl;
        cout << "seedShift: " << seedShift << endl;
        cout << "nentries: " << nentries << endl; 
        cout << "nRndRepSubset: " << nRndRepSubset << endl;
        cout << "nBS: " << nBS << endl;
        cout << "saveBShistsAlso: " << saveBShistsAlso << endl;
        cout << "override_nentries: " << override_nentries << endl;
	cout << "varString: " << varString << endl; 
        cout << "standardizationType: " << standardizationType << endl;
        cout << "fitLocation: " << fitLocation << endl;
        cout << "verbose: " << verbose  << endl; 
        cout << "redistributeBkgSigFits: " << redistributeBkgSigFits << endl;
        cout << "doKRandomNeighbors: " << doKRandomNeighbors << endl;
        cout << "----------------------------" << endl;
        //cout << "Sleeping for 10 seconds so you can look at these settings" << endl;
        //sleep(10);
    
	QFactorAnalysis analysisControl(kDim, varString, standardizationType, redistributeBkgSigFits, doKRandomNeighbors, 
                numberEventsToSavePerProcess, nProcess, seedShift, nentries, nRndRepSubset, nBS, saveBShistsAlso, override_nentries, verbose);
	analysisControl.loadTree(rootFileLoc, rootTreeName);
	analysisControl.loadFitParameters(fitLocation);
	analysisControl.loadData();
	analysisControl.runQFactorThreaded();
        return 0;
}


