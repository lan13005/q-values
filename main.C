#include "main.h"
#include "auxilliary/drawPlots/drawPlots.C"

using namespace RooFit;
int extra=0;

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
        
        cout << "Reserving space in vectors..." << endl;
	// Use these vectors to import all the data to RAM instead of reading from root file
        // First we should reserve the space so there is no resizing of the vectors
        std::vector<double> emptyVec;
        for (auto it=0; it<dim; ++it){
            // copy over emptyVec and then expand
            phaseSpaceVars.push_back(emptyVec);
            phaseSpaceVars[it].reserve(nentries);
        }
	parseDiscrimVars.parseString(s_discrimVar);
        std::vector<double> temp;
	for (int iVar=0; iVar<parseDiscrimVars.varStringSet.size(); ++iVar){
            cout << "Reserving space for discrimVar: " << parseDiscrimVars.varStringSet[iVar] << endl;
            discrimVars.push_back(temp);
            discrimVar.push_back(0);
        }
        // reserve space for the discriminating vars
	for (int iVar=0; iVar<parseDiscrimVars.varStringSet.size(); ++iVar){
            discrimVars[iVar].reserve(nentries);
        }
	AccWeights.reserve(nentries);
        mcprocesses.reserve(nentries);
        utWeights.reserve(nentries);
	// will hold all the ids of the unique combos
	phasePoint2PotentialNeighbor.reserve(nentries);
}

void QFactorAnalysis::loadFitParameters(string fitLocation,string cwd){
	cout << "Loading the fit parameters" << endl;
	double eventRatioSigToBkg; // We try 3 different initializations: {100bkg, 100%sig, 50/50 bkg/sig}. We will use eventRatioSigToBkg to scale the amplitude parameters
	// -----------------------------------------------------
	// -----------------------------------------------------
	//                         LOAD IN THE FITTED PARAMETERS TO THE FULL DISTRIBUTION
	// -----------------------------------------------------
	// -----------------------------------------------------
        // Import all the fitted values from getInitParams
	string varName;
	double varVal;
	ifstream inFile;
	inFile.open((cwd+"/"+fitLocation).c_str());
        std::vector<string> initParNames;
        cout << "RootFile(treeName)(fileTag): " << rootFileLoc << "(" << rootTreeName << ")(" << fileTag << ")" << endl;
	while (inFile >> varName >> varVal){
            if (varName.at(0) != '#'){ //not all output from getInitParams is useful. They start with a # 
                initializationParMap[varName.c_str()] = varVal;
                initParNames.push_back(varName);
                cout << varName << ": " << varVal << endl;
            }
            if (varName == "#eventRatioSigToBkg"){
                eventRatioSigToBkg = varVal;
            }
	}
	cout << "eventRatioSigToBkg: " << eventRatioSigToBkg << endl;
	inFile.close();

        // We will do 3 iterations. Not sure if I am doing this correctly
        // The goal would be to use 100 bkg, 50/50, 100% signal. We can scale the amplitudes by a certain factor related to eventRatioSigToBkg
        // If not redistributing we will use what you want to set it at
        if (redistributeBkgSigFits) { sigFracs={0,0.5,1}; }
        else { sigFracs={eventRatioSigToBkg}; }
}


void QFactorAnalysis::loadData(){
	cout << "Loading the data into arrays" << endl;
	// -----------------------------------------------------
	// -----------------------------------------------------
	//                         LOAD IN THE DATA
	// -----------------------------------------------------
	// -----------------------------------------------------
	// Create variables to hold the data as we read in the data from the tree
        double phaseSpaceVar[dim];
	double AccWeight;

        // vars we will use to fill but not use directly
	ULong64_t eventNumber;
        double utWeight;

	// Set branch addresses so we can read in the data
	parsePhaseSpace.parseString(varString);
        if ( parsePhaseSpace.varStringSet.size() != dim ) { cout << "Uh-oh something went wrong. varString size not the same as dim" << endl; }
	for (int iVar=0; iVar<dim; ++iVar){
            cout << "Setting " << parsePhaseSpace.varStringSet[iVar] << " to phaseSpaceVar index " << iVar << endl; 
            dataTree->SetBranchAddress(parsePhaseSpace.varStringSet[iVar].c_str(),&phaseSpaceVar[iVar]);
        }

	for (int iVar=0; iVar<parseDiscrimVars.varStringSet.size(); ++iVar){
	    dataTree->SetBranchAddress(parseDiscrimVars.varStringSet[iVar].c_str(), &discrimVar[iVar]);
        }
        
        // vars we will use to fill but not use directly
	dataTree->SetBranchAddress("event",&eventNumber);

        //////////////////////////////////////////////////////////
        // SET UP WEIGHTS
        //////////////////////////////////////////////////////////
        cout << "\n" << endl;
        // Setting up tracking of accidental weights
        if (!s_accWeight.empty()){ // if string is not empty we will set the branch address
	    dataTree->SetBranchAddress(s_accWeight.c_str(),&AccWeight);
            cout << "Using accidental weights in branch: "+s_accWeight << endl;
        }
        else{
            AccWeight=1;
            cout << "No accidental weights used" << endl;
        }
        // Setting up tracking of uniqueness tracking weights
        if (!s_utBranch.empty()){ // if string is not empty we will set the branch address
            dataTree->SetBranchAddress(s_utBranch.c_str(),&utWeight);
            cout << "Using uniqueness tracking weights in branch: "+s_utBranch << endl;
        }
        else{
            utWeight=1;
            cout << "No uniqueness tracking weights used" << endl;
        }



        ////////////////////////////////////////////////////
        // CHECKING TO SEE IF MCPROCESSES BRANCH EXIST
        ////////////////////////////////////////////////////
        int mcprocess; 
        bool includeMCprocessInfo;
        if (!s_mcprocessBranch.empty()){ // if string is not empty we will set the branch address
            dataTree->SetBranchAddress(s_mcprocessBranch.c_str(),&mcprocess);
            includeMCprocessInfo=true;
            cout << "mcprocess branch exists at branch: "+s_mcprocessBranch << endl;
        }
        else{
            includeMCprocessInfo=false;
            cout << "no mcprocess branch" << endl;
        }

	
        /////////////////////////////////////////////////////////////
        // LOAD THE DATA
        /////////////////////////////////////////////////////////////
	// We will use a ientry to keep track of which entries we will get from the tree. We will simply use ientry when filling the arrays.  
	for (Long64_t ientry=0; ientry<nentries; ientry++)
	{
		dataTree->GetEntry(ientry);
	        for (int iVar=0; iVar<parsePhaseSpace.varStringSet.size(); ++iVar){
                    phaseSpaceVars[iVar].push_back(phaseSpaceVar[iVar]);
                }
	        for (int iVar=0; iVar<parseDiscrimVars.varStringSet.size(); ++iVar){
                    discrimVars[iVar].push_back(discrimVar[iVar]);
                }

	        AccWeights.push_back(AccWeight);
                utWeights.push_back(utWeight);
                if(includeMCprocessInfo){
                   mcprocesses.push_back(mcprocess); 
                }
	}

	if ( verbose_outputDistCalc ) {
	    cout << "Before standarization the first nentries of phaseSpaceVar[0]" << endl;
	    for ( int ientry=0 ; ientry < nentries; ientry++){
	        cout << phaseSpaceVars[0][ientry] << endl;
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
	    }
	}

	
	// phasePoint1 will consider all events from lowest to largest since these will be our attached q values. 
        // phasePoint2 can use a random subset (hopefully representative of the entire dataset), or accept/reject an entry given a certain criteria, i.e. if some unique indentifier has been seen before
        //      Currently just considers all other entries as potential neighbors
        if(nRndRepSubset!=0 && nRndRepSubset<nentries){
            std::set<Int_t> rndRepSubset; //probably should make the subset have unique elements. So sampling without replacement
            srand (time(NULL)); // set a random seed
            while(rndRepSubset.size() < nRndRepSubset){
                rndRepSubset.insert(rand()%nentries);
            }
            phasePoint2PotentialNeighbor.assign(rndRepSubset.begin(),rndRepSubset.end());
            rndRepSubset.clear();
        }
        else {
	    set<Int_t> setUsedSpectroscopicIDs;
	    for (Int_t ientry=0; ientry<nentries; ientry++){ 
	        phasePoint2PotentialNeighbor.push_back(ientry);
                
	    }
        }
        cout << "\n\n--------------------------------------" << endl;
	cout << phasePoint2PotentialNeighbor.size() << "/" << nentries << " are used as potential neighbors" << endl;
        cout << "Some indicies that will be used as potential neighbors: " << endl;
        for (auto iPhasePoint2=0; iPhasePoint2 < 20; ++iPhasePoint2){
            cout << phasePoint2PotentialNeighbor[iPhasePoint2] << " ";
        }
        cout << endl;

}

// This method was originally designed for multithreading. Turns out RooFit is not thread safe we we had to resort back to spawning multiple root processes. The orignially uses a lambda functon which contains the code to extract the q-values in batches. Some extra code to spawn the threads and waits for all of them to execute is there also
void QFactorAnalysis::runQFactorThreaded(int iProcess){
	// make sure the global variables are read in correctly
	cout << "kDim: " << kDim << endl;
	cout << "numberEventsToSavePerProcess: " << numberEventsToSavePerProcess << endl;
	cout << "nProcess: " << nProcess << endl;
	cout << "seedShift: " << seedShift << endl;
	cout << "nentries: " << nentries << endl;
	cout << "override_nentries: " << override_nentries << endl;
	cout << "verbose: " << verbose << endl;


	// [=] refers to a capture list which is used by this lambda expression. The lambda gets a copy of all the local variables that it uses when it is created. If we
	// just use [] we will get an error since the lambda will have no idea what these variables are	
	//ROOT::EnableThreadSafety();
	//auto f = [=](int iProcess){
	// ----------------------------------------
	// Open up a root file to save the q-factors and other diagnostics to it
	// ---------------------------------------
        // Initializing some variables we can track during the q-value extraction
        double NLL;
        ULong64_t flatEntryNumber;
        double bestNLL;
	double worstNLL;
        double qvalue;
        double best_qvalue;
        double worst_qvalue;
        double eff_nentries;
        double qvalueBS_std=0;
        double best_nsig;
        double best_nbkg;
        double best_ntot;
        double effNentriesMinusTotal;
        int neighbors[kDim];

        // Saving the results along with some diagnostics
        TFile *resultsFile = new TFile((cwd+"/logs"+runTag+"/"+fileTag+"/results"+to_string(iProcess)+".root").c_str(),"RECREATE");
        TTree* resultsTree = new TTree("resultsTree","results");
        resultsTree->Branch("flatEntryNumber",&flatEntryNumber,"flatEntryNumber/l");
        resultsTree->Branch("qvalue",&best_qvalue,"qvalue/D");
        resultsTree->Branch("worst_qvalue",&worst_qvalue,"worst_qvalue/D");
        resultsTree->Branch("qvalueBS_std",&qvalueBS_std,"qvalueBS_std/D");
        resultsTree->Branch("bestNLL",&bestNLL,"bestNLL/D");
        resultsTree->Branch("worstNLL",&worstNLL,"worstNLL/D");
        resultsTree->Branch("best_nsig",&best_nsig,"best_nsig/D");
        resultsTree->Branch("best_nbkg",&best_nbkg,"best_nbkg/D");
        resultsTree->Branch("best_ntot",&best_ntot,"best_ntot/D");
        resultsTree->Branch("eff_nentries",&eff_nentries,"eff_nentries/D");
        resultsTree->Branch("effNentriesMinusTotal",&effNentriesMinusTotal,"effNentriesMinusTotal/D");
        if (saveBranchOfNeighbors){
            resultsTree->Branch("kDim",&kDim,"kDim/I"); // 32 bit integer. Could have used unsigned but not really worth the change...
            resultsTree->Branch("neighbors",neighbors,"neighbors[kDim]/I"); // I = 32 bit integer
        }
        cout << "Set up branch addresses" << endl;

	// Define some needed variables like canvases, histograms, and legends
	cout << "Creating canvas " << iProcess << endl;
    	TCanvas *allCanvases = new TCanvas(("anyHists"+to_string(iProcess)).c_str(),"",1440,900);
        auto legend_fit = new TLegend(0.1,0.7,0.4,0.9);
        auto legend_qVal = new TLegend(0.1,0.7,0.4,0.9);
        TLine* qValLine;
        TH1F* dHist_qvaluesBS = new TH1F(("qvaluesBS"+to_string(iProcess)).c_str(),"Bootstrapped Q-Factors",100,0,1);
        TH1F* dHist_mcprocess = new TH1F(("mcprocess"+to_string(iProcess)).c_str(),"mcprocess",10,0,10);

        // Variables to keep track of phase space neighbors
	double phasePoint1[dim];
	double phasePoint2[dim];
	double distance;
        distSort_kNN distKNN(kDim);
        pair<double,int> newPair;

	// opening a file to write my log data to
    	ofstream logFile;
    	logFile.open((cwd+"/logs"+runTag+"/"+fileTag+"/processLog"+to_string(iProcess)+".txt").c_str());
	
	// Determine what events each thread should run
	int batchEntries = (int)nentries/nProcess; // batchEntries the size of the batch
	int lowest_nentry = iProcess*batchEntries;
	int largest_nentry;
        if (iProcess!=(nProcess-1)) {
            largest_nentry  = (iProcess+1)*batchEntries;
        }
        else {
            largest_nentry = nentries; 
        }
	cout << "nentries we will use for this process: " << lowest_nentry << ", " << largest_nentry << endl;

	// randomly select some events (but deterministically since we migth want to double check) to write histograms for 
        parseEventsToSave.parseString(alwaysSaveTheseEvents);
        TFile *qHistsFile;
	set<int> selectRandomIdxToSave;
        bool saveAllHistograms=false;
        if (numberEventsToSavePerProcess>=0){
	    int randomEvent;
	    srand(iProcess+seedShift);
	    for (int i=0; i<numberEventsToSavePerProcess; i++){
	    	// basically randomly sample a uniform number between (lowest_nentry, largest_nentry)
	    	randomEvent = rand() % (int)batchEntries; // batchEntries is the size of the batch
	    	randomEvent += lowest_nentry; // shift by the lowest entry of the batch
	    	selectRandomIdxToSave.insert( randomEvent );
	    }
            cout << "randomly selected some events to save" << endl;
            cout << "now gathering the events we are told to always save:";
            for (string saveEvent : parseEventsToSave.varStringSet){
                selectRandomIdxToSave.insert(stoi(saveEvent));
                cout << " " << saveEvent;
            }
            cout << endl;
        }
        else {
            saveAllHistograms=true;
        }

        // Determining how many fits we do, depending on if we want to redistribute signal and bkg ratios
        cout << "Eta hist range: " << binRangeEta[0] << ", " << binRangeEta[1] << ", " << binRangeEta[2] << endl;
        cout << "Eta fit range: " << fitRangeEta2[0] << ", " << fitRangeEta2[1] << endl;
        cout << "Pi0 hist range: " << binRangePi0[0] << ", " << binRangePi0[1] << ", " << binRangePi0[2] << endl;
        cout << "Pi0 fit range: " << fitRangePi02[0] << ", " << fitRangePi02[1] << endl;
        
        // ---------------------------
        // RooFit Variables
        // ---------------------------
        // Some vars to initialize but not declare yet
        double sigPdfVal;
        double bkgPdfVal;
        double totPdfVal;
        double weight;

        // Fit parameters 
        double fittedMassX = initializationParMap["massx"];
        double fittedMassY = initializationParMap["massy"];
        double fittedSigmaX = initializationParMap["sigmax"];
        double fittedSigmaY = initializationParMap["sigmay"];
        double fittedBernA = 0.5;//initializationParMap["bernx01"];
        double fittedBernB = 0.5;//initializationParMap["bernx11"];
        double fittedBernC = 0.5;//initializationParMap["berny01"];
        double fittedBernD = 0.5;//initializationParMap["berny11"];
        double fittedBernE = 0.5;//initializationParMap["berny11"];
        cout << "fittedMassX " << fittedMassX << endl;
        cout << "fittedMassY " << fittedMassY << endl;
        cout << "fittedSigmaX " << fittedSigmaX << endl;
        cout << "fittedSigmaY " << fittedSigmaY << endl;
        cout << "fittedBernA " << fittedBernA << endl;
        cout << "fittedBernB " << fittedBernB << endl;
        cout << "fittedBernC " << fittedBernC << endl;
        cout << "fittedBernD " << fittedBernD << endl;

        // For loading the data
        string sThread = to_string(iProcess);
        RooRealVar roo_Meta(("roo_Meta"+sThread).c_str(),"Mass GeV",fitRangeEta2[0],fitRangeEta2[1]);
        roo_Meta.setRange(("roo_fitRangeMeta"+sThread).c_str(),fitRangeEta2[0], fitRangeEta2[1]);
        RooRealVar roo_Mpi0(("roo_Mpi0"+sThread).c_str(),"Mass GeV",fitRangePi02[0],fitRangePi02[1]);
        roo_Mpi0.setRange(("roo_fitRangeMpi0"+sThread).c_str(),fitRangePi02[0], fitRangePi02[1]);
        RooRealVar roo_Weight(("roo_Weight"+sThread).c_str(), "Weight", -10, 10); // Weights can take a wide range
        roo_Meta.setBins(100);
        roo_Mpi0.setBins(100);

        // We need to declare the weight variable with WeightVar and when we fill the weights we have to include the weight var in the argset AND include the weight
        // So when we fill the data we use: rooData.add(RooArgSet(roo_Mpi0,roo_Meta,roo_Weight),weight);
        RooDataSet rooData(("rooData"+sThread).c_str(),"rooData",RooArgSet(roo_Mpi0,roo_Meta,roo_Weight),WeightVar(roo_Weight));
        RooRealVar peak_pi0(("peak_pi0"+sThread).c_str(),"peak_pi0",fittedMassX);//,fittedMassX,fittedMassX);
        RooRealVar width_pi0(("width_pi0"+sThread).c_str(),"width_pi0",fittedSigmaX,fittedSigmaX*0.7,fittedSigmaX*1.3);
        RooRealVar peak_eta(("peak_eta"+sThread).c_str(),"peak_eta",fittedMassY);//,fittedMassY,fittedMassY);
        RooRealVar width_eta(("width_eta"+sThread).c_str(),"width_eta",fittedSigmaY,fittedSigmaY*0.7,fittedSigmaY*1.3);

        RooRealVar bern_parA(("bern_parA"+sThread).c_str(),"bern_parA",fittedBernA,0,1);
        RooRealVar bern_parB(("bern_parB"+sThread).c_str(),"bern_parB",fittedBernB,0,1);
        RooRealVar bern_parC(("bern_parC"+sThread).c_str(),"bern_parC",fittedBernC,0,1);
        RooRealVar bern_parD(("bern_parD"+sThread).c_str(),"bern_parD",fittedBernD,0,1);
        //RooRealVar bern_parE(("bern_parE"+sThread).c_str(),"bern_parE",fittedBernE,0,1);///,fittedBernE,0,1);
        RooGenericPdf rooBkgX(("rooBkgX"+sThread).c_str(), "rooBkgX", ("bern_parA"+sThread+"*roo_Mpi0"+sThread+"+bern_parB"+sThread+"*(1-roo_Mpi0"+sThread+")").c_str(),RooArgSet(bern_parA,bern_parB,roo_Mpi0));
        //RooGenericPdf rooBkgY(("rooBkgY"+sThread).c_str(), "rooBkgY", ("bern_parC"+sThread+"*(1-roo_Meta"+sThread+")**2+bern_parD"+sThread+"*2*roo_Meta"+sThread+"*(1-roo_Meta"+    sThread+")+bern_parE"+sThread+"*(roo_Meta"+sThread+")**2").c_str(),RooArgSet(bern_parC,bern_parD,bern_parE,roo_Meta));
        RooGenericPdf rooBkgY(("rooBkgY"+sThread).c_str(), "rooBkgY", ("bern_parC"+sThread+"*roo_Meta"+sThread+"+bern_parD"+sThread+"*(1-roo_Meta"+sThread+")").c_str(),RooArgSet(bern_parC,bern_parD,roo_Meta));
        RooGaussian rooGausPi0_bkg(("rooGausPi0_bkg"+sThread).c_str(), "rooGausPi0_bkg", roo_Mpi0, peak_pi0, width_pi0);
        RooRealVar bkgPeakFrac(("bkgPeakFrac"+sThread).c_str(),"bkgPeakFrac",1,0,1);
        RooAddPdf rooBkgXplusPi0Peak(("rooBkgXplusPi0Peak"+sThread).c_str(), "rooBkgXplusPi0Peak", RooArgList(rooGausPi0_bkg,rooBkgX),RooArgSet(bkgPeakFrac));
        RooProdPdf rooBkg(("rooBkg"+sThread).c_str(),"rooBkg",RooArgList(rooBkgXplusPi0Peak,rooBkgY));

        RooGaussian rooGausPi0(("rooGausPi0"+sThread).c_str(), "rooGausPi0", roo_Mpi0, peak_pi0, width_pi0);
        RooGaussian rooGausEta(("rooGausEta"+sThread).c_str(), "rooGausEta", roo_Meta, peak_eta, width_eta);
        RooProdPdf rooGaus2D(("rooGaus2D"+sThread).c_str(), "rooGaus2D", RooArgSet(rooGausPi0,rooGausEta));

        //RooGenericPdf rooZero(("rooZero"+sThread).c_str(), "rooZero", "0",RooArgSet(roo_Mpi0,roo_Meta)); // used just to do some tests
        
        RooRealVar nsig(("nsig"+sThread).c_str(),"nsig",0,kDim);
        RooRealVar nbkg(("nbkg"+sThread).c_str(),"nbkg",0,kDim);
        RooAddPdf rooSigPlusBkg(("rooSumPdf"+sThread).c_str(), "rooSumPdf", RooArgList(rooGaus2D,rooBkg),RooArgSet(nsig,nbkg));

        RooArgSet* savedParams; 
        RooArgSet* params; // intermedate parameter values for the pdfs
        // ---------------------------

        // Saving bootstrap results
        vector<double> qvalues; qvalues.reserve(nBS);

	// the main loop where we loop through all events in a double for loop to calculate dij. Iterating through all j we can find the k nearest neighbors to event i.
	// Then we can plot the k nearest neighbors in the discriminating distribution
        // Finally calculate the q-value by fitting and getting signal fraction
	//logFile << std::fixed << std::setprecision(6);
	int randomEntry;
        int counter=-1;
        RooTrace::mark();
        int skipInitial=0;
        //if (iProcess==0){
        //    skipInitial=229079;
        //}
        for (int ientry=lowest_nentry+skipInitial; ientry<largest_nentry+skipInitial; ientry++){
                //ientry = ientry2+38145;
                //if (ientry > largest_nentry){ ientry = lowest_nentry; } 
                //if ( (discrimVars[1][ientry] < 0.548625-3*0.0191) || (discrimVars[1][ientry] > 0.548625+3*0.0191) 
                //        || (discrimVars[0][ientry] < 0.135881-3*0.0076) || (discrimVars[0][ientry] > 0.135881+3*0.0076) 
                //        || (phaseSpaceVars[1][ientry] < 0.75) || (phaseSpaceVars[1][ientry] > 0.85) 
                //){
                //    continue;
                //}
                //++counter;
                //if (counter > extra){
                //    exit(0);
                //}
                dHist_qvaluesBS->Reset();
                dHist_mcprocess->Reset();
                if (saveBranchOfNeighbors){
                    memset(neighbors,-1,sizeof(neighbors)); // last argument sets that number of bytes to the specified value. Dont want just kDim here since it is 4 Bytes per
                }
                qvalues.clear();
        	if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end() || saveAllHistograms) {
                    qHistsFile = new TFile((cwd+"/histograms"+runTag+"/"+fileTag+"/qValueHists_"+to_string(ientry)+".root").c_str(),"RECREATE");
                }

                // Outputting the progress of each thread
		flatEntryNumber=ientry;
		auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
		auto duration_beginEvent = std::chrono::high_resolution_clock::now();
		if(verbose) { logFile << "Starting event " << ientry << "/" << largest_nentry << " ---- Global Time: " << duration2 << "ms" << endl; }
		cout << "Starting event " << ientry << "/" << largest_nentry << " ---- Global Time: " << duration2 << "ms" << endl; 
		
                // Phase space Definitions
		for ( int iVar=0; iVar<dim; ++iVar ){
			phasePoint1[iVar] = phaseSpaceVars[iVar][ientry];
		}
                int nPotentialNeighbors=(int)phasePoint2PotentialNeighbor.size();
                vector<int> phasePoint2PotentailNeighbor_BS;
                phasePoint2PotentailNeighbor_BS.reserve(nPotentialNeighbors);

                // --------------------------------------------
                // Random generator for resampling of the input data, in this case the input data will be the set of potential neighbors for phasePoint2.
                // --------------------------------------------
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); 
                // RNG using distributions are really cheap. It is also threadsafe also whereas rand() is not. You will see all the threads are not operating at max potential if you use rand
                // https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers
                // Mersenne twistor requires large state storage but has the highest quality + period. 
                // Lagged Fibonoci Generator liek ranlux24 has a much smaller period but should be much faster
                // https://books.google.com/books?id=zi08DQAAQBAJ&pg=PA237&lpg=PA237&dq=ranlux24+%22period%22&source=bl&ots=N4HXolKGDr&sig=ACfU3U1JfZcQczw-xoAv-lyHTOCOf4Xpjw&hl=en&sa=X&ved=2ahUKEwiosZP8-_npAhW2RTABHatDDe8Q6AEwAnoECAsQAQ#v=onepage&q=ranlux24%20%22period%22&f=false
                // mt19937 has a period around 10^6000 whereas ranlux24 has a period around 10^171. ranlux24 produces 24 bit integers. 24 bits ~ 17M which is much larger than the number of potential neighbors
                //      so in our case ranlux is probably just fine
                std::uniform_int_distribution<int> distribution(0,nPotentialNeighbors-1); //if we ask for a range of 0 to 10 it will include 10. phasePoint2PotentialNeighbor has a zero index. Need to subtract by 1
                //static thread_local std::mt19937 generator(seed); 
                static thread_local std::ranlux24_base generator(seed); 

                // ---------------------
                // NOW FIND NEIGHBORS AND EXTRACT Q-FACTORS.
                // - The last iteartion is also the full data. This is because we will draw the histograms on the last iteration to 
                //   skip intermediate saving of the histograms
                //  if nBS > 0 then we will rerun and resample to get the bootstrapped q-factors
                // ---------------------
                for (int iBS=0; iBS<nBS+1; ++iBS){ 
                    // clean up and reserve for next entry
                    rooData.reset();
                    // resetting some variables
		    bestNLL=DBL_MAX;
		    worstNLL=-1*DBL_MAX; // DBL_MIN is basically 0. We want a very negative number
		    duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		    if(verbose){logFile << "\tBegin bootstrapping potential neighbors: " << duration2 << "ms" << endl; }
                    phasePoint2PotentailNeighbor_BS.clear();
                    if (iBS!=nBS){ // resampling neighbors with replacement if we want to do bootstrapping
                        for(int iNeighbor=0; iNeighbor<nPotentialNeighbors; ++iNeighbor){
                            phasePoint2PotentailNeighbor_BS.push_back(phasePoint2PotentialNeighbor[distribution(generator)]);
                        }
                    }
                    else{
                        phasePoint2PotentailNeighbor_BS = phasePoint2PotentialNeighbor;
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
                            if (jentry == ientry){ continue; } 
		            if ( verbose_outputDistCalc ) { cout << "event i,j = " << ientry << "," << jentry << endl;} 
		            for ( int iVar=0; iVar<dim; ++iVar ){
		               	phasePoint2[iVar] = phaseSpaceVars[iVar][jentry];
		            }
		            distance = calc_distance(dim,phasePoint1,phasePoint2,verbose_outputDistCalc);
		            distKNN.insertPair(make_pair(distance,jentry));
		        }
                    }
		    duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		    if(verbose){logFile << "\tFound neighbors: " << duration2 << "ms" << endl; }
		    if (distKNN.kNN.size() != kDim){ cout << "size of distKNN is not equal to kDim! size,kDim="<< distKNN.kNN.size() << "," << kDim 
		        << "\n    -- if size is 1 less than kDim it is probably because kDim=nentries and event i cannot be a neighbor to itself" << 
		        "\n    -- if size != kDim it could also mean that the number of spectroscopically unique neighbors reduces the number of poential neighbors below kDim" << endl;}
                    if(verbose_outputDistCalc){ cout << "These are our neighbors" << endl; }


                    int iNeighbor=-1;
		    while ( distKNN.kNN.empty() == false ){
		            newPair = distKNN.kNN.top();
		            distKNN.kNN.pop();
                            if (weightingScheme==""){ weight=1; }
                            if (weightingScheme=="as"){ weight=AccWeights[newPair.second]; }
                            else { weight=1; } 
                            weight=weight*utWeights[newPair.second];

                            roo_Meta = discrimVars[1][newPair.second];
                            roo_Mpi0 = discrimVars[0][newPair.second];
                            // with 600 neighbors it seems like the fitTo command takes ~2x longer when using Range() argument which selects the fit range.
                            // This is equivalent to shrinking the data set range and fitting over the full range which will save time.
                            // It might be useful to think of setting a fit range as to lower the effective number of nearest neighbors. Another thing that lowers
                            // the effective number of neighbors is any weights we apply to the filling of the histograms
                            if ( discrimVars[1][newPair.second] > fitRangeEta2[0] && 
                                    discrimVars[1][newPair.second] < fitRangeEta2[1] &&
                                    discrimVars[0][newPair.second] > fitRangePi02[0] &&
                                    discrimVars[0][newPair.second] < fitRangePi02[1] ){ 
                                roo_Weight.setVal(weight);
                                // roo_Weight will get overwritten here when adding to RooDataSet but actually does not pick up the value. So we cannot use 
                                // roo_Weight.getVal() but the dataset is weighted: https://root-forum.cern.ch/t/fit-to-a-weighted-unbinned-data-set/33495
                                rooData.add(RooArgSet(roo_Mpi0,roo_Meta,roo_Weight),weight);
                                if (saveBranchOfNeighbors){
                                    neighbors[++iNeighbor]=newPair.second;
                                }
                            }

                            if(verbose_outputDistCalc){
		                cout << "(" << newPair.first << ", " << newPair.second << ", " << discrimVars[0][newPair.second] << ", " << discrimVars[1][newPair.second] << ")" << endl; 
                            }

                            dHist_mcprocess->Fill(mcprocesses[newPair.second]);
		    }
                    eff_nentries = rooData.sumEntries();
                    //rooData.Print();
		    
		    duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		    if(verbose){logFile <<	"\tFilled neighbors: " << duration2 << "ms" << endl;}
		    
		    // /////////////////////////////////////////
		    // Calcuclate q-value
		    // /////////////////////////////////////////
		    //// We use a normalized gaussian and a flat function. 
		    for ( auto initSigFrac : sigFracs ){
                        // intitialize the variables for the fit PDF
		        duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		        if(verbose){logFile <<	"\tPrepping 2D fits:	 " << duration2 << "ms" << endl;}
                        peak_pi0.setVal(fittedMassX);
                        peak_eta.setVal(fittedMassY);
                        width_pi0.setVal(fittedSigmaX);
                        width_eta.setVal(fittedSigmaY);

                        bern_parA.setVal(fittedBernA);
                        bern_parB.setVal(fittedBernB);
                        bern_parC.setVal(fittedBernC);
                        bern_parD.setVal(fittedBernD);
                        bkgPeakFrac.setVal(1);

                        nsig.setVal(initSigFrac*eff_nentries);
                        nbkg.setVal((1-initSigFrac)*eff_nentries);

                        // need RooFit::SumW2Error(true) since we are using a weighted dataset
		        duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		        if(verbose){logFile <<	"\tBeginning Fit:	 " << duration2 << "ms" << endl;}
                        RooFitResult* roo_result;
                        ////////////////////////////////////////////////////////////////////////////////////////////
                        //  BEGIN 2D GAUS FIT
                        ////////////////////////////////////////////////////////////////////////////////////////////
                        roo_result = rooSigPlusBkg.fitTo(rooData, Minos(kTRUE), Save(), RooFit::SumW2Error(true), PrintLevel(-1), BatchMode(kTRUE));//, Range(("roo_fitRangeMpi0"+sThread+",roo_fitRangeMeta"+sThread).c_str()));// Hesse(kFALSE));
                        
		        duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		        if(verbose){logFile <<	"\tFitted hist with some initialization: " << duration2 << "ms" << endl;}
                        
                        // setting parameters for bkg/sig and extracting q-value
                        double sigFrac = nsig.getVal()/(nsig.getVal()+nbkg.getVal());
                        roo_Mpi0.setVal(discrimVars[0][ientry]);
                        roo_Meta.setVal(discrimVars[1][ientry]);
                        sigPdfVal = sigFrac*rooGaus2D.getVal(RooArgSet(roo_Mpi0,roo_Meta));
                        bkgPdfVal = (1-sigFrac)*rooBkg.getVal(RooArgSet(roo_Mpi0,roo_Meta));
                        totPdfVal = rooSigPlusBkg.getVal(RooArgSet(roo_Mpi0,roo_Meta));
                        qvalue = sigPdfVal/(sigPdfVal+bkgPdfVal);
                        //if (width_eta.getVal() <= width_eta.getError()){
                        //    cout << "event " << ientry << " has width_eta="<<width_eta.getVal()<< " with errors="<<width_eta.getError()<<" within 0! Set Q to zero" << endl;
                        //    qvalue=0;
                        //}
		        //if(verbose){logFile <<	"\tExtracted q-value (" << qvalue << "): " << duration2 << "ms" << endl;}
                        
                        
		    	NLL = roo_result->minNll();
		    	if (NLL < bestNLL){
		    		best_qvalue = qvalue;
                                best_nsig = nsig.getVal();
                                best_nbkg = nbkg.getVal();
                                best_ntot = best_nsig + best_nbkg;
                                effNentriesMinusTotal = eff_nentries-best_ntot;
		    		bestNLL=NLL;
                                params=rooSigPlusBkg.getParameters(RooArgList(roo_Mpi0,roo_Meta));
                                savedParams = (RooArgSet*)params->snapshot();
		    	} 
		    	if (NLL > worstNLL){
		    		worstNLL = NLL;
                                worst_qvalue = qvalue;
		    	}
                        
                        
                        // Different ways to calculate chiSq: https://nbviewer.jupyter.org/gist/wiso/443934add13fd7226e4b
                        // We can calculate chiSq by binning the data and using RooChi2Var. Example calculation is at:
                        // https://hep.lancs.ac.uk/~ajf/root/RooChi2MCSModule_8cxx_source.html
                        // Think we want to use nDataPts - nConstraints as in https://ned.ipac.caltech.edu/level5/Leo/Stats7_2.html
                        //RooDataHist* binnedData = rooData.binnedClone();
                        //RooChi2Var chiSq("chiSq","chiSq",rooSigPlusBkg,*binnedData);//,DataError(RooAbsData::SumW2));
                        //RooPlot* roo_Mpi0_Meta_frame = new RooPlot(roo_Mpi0,roo_Meta);
                        //rooData.plotOn(roo_Mpi0_Meta_frame);
                        //rooSigPlusBkg.plotOn(roo_Mpi0_Meta_frame);
                        //double chiSq = roo_Mpi0_Meta_frame->chiSquare();
                        //RooArgSet* floatPars = (RooArgSet*)rooSigPlusBkg.getParameters(rooData)->selectByAttrib("Constant",kFALSE);
                        //int nParams = floatPars->getSize();
                        //int nBins = 900;
                        //double chiSqPerDOF = chiSq/(nBins-nParams);
		    	// now that the q-value is found we can get the NLL and save the parameters with the best NLL
                        // for more information look at RooFitResult class https://root.cern.ch/doc/master/classRooFitResult.html
                        //cout << "NLL, chiSq, reduced ChiSq, nParams: " << NLL << ", " << chiSq << ", " << chiSqPerDOF << ", " << nParams << endl;
		    } // close redistribute initialization fit loop


                    // Save the qvalues for this iteration
                    if(iBS<nBS){  
                        qvalues.push_back(best_qvalue);
                        dHist_qvaluesBS->Fill(best_qvalue);
                    }

                    // EITHER WE POTENTIALLY SAVE THE HISTOGRAM AFTER ALL BS IS DONE OR WE SAVE IT EVERY ITERATION DURING THE BS. ALLOWS US TO SEE HOW THINGS DEVELOP
                    // THERE IS STILL AN INNER CONDITION THAT ONLY SELECTS A SUBSET OF THESE. SINCE WE DONT WANT 100K+ HISTOGRAMS
                    if (iBS==nBS || saveBShistsAlso){
        	        // /////////////////////////////////////////
        	        // Drawing histogram 
        	        // /////////////////////////////////////////
        	        if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end() || saveAllHistograms) {
        	                // Here we draw the histograms that were randomly selected
        	                allCanvases->Clear();
        	                legend_fit->Clear();
                                legend_qVal->Clear();

                                TLine* varLine = new TLine(0,0,0,0);
        	        	varLine->SetLineColor(kOrange);
                                varLine->SetLineWidth(2);

		                cout <<	"\tSaving diagnostic histogram: " << duration2 << "ms" << endl;
		                if(verbose){logFile <<	"\tSaving 2D histogram: " << duration2 << "ms" << endl;}
                                allCanvases->Divide(3,2);
        
                                allCanvases->cd(1);
                                /// -- Might need to use this code section if RooFit ever becomes thread safe. If we draw rooSigPlusBkg with the various Components we get an error when using plotOn
                                /// -- when using multi threads. Projecting each sub PDF makes it work fine. 
                                //RooAbsPdf* rooGaus2DProjMpi0 = rooGaus2D.createProjection(roo_Mpi0);
                                //rooGaus2DProjMpi0->setNormRange(("roo_fitRangeMpi0"+sThread).c_str());
                                //rooGaus2DProjMpi0->plotOn(roo_Meta_frame);
                                //RooAbsPdf* rooBkgProjMpi0 = rooBkg.createProjection(roo_Mpi0);
                                //rooBkgProjMpi0->setNormRange(("roo_fitRangeMpi0"+sThread).c_str());
                                //rooBkgProjMpi0->plotOn(roo_Meta_frame, LineStyle(kDashed),LineColor(kOrange));
                                //RooAbsPdf* rooSigPlusBkgProjMpi0 = rooSigPlusBkg.createProjection(roo_Mpi0);
                                //rooSigPlusBkgProjMpi0->setNormRange(("roo_fitRangeMpi0"+sThread).c_str());
                                //rooSigPlusBkgProjMpi0->plotOn(roo_Meta_frame);
                                //rooSigPlusBkgProjMpi0->plotOn(roo_Meta_frame, Components("rooBkg*"),LineStyle(kDashed),LineColor(kOrange));
                                //rooSigPlusBkgProjMpi0->paramOn(roo_Meta_frame);
                                
                                // OMG FOUND HOW TO FIX THE BUG! WE NEED TO SET THE NORMRANGE AND RANGE HERE IF WE ARE GOING TO FIT USING THE SAME PDF ON THE SAME DATA MULTIPLE TIMES. 
                                // THERE IS AN ERROR ABOUT THE NORMALIZATIONS AND SOMEHOW MORE FUNCTIONS ARE CREATED AND WE GET A BUNCH OF SAME OF THE SAME OBJECTS THAT GETS ADDED
                                // ONTO OUR ARGLIST
                                // https://root-forum.cern.ch/t/roofit-normalization/23644/2
                                //RooPlot* roo_Meta_frame = roo_Meta.frame();
                                //rooData.plotOn(roo_Meta_frame);

                                // reload the best params
                                TH1* model_hist;
                                TH1* model_sig;
                                TH1* model_bkg; // need these pointers to pass to drawPlots so we can properly clean these objects afterwards
                                RooArgSet* params=rooSigPlusBkg.getParameters(RooArgList(roo_Mpi0,roo_Meta));
                                *params = *savedParams;
                                drawPlots(&roo_Mpi0, &roo_Meta, discrimVars[0][ientry], discrimVars[1][ientry], eff_nentries, &rooSigPlusBkg, &rooBkg, &rooGaus2D, &rooData, &nsig, &nbkg, allCanvases, model_hist, model_sig, model_bkg);
                                
                                //rooSigPlusBkg.plotOn(roo_Meta_frame, NormRange(("roo_fitRangeMeta"+sThread).c_str()),Range(("roo_fitRangeMeta"+sThread).c_str()));
                                //rooSigPlusBkg.plotOn(roo_Meta_frame, NormRange(("roo_fitRangeMeta"+sThread).c_str()),Range(("roo_fitRangeMeta"+sThread).c_str()), Components(rooBkg),LineStyle(kDashed),LineColor(kOrange));
                                //rooSigPlusBkg.plotOn(roo_Meta_frame, NormRange(("roo_fitRangeMeta"+sThread).c_str()),Range(("roo_fitRangeMeta"+sThread).c_str()), Components(rooGaus2D),LineStyle(kDashed),LineColor(kMagenta));
                                //roo_Meta_frame->Draw();
                                //roo_Meta_frame->GetXaxis()->SetTitle("M_{#eta}");
                                //roo_Meta_frame->SetTitle(("BEST VALS:  QValue="+to_string(best_qvalue)).c_str());
        	          	//varLine->DrawLine(discrimVars[1][ientry],0,discrimVars[1][ientry],roo_Meta_frame->GetMaximum());
		                //if(verbose){logFile <<	"\tCompleted drawing on pad 1: " << duration2 << "ms" << endl;}

                                //allCanvases->cd(2);
                                //RooPlot* roo_Mpi0_frame = roo_Mpi0.frame();
                                //rooData.plotOn(roo_Mpi0_frame);
                                //rooSigPlusBkg.plotOn(roo_Mpi0_frame, NormRange(("roo_fitRangeMpi0"+sThread).c_str()),Range(("roo_fitRangeMpi0"+sThread).c_str()));
                                //rooSigPlusBkg.plotOn(roo_Mpi0_frame, NormRange(("roo_fitRangeMpi0"+sThread).c_str()),Range(("roo_fitRangeMpi0"+sThread).c_str()), Components(rooBkg),LineStyle(kDashed),LineColor(kOrange));
                                //rooSigPlusBkg.plotOn(roo_Mpi0_frame, NormRange(("roo_fitRangeMpi0"+sThread).c_str()),Range(("roo_fitRangeMpi0"+sThread).c_str()), Components(rooGaus2D),LineStyle(kDashed),LineColor(kMagenta));
                                //roo_Mpi0_frame->Draw();
                                //roo_Mpi0_frame->GetXaxis()->SetTitle("M_{#pi}");
        	          	//varLine->DrawLine(discrimVars[0][ientry],0,discrimVars[0][ientry],roo_Mpi0_frame->GetMaximum());
		                //if(verbose){logFile <<	"\tCompleted drawing on pad 2: " << duration2 << "ms" << endl;}
        
                                //allCanvases->cd(3);
                                //TH1* rooDataHist = rooData.createHistogram(("roo_Mpi0"+sThread+",roo_Meta"+sThread).c_str(),50,50);
                                //rooDataHist->SetFillColor(kYellow);
                                //rooDataHist->Draw("lego1 0"); 
                                //TPolyLine3D *pl3d = new TPolyLine3D(2);
                                //pl3d->SetLineColor(kRed);
                                //pl3d->SetLineWidth(3);
                                //pl3d->SetPoint(0,discrimVars[0][ientry],discrimVars[1][ientry],0);
                                //pl3d->SetPoint(1,discrimVars[0][ientry],discrimVars[1][ientry],rooDataHist->GetMaximum());
                                //pl3d->Draw("same");
                                //gPad->SetTheta(75);
                                //gPad->Update();
		                //if(verbose){logFile <<	"\tCompleted drawing on pad 3: " << duration2 << "ms" << endl;}

                                //allCanvases->cd(4);
                                //rooDataHist->Draw("COLZ"); 
        	          	//varLine->DrawLine(fitRangePi02[0],discrimVars[1][ientry],fitRangePi02[1],discrimVars[1][ientry]);
        	          	//varLine->DrawLine(discrimVars[0][ientry],fitRangeEta2[0],discrimVars[0][ientry],fitRangeEta2[1]);
                                //gPad->Update();
		                //if(verbose){logFile <<	"\tCompleted drawing on pad 4: " << duration2 << "ms" << endl;}

                                allCanvases->cd(6);
                                dHist_mcprocess->SetTitle(("Current mcprocess: "+std::to_string(mcprocesses[ientry])).c_str());
                                dHist_mcprocess->Scale(1.0/kDim);
                                dHist_mcprocess->Draw("HIST");
                                dHist_mcprocess->GetXaxis()->SetTitle("mcprocess");
                                dHist_mcprocess->GetYaxis()->SetTitle("percentage of neighbors");


                                /////////////////////////////////////////////////////////////////////////
                                ////////////////////////////////////// BOOTSTRAP HISTOGRAM OF Q-FACTORS
                                //if (iBS==nBS){
                                //    qvalueBS_std=calculateStd(nBS,&qvalues[0]);
                                //    dHist_qvaluesBS->SetTitle(("STD: "+to_string(qvalueBS_std)).c_str());
                                //}
                                //else{
                                //    dHist_qvaluesBS->SetTitle("Bootstrapped Q-Factors");

                                //}
                                //dHist_qvaluesBS->Draw();
        	                //qValLine = new TLine(best_qvalue,0,best_qvalue,nBS);
        	                //qValLine->SetLineColor(kOrange);
                                //legend_qVal->AddEntry(qValLine,"True Q Value");
                                //qValLine->Draw("SAME");
                                //legend_qVal->AddEntry(dHist_qvaluesBS,"Bootstrapped Q Values");
                                //legend_qVal->Draw();
		                //if(verbose){logFile <<	"\tCompleted drawing on pad 6: " << duration2 << "ms" << endl;}
                                /////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////


        	        	// need to save as a root file first then convert to pngs or whatever. Seems like saveas doesnt like threaded since the processes might make only one png converter
        	        	// or whatever and maybe if multiple threads calls it then a blocking effect can happen
                                qHistsFile->cd();
                                if(iBS==nBS){
                                    allCanvases->Write(("Mass-event"+std::to_string(ientry)).c_str());
                                    qHistsFile->Close();
                                }
                                else {
                                        allCanvases->Write(("Mass-event"+std::to_string(ientry)+"BS"+to_string(iBS)).c_str());
                                }
                                
        	        	    
        	        	if(verbose){
        	        	     duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
        	        	     logFile << "\tSaved this histogram since it was randomly selected: " << duration2 <<  "ms" << endl;
        	        	}
		        } 
                    } // closes condition to save histograms

                } // finishes nBS loop
		if(verbose){logFile << "\tCurrent Best NLL = " << to_string(bestNLL) << ": " << duration2 << "ms" << endl; }
		resultsTree->Fill();
	}
        resultsFile->cd();
        resultsTree->Write();
        cout << "nentries: " << nentries << endl;
        
	// Finish the log files by including an elapsed time and finally closing the file
	if (verbose){
		auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
		logFile << "Total time: " << duration2 << " ms" << endl;
		logFile << "Time Per Event: " << duration2/(largest_nentry-lowest_nentry) << " ms" << endl;
	}
	logFile.close();

        // OLD CODE FOR SPAWNING THREADS INSTEAD OF WHAT WE USE NOW, WHICH IS JUST PROCESSES
	//};
        //// Now that we have the lambda function we can start to spawn threads
	//cout << "Launching " << nProcess << " threads 1 second apart!" << endl;
	//vector<thread> threads;
	//for ( int iProcess=0; iProcess<nProcess; ++iProcess){
	//	cout << "(Thread " << iProcess << ") is starting" << endl;
	//	threads.emplace_back( [f, iProcess] { f(iProcess); } );
        //        sleep(1);
	//	//threads[iProcess] = std::thread(QFactorAnalysis::staticEntryPoint, this, iProcess);
	//}
	//for (auto&& t : threads) t.join(); // join waits for completion
	//cout << "Threads have completed running!" << endl;
}

int main( int argc, char* argv[] ){
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); // For thread safety we need this
	TH1::AddDirectory(kFALSE);

	// This suppresses all the "info" messages from roofit. Like saying that it will use numeric integrator for the generic pdf we define
	// https://root-forum.cern.ch/t/suppressing-info-messages/14642/6
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
        RooTrace::active(kTRUE);
    	gStyle->SetOptFit(1111);
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
        int iProcess=atoi(argv[8]);
	int nProcess=atoi(argv[9]);
        int seedShift= atoi(argv[10]);
        Long64_t nentries= atoi(argv[11]);
        int nRndRepSubset=atoi(argv[12]);
        int nBS=atoi(argv[13]);
        bool saveBShistsAlso;
        if ( atoi(argv[14]) ==1 ){ saveBShistsAlso=true; }
        else { saveBShistsAlso=false; }
        bool override_nentries;
        cout << "Num Vars="<< dim << endl;
        if ( atoi(argv[15]) ==1 ){ override_nentries=true; }
        else{ override_nentries=false; }
        bool verbose;
        if ( atoi(argv[16])==1 ){ verbose=true;}
        else{ verbose=false;}
        std::string cwd=argv[17];
        std::string alwaysSaveTheseEvents=argv[18];
        bool saveBranchOfNeighbors;
        if ( atoi(argv[19])==1) { saveBranchOfNeighbors=true; }
        else { saveBranchOfNeighbors=false; } 
        cout << "----------------------------" << endl;
        cout << "kDim: " << kDim << endl;
        cout << "iProcess: " << iProcess << endl;
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
        cout << "verbose: " << verbose  << endl; 
        cout << "redistributeBkgSigFits: " << redistributeBkgSigFits << endl;
        cout << "doKRandomNeighbors: " << doKRandomNeighbors << endl;
        cout << "cwd: " << cwd << endl;
        cout << "alwaysSaveTheseEvents: " << alwaysSaveTheseEvents << endl;
        cout << "saveBranchOfNeighbors: " << saveBranchOfNeighbors << endl;
        cout << "fitLocation: " << fitLocation << endl;
        cout << "----------------------------" << endl;

        //cout << "Sleeping for 10 seconds so you can look at these settings" << endl;
        //sleep(10);
    
	QFactorAnalysis analysisControl(kDim, varString, cwd, standardizationType, redistributeBkgSigFits, doKRandomNeighbors, 
                numberEventsToSavePerProcess, nProcess, seedShift, nentries, nRndRepSubset, nBS, saveBShistsAlso, override_nentries, verbose, alwaysSaveTheseEvents,
                saveBranchOfNeighbors);
	analysisControl.loadTree(rootFileLoc, rootTreeName);
	analysisControl.loadFitParameters(fitLocation,cwd);
	analysisControl.loadData();
        analysisControl.runQFactorThreaded(iProcess);
        return 0;
}


