#ifndef MAIN_H
#define MAIN_H

#include "helperFuncs.h"
#include <thread>
#include "Math/MinimizerOptions.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <TCanvas.h>
#include <TRandom.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLegend.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGenericPdf.h>
#include <RooConstVar.h>
#include <RooArgList.h>
#include <RooPlot.h>
#include <RooMinuit.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooChi2Var.h>
//#include "RooAddPdf.h"
//#include "RooFormulaVar.h"

//#include "Math/MinimizerOptions.h"
//using namespace RooFit;


const int dim=3; // will get replaced by run.py
bool verbose_outputDistCalc=false;
TRandom rgen;

using namespace std;
// NO SPACES BETWEEN THE = SIGNS. I USE SED TO REPLACE
string rootFileLoc="/d/grid15/ln16/pi0eta/q-values/degALL_bcal_treeFlat_DSelector.root";
string rootTreeName="degALL_bcal_tree_flat";
string fileTag="bcal";
string weightingScheme="as"; // "" or "as*bs"
string s_accWeight="AccWeight";
string s_discrimVar="Meta";
string s_sideBandVar="Mpi0";


// OUT OF DATED CODE THAT USES ROOFIT TO DO UNBINNED MAX LIKELIHOOD FIT. WILL PROBABLY NEED TO REIMPLEMENT THIS
//class rooFitML{
//	private:
//        	double peakWidtheta[2] = {0.545928, 0.0196892};
//        	double peakWidthpi0[2] = {0.135399, 0.00760648};
//        	double par0eta[3] = {3.4528, 1.7264, 0};
//        	double par1eta[3] = {5.49805, 2.74902, 0}; 
//        	double par2eta[3] = {0, 0.9, 1.18}; 
//        	double par0pi0[3] = {7.30661, 3.6533, 0}; 
//        	double par1pi0[3] = {30.05331, 15.2665, 0}; 
//        	double par2pi0[3] = {0, 0.400002, 0.800003}; 
//
//		/// we had to make everything deal with pointers since we couldn't define and initialize at the same time here. It can only be done with static const primitive types
//		RooRealVar par0_eta("par0_eta","par0_eta", 0, 1, 5.25);
//		RooRealVar par1_eta("par1_eta","par1_eta", 0, -8.25, 8.25);
//		RooRealVar par2_eta("par2_eta","par2_eta", 0, 0, 1.77);
//		RooRealVar peak_eta("peak_eta","peak_eta",0.545928,0.52,0.58);
//		RooRealVar x("x","x",0, 1);
//		RooRealVar width_eta("width_eta","width_eta",0.0196892,0.017,0.027);
//		RooDataSet data("data","data",RooArgSet(x));
//
//		RooPlot xframe = x->frame();
//	public:
//		//RooGaussian gauss("gauss","gauss(x,m,s)", x, peak_eta, width_eta) ;
//		//RooFormulaVar poly("poly","polynomial",x,RooArgList(par0_eta,par1_eta)) ;
//		//RooAddPdf sigWithBkg("sigWithBkg","gaus+polyBkg",RooArgList(gauss,poly),RooArgList(scale
//		RooGenericPdf sigWithBkg("sigWithBkg","sigWithBkg","par0_eta+par1_eta*x+par2_eta*exp(-0.5*((x-peak_eta)/width_eta)**2)", RooArgSet(x,par0_eta,par1_eta,par2_eta,peak_eta,width_eta));
//
//
//		void fillValues( double mass ){
//			data.add(RooArgSet(mass));
//		}
//
//		void fitToData( UInt_t iFit ){
//			par0_eta.setVal(par0eta[iFit]);
//			par1_eta.setVal(par1eta[iFit]);
//			par2_eta.setVal(par2eta[iFit]); 
//			// not sure what these flags do
//			sigWithBkg.fitTo(data,"mhvr");
//		}
//
//		void drawFitData(){
//			data.plotOn(xframe);
//			sigWithBkg.plotOn(xframe);
//			xframe.Draw();
//		}
//
//		void clear(){
//			// think this deletes the pool of data
//			data.Clear();
//			// we might also have to reinitialize the RooRealVars if the fitTo function decides to change things, idky it would though
//		}
//
//};


class QFactorAnalysis{
	private:
		// variables to help load the data
		TTree *dataTree;
		Long64_t total_nentries;
	
		// we define the variables that we need in the thread functions
		int kDim;
		int numberEventsToSavePerProcess;
		int nProcess;
		int seedShift;
		Long64_t nentries;
		bool override_nentries;
		bool verbose;
		string varString;
                string standardizationType;
                bool redistributeBkgSigFits;
                bool doKRandomNeighbors;
		std::chrono::time_point<std::chrono::high_resolution_clock> start2;
                
                // These block of variables will be used to hold the initialization parameters. In the Q-factor paper they use 3 different initializations which
                // correspond to 100% bkg, 50/50, and 100% sig. If we want to do this here, the yields in the bkg and signal need to be modified. These vectors
                // will hold that information
		std::vector<double> initPars;
                std::vector<double> redistributeFactorBkg;
                std::vector<double> redistributeFactorSig;
                parameterLimits parLimits;

                // initialize vectors to hold the discriminating and phase space variables
		std::vector<double> discrimVars; 
		std::vector<double> sideBandVars; 
                std::vector<std::vector<double>> phaseSpaceVars;
                std::vector<double> phaseSpaceVars0;
		std::vector<double> AccWeights; 
		std::vector<double> sbWeights; 
		std::vector<ULong64_t> spectroscopicComboIDs; 
                // Not all combinations will be a valid pairing. Suppose we only care about spectroscopically unique pairs, then we can fill phasePoint2PotentailNeighbor with
                // only unique combos.
		std::vector<int> phasePoint2PotentailNeighbor; 
	
	public:
		QFactorAnalysis(int kDim1, string varString1, string standardizationType1, bool redistributeBkgSigFits1, bool doKRandomNeighbors1, int numberEventsToSavePerProcess1, int nProcess1, int seedShift1, Long64_t nentries1, bool override_nentries1, bool verbose1){ 
			cout << "Constructed QFactorAnalysis class..." << endl;
			kDim=kDim1;
			numberEventsToSavePerProcess=numberEventsToSavePerProcess1;
			nProcess = nProcess1;
			seedShift=seedShift1;
			nentries=nentries1;
			override_nentries=override_nentries1;
			verbose=verbose1;
			start2 = std::chrono::high_resolution_clock::now();
			varString=varString1;
                        standardizationType=standardizationType1;
                        redistributeBkgSigFits=redistributeBkgSigFits1;
                        doKRandomNeighbors=doKRandomNeighbors1;
	
		        // Use these vectors to import all the data to RAM instead of reading from root file
                        // First we should reserve the space so there is no resizing of the vectors
			discrimVars.reserve(nentries);
                        std::vector<double> emptyVec;
                        for (auto it=0; it<dim; ++it){
                            // copy over emptyVec and then expand
                            phaseSpaceVars.push_back(emptyVec);
                            phaseSpaceVars[it].reserve(nentries);
                        }
			sideBandVars.reserve(nentries);
			AccWeights.reserve(nentries);
			sbWeights.reserve(nentries);
			spectroscopicComboIDs.reserve(nentries);

			// will hold all the ids of the unique combos
			phasePoint2PotentailNeighbor.reserve(nentries);
		}
		void loadTree(string rootFileLoc, string rootTreeName);
		void loadFitParameters(string fitLocation);
		void loadData();
		void runQFactorThreaded();

};


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
            if (varName.at(0) != '#'){ 
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
	set<Int_t> setUsedSpectroscopicIDs;
	for (Int_t ientry=0; ientry<nentries; ientry++){ 
	    if ( setUsedSpectroscopicIDs.find( spectroscopicComboIDs[ientry] ) == setUsedSpectroscopicIDs.end() ) {
	        setUsedSpectroscopicIDs.insert( spectroscopicComboIDs[ientry] );
	        phasePoint2PotentailNeighbor.push_back(ientry);
	    }
	    ///phasePoint2PotentailNeighbor.push_back(ientry);
	}
        cout << "\n\n--------------------------------------" << endl;
	cout << phasePoint2PotentailNeighbor.size() << "/" << nentries << " are used as potential neighbors" << endl;
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

	// [=] refers to a capture list which is used by this lambda expression. The lambda gets a copy of all the local variables that it uses when it is created. If we
	// just use [] we will get an error since the lambda will have no idea what these variables are	
	auto f = [=](int iProcess){
		int iFit = 0;
		int iThread = iProcess;
		// ----------------------------------------
		// Open up a root file to save the q-factors and other diagnostics to it
		// ---------------------------------------
                // Initializing some variables we can track during the q-value extraction
        	//double comboStd; 
        	double chiSq;
        	double chiSq_pi0;
        	//double comboStd2; 
        	ULong64_t flatEntryNumber;
        	double bestChiSq;
		double worstChiSq;
        	double best_qvalue;

                // Saving the results along with some diagnostics
		TBranch* b_sbWeight;
		TBranch* b_flatEntryNumber;
        	TFile *resultsFile = new TFile(("logs/"+fileTag+"/results"+to_string(iThread)+".root").c_str(),"RECREATE");
        	TTree* resultsTree = new TTree("resultsTree","results");
        	resultsTree->Branch("flatEntryNumber",&flatEntryNumber,"flatEntryNumber/l");
        	resultsTree->Branch("qvalue",&best_qvalue,"qvalue/D");
        	resultsTree->Branch("bestChiSq",&bestChiSq,"bestChiSq/D");
        	resultsTree->Branch("worstChiSq",&worstChiSq,"worstChiSq/D");
        	//resultsTree->Branch("combostd",&comboStd,"combostd/D");
        	//resultsTree->Branch("combostd2",&comboStd2,"combostd2/D");
        	cout << "Set up branch addresses" << endl;

		// Define some needed variables like canvases, histograms, and legends
        	auto legend_init = new TLegend(0.1,0.7,0.4,0.9);
        	auto legend_fit = new TLegend(0.1,0.7,0.4,0.9);

    		TCanvas *allCanvases = new TCanvas(("anyHists"+to_string(iThread)).c_str(),"",1440,900);
    		TCanvas *allCanvases_badFit = new TCanvas(("anyHists_badFit"+to_string(iThread)).c_str(),"",1440,900);
		cout << "Creating canvas " << iThread << endl;
        	TLine* discrimVarLine;
		TLine* qSigLine;
		TLine* qBkgLine;
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

                double _SET_DiscrimBinLower = 0.25;
                double _SET_DiscrimBinUpper = 0.85;
                double _SET_DiscrimBinNum = 200;
                double _SET_DiscrimFitLower = 0.3;
                double _SET_DiscrimFitUpper = 0.8;

		binRangeEta={_SET_DiscrimBinNum,_SET_DiscrimBinLower,_SET_DiscrimBinUpper};
		fitRangeEta={_SET_DiscrimFitLower,_SET_DiscrimFitUpper};

        	// Going to keep a track of the fits to see how they move around
        	// TF1 *showInit[3];
        	// TF1 *showConv[3];
        	// int colors[3] = {kRed, kBlue, kGreen};
        	//string initNames[3] = {"100% Bkg", "50/50% Bkg Sig", "100% Sig"};
        	// for (int iFit=0; iFit<3; ++iFit){
		//     showInit[iFit] = new TF1(("initFit"+to_string(iFit)).c_str(),fitFunc,fitRangeEta[0],fitRangeEta[1],numDOFbkg+numDOFsig);
		//     showConv[iFit] = new TF1(("convFit"+to_string(iFit)).c_str(),fitFunc,fitRangeEta[0],fitRangeEta[1],numDOFbkg+numDOFsig);
        	// }
                

                // Saving the chiSqs to compare the 3 different initializations
		double chiSqs[3];
        	int nBest100Bkg=0;
        	int nBest100Sig=0;
        	int nBest50Bkg50Sig=0;
        	// if(verbose) {cout << "Made array of fits" << endl;}

        	// Getting the scaling of the flat bkg is a bit tricky. Have to count how much bins we have in our fit range. kDim divided by the bins in fit range is the height of the flat function.
        	double numBinsInFit = (fitRangeEta[1]-fitRangeEta[0])/((binRangeEta[2]-binRangeEta[1])/binRangeEta[0]);
        	double binSize=((binRangeEta[2]-binRangeEta[1])/binRangeEta[0]);

        	if (verbose) {cout << "Eta hist range: " << binRangeEta[0] << ", " << binRangeEta[1] << ", " << binRangeEta[2] << endl;}
        	if (verbose) {cout << "Eta fit range: " << fitRangeEta[0] << ", " << fitRangeEta[1] << endl;}
        	

		TF1 *fit;
        	TF1 *bkgFit;
        	TF1 *sigFit;
		// the main loop where we loop through all events in a double for loop to calculate dij. Iterating through all j we can find the k nearest neighbors to event i.
		// Then we can plot the k nearest neighbors in the discriminating distribution which happens to be a double gaussian and flat bkg. Calculate Q-Value from event
		// i's discriminating variable's value. This value will be plugged into the signal PDF and the total PDF. The ratio of these two values are taken which is the Q-Value.
		//logFile << std::fixed << std::setprecision(6);
		int randomEntry;
		int saveN_badEvents=1;
		int savedN_badEvents=0;
        	for (int ientry=lowest_nentry; ientry<largest_nentry; ientry++){ 
			double unshiftedEntry = (double)(ientry-lowest_nentry);
                        if(batchEntries<20){ cout << "Due to how we output the progress of the threads, we must have > 20 nentries PER nProcess\nEXITING\nEXITING" << endl; exit(0);}
			int percentOfBatchEntries = (int)(batchEntries/20);
			if ( (ientry-lowest_nentry) % percentOfBatchEntries == 0) { 
				cout << "(Process " << iProcess << ") Percent done: " << (int)round( unshiftedEntry/(largest_nentry-lowest_nentry)*100 ) << "%" << endl;
		       	}
			flatEntryNumber=ientry;
			//cumulativeStd stdCalc(kDim);
			//cumulativeStd stdCalc2(kDim);
			
			auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
			auto duration_beginEvent = std::chrono::high_resolution_clock::now();
			if(verbose) { logFile << "Starting event " << ientry << "/" << largest_nentry << " ---- Global Time: " << duration2 << "ms" << endl; }
			
			// clear data from previous events
			mapDistToEntry.clear();
			distances.clear();
			discriminatorHist = new TH1F(("discriminatorHist"+to_string(iThread)).c_str(),"",binRangeEta[0],binRangeEta[1],binRangeEta[2]);
			
			for ( int iVar=0; iVar<dim; ++iVar ){
				phasePoint1[iVar] = phaseSpaceVars[iVar][ientry];
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
			    for (int jentry : phasePoint2PotentailNeighbor) {  
			    	if ( verbose_outputDistCalc ) { cout << "event i,j = " << ientry << "," << jentry << endl;} 
			    	
			    	for ( int iVar=0; iVar<dim; ++iVar ){
			    	   	phasePoint2[iVar] = phaseSpaceVars[iVar][jentry];
			    	}
			    	if (spectroscopicComboIDs[jentry] != spectroscopicComboIDs[ientry]){
			    	        distance = calc_distance(dim,phasePoint1,phasePoint2,verbose_outputDistCalc);
			    	        distKNN.insertPair(make_pair(distance,jentry));
			    	}
			    	//if ( verbose) { 
			    	//	cout << "CURRENT SET: " << endl;
			    	//	for(auto elem : mapDistToEntry){
			    	//		std::cout << elem.first << " " << elem.second << "\n";
			    	//	}
			    	//}
			    }
                        }
			duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
			if(verbose){logFile << "\tFound neighbors: " << duration2 << "ms" << endl; }
			if (distKNN.kNN.size() != kDim){ cout << "size of distKNN is not equal to kDim! size,kDim="<< distKNN.kNN.size() << "," << kDim 
			    << "\n    -- if size is 1 less than kDim it is probably because kDim=nentries and event i cannot be a neighbor to itself" << 
			    "\n    -- if size != kDim it could also mean that the number of spectroscopically unique neighbors reduces the number of poential neighbors below kDim" << endl;}
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
			        discriminatorHist->Fill(discrimVars[newPair.second],weight);

			        //stdCalc.insertValue(discrimVars[newPair.second]);
			        //cout << "(" << newPair.first << ", " << newPair.second << ")"; 
			        //cout << endl; 
			}
			// We will not calculate std for now, but lets just leave the structure here
			//comboStd = 1;// stdCalc.calcStd();
			//comboStd2 = 1;//stdCalc2.calcStd();
			
			duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
			if(verbose){logFile <<	"\tFilled neighbors: " << duration2 << "ms" << endl;}
			
			// Building the fit functions. We have to set some parameter limits to try to guide Minuit. We choose a double gaussian since for whatever reason the discrimVar has a asymmetry
			// such that there are more events to the left side of the peak. The second gaussian will have lower amplitude and a large sigma with a left shifted mean. We also want to 
			// fix the amplitudes and constnat background to be strictly positive and choose kDim as the max value ( which is reached only if all values are filled into a single bin)
			
			// We will always do 3 fits per event. In the original paper on Q-factors they use 100% bkg, %100 signal, 50% bkg and 50% signal.
			//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
			
			bestChiSq=DBL_MAX;
			worstChiSq=DBL_MIN;
			int best_iFit=0;
			int countSig;
			int countBkg;

		        Double_t par[numDOFbkg+numDOFsig]; // needed to calculate the qvalue
		        Double_t parBest[numDOFbkg+numDOFsig]; // needed to calculate the qvalue
		        Double_t parFlat[numDOFbkg]; // needed to calculate the qvalue
			
			// /////////////////////////////////////////
			// Calcuclate q-value
			// /////////////////////////////////////////
			// We use a normalized gaussian and a flat function. 
			{
				// have to use a mutex here or else we get some weird error
				//R__LOCKGUARD(gGlobalMutex);
				fit = new TF1(("fit"+to_string(iThread)).c_str(),fitFunc,fitRangeEta[0],fitRangeEta[1],numDOFbkg+numDOFsig);
				bkgFit = new TF1(("bkgFit"+to_string(iThread)).c_str(),background,fitRangeEta[0],fitRangeEta[1],numDOFbkg);
				sigFit = new TF1(("sigFit"+to_string(iThread)).c_str(),signal,fitRangeEta[0],fitRangeEta[1],numDOFsig);
			}
                        // saving the initializations to see if they sort of make sense
                        int numFits;
                        if(redistributeBkgSigFits){ 
                            numFits=3;
                        }
                        else{
                            numFits=1;
                        }
                        std::vector<TF1*> initFits;
                        TF1* initFit;

			duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
			if(verbose){logFile <<	"\tInitialized fits: " << duration2 << "ms" << endl;}

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
				initFit = new TF1(("initFit"+to_string(iFit)+to_string(iThread)).c_str(),fitFunc,fitRangeEta[0],fitRangeEta[1],numDOFbkg+numDOFsig);
				initFit->SetParameters(&initParsVaryPercentSig[0]);
                                initFits.push_back(initFit);
                                fit->SetParameters(&initParsVaryPercentSig[0]);


                                for (int iPar=0; iPar<parLimits.numDOF; ++iPar){
                                    fit->SetParLimits(iPar,parLimits.lowerParLimits[iPar], parLimits.upperParLimits[iPar]);
                                }
                                

				{
					// have to use a mutex here or else we get some weird error
					R__LOCKGUARD(gGlobalMutex);
					discriminatorHist->Fit(("fit"+to_string(iThread)).c_str(),"RQBNL"); // B will enforce the bounds, N will be no draw
				}
				fit->GetParameters(par);
				bkgFit->SetParameters(par);
				sigFit->SetParameters(&par[numDOFbkg]);
				qvalue=sigFit->Eval(discrimVars[ientry])/fit->Eval(discrimVars[ientry]);
				
				// Actually a lot of the times the value -> -0 so we should just make it zero at this point rather than refitting....
				if ( abs(qvalue) < 0.00001 ){	
					qvalue = 0;
				}	

				// /////////////////////////////////////////
				// need to calcuclate new q-value since it is out of bounds
				// /////////////////////////////////////////
				if (qvalue>1 || qvalue<0){
					cout << "Using less complex fit instead on event: " << ientry << " -- QFactor = " << qvalue << endl;
					// first we will save the bad event to get a sample then fix the linear component of the bkg
					if ( savedN_badEvents < saveN_badEvents ) {
						allCanvases_badFit->cd();
        	        			fit->SetLineColor(kRed+2);
  		        			bkgFit->SetFillColor(kMagenta+2);
        	        			bkgFit->SetLineColor(kMagenta+2);
  		        			bkgFit->SetFillStyle(3004);
  		        			sigFit->SetFillColor(kBlue+2);
        	        			sigFit->SetLineColor(kBlue+2);
  		        			sigFit->SetFillStyle(3005);
		        			discriminatorHist->Draw();
        	        			fit->Draw("SAME");
  		        			bkgFit->Draw("SAME FC");
  		        			sigFit->Draw("SAME FC");
						discrimVarLine = new TLine(discrimVars[ientry],0,discrimVars[ientry],kDim);
						discrimVarLine->SetLineColor(kOrange);
		        			discrimVarLine->Draw("same");
						discriminatorHist->SetTitle(("q-value: "+to_string(qvalue)).c_str());
				                {
				                    R__LOCKGUARD(gGlobalMutex);
				                    allCanvases_badFit->SaveAs(("histograms/"+fileTag+"/bad-Mass-event"+std::to_string(ientry)+".root").c_str());
                                                }
						++savedN_badEvents;
					}
					
                                        fit->SetParameters(&initParsVaryPercentSig[0]);
                                        for (auto iPar : parLimits.zeroTheseParsOnFail){
                                            cout << "Fixing to 0 Par" << iPar << endl;
                                            fit->SetParameter(iPar,0);
                                            fit->FixParameter(iPar,0);
                                        } 
					discriminatorHist->Fit(("fit"+to_string(iThread)).c_str(),"RQBNL"); // B will enforce the bounds, N will be no draw
					fit->GetParameters(par);
					bkgFit->SetParameters(par);
					sigFit->SetParameters(&par[numDOFbkg]);

					qvalue=sigFit->Eval(discrimVars[ientry])/fit->Eval(discrimVars[ientry]);
					if (qvalue>1 || qvalue<0){
					    cout << "Not sure why qvalue is still >1 or <0. Need to fix this!" << endl;
					    cout << "These are the parameters:"<<endl;
					    for ( double parVal : par ){
							cout << " " << parVal << endl;
					    }
                                            cout << " Is your kDim very small?" << endl;
					    cout << " **************** BREAKING (q-value out of bounds) ****************** " << endl;
					    exit(0);
					}
			                duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
			                if(verbose){logFile <<	"\tFinished Flat fit: " << duration2 << "ms" << endl;}
				}

				// now that the q-value is found we can get the chiSq and save the parameters with the best chiSq
				chiSq = fit->GetChisquare()/(fit->GetNDF());
				
				// save some fit and the chiSq for all the fits
				//showConv[iFit]->SetParameters(par); // save all the fits
				//chiSqs[iFit]=chiSq;
			
				//if (verbose) { logFile << "\tcurrent ChiSq, best ChiSq: " << chiSq << ", " << bestChiSq << endl; }
				
				if (chiSq < bestChiSq){
					best_qvalue = qvalue;
					bestChiSq=chiSq;
					best_iFit=iFit;
					for (int i=0; i < sizeof(par)/sizeof(Double_t); ++i){
						parBest[i]=par[i];
					}
					if(qvalue>1 || qvalue<0) {
					      cout << "qvalue out of bounds!\n-------------" << endl;
					      for (double parVal : par){
					          cout << parVal << endl;
					      } 
					      //exit(0);
					}
				} 
				if (chiSq > worstChiSq){
					worstChiSq = chiSq;
				}
				duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
				if(verbose){logFile << "\t("+to_string(iFit+1)+"st init config) Fitted the reference distribution: " << duration2 << "ms" << endl; }
			}

			if(verbose){logFile << "\tDelta b/w best and worst chiSq = " << to_string(bestChiSq-worstChiSq) << ": " << duration2 << "ms" << endl; }
				

			// /////////////////////////////////////////
			// Calculating some chiSq differences 
			// /////////////////////////////////////////
			 
			//if ( best_iFit == 0){ ++nBest100Bkg; }
			//else if ( best_iFit == 1){ ++nBest100Sig; }
			//else if ( best_iFit == 2){ ++nBest50Bkg50Sig; }
			
			// Here we draw the histograms that were randomly selected
			allCanvases->Clear();
                        // LEFT HISTOGRAM IS FOR THE FITTED AND THE RIGHT HISTOGRAM IS FOR THE INITIALIZATIONS
			allCanvases->Divide(2,1);
			if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end()) {
			        legend_init->Clear();
			        legend_fit->Clear();

				if(verbose) { cout << "Keeping this event" << endl; }
				discriminatorHist->SetTitle(("BEST VALS:  QValue="+to_string(best_qvalue)+"  ChiSq="+to_string(bestChiSq) + " (#Delta="+to_string(bestChiSq-worstChiSq)+ ") iFit="+to_string(best_iFit)).c_str() );
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



        	          	//allCanvases->cd();
        	          	allCanvases->cd(1);
		          	discriminatorHist->Draw();
        	          	drawText(parBest,numDOFbkg+numDOFsig,"par",sigFit->Eval(discrimVars[ientry]),bkgFit->Eval(discrimVars[ientry]),fit->Eval(discrimVars[ientry]));
		          	discrimVarLine->Draw("same");

        	          	fit->Draw("SAME");
  		          	bkgFit->Draw("SAME FC");
  		          	sigFit->Draw("SAME FC");
				qBkgLine->Draw("SAME");
				qSigLine->Draw("SAME");
                                legend_fit->Draw();
				// INTERESTING, IF I WERE TO SAVE THE CANVAS AS A ROOT FILE I GET AN ERROR IF I PUT THE SAME HISTOGRAM ON TWO DIFFERENT PADS. SEEMS LIKE THE CANVAS
				// SAVES A TList OF HISTS+TF1'S AND IF THERE ARE MULTIPLE CALLS TO THE SAME HISTOGRAM IT MIGHT DELETE THE HISTOGRAM AFTER SEEING IT FOR THE FIRST TIME AND
				// THEN IT WOULD NOT BE ABLE TO FIND THE HISTOGRAM AGAIN THE SECOND TIME AROUND. WE HAVE TO CLONE THE HISTOGRAM FIRST AND THEN SAVE THE ROOT FILE SO THE CANVAS
				// ARE TWO DIFFERENT ELEMENTS.
        	          	allCanvases->cd(2);
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
                                }
				TF1* initFit4 = new TF1(("initFit4"+to_string(iThread)).c_str(),fitFunc,fitRangeEta[0],fitRangeEta[1],numDOFbkg+numDOFsig);
				initFit4->SetParameters(initPars[0],initPars[1],initPars[2],initPars[3],initPars[4]);
				initFit4->SetLineColor(kOrange+1);
				initFit4->Draw("SAME");

                                legend_init->AddEntry(initFits[0],"100% bkg");
                                legend_init->AddEntry(initFits[1],"50/50 bkg/sig");
                                legend_init->AddEntry(initFits[2],"100% sig");
                                legend_init->AddEntry(initFit4,"Scaled Full Fit");
                                legend_init->Draw();
				// need to save as a root file first then convert to pngs or whatever. Seems like saveas doesnt like threaded since the processes might make only one pdf converter
				// or whatever and maybe if multiple threads calls it then a blocking effect can happen
				cout << "Choosing to save event " << ientry << endl;
                                {
				    R__LOCKGUARD(gGlobalMutex);
				    allCanvases->SaveAs(("histograms/"+fileTag+"/Mass-event"+std::to_string(ientry)+".root").c_str());
                                }
				

				if(verbose){
				     duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
				     logFile << "\tSaved this histogram since it was randomly selected: " << duration2 <<  "ms" << endl;
				}
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
        	//cout << "Number of times 100% bkg initialization was the best: " << nBest100Bkg << endl;
        	//cout << "Number of times 100% sig initialization was the best: " << nBest100Sig << endl;
        	//cout << "Number of times 50/50 bkg/sig initialization was the best: " << nBest50Bkg50Sig << endl;
        	
        	
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


#endif
