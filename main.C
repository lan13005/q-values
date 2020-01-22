#include "main.h"

//bool verbose = false;
//int kDim=200;
//int numberEventsToSavePerProcess=5;
//bool override_nentries=true;
//Long64_t nentries=100;
//int nProcess=5;
//int seedShift=123125;

bool useEta=true;
using namespace RooFit;
string detector="bcal";

//void main(int iProcess, int kDim, int numberEventsToSavePerProcess, int nProcess, int seedShift, Long64_t nentries, bool override_nentries, bool verbose){
int main( int argc, char* argv[] ){
	// This suppresses all the "info" messages from roofit. Like saying that it will use numeric integrator for the generic pdf we define
	// https://root-forum.cern.ch/t/suppressing-info-messages/14642/6
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
        cout << "Num Vars="<< dim << endl;

    	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
        //int iProcess = std::stoi(argv[0]);
        //int kDim= std::stoi(argv[1]);
        //int numberEventsToSavePerProcess= std::stoi(argv[2]);
        //int nProcess= std::stoi(argv[3]);
        //int seedShift= std::stoi(argv[4]);
        //Long64_t nentries= std::stoll(argv[5]);
        int iProcess = atoi(argv[1]);
        int kDim= atoi(argv[2]);
        int numberEventsToSavePerProcess= atoi(argv[3]);
        int nProcess= atoi(argv[4]);
        int seedShift= atoi(argv[5]);
        Long64_t nentries= atoi(argv[6]);
        bool override_nentries;
        if ( atoi(argv[7]) ==1 ){ 
            override_nentries=true;
        }
        else{ 
            override_nentries=false;
        }
        bool verbose;
        if ( atoi(argv[8])==1 ){ verbose=true;}
        else{ verbose=false;}
        std::string varString=argv[9];
        cout << "----------------------------" << endl;
        cout << "iProcess: " << iProcess << endl;
        cout << "kDim: " << kDim << endl;
        cout << "numberEventsToSavePerProcess: " << numberEventsToSavePerProcess << endl;
        cout << "nProcess: " << nProcess << endl;
        cout << "seedShift: " << seedShift << endl;
        cout << "nentries: " << nentries << endl; 
        cout << "override_nentries: " << override_nentries << endl;
        cout << "verbose: " << verbose  << endl; 
    
        parseVarString parse(varString);
        parse.parseString();
        cout << "----------------------------" << endl;
        cout << "----------------------------" << endl;
        for (int iVar=0; iVar<parse.varStringSet.size(); ++iVar){
            cout << "var" << std::to_string(iVar) << ": " << parse.varStringSet[iVar] << endl;
        }
        cout << "----------------------------" << endl;
        cout << "----------------------------" << endl;


	// Starting timing
	//clock_t start;
	//double duration;
	auto start2 = std::chrono::high_resolution_clock::now();
	//start = clock();
	
	
	///////////////////// INITIALIZING THE PARAMETERS ////////////////////
	///////// Pi0 --- NOT SB subtracted
	double fittedConst_pi0 = 1;
	double fittedLinear_pi0 = 1;
	double fittedAmp_pi0 = 1;
	double peak_pi0 = 1;
	double sigma_pi0 = 1;
	double ampRatio_pi0 = 1;
	double sigmaRatio_pi0 = 1;
	// ------------- tLT1 --- SB subtracted --- AccSub
	double fittedConst_eta;
	double fittedLinear_eta;
	double fittedAmp_eta;
	double peak_eta;
	double sigma_eta;
	double ampRatio_eta;
	double sigmaRatio_eta;

	string varName;
	double varVal;
	ifstream inFile;
	inFile.open(("fitResults/etaFitNoAccSub_"+detector+".txt").c_str());
	while (inFile >> varName >> varVal){
		if (varName == "const"){ fittedConst_eta = varVal; }
		if (varName == "linear"){ fittedLinear_eta = varVal;}
		if (varName == "amp1"){ fittedAmp_eta = varVal;}
		if (varName == "mass"){ peak_eta = varVal;}
		if (varName == "sigma1"){ sigma_eta = varVal;}
		if (varName == "ampRatio"){ ampRatio_eta = varVal;}
		if (varName == "sigmaRatio"){ sigmaRatio_eta = varVal;}
	}
	cout << detector << endl;
	cout << "const: " << fittedConst_eta << endl;
	cout << "linear: " << fittedLinear_eta << endl;
	cout << "amp1: " << fittedAmp_eta << endl;
	cout << "mass: " << peak_eta << endl;
	cout << "sigma1: " << sigma_eta << endl;
	cout << "ampRatio: " << ampRatio_eta << endl;
	cout << "sigmaRatio: " << sigmaRatio_eta << endl;
	inFile.close();
	
	// setting up some basic root stuff and getting the file and tree
	//TFile* dataFile=new TFile("pi0eta_a0_recotreeFlat_DSelector.root");
	TFile* dataFile=new TFile(("pi0eta_"+detector+"_tLT1treeFlat_DSelector.root").c_str());
	TTree *dataTree;
	dataFile->GetObject(("pi0eta_"+detector+"_tLT1tree_flat").c_str(),dataTree);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
    	TCanvas *allCanvases_badFit = new TCanvas("anyHists_badFit","",1440,900);
        auto legend_init = new TLegend(0.1,0.8,0.3,0.9);
        auto legend_conv = new TLegend(0.1,0.8,0.3,0.9);
        TLine* etaLine;
        TLine* pi0Line;
	TLine* qSigLine;
	TLine* qBkgLine;
	TH1F* discriminatorHist;
	TH1F* discriminatorHist2;
	double Meta;
	double Mpi0;
	double Mpi0eta;
	double cosTheta_X_cm;
	double phi_X_cm;
	double cosTheta_eta_gj;
	double phi_eta_gj;
	double cosThetaHighestEphotonIneta_gj;
	double cosThetaHighestEphotonInpi0_cm;
        double vanHove_x;
        double vanHove_y;
        double vanHove_omega;
        double pi0_energy;
        double mandelstam_tp;
        ULong64_t eventNumber;
        double uniqueComboID;
        double AccWeight;
        bool isUniqueEtaB;
        bool isUniquePi0B;
        bool isUniquePi0EtaB;
        ULong64_t  spectroscopicComboID;

	// tree contains measured and kinfit values for these, except phi_X_cm
	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("Mpi0",&Mpi0);
	dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
        dataTree->SetBranchAddress("cosTheta_X_cm", &cosTheta_X_cm); 
        dataTree->SetBranchAddress("phi_X_cm",&phi_X_cm); 
        dataTree->SetBranchAddress("cosTheta_eta_gj",&cosTheta_eta_gj);
        dataTree->SetBranchAddress("phi_eta_gj",&phi_eta_gj); 

	// doesn't contain measured values I think
        dataTree->SetBranchAddress("cosThetaHighestEphotonIneta_gj",&cosThetaHighestEphotonIneta_gj);
        dataTree->SetBranchAddress("cosThetaHighestEphotonInpi0_cm",&cosThetaHighestEphotonInpi0_cm);
        dataTree->SetBranchAddress("vanHove_x",&vanHove_x);
        dataTree->SetBranchAddress("vanHove_y",&vanHove_y);
        dataTree->SetBranchAddress("vanHove_omega",&vanHove_omega);
        dataTree->SetBranchAddress("pi0_energy", &pi0_energy);
        dataTree->SetBranchAddress("mandelstam_tp", &mandelstam_tp);
        dataTree->SetBranchAddress("uniqueComboID",&uniqueComboID);
        dataTree->SetBranchAddress("event",&eventNumber);
        dataTree->SetBranchAddress("AccWeight",&AccWeight);
        dataTree->SetBranchAddress("isNotRepeated_eta",&isUniqueEtaB);
        dataTree->SetBranchAddress("isNotRepeated_pi0",&isUniquePi0B);
        dataTree->SetBranchAddress("isNotRepeated_pi0eta",&isUniquePi0EtaB);
        dataTree->SetBranchAddress("spectroscopicComboID",&spectroscopicComboID);



	const Long64_t total_nentries = (const Long64_t)dataTree->GetEntries();
	if (!override_nentries){
		nentries=dataTree->GetEntries();
	}
	cout << "Chosen Total Entries: " << nentries << endl;


	// opening a file to write my log data to
    	ofstream logFile;
    	logFile.open(("logs/processLog"+detector+"_"+to_string(iProcess)+".txt").c_str());
	//logFile << "Event\tQ-Value\tChiSq\tMpi0" << endl;
	


	// importing all the data to RAM instead of reading from root file
	std::vector<double> Metas; Metas.reserve(nentries);
        std::vector<double> Mpi0s; Mpi0s.reserve(nentries);
        std::vector<double> Mpi0etas; Mpi0etas.reserve(nentries);
        std::vector<double> cosTheta_X_cms; cosTheta_X_cms.reserve(nentries);
        //std::vector<double> phi_X_cms; phi_X_cms.reserve(nentries);
        std::vector<double> cosTheta_eta_gjs; cosTheta_eta_gjs.reserve(nentries);
        std::vector<double> phi_eta_gjs; phi_eta_gjs.reserve(nentries);
        //std::vector<double> cosThetaHighestEphotonIneta_gjs; cosThetaHighestEphotonIneta_gjs.reserve(nentries);
        //std::vector<double> cosThetaHighestEphotonInpi0_cms; cosThetaHighestEphotonInpi0_cms.reserve(nentries);
        //std::vector<double> pi0_energies; pi0_energies.reserve(nentries);
        //std::vector<double> mandelstam_tps; mandelstam_tps.reserve(nentries);
	//std::vector<double> vanHove_xs; vanHove_xs.reserve(nentries);
	//std::vector<double> vanHove_ys; vanHove_ys.reserve(nentries);
	//std::vector<double> vanHove_omegas; vanHove_omegas.reserve(nentries);
        std::vector<double> AccWeights; AccWeights.reserve(nentries);
        std::vector<double> sbWeights; sbWeights.reserve(nentries);
        std::vector<ULong64_t> spectroscopicComboIDs; spectroscopicComboIDs.reserve(nentries);

	double sbRL = 0.09; // Right sideband left line
	double sbRR = 0.105; // Right sideband right line
	double sigL = 0.12;
	double sigR = 0.15;
	double sbLL = 0.165;
	double sbLR = 0.18;
	double sbWeight;
        
	// We will use a ientry to keep track of which entries we will get from the tree. We will simply use ientry when filling the arrays.  
	for (Long64_t ientry=0; ientry<nentries; ientry++)
	{
        	dataTree->GetEntry(ientry);
                //if ( ientry != eventNumber){ cout << "ientry != eventNumber. The events are in the root tree are out of order or missing!" <<
                //        "\n ientry,eventNumber: " << ientry << ", " << eventNumber << endl; break; }
                //                ***** THE COMBINATIONS COME IN AT RANDOM ORDERS! ***** DOESNT MATTER I THINK
		Metas.push_back(Meta);
		Mpi0s.push_back(Mpi0);
		Mpi0etas.push_back(Mpi0eta);
		cosTheta_X_cms.push_back(cosTheta_X_cm);
		//phi_X_cms.push_back(phi_X_cm);
		cosTheta_eta_gjs.push_back(cosTheta_eta_gj);
		phi_eta_gjs.push_back(phi_eta_gj);
		//cosThetaHighestEphotonIneta_gjs.push_back(cosThetaHighestEphotonIneta_gj);	 
		//cosThetaHighestEphotonInpi0_cms.push_back(cosThetaHighestEphotonInpi0_cm);	 
                //vanHove_xs.push_back(vanHove_x);
                //vanHove_ys.push_back(vanHove_y);
                //vanHove_omegas.push_back(vanHove_omega);
                //pi0_energies.push_back(pi0_energy);
                //mandelstam_tps.push_back(mandelstam_tp);
                AccWeights.push_back(AccWeight);
		if ( Mpi0 > sbRL && Mpi0 < sbRR ) { sbWeight = -1; } 
		else if ( Mpi0 > sbLL && Mpi0 < sbLR ) { sbWeight = -1; } 
		else if ( Mpi0 > sigL && Mpi0 < sigR ) { sbWeight = 1; } 
		else { sbWeight = 0; }
		sbWeights.push_back(sbWeight);

                spectroscopicComboIDs.push_back(spectroscopicComboID);
	}

	int batchEntries = (int)nentries/nProcess;
	int lowest_nentry = iProcess*batchEntries;
	int largest_nentry;
        if (iProcess!=(nProcess-1)) {
            largest_nentry  = (iProcess+1)*batchEntries;
        }
        else {
            largest_nentry = nentries; 
        }

	cout << "nentries we will use for this process: " << lowest_nentry << ", " << largest_nentry << endl;

	// randomly select some events to write histograms for 
	set<int> selectRandomIdxToSave;
	int randomEvent;
	srand(iProcess+seedShift);
	for (int i=0; i<numberEventsToSavePerProcess; i++){
		randomEvent = rand() % (int)batchEntries;
		randomEvent += lowest_nentry;
		selectRandomIdxToSave.insert( randomEvent );
	}
        if(verbose){cout << "randomly selected some events to save" << endl; }


	
        if ( verbose_outputDistCalc ) {
            cout << "Before standarization" << endl;
            for ( int ientry=0 ; ientry < nentries; ientry++){
                cout << cosTheta_X_cms[ientry] << endl;//"," << phi_X_cms[ientry] << endl;
                //cout << phi_X_cms[ientry] <<endl;// "," << cosTheta_eta_gjs[ientry] << endl;
            }
        }
	//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
	//logFile << "Time marker to load data: " << duration << "ms" << endl;

	// outputting the results before and after standardizeArray will show that it works
	// for(auto& cosTheta_X_cm1 : cosTheta_X_cms){ cout << cosTheta_X_cm1 << endl; }
	// **** WE FIRST IMPORT THE DATA INTO A CLASS, SAVES AN INTERNAL COPY
        standardizeArray class_cosTheta_X_cms(cosTheta_X_cms,nentries);
        //standardizeArray class_phi_X_cms(phi_X_cms,nentries);
        standardizeArray class_cosTheta_eta_gjs(cosTheta_eta_gjs,nentries);
        standardizeArray class_phi_eta_gjs(phi_eta_gjs,nentries);
        //standardizeArray class_cosThetaHighestEphotonIneta_gjs(cosThetaHighestEphotonIneta_gjs,nentries);
        //standardizeArray class_cosThetaHighestEphotonInpi0_cms(cosThetaHighestEphotonInpi0_cms,nentries);
        //standardizeArray class_vanHove_xs(vanHove_xs,nentries);
        //standardizeArray class_vanHove_ys(vanHove_ys,nentries);
        //standardizeArray class_vanHove_omegas(vanHove_omegas,nentries);
        //standardizeArray class_pi0_energies(pi0_energies, nentries);
        //standardizeArray class_mandelstam_tps(mandelstam_tps, nentries);
        //standardizeArray class_Mpi0s(Mpi0s,nentries);
        //standardizeArray class_Metas(Metas,nentries);

        // ***** STANDARDIZES THE INTERAL COPY 
        class_cosTheta_X_cms.rangeStandardization();
        //class_phi_X_cms.rangeStandardization();
        class_cosTheta_eta_gjs.rangeStandardization();
        class_phi_eta_gjs.rangeStandardization();
        //class_cosThetaHighestEphotonIneta_gjs.rangeStandardization();
        //class_cosThetaHighestEphotonInpi0_cms.rangeStandardization();
        //class_vanHove_xs.rangeStandardization();
        //class_vanHove_ys.rangeStandardization();
        //class_vanHove_omegas.rangeStandardization();
        //class_pi0_energies.rangeStandardization();
        //class_mandelstam_tps.rangeStandardization();
        //class_cosTheta_X_cms.stdevStandardization();
        //class_phi_X_cms.stdevStandardization();
        //class_cosTheta_eta_gjs.stdevStandardization();
        //class_phi_eta_gjs.stdevStandardization();
        //class_cosThetaHighestEphotonIneta_gjs.stdevStandardization();
        //class_cosThetaHighestEphotonInpi0_cms.stdevStandardization();
        //class_vanHove_xs.stdevStandardization();
        //class_vanHove_ys.stdevStandardization();
        //class_vanHove_omegas.stdevStandardization();
        //class_pi0_energies.stdevStandardization();
        //class_mandelstam_tps.stdevStandardization();
        //class_Mpi0s.stdevStandardization();
        //class_Metas.stdevStandardization();
    
        // ******* REPLACES THE ORIGNAL WITH THE STANDARDIZED COPY INSIDE THE ARRAY
        cosTheta_X_cms = class_cosTheta_X_cms.getVector();
        //phi_X_cms=class_phi_X_cms.getVector();
        cosTheta_eta_gjs=class_cosTheta_eta_gjs.getVector();
        phi_eta_gjs=class_phi_eta_gjs.getVector();
        //cosThetaHighestEphotonIneta_gjs=class_cosThetaHighestEphotonIneta_gjs.getVector();
        //cosThetaHighestEphotonInpi0_cms=class_cosThetaHighestEphotonInpi0_cms.getVector();
        //vanHove_xs=class_vanHove_xs.getVector();
        //vanHove_ys=class_vanHove_xs.getVector();
        //vanHove_omegas=class_vanHove_omegas.getVector();
        //pi0_energies=class_pi0_energies.getVector();
        //mandelstam_tps=class_mandelstam_tps.getVector();
        //Mpi0s=class_Mpi0s.getVector();
        //Metas=class_Metas.getVector();

        map<std::string, std::vector<double>> nameToVec;
        nameToVec["cosTheta_X_cms"] = cosTheta_X_cms;
        //nameToVec["phi_X_cms"] = phi_X_cms; 
        nameToVec["cosTheta_eta_gjs"] = cosTheta_eta_gjs; 
        nameToVec["phi_eta_gjs"] = phi_eta_gjs; 
        //nameToVec["cosThetaHighestEphotonIneta_gjs"] = cosThetaHighestEphotonIneta_gjs; 
        //nameToVec["cosThetaHighestEphotonInpi0_cms"] = cosThetaHighestEphotonInpi0_cms; 
        //nameToVec["vanHove_omegas"] = vanHove_omegas; 
        //nameToVec["vanHove_xs"] = vanHove_xs; 
        //nameToVec["vanHove_ys"] = vanHove_ys; 
        //nameToVec["pi0_energies"] = pi0_energies; 
        //nameToVec["mandelstam_tps"] = mandelstam_tps; 
        //nameToVec["Mpi0s"] = Mpi0s; 
        //nameToVec["Metas"] = Metas; 

        //int nameMapCounter=0;
	//for(auto elem : nameToVec){
        //    phasePoint1[nameMapCounter] = elem.second[ientry];
        //    ++nameMapCounter;
	//}

        if ( verbose_outputDistCalc ) {
	    cout << "After standardization" << endl;
            for ( int ientry=0 ; ientry < nentries; ientry++){
                cout << cosTheta_X_cms[ientry] << endl;//"," << phi_X_cms[ientry] << endl;
                //cout << phi_X_cms[ientry] << endl;//"," << cosTheta_eta_gjs[ientry] << endl;
            }
        }


	// defining some variables we will use in the main loop to get the distances and then the q-values
	map<double, int> mapDistToJ;
	set<double> distances;
	double phasePoint1[dim];
	double phasePoint2[dim];
	double distance;
        double qvalue;
        double best_qvalue;

        distSort_kNN distKNN(kDim);
        pair<double,int> newPair;


        // It is much slower to constantly read a map to get the vector rather than just importing it all into a vector first.
        std::vector< std::vector< double > > varVector;
        int numVars=parse.varStringSet.size();
        //cout << numVars << ", " << dim << endl;
        if (numVars!=dim){ cout << " **** ERROR THE NUMBER OF VARIABLES USED IS NOT THE SAME HERE! *******" << endl; exit(0); }
        for (int iVar=0; iVar<numVars; ++iVar){
            varVector.push_back(nameToVec[parse.varStringSet[iVar]]);
        }

        double comboStd; 
        double chiSq;
        double chiSq_pi0;
        double comboStd2; 
        ULong64_t flatEntryNumber;
        double bestChiSq;
	double chiSq_eta_01;
	double chiSq_eta_02;


	TBranch* b_sbWeight;
	TBranch* b_flatEntryNumber;
        TFile *resultsFile = new TFile(("logs/results"+to_string(iProcess)+".root").c_str(),"RECREATE");
        TTree* resultsTree = new TTree("resultsTree","results");
        resultsTree->Branch("flatEntryNumber",&flatEntryNumber,"flatEntryNumber/l");
        resultsTree->Branch("qvalue",&best_qvalue,"qvalue/D");
        resultsTree->Branch("chisq_eta",&bestChiSq,"chisq_eta/D");
        resultsTree->Branch("chisq_eta_01",&chiSq_eta_01,"chisq_eta/D");
        resultsTree->Branch("chisq_eta_02",&chiSq_eta_02,"chisq_eta/D");
        resultsTree->Branch("combostd",&comboStd,"combostd/D");
        resultsTree->Branch("combostd2",&comboStd2,"combostd2/D");
        if(verbose) {cout << "Set up branch addresses" << endl;}

        std::vector<double> binRange1;
        std::vector<double> fitRange1;
        std::vector<double> binRange2;
        std::vector<double> fitRange2;
	TF1 *fit,*fit2;
        TF1 *bkgFit,*bkgFit2;
        TF1 *sigFit,*sigFit2;
        binRange1={100,0.25,0.8};
        //fitRange1={0.42,0.68};
        fitRange1={0.38,0.65};

        binRange2={100,0.05,0.25};
        //fitRange2={0.05,0.25};
        fitRange2={0.1,0.17};

        // Going to keep a track of the fits to see how they move around
        TF1 *showInit[3];
        TF1 *showConv[3];
        int colors[3] = {kRed, kBlue, kGreen};
	double chiSqs[3];
        string initNames[3] = {"100% Bkg", "50/50% Bkg Sig", "100% Sig"};
        for (int iFit=0; iFit<3; ++iFit){
	    showInit[iFit] = new TF1(("initFit"+to_string(iFit)).c_str(),fitFunc,fitRange1[0],fitRange1[1],numDOFbkg+numDOFsig);
	    showConv[iFit] = new TF1(("convFit"+to_string(iFit)).c_str(),fitFunc,fitRange1[0],fitRange1[1],numDOFbkg+numDOFsig);
        }
        if(verbose) {cout << "Made array of fits" << endl;}

        int nBest100Bkg=0;
        int nBest100Sig=0;
        int nBest50Bkg50Sig=0;
        int nBest100Bkg2=0;
        int nBest100Sig2=0;
        int nBest50Bkg50Sig2=0;
        // Getting the scaling of the flat bkg is a bit tricky. Have to count how much bins we have in our fit range. kDim divided by the bins in fit range is the height of the flat function.
        double numBinsInFit = (fitRange1[1]-fitRange1[0])/((binRange1[2]-binRange1[1])/binRange1[0]);
        double binSize=((binRange1[2]-binRange1[1])/binRange1[0]);
        double numBinsInFit2 = (fitRange2[1]-fitRange2[0])/((binRange2[2]-binRange2[1])/binRange2[0]);
        double binSize2=((binRange2[2]-binRange2[1])/binRange2[0]);

        if (verbose) {cout << "Eta hist range: " << binRange1[0] << ", " << binRange1[1] << ", " << binRange1[2] << endl;}
        if (verbose) {cout << "Eta fit range: " << fitRange1[0] << ", " << fitRange1[1] << endl;}
        if (verbose) {cout << "Pi0 hist range: " << binRange2[0] << ", " << binRange2[1] << ", " << binRange2[2] << endl;}
        if (verbose) {cout << "Pi0 fit range: " << fitRange2[0] << ", " << fitRange2[1] << endl;}
        

	// KinFit values
	
	
	// Measured values
	double scaleFactor = (double)kDim/total_nentries;
	cout << "scaleFactor: " << scaleFactor << endl;
	double fittedConst;
	double fittedLinear;
	double fittedAmp;
	double ampRatio;
	double widthRatio;
	double peakLoc;
	double sigValue;

	if (useEta) { 
		fittedConst = fittedConst_eta;
		fittedLinear = fittedLinear_eta;
		fittedAmp = fittedAmp_eta;
		ampRatio = ampRatio_eta;
		widthRatio = sigmaRatio_eta;
		peakLoc = peak_eta;
		sigValue = sigma_eta;
	}
	else { 
		fittedConst = fittedConst_pi0;
		fittedLinear = fittedLinear_pi0;
		fittedAmp = fittedAmp_pi0;
		ampRatio = ampRatio_pi0;
		widthRatio = sigmaRatio_pi0;
		peakLoc = peak_pi0;
		sigValue = sigma_pi0;
	}
	cout << "Scaled Const: " << scaleFactor*fittedConst << endl;
	cout << "Scaled Linear: " << scaleFactor*fittedLinear << endl;
	cout << "Scaled Amp: " << scaleFactor*fittedAmp << endl;

	std::vector<double> par0 = { scaleFactor*fittedConst, scaleFactor/2*fittedConst, 0 }; 
	std::vector<double> par1 = { scaleFactor*fittedLinear, scaleFactor/2*fittedLinear, 0 }; 
	std::vector<double> par2 = { 0, scaleFactor/2*fittedAmp, scaleFactor*fittedAmp }; 

        showInit[0]->SetParameters(par0[0],par1[0],par2[0],peakLoc,sigValue, ampRatio, widthRatio);
        showInit[1]->SetParameters(par0[1],par1[1],par2[1],peakLoc,sigValue, ampRatio, widthRatio);
        showInit[2]->SetParameters(par0[2],par1[2],par2[2],peakLoc,sigValue, ampRatio, widthRatio);

        // phasePoint1 will consider all events from lowest to largest since these will be our attached q values. phasePoint2 on the other hand will only look at a subset of the events where the
        // elements must be spectroscopically distinct, i.e. the 4 photons in consideration are different. Not sure if this is the value I should consider or is it better to do a pair of maps
        // where I am tracking the two photon pairs that make up the eta and pi0.         
        set<Int_t> setUsedSpectroscopicIDs;
        std::vector<int> phasePoint2PotentailNeighbor; phasePoint2PotentailNeighbor.reserve(nentries);
        for (Int_t ientry=0; ientry<nentries; ientry++){ 
            if ( setUsedSpectroscopicIDs.find( spectroscopicComboIDs[ientry] ) == setUsedSpectroscopicIDs.end() ) {
                setUsedSpectroscopicIDs.insert( spectroscopicComboIDs[ientry] );
                phasePoint2PotentailNeighbor.push_back(ientry);
            }
            ///phasePoint2PotentailNeighbor.push_back(ientry);
        }
        cout << phasePoint2PotentailNeighbor.size() << "/" << nentries << " are used as potential neighbors" << endl;

	// the main loop where we loop through all events in a double for loop to calculate dij. Iterating through all j we can find the k nearest neighbors to event i.
	// Then we can plot the k nearest neighbors in the discriminating distribution which happens to be a double gaussian and flat bkg. Calculate Q-Value from event
	// i's discriminating variable's value. This value will be plugged into the signal PDF and the total PDF. The ratio of these two values are taken which is the Q-Value.
	//logFile << std::fixed << std::setprecision(6);
	int randomEntry;
	int saveN_badEvents=1;
	int savedN_badEvents=0;
        for (int ientry=lowest_nentry; ientry<largest_nentry; ientry++){ 
		//if (ientry % 100 == 0) { logFile << "Current event: " << ientry << "/" << largest_nentry <<  endl; }
		legend_init->Clear();
		legend_conv->Clear();
		flatEntryNumber=ientry;
		cumulativeStd stdCalc(kDim);
		cumulativeStd stdCalc2(kDim);
		
		if(verbose) { cout << "Getting next event!\n--------------------------------\n" << endl;  }
		auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
		auto duration_beginEvent = std::chrono::high_resolution_clock::now();
		if(verbose2) { logFile << "Starting event " << ientry << "/" << largest_nentry << " ---- Time: " << duration2 << "ms" << endl; }
		
		// clear data from previous events
		mapDistToJ.clear();
		distances.clear();
		allCanvases->Clear();
		allCanvases_badFit->Clear();
		allCanvases->Divide(2,1);
		if(useEta){
			discriminatorHist = new TH1F("","",binRange1[0],binRange1[1],binRange1[2]);
			discriminatorHist2 = new TH1F("","",binRange2[0],binRange2[1],binRange2[2]);
		}
		else{
			discriminatorHist = new TH1F("","",binRange2[0],binRange2[1],binRange2[2]);
			discriminatorHist2 = new TH1F("","",binRange1[0],binRange1[1],binRange1[2]);
		}
		
		for ( int iVar=0; iVar<numVars; ++iVar ){
			phasePoint1[iVar] = varVector[iVar][ientry];
		}
		//for (int jentry=0; jentry<200;++jentry) {  
		//      randomEntry = rand() % nentries;
		//      distKNN.insertPair(make_pair(1, randomEntry) );
		//}
		for (int jentry : phasePoint2PotentailNeighbor) {  
			if ( verbose_outputDistCalc ) { cout << "event i,j = " << ientry << "," << jentry << endl;} 
			
			for ( int iVar=0; iVar<numVars; ++iVar ){
			   	phasePoint2[iVar] = varVector[iVar][jentry];
			}
			if (spectroscopicComboIDs[jentry] != spectroscopicComboIDs[ientry]){
			        distance = calc_distance(phasePoint1,phasePoint2);
			        //distance = rgen.Uniform(nentries);
			        distKNN.insertPair(make_pair(distance,jentry));
			}
			//if ( verbose) { 
			//	cout << "CURRENT SET: " << endl;
			//	for(auto elem : mapDistToJ){
			//		std::cout << elem.first << " " << elem.second << "\n";
			//	}
			//}
		}
		duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
		if(verbose2){logFile << "	Found neighbors: " << duration2 << "ms" << endl; }
		if (distKNN.kNN.size() != kDim){ cout << "size of distKNN is not equal to kDim! size,kDim="<< distKNN.kNN.size() << "," << kDim 
		    << "\n    -- if size is 1 less than kDim it is probably because kDim=nentries and event i cannot be a neighbor to itself" << 
		    "\n    -- if size != kDim it could also mean that the number of spectroscopically unique neighbors reduces the number of poential neighbors below kDim" << endl;}
		while ( distKNN.kNN.empty() == false ){
		        newPair = distKNN.kNN.top();
		        distKNN.kNN.pop();
			if(useEta){
		        	discriminatorHist->Fill(Metas[newPair.second],AccWeights[newPair.second]);//*sbWeights[newPair.second]);
		        	discriminatorHist2->Fill(Mpi0s[newPair.second],AccWeights[newPair.second]);//*sbWeights[newPair.second]);
		        	stdCalc.insertValue(Metas[newPair.second]);
		        	stdCalc2.insertValue(Mpi0s[newPair.second]);
			}
			else{
		        	discriminatorHist->Fill(Mpi0s[newPair.second],AccWeights[newPair.second]);//*sbWeights[newPair.second]);
		        	discriminatorHist2->Fill(Metas[newPair.second],AccWeights[newPair.second]);//*sbWeights[newPair.second]);
		        	stdCalc.insertValue(Mpi0s[newPair.second]);
		        	stdCalc2.insertValue(Metas[newPair.second]);
			}
		        //cout << "(" << newPair.first << ", " << newPair.second << ")"; 
		        //cout << endl; 
		}
		comboStd = stdCalc.calcStd();
		comboStd2 = stdCalc2.calcStd();
		
		//duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
		//if(verbose2){logFile << "	Filled neighbors: " << duration2 << "ms" << endl;}
		
		// Building the fit functions. We have to set some parameter limits to try to guide Minuit. We choose a double gaussian since for whatever reason the Meta has a asymmetry
		// such that there are more events to the left side of the peak. The second gaussian will have lower amplitude and a large sigma with a left shifted mean. We also want to 
		// fix the amplitudes and constnat background to be strictly positive and choose kDim as the max value ( which is reached only if all values are filled into a single bin)
		
		// We will always do 3 fits per event. In the original paper on Q-factors they use 100% bkg, %100 signal, 50% bkg and 50% signal.
		//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
		
		bestChiSq=DBL_MAX;
		int best_iFit=0;
		int countSig;
		int countBkg;

	        Double_t par[numDOFbkg+numDOFsig]; // needed to calculate the qvalue
	        Double_t parBest[numDOFbkg+numDOFsig]; // needed to calculate the qvalue
	        Double_t par2[numDOFbkg+numDOFsig]; // needed to calculate the qvalue
	        Double_t parFlat[numDOFbkg]; // needed to calculate the qvalue
		
		// /////////////////////////////////////////
		// Calcuclate q-value
		// /////////////////////////////////////////
		for (UInt_t iFit=0; iFit<3; ++iFit){
			if ( iFit != 0 ) { continue; }
			// We use a normalized gaussian and a flat function. 
			if(useEta){
				fit = new TF1("fit",fitFunc,fitRange1[0],fitRange1[1],numDOFbkg+numDOFsig);
				bkgFit = new TF1("bkgFit",background,fitRange1[0],fitRange1[1],numDOFbkg);
				sigFit = new TF1("sigFit",signalDG,fitRange1[0],fitRange1[1],numDOFsig);
			}
			else{
				fit = new TF1("fit",fitFunc,fitRange2[0],fitRange2[1],numDOFbkg+numDOFsig);
				bkgFit = new TF1("bkgFit",background,fitRange2[0],fitRange2[1],numDOFbkg);
				sigFit = new TF1("sigFit",signalDG,fitRange2[0],fitRange2[1],numDOFsig);
			}
			// Should use getInitParams.C whenever we get a new dataset to initialize the peak and width of the pi0 and eta
			fit->SetParameters(par0[iFit],0,par2[iFit],peakLoc,sigValue,ampRatio,widthRatio);
			fit->SetParLimits(0,0,kDim);
			fit->FixParameter(1,0); 
			fit->SetParLimits(2,0,kDim);
			fit->SetParLimits(3,peakLoc*0.95, peakLoc*1.05);
			fit->SetParLimits(4,sigValue*0.95, sigValue*1.05); 
			fit->SetParLimits(5,ampRatio*0.95,ampRatio*1.05 );
			fit->SetParLimits(6,widthRatio*0.95, widthRatio*1.05 );

			// we have to calculate the q-value for the eta distribution to check if it is between 0 and 1 
			discriminatorHist->Fit("fit","RQBNL"); // B will enforce the bounds, N will be no draw
			fit->GetParameters(par);
			bkgFit->SetParameters(par);
			sigFit->SetParameters(&par[numDOFbkg]);
			if(useEta){
				qvalue=sigFit->Eval(Metas[ientry])/fit->Eval(Metas[ientry]);
			}
			else{
				qvalue=sigFit->Eval(Mpi0s[ientry])/fit->Eval(Mpi0s[ientry]);
			}
			

			// /////////////////////////////////////////
			// need to calcuclate new q-value since it is out of bounds
			// /////////////////////////////////////////
			if (qvalue>1 || qvalue<0){
				cout << "Using flat fit instead of linear on event:\n" << ientry << endl;
				// first we will save the bad event then fix the linear component of the bkg
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
					etaLine = new TLine(Metas[ientry],0,Metas[ientry],kDim);
					etaLine->SetLineColor(kOrange);
	        			etaLine->Draw("same");
					allCanvases_badFit->SaveAs(("histograms/bad-Mass-event"+std::to_string(ientry)+".pdf").c_str());
					++savedN_badEvents;
				}
				
			        fit->SetParameters(par0[iFit],0,par2[iFit],peakLoc,sigValue,ampRatio,widthRatio);
				fit->FixParameter(1,0); 
				discriminatorHist->Fit("fit","RQBNL"); // B will enforce the bounds, N will be no draw
				fit->GetParameters(par);
				bkgFit->SetParameters(par);
				sigFit->SetParameters(&par[numDOFbkg]);

				qvalue=sigFit->Eval(Metas[ientry])/fit->Eval(Metas[ientry]);
				if (qvalue>1 || qvalue<0){
				    cout << "Not sure why qvalue is still >1 or <0. Need to fix this!" << endl;
				    cout << "These are the parameters:"<<endl;
				    for ( double parVal : par ){
						cout << " " << parVal << endl;
				    }
				    cout << " **************** BREAKING ****************** " << endl;
				    //exit(0);
				}
			}


			// now that the q-value is found we can get the chiSq and save the parameters with the best chiSq
			chiSq = fit->GetChisquare()/(fit->GetNDF());
			
			// save some fit and the chiSq for all the fits
			showConv[iFit]->SetParameters(par); // save all the fits
			chiSqs[iFit]=chiSq;
			
			if (verbose2) { logFile << "current ChiSq, best ChiSq: " << chiSq << ", " << bestChiSq << endl; }
			
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
		}


		// /////////////////////////////////////////
		// Calculating some chiSq differences 
		// /////////////////////////////////////////
		chiSq_eta_01 = chiSqs[1]-chiSqs[0];
		chiSq_eta_02 = chiSqs[2]-chiSqs[0];
		 
		if ( best_iFit == 0){ ++nBest100Bkg; }
		else if ( best_iFit == 1){ ++nBest100Sig; }
		else if ( best_iFit == 2){ ++nBest50Bkg50Sig; }
		
		// Here we draw the histograms that were randomly selected
		if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end()) {
			if(verbose) { cout << "Keeping this event" << endl; }
			discriminatorHist->SetTitle(("BEST VALS:  QValue="+to_string(best_qvalue)+"  ChiSq="+to_string(bestChiSq)+"     Std="+to_string(comboStd)+"    iFit="+to_string(best_iFit)).c_str() );
			etaLine = new TLine(Metas[ientry],0,Metas[ientry],kDim);
			pi0Line = new TLine(Mpi0s[ientry],0,Mpi0s[ientry],kDim);
			etaLine->SetLineColor(kOrange);
			pi0Line->SetLineColor(kOrange);
			
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

			if(useEta){
				qSigLine = new TLine(0,sigFit->Eval(Metas[ientry]),binRange1[2],sigFit->Eval(Metas[ientry]));
				qBkgLine = new TLine(0,bkgFit->Eval(Metas[ientry]),binRange1[2],bkgFit->Eval(Metas[ientry]));
			}
			else{
				qSigLine = new TLine(0,sigFit->Eval(Mpi0s[ientry]),binRange2[2],sigFit->Eval(Mpi0s[ientry]));
				qBkgLine = new TLine(0,bkgFit->Eval(Mpi0s[ientry]),binRange2[2],bkgFit->Eval(Mpi0s[ientry]));
			}
                  	qSigLine->SetLineColor(kBlue);
			qSigLine->SetLineStyle(9);
                  	qBkgLine->SetLineColor(kMagenta);
			qBkgLine->SetLineStyle(9);


                  	allCanvases->cd(1);
	          	discriminatorHist->Draw();
                  	drawText(parBest,numDOFbkg+numDOFsig,"par",sigFit->Eval(Metas[ientry]),bkgFit->Eval(Metas[ientry]),fit->Eval(Metas[ientry]));
			if(useEta){
	          		etaLine->Draw("same");
			}
                  	fit->Draw("SAME");
  	          	bkgFit->Draw("SAME FC");
  	          	sigFit->Draw("SAME FC");
			qBkgLine->Draw("SAME");
			qSigLine->Draw("SAME");
                  	allCanvases->cd(2);
	          	discriminatorHist2->Draw();
			if(!useEta){
	          		pi0Line->Draw("same");
			}

			allCanvases->SaveAs(("histograms/Mass-event"+std::to_string(ientry)+".pdf").c_str());

        		allCanvases->Clear();
        		allCanvases->Divide(2,1);
        		allCanvases->cd(1);
			discriminatorHist->Draw();
        		discriminatorHist->SetTitle("initialization");
        		allCanvases->cd(2);
			discriminatorHist->Draw();
        		discriminatorHist->SetTitle("converged");
        		for (int iFit=0; iFit<3; ++iFit){
        			allCanvases->cd(1);
        			showInit[iFit]->SetLineColor(colors[iFit]);
        			showInit[iFit]->Draw("SAME");
        			legend_init->AddEntry(showInit[iFit],initNames[iFit].c_str());
        			allCanvases->cd(2);
        			showConv[iFit]->SetLineColorAlpha(colors[iFit],0.2);
        			showConv[iFit]->Draw("SAME");
        			legend_conv->AddEntry(showConv[iFit],initNames[iFit].c_str());
        		}
        		allCanvases->cd(1);
        		legend_init->Draw();
        		allCanvases->cd(2);
			drawText(chiSqs,3,"chiSq",0,0,0);
        		legend_conv->Draw();
			//allCanvases->SaveAs(("histograms/fitCheck-event"+std::to_string(ientry)+".pdf").c_str());

			//if(verbose2){
			//     duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
			//     logFile << "	Best ChiSq (" << bestChiSq << ") with Q-value: " << qvalue << " -- time: " << duration2 <<  "ms" << endl;
			//     logFile << "	Delta T to finish event: " << duration2 <<  "ms" << endl;
			//}
		} 
		sbWeight = sbWeights[ientry];
		resultsTree->Fill();
	}

        resultsFile->cd();
        resultsTree->Write();
        //resultsFile->Write();
        if(verbose){cout << "nentries: " << nentries << endl;}
        if(verbose){cout << "Number of times 100% bkg initialization was the best: " << nBest100Bkg << endl;}
        if(verbose){cout << "Number of times 100% sig initialization was the best: " << nBest100Sig << endl;}
        if(verbose){cout << "Number of times 50/50 bkg/sig initialization was the best: " << nBest50Bkg50Sig << endl;}
        if(verbose){cout << "Number of times 100% bkg initialization for pi0 was the best: " << nBest100Bkg2 << endl;}
        if(verbose){cout << "Number of times 100% sig initialization for pi0 was the best: " << nBest100Sig2 << endl;}
        if(verbose){cout << "Number of times 50/50 bkg/sig initialization for pi0 was the best: " << nBest50Bkg50Sig2 << endl;}
        
        
	// Finish the log files by including an elapsed time and finally closing the file
	auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
	//logFile << "Total time: " << duration2 << " ms" << endl;
	//logFile << "Time Per Event: " << ( std::clock() - start ) / (double) CLOCKS_PER_SEC / nentries << "ns" << endl;
	logFile.close();
        return 0;
}

