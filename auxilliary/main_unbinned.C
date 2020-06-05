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

//void main(int iProcess, int kDim, int numberEventsToSavePerProcess, int nProcess, int seedShift, Long64_t nentries, bool override_nentries, bool verbose){
int main( int argc, char* argv[] ){
	// This suppresses all the "info" messages from roofit. Like saying that it will use numeric integrator for the generic pdf we define
	// https://root-forum.cern.ch/t/suppressing-info-messages/14642/6
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
        cout << "Num Vars="<< dim << endl;

    	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
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
        if(verbose){cout << "----------------------------" << endl;}
        if(verbose){cout << "iProcess: " << iProcess << endl;}
        if(verbose){cout << "kDim: " << kDim << endl;}
        if(verbose){cout << "numberEventsToSavePerProcess: " << numberEventsToSavePerProcess << endl;}
        if(verbose){cout << "nProcess: " << nProcess << endl;}
        if(verbose){cout << "seedShift: " << seedShift << endl;}
        if(verbose){cout << "nentries: " << nentries << endl; }
        if(verbose){cout << "override_nentries: " << override_nentries << endl;}
        if(verbose){cout << "verbose: " << verbose  << endl; }
    
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
	
	// setting up some basic root stuff and getting the file and tree
	//TFile* dataFile=new TFile("pi0eta_a0_recotreeFlat_DSelector.root");
	TFile* dataFile=new TFile("pi0eta_datatreeFlat_DSelector.root");
	TTree *dataTree;
	dataFile->GetObject("pi0eta_datatree_flat",dataTree);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
        //auto legend_init = new TLegend(0.1,0.8,0.3,0.9);
        //auto legend_conv = new TLegend(0.1,0.8,0.3,0.9);
        TLine* etaLine;
        TLine* pi0Line;
	//TH1F* discriminatorHist;
	//TH1F* discriminatorHist2;
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
        Int_t uniqueSpectroscopicPi0EtaID;

	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("Mpi0",&Mpi0);
	dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
        dataTree->SetBranchAddress("cosTheta_X_cm", &cosTheta_X_cm); 
        dataTree->SetBranchAddress("phi_X_cm",&phi_X_cm); 
        dataTree->SetBranchAddress("cosTheta_eta_gj",&cosTheta_eta_gj);
        dataTree->SetBranchAddress("phi_eta_gj",&phi_eta_gj); 
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
        dataTree->SetBranchAddress("uniqueSpectroscopicPi0EtaID",&uniqueSpectroscopicPi0EtaID);

	if (!override_nentries){
		nentries=dataTree->GetEntries();
	}
	cout << "Chosen Total Entries: " << nentries << endl;


	// opening a file to write my log data to
    	ofstream logFile;
    	logFile.open(("logs/logEventChiSqQValue_process"+to_string(iProcess)+".txt").c_str());
	//logFile << "Event\tQ-Value\tChiSq\tMpi0" << endl;
	

	const Long64_t c_nentries = (const Long64_t)nentries;

	// importing all the data to RAM instead of reading from root file
	std::vector<double> Metas; Metas.reserve(c_nentries);
        std::vector<double> Mpi0s; Mpi0s.reserve(c_nentries);
        std::vector<double> Mpi0etas; Mpi0etas.reserve(c_nentries);
        std::vector<double> cosTheta_X_cms; cosTheta_X_cms.reserve(c_nentries);
        std::vector<double> phi_X_cms; phi_X_cms.reserve(c_nentries);
        std::vector<double> cosTheta_eta_gjs; cosTheta_eta_gjs.reserve(c_nentries);
        std::vector<double> phi_eta_gjs; phi_eta_gjs.reserve(c_nentries);
        std::vector<double> cosThetaHighestEphotonIneta_gjs; cosThetaHighestEphotonIneta_gjs.reserve(c_nentries);
        std::vector<double> cosThetaHighestEphotonInpi0_cms; cosThetaHighestEphotonInpi0_cms.reserve(c_nentries);
        std::vector<double> pi0_energies; pi0_energies.reserve(c_nentries);
        std::vector<double> mandelstam_tps; mandelstam_tps.reserve(c_nentries);
	std::vector<double> vanHove_xs; vanHove_xs.reserve(c_nentries);
	std::vector<double> vanHove_ys; vanHove_ys.reserve(c_nentries);
	std::vector<double> vanHove_omegas; vanHove_omegas.reserve(c_nentries);
        std::vector<double> AccWeights; AccWeights.reserve(c_nentries);
        std::vector<Int_t> uniqueSpectroscopicPi0EtaIDs; uniqueSpectroscopicPi0EtaIDs.reserve(c_nentries);
        
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
		phi_X_cms.push_back(phi_X_cm);
		cosTheta_eta_gjs.push_back(cosTheta_eta_gj);
		phi_eta_gjs.push_back(phi_eta_gj);
		cosThetaHighestEphotonIneta_gjs.push_back(cosThetaHighestEphotonIneta_gj);	 
		cosThetaHighestEphotonInpi0_cms.push_back(cosThetaHighestEphotonInpi0_cm);	 
                vanHove_xs.push_back(vanHove_x);
                vanHove_ys.push_back(vanHove_y);
                vanHove_omegas.push_back(vanHove_omega);
                pi0_energies.push_back(pi0_energy);
                mandelstam_tps.push_back(mandelstam_tp);
                AccWeights.push_back(AccWeight);
                uniqueSpectroscopicPi0EtaIDs.push_back(uniqueSpectroscopicPi0EtaID);
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
        standardizeArray class_phi_X_cms(phi_X_cms,nentries);
        standardizeArray class_cosTheta_eta_gjs(cosTheta_eta_gjs,nentries);
        standardizeArray class_phi_eta_gjs(phi_eta_gjs,nentries);
        standardizeArray class_cosThetaHighestEphotonIneta_gjs(cosThetaHighestEphotonIneta_gjs,nentries);
        standardizeArray class_cosThetaHighestEphotonInpi0_cms(cosThetaHighestEphotonInpi0_cms,nentries);
        standardizeArray class_vanHove_xs(vanHove_xs,nentries);
        standardizeArray class_vanHove_ys(vanHove_ys,nentries);
        standardizeArray class_vanHove_omegas(vanHove_omegas,nentries);
        standardizeArray class_pi0_energies(pi0_energies, nentries);
        standardizeArray class_mandelstam_tps(mandelstam_tps, nentries);
        //standardizeArray class_Mpi0s(Mpi0s,nentries);
        //standardizeArray class_Metas(Metas,nentries);

        // ***** STANDARDIZES THE INTERAL COPY 
        class_cosTheta_X_cms.rangeStandardization();
        class_phi_X_cms.rangeStandardization();
        class_cosTheta_eta_gjs.rangeStandardization();
        class_phi_eta_gjs.rangeStandardization();
        class_cosThetaHighestEphotonIneta_gjs.rangeStandardization();
        class_cosThetaHighestEphotonInpi0_cms.rangeStandardization();
        class_vanHove_xs.rangeStandardization();
        class_vanHove_ys.rangeStandardization();
        class_vanHove_omegas.rangeStandardization();
        class_pi0_energies.rangeStandardization();
        class_mandelstam_tps.rangeStandardization();
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
        phi_X_cms=class_phi_X_cms.getVector();
        cosTheta_eta_gjs=class_cosTheta_eta_gjs.getVector();
        phi_eta_gjs=class_phi_eta_gjs.getVector();
        cosThetaHighestEphotonIneta_gjs=class_cosThetaHighestEphotonIneta_gjs.getVector();
        cosThetaHighestEphotonInpi0_cms=class_cosThetaHighestEphotonInpi0_cms.getVector();
        vanHove_xs=class_vanHove_xs.getVector();
        vanHove_ys=class_vanHove_xs.getVector();
        vanHove_omegas=class_vanHove_omegas.getVector();
        pi0_energies=class_pi0_energies.getVector();
        mandelstam_tps=class_mandelstam_tps.getVector();
        //Mpi0s=class_Mpi0s.getVector();
        //Metas=class_Metas.getVector();

        map<std::string, std::vector<double>> nameToVec;
        nameToVec["cosTheta_X_cms"] = cosTheta_X_cms;
        nameToVec["phi_X_cms"] = phi_X_cms; 
        nameToVec["cosTheta_eta_gjs"] = cosTheta_eta_gjs; 
        nameToVec["phi_eta_gjs"] = phi_eta_gjs; 
        nameToVec["cosThetaHighestEphotonIneta_gjs"] = cosThetaHighestEphotonIneta_gjs; 
        nameToVec["cosThetaHighestEphotonInpi0_cms"] = cosThetaHighestEphotonInpi0_cms; 
        nameToVec["vanHove_omegas"] = vanHove_omegas; 
        nameToVec["vanHove_xs"] = vanHove_xs; 
        nameToVec["vanHove_ys"] = vanHove_ys; 
        nameToVec["pi0_energies"] = pi0_energies; 
        nameToVec["mandelstam_tps"] = mandelstam_tps; 
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
        double qvalue_eta;
        double qvalue_pi0;
        double best_qvalue_eta=1;
        double best_qvalue_pi0=1;

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
	double allChiSq_eta[3];
	double chiSq_eta, chiSq_pi0;
	double chiSq_eta_01, chiSq_eta_02;
        double comboStd2; 
        ULong64_t flatEntryNumber;
        double bestChiSq_eta=1;
        double bestChiSq_pi0=1;
	bool useFlat;

        TFile *resultsFile = new TFile(("logs/results"+to_string(iProcess)+".root").c_str(),"RECREATE");
        TTree* resultsTree = new TTree("resultsTree","results");
        resultsTree->Branch("flatEntryNumber",&flatEntryNumber,"flatEntryNumber/l");
        resultsTree->Branch("qvalue",&best_qvalue_eta,"qvalue/D");
        resultsTree->Branch("chisq_eta",&bestChiSq_eta,"chisq_eta/D");
        resultsTree->Branch("chisq_eta_01",&chiSq_eta_01,"chisq_eta/D");
        resultsTree->Branch("chisq_eta_02",&chiSq_eta_02,"chisq_eta/D");
        resultsTree->Branch("combostd",&comboStd,"combostd/D");
        resultsTree->Branch("chisq_pi0",&bestChiSq_pi0,"chisq_pi0/D");
        resultsTree->Branch("combostd2",&comboStd2,"combostd2/D");
        resultsTree->Branch("useFlat",&useFlat,"useFlat/B");
        if(verbose) {cout << "Set up branch addresses" << endl;}

        std::vector<double> binRange;
        std::vector<double> fitRange;
        std::vector<double> binRange2;
        std::vector<double> fitRange2;
	TF1 *fit,*fit2;
        TF1 *bkgFit,*bkgFit2;
        TF1 *sigFit,*sigFit2;
        binRange={50,0.35,0.8};
        fitRange={0.45,0.7};
        //fitRange={0.35,0.8};
        binRange2={50,0.05,0.25};
        //fitRange2={0.05,0.25};
        fitRange2={0.1,0.17};

        // Going to keep a track of the fits to see how they move around
        TF1 *showInit[3];
        TF1 *showConv[3];
        int colors[3] = {kRed, kBlue, kGreen};
        string initNames[3] = {"100% Bkg", "50/50% Bkg Sig", "100% Sig"};
        for (int iFit=0; iFit<3; ++iFit){
	    showInit[iFit] = new TF1(("initFit"+to_string(iFit)).c_str(),fitFunc,fitRange[0],fitRange[1],numDOFbkg+numDOFsig);
	    showConv[iFit] = new TF1(("convFit"+to_string(iFit)).c_str(),fitFunc,fitRange[0],fitRange[1],numDOFbkg+numDOFsig);
        }
        if(verbose) {cout << "Made array of fits" << endl;}

        int nBest100Bkg=0;
        int nBest100Sig=0;
        int nBest50Bkg50Sig=0;
        int nBest100Bkg2=0;
        int nBest100Sig2=0;
        int nBest50Bkg50Sig2=0;
        // Getting the scaling of the flat bkg is a bit tricky. Have to count how much bins we have in our fit range. kDim divided by the bins in fit range is the height of the flat function.
        double numBinsInFit = (fitRange[1]-fitRange[0])/((binRange[2]-binRange[1])/binRange[0]);
        double binSize=((binRange[2]-binRange[1])/binRange[0]);
        double numBinsInFit2 = (fitRange2[1]-fitRange2[0])/((binRange2[2]-binRange2[1])/binRange2[0]);
        double binSize2=((binRange2[2]-binRange2[1])/binRange2[0]);

        if (verbose) {cout << "Eta hist range: " << binRange[0] << ", " << binRange[1] << ", " << binRange[2] << endl;}
        if (verbose) {cout << "Eta fit range: " << fitRange[0] << ", " << fitRange[1] << endl;}
        if (verbose) {cout << "Pi0 hist range: " << binRange2[0] << ", " << binRange2[1] << ", " << binRange2[2] << endl;}
        if (verbose) {cout << "Pi0 fit range: " << fitRange2[0] << ", " << fitRange2[1] << endl;}
        
        // These initializations with 50 bins are related to the fitted gaussian on all the data
        double peakWidtheta[2] = {0.545928, 0.0196892};
        double peakWidthpi0[2] = {0.135399, 0.00760648};
        double par0eta[3] = {3.4528, 1.7264, 0};
        double par1eta[3] = {5.49805, 2.74902, 0}; 
        double par2eta[3] = {0, 0.9, 1.18}; 
        double par0pi0[3] = {7.30661, 3.6533, 0}; 
        double par1pi0[3] = {30.05331, 15.2665, 0}; 
        double par2pi0[3] = {0, 0.400002, 0.800003}; 

        showInit[0]->SetParameters(par0eta[0],par1eta[0],par2eta[0],peakWidtheta[0],peakWidtheta[1]);
        showInit[1]->SetParameters(par0eta[1],par1eta[1],par2eta[1],peakWidtheta[0],peakWidtheta[1]);
        showInit[2]->SetParameters(par0eta[2],par1eta[2],par2eta[2],peakWidtheta[0],peakWidtheta[1]);

	double sigfrac, sigPdfAtMass_i_eta, bkgPdfAtMass_i_eta, sigPdfAtMass_i_pi0, bkgPdfAtMass_i_pi0;
	int numDOF;

        // phasePoint1 will consider all events from lowest to largest since these will be our attached q values. phasePoint2 on the other hand will only look at a subset of the events where the
        // elements must be spectroscopically distinct, i.e. the 4 photons in consideration are different. Not sure if this is the value I should consider or is it better to do a pair of maps
        // where I am tracking the two photon pairs that make up the eta and pi0.         
        set<Int_t> setUsedSpectroscopicIDs;
        std::vector<int> phasePoint2PotentailNeighbor; phasePoint2PotentailNeighbor.reserve(nentries);
        for (Int_t ientry=0; ientry<nentries; ientry++){ 
            if ( setUsedSpectroscopicIDs.find( uniqueSpectroscopicPi0EtaIDs[ientry] ) == setUsedSpectroscopicIDs.end() ) {
                setUsedSpectroscopicIDs.insert( uniqueSpectroscopicPi0EtaIDs[ientry] );
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
        for (int ientry=lowest_nentry; ientry<largest_nentry; ientry++){ 
		etaLine = new TLine(Metas[ientry],0,Metas[ientry],kDim);
		pi0Line = new TLine(Mpi0s[ientry],0,Mpi0s[ientry],kDim);
		etaLine->SetLineColor(kOrange);
		pi0Line->SetLineColor(kOrange);
		RooRealVar x_eta("x_eta","Mass GeV",0.35 , 0.8);
		x_eta.setRange("fitRange",0.4,0.7);
		RooRealVar x_pi0("x_pi0","Mass GeV",0.05, 0.25);
		RooDataSet data_eta("data_eta","data_eta",RooArgSet(x_eta));
		RooDataSet data_pi0("data_pi0","data_pi0",RooArgSet(x_pi0));

		flatEntryNumber=ientry;
		cumulativeStd stdCalc(kDim);
		cumulativeStd stdCalc2(kDim);
		
		auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
		auto duration_beginEvent = std::chrono::high_resolution_clock::now();
		if(verbose2){cout << "Starting event " << ientry << "/" << largest_nentry << " ---- Time: " << duration2 << "ms" << endl; }
		
		// clear data from previous events
		mapDistToJ.clear();
		distances.clear();
		
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
			if (uniqueSpectroscopicPi0EtaIDs[jentry] != uniqueSpectroscopicPi0EtaIDs[ientry]){
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
		if(verbose2){cout << "	Found neighbors: " << duration2 << "ms" << endl; }
		if (distKNN.kNN.size() != kDim){ cout << "size of distKNN is not equal to kDim! size,kDim="<< distKNN.kNN.size() << "," << kDim 
		    << "\n    -- if size is 1 less than kDim it is probably because kDim=nentries and event i cannot be a neighbor to itself" << 
		    "\n    -- if size != kDim it could also mean that the number of spectroscopically unique neighbors reduces the number of poential neighbors below kDim" << endl;}
		while ( distKNN.kNN.empty() == false ){
		        newPair = distKNN.kNN.top();
		        distKNN.kNN.pop();
		        //discriminatorHist->Fill(Metas[newPair.second]);//,AccWeights[newPair.second]);
		        //discriminatorHist2->Fill(Mpi0s[newPair.second]);//,AccWeights[newPair.second]);
		        stdCalc.insertValue(Metas[newPair.second]);
		        stdCalc2.insertValue(Mpi0s[newPair.second]);
			x_eta = Metas[newPair.second];
			x_pi0 = Mpi0s[newPair.second];
			data_eta.add(RooArgSet(x_eta));
			data_pi0.add(RooArgSet(x_pi0));
		        //else {
		        //    discriminatorHist->Fill(Mpi0s[newPair.second]);//,AccWeights[newPair.second]);
		        //    stdCalc.insertValue(Mpi0s[newPair.second]);
		        //}
		        //cout << "(" << newPair.first << ", " << newPair.second << ")"; 
		        //cout << endl; 
		}
		if (verbose2){ cout << "Loaded neighbors into dataset" << endl; }
		comboStd = stdCalc.calcStd();
		comboStd2 = stdCalc2.calcStd();
		
		//duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
		//if(verbose2){logFile << "	Filled neighbors: " << duration2 << "ms" << endl;}
		
		// Building the fit functions. We have to set some parameter limits to try to guide Minuit. We choose a double gaussian since for whatever reason the Meta has a asymmetry
		// such that there are more events to the left side of the peak. The second gaussian will have lower amplitude and a large sigma with a left shifted mean. We also want to 
		// fix the amplitudes and constnat background to be strictly positive and choose kDim as the max value ( which is reached only if all values are filled into a single bin)
		
		
		// We will always do 3 fits per event. In the original paper on Q-factors they use 100% bkg, %100 signal, 50% bkg and 50% signal.
		//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
		
		bestChiSq_eta=DBL_MAX;
		int best_iFit=0;
		bestChiSq_pi0=DBL_MAX;
		int best_iFit2=0;
		int countSig;
		int countBkg;
		
		// So even though we assign pointers to this object it doesnt act like we are assigning to the same underlying object if we call frame twice. Can 
		// check this by loading the varaible and plotOn one of the frames to see that numItems() only increments one of them
		RooPlot * xframe_eta = x_eta.frame();
		//RooPlot * xframe_pi0 = x_pi0.frame();
		RooPlot * xframe_best_eta = x_eta.frame();
		RooPlot * xframe_best_pi0 = x_pi0.frame();

		if (verbose2) { cout << "made roo frames" << endl; }

		for (UInt_t iFit=0; iFit<3; ++iFit){
			// ==================
			// ****** ETA *******
			// ==================
			useFlat=false;
			RooRealVar peak_eta("peak_eta","peak_eta",0.545928,0.53,0.57);
			RooConstVar width_eta("width_eta","width_eta", 0.017712);//0.0196892);
			RooRealVar par0_eta("par0_eta","par0_eta", par0eta[iFit], 0, 5.25);
			RooRealVar par1_eta("par1_eta","par1_eta", par1eta[iFit], 0, 8.25);
			RooRealVar par2_eta("par2_eta","par2_eta", par2eta[iFit], 0, 1.77);
			RooGenericPdf rooBkgFit_eta("rooBkgFit_eta","rooBkgFit_eta","par0_eta+par1_eta*x_eta",RooArgSet(x_eta,par0_eta,par1_eta));
			RooGenericPdf rooSigFit_eta("rooSigFit_eta","rooSigFit_eta","par2_eta/sqrt(2*pi)/width_eta*exp(-0.5*((x_eta-peak_eta)/width_eta)**2)", RooArgSet(x_eta,par0_eta,par1_eta,par2_eta,peak_eta,width_eta));
  			RooRealVar signalfraction_eta("signalFraction_eta","signalFraction_eta",0.36,0,1); // Sung had fixed signal fraction to 0.36 for his method.
			RooAddPdf rooSigPlusBkg_eta("sum_eta","g+a",RooArgList(rooSigFit_eta,rooBkgFit_eta),signalfraction_eta) ;

			// Calculating the qvalue
			RooFitResult* r = rooSigPlusBkg_eta.fitTo(data_eta, Range("fitRange"), PrintLevel(-1),Save());
			//r->Print();
  			RooArgSet* value =  rooSigPlusBkg_eta.getObservables(x_eta); // if we cout << value  we get (x). 
    			value->setRealValue(x_eta.GetName(),Metas[ientry]);
    			sigfrac = signalfraction_eta.getVal();
    			sigPdfAtMass_i_eta = rooSigFit_eta.getVal(value);
    			bkgPdfAtMass_i_eta = rooBkgFit_eta.getVal(value);
    			bkgPdfAtMass_i_eta = (1.- sigfrac) * bkgPdfAtMass_i_eta;
    			sigPdfAtMass_i_eta =  sigfrac * sigPdfAtMass_i_eta;
    			qvalue_eta = sigPdfAtMass_i_eta /(sigPdfAtMass_i_eta+bkgPdfAtMass_i_eta);

			if (qvalue_eta>1 || qvalue_eta<0){
				useFlat=true;
				RooRealVar peak_eta("peak_eta","peak_eta",0.545928,0.53,0.57);
				RooConstVar width_eta("width_eta","width_eta",0.0196892);
				RooRealVar par0_eta("par0_eta","par0_eta", par0eta[iFit], 0, 5.25);
				RooRealVar par2_eta("par2_eta","par2_eta", par2eta[iFit], 0, 1.77);
				RooGenericPdf rooBkgFit_eta("rooBkgFit_eta","rooBkgFit_eta","par0_eta",RooArgSet(x_eta,par0_eta));
				RooGenericPdf rooSigFit_eta("rooSigFit_eta","rooSigFit_eta","par2_eta/sqrt(2*pi)/width_eta*exp(-0.5*((x_eta-peak_eta)/width_eta)**2)", RooArgSet(x_eta,par0_eta,par1_eta,par2_eta,peak_eta,width_eta));
  				RooRealVar signalfraction_eta("signalFraction_eta","signalFraction_eta",0.36,0,1); // Sung had fixed signal fraction to 0.36 for his method.
				RooAddPdf rooSigPlusBkg_eta("sum_eta","g+a",RooArgList(rooSigFit_eta,rooBkgFit_eta),signalfraction_eta) ;

				// Calculating the qvalue
				RooFitResult* r = rooSigPlusBkg_eta.fitTo(data_eta,PrintLevel(-1),Save());
				//r->Print();
  				RooArgSet* value =  rooSigPlusBkg_eta.getObservables(x_eta); // if we cout << value  we get (x). 
    				value->setRealValue(x_eta.GetName(),Metas[ientry]);
    				sigfrac = signalfraction_eta.getVal();
    				sigPdfAtMass_i_eta = rooSigFit_eta.getVal(value);
    				bkgPdfAtMass_i_eta = rooBkgFit_eta.getVal(value);
    				bkgPdfAtMass_i_eta = (1.- sigfrac) * bkgPdfAtMass_i_eta;
    				sigPdfAtMass_i_eta =  sigfrac * sigPdfAtMass_i_eta;
    				qvalue_eta = sigPdfAtMass_i_eta /(sigPdfAtMass_i_eta+bkgPdfAtMass_i_eta);
			}

			if (qvalue_eta>1 || qvalue_eta<0){
				cout << "STILL QVALUE IS OUT OF BOUNDS! FIX ME!" << endl;
				exit(0);
			}

			if (verbose2) { cout << "Fitted the eta mass distribution - got qvalue" << endl; }

			RooAbsCollection* flparams = rooSigPlusBkg_eta.getParameters(value)->selectByAttrib("Constant",kFALSE) ;
			numDOF = flparams->getSize();
			data_eta.plotOn(xframe_eta); 
			rooSigPlusBkg_eta.plotOn(xframe_eta) ;
			//rooSigPlusBkg_eta.paramOn(xframe_eta,Layout(0.7,0.9,0.9));
			//xframe_eta->getAttText()->SetTextSize(0.02) ;
			//rooSigPlusBkg_eta.plotOn(xframe_eta,Components(rooBkgFit_eta),LineStyle(kDashed),LineColor(kMagenta)) ;
			//rooSigPlusBkg_eta.plotOn(xframe_eta,Components(rooSigFit_eta),LineStyle(kDashed),LineColor(kGreen)) ;
			xframe_eta->Draw();
			// need to plot the data and the fit first before calculating the chiSq
			chiSq_eta = xframe_eta->chiSquare(flparams->getSize()) ;
			allChiSq_eta[iFit] = chiSq_eta;
			if (verbose2) { cout << "plotOn data and fit to get the chiSq" << endl; }
			// We should check if the q-value is in bounds. If it is we can skip the next loop (loop with flat bkg).
					
			// ==================
			// ****** PI0 ******
			// ==================
	//		RooRealVar peak_pi0("peak_pi0","peak_pi0",0.135399,0.125,0.15);
	//		RooConstVar width_pi0("width_pi0","width_pi0",0.00760648);
	//		//RooRealVar width_pi0("width_pi0","width_pi0",0.00760648,0.005,0.025);
	//		RooRealVar par0_pi0("par0_pi0","par0_pi0", par0pi0[iFit], 0, 10.95);
	//		RooRealVar par2_pi0("par2_pi0","par2_pi0", par2pi0[iFit], 0, 1.2);
	//		RooGenericPdf rooBkgFit_pi0("rooBkgFit_pi0","rooBkgFit_pi0","par0_pi0",RooArgSet(x_pi0,par0_pi0));
	//		RooGenericPdf rooSigFit_pi0("rooSigFit_pi0","rooSigFit_pi0","par2_pi0/sqrt(2*pi)/width_pi0*exp(-0.5*((x_pi0-peak_pi0)/width_pi0)**2)", RooArgSet(x_pi0,par0_pi0,par2_pi0,peak_pi0,width_pi0));
  	//		RooRealVar signalfraction_pi0("signalFraction_pi0","signalFraction_pi0",0.36,0,1); // Sung had fixed signal fraction to 0.36 for his method.
	//		RooAddPdf rooSigPlusBkg_pi0("sum_pi0","g+a",RooArgList(rooSigFit_pi0,rooBkgFit_pi0),signalfraction_pi0) ;

	//		r = rooSigPlusBkg_pi0.fitTo(data_pi0,PrintLevel(-1),Save());
	//		//r->Print();
  	//		value =  rooSigPlusBkg_pi0.getObservables(x_pi0); // if we cout << value  we get (x). 
    	//		value->setRealValue(x_pi0.GetName(),Mpi0s[ientry]);
    	//		sigfrac = signalfraction_pi0.getVal();
    	//		sigPdfAtMass_i_pi0 = rooSigFit_pi0.getVal(value);
    	//		bkgPdfAtMass_i_pi0 = rooBkgFit_pi0.getVal(value);
    	//		bkgPdfAtMass_i_pi0 = (1.- sigfrac) * bkgPdfAtMass_i_pi0;
    	//		sigPdfAtMass_i_pi0 =  sigfrac * sigPdfAtMass_i_pi0;
    	//		qvalue_pi0 = sigPdfAtMass_i_pi0 /(sigPdfAtMass_i_pi0+bkgPdfAtMass_i_pi0);
	//		flparams = rooSigPlusBkg_pi0.getParameters(value)->selectByAttrib("Constant",kFALSE) ;
	//		if( iFit==0) { data_pi0.plotOn(xframe_pi0); }
	//		//rooSigPlusBkg_pi0.paramOn(xframe_pi0,Layout(0.7,0.9,0.9));
	//		//xframe_pi0->getAttText()->SetTextSize(0.02) ;
	//		//rooSigPlusBkg_pi0.plotOn(xframe_pi0,Components(rooBkgFit_pi0),LineStyle(kDashed),LineColor(kMagenta)) ;
	//		//rooSigPlusBkg_pi0.plotOn(xframe_pi0,Components(rooSigFit_pi0),LineStyle(kDashed),LineColor(kGreen)) ;
	//		rooSigPlusBkg_pi0.plotOn(xframe_pi0) ;
	//		chiSq_pi0 = xframe_pi0->chiSquare(flparams->getSize()) ;

	//		if (qvalue_pi0>1 || qvalue_pi0<0){
	//			cout << "STILL QVALUE IS OUT OF BOUNDS! FIX ME!" << endl;
	//			exit(0);
	//		}


			// check out this code. It shows that when we use plotOn it seems to take a snapshot of the current pdf and its pars
			// this allows us to update the xframe as the pdf progresses. Useful in the case of cases where we are trying to find the best
			// /d/grid15/ln16/pi0eta/q-values/testCode/testRooFit.C
			if (chiSq_eta < bestChiSq_eta){
			        xframe_best_eta  = x_eta.frame();
				bestChiSq_eta = chiSq_eta;
			        best_qvalue_eta = qvalue_eta;
			        best_iFit = iFit;
			        data_eta.plotOn(xframe_best_eta);
			        rooSigPlusBkg_eta.paramOn(xframe_best_eta,Layout(0.7,0.9,0.9));
			        xframe_best_eta->getAttText()->SetTextSize(0.02) ;
			        rooSigPlusBkg_eta.plotOn(xframe_best_eta,Components(rooBkgFit_eta),LineStyle(kDashed),LineColor(kMagenta)) ;
			        rooSigPlusBkg_eta.plotOn(xframe_best_eta,Components(rooSigFit_eta),LineStyle(kDashed),LineColor(kGreen)) ;
			        rooSigPlusBkg_eta.plotOn(xframe_best_eta) ;
			}
			if (verbose2) { cout << "Saving the current best fit" << endl; }

	//		if (chiSq_pi0 < bestChiSq_pi0){
				//xframe_best_pi0  = x_pi0.frame();
	//			best_qvalue_pi0 = qvalue_pi0;
	//			best_iFit2 = iFit;
				//data_pi0.plotOn(xframe_best_pi0);
	//			rooSigPlusBkg_pi0.paramOn(xframe_best_pi0,Layout(0.7,0.9,0.9));
	//			xframe_best_pi0->getAttText()->SetTextSize(0.02) ;
	//			rooSigPlusBkg_pi0.plotOn(xframe_best_pi0,Components(rooBkgFit_pi0),LineStyle(kDashed),LineColor(kMagenta)) ;
	//			rooSigPlusBkg_pi0.plotOn(xframe_best_pi0,Components(rooSigFit_pi0),LineStyle(kDashed),LineColor(kGreen)) ;
	//			rooSigPlusBkg_pi0.plotOn(xframe_best_pi0) ;
	//		}
		}

		if ( best_iFit == 0){ ++nBest100Bkg; }
		else if ( best_iFit == 1){ ++nBest100Sig; }
		else if ( best_iFit == 2){ ++nBest50Bkg50Sig; }
		if ( best_iFit2 == 0){ ++nBest100Bkg2; }
		else if ( best_iFit2 == 1){ ++nBest100Sig2; }
		else if ( best_iFit2 == 2){ ++nBest50Bkg50Sig2; }

		chiSq_eta_01 = allChiSq_eta[0]-allChiSq_eta[1];
		chiSq_eta_02 = allChiSq_eta[0]-allChiSq_eta[2];
		
		// Here we draw the histograms that were randomly selected
		if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end()) {
			allCanvases->Clear();
			allCanvases->Divide(2,1);

			allCanvases->cd(1);
			xframe_best_eta->SetTitle(("best_qvalue_eta: "+std::to_string(best_qvalue_eta)).c_str());
			xframe_best_eta->Draw();
			etaLine->Draw("same");
			cout << "===== ETA =====" << endl;
			cout << "NDOF: " << numDOF << endl;
			cout << "chiSq_eta: " <<  bestChiSq_eta << endl ;
			cout << "Mass at event i:" << Metas[ientry] << endl;
			cout << "qvalue_eta: " << best_qvalue_eta << endl;

			allCanvases->cd(2);

			xframe_best_pi0->SetTitle(("best qvalue_pi0: "+std::to_string(best_qvalue_pi0)).c_str());
			data_pi0.plotOn(xframe_best_pi0);
			xframe_best_pi0->Draw();
			pi0Line->Draw("same");
			cout << "===== PI0 =====" << endl;
			//cout << "NDOF: " << numDOF << endl;
			//cout << "chiSq_pi0: " <<  chiSq_pi0 << endl ;
			cout << "Mass at event i:" << Metas[ientry] << endl;
			//cout << "qvalue_pi0: " << qvalue_pi0 << endl;
			
			allCanvases->SaveAs(("histograms/Mass-event"+std::to_string(ientry)+".png").c_str());

			//if(verbose2){
			//     duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
			//     logFile << "	Best ChiSq (" << bestChiSq << ") with Q-value: " << qvalueEta << " -- time: " << duration2 <<  "ms" << endl;
			//     logFile << "	Delta T to finish event: " << duration2 <<  "ms" << endl;
			//}
			
		} 
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
	// duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
	//logFile << "Total time: " << duration2 << " ms" << endl;
	//logFile << "Time Per Event: " << ( std::clock() - start ) / (double) CLOCKS_PER_SEC / nentries << "ns" << endl;
	logFile.close();
        return 0;
}

