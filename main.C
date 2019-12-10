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
	TFile* dataFile=new TFile("pi0eta_fcal_treeFlat_DSelector.root");
	TTree *dataTree;
	dataFile->GetObject("pi0eta_datatree_flat",dataTree);
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
        Int_t uniqueSpectroscopicPi0EtaID;

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
        //std::vector<double> phi_X_cms; phi_X_cms.reserve(c_nentries);
        std::vector<double> cosTheta_eta_gjs; cosTheta_eta_gjs.reserve(c_nentries);
        std::vector<double> phi_eta_gjs; phi_eta_gjs.reserve(c_nentries);
        //std::vector<double> cosThetaHighestEphotonIneta_gjs; cosThetaHighestEphotonIneta_gjs.reserve(c_nentries);
        //std::vector<double> cosThetaHighestEphotonInpi0_cms; cosThetaHighestEphotonInpi0_cms.reserve(c_nentries);
        //std::vector<double> pi0_energies; pi0_energies.reserve(c_nentries);
        //std::vector<double> mandelstam_tps; mandelstam_tps.reserve(c_nentries);
	//std::vector<double> vanHove_xs; vanHove_xs.reserve(c_nentries);
	//std::vector<double> vanHove_ys; vanHove_ys.reserve(c_nentries);
	//std::vector<double> vanHove_omegas; vanHove_omegas.reserve(c_nentries);
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
        double qvalue_eta;
        double qvalue_pi0;
        double best_qvalue_eta;
        double best_qvalue_pi0;

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
        double chiSq_eta;
        double chiSq_pi0;
        double comboStd2; 
        ULong64_t flatEntryNumber;
        double bestChiSq_eta;
        double bestChiSq_pi0;
	double chiSq_eta_01;
	double chiSq_eta_02;

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
        if(verbose) {cout << "Set up branch addresses" << endl;}

        std::vector<double> binRange;
        std::vector<double> fitRange;
        std::vector<double> binRange2;
        std::vector<double> fitRange2;
	TF1 *fit,*fit2;
        TF1 *bkgFit,*bkgFit2;
        TF1 *sigFit,*sigFit2;
        binRange={50,0.35,0.8};
        fitRange={0.42,0.68};
        //fitRange={0.35,0.8};
        binRange2={50,0.05,0.25};
        //fitRange2={0.05,0.25};
        fitRange2={0.1,0.17};

        // Going to keep a track of the fits to see how they move around
        TF1 *showInit[3];
        TF1 *showConv[3];
        int colors[3] = {kRed, kBlue, kGreen};
	double chiSqs[3];
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
        //double peakWidtheta[2] = {0.545928, 0.0196892};
        //double peakWidthpi0[2] = {0.135399, 0.00760648};
        //double par0eta[3] = {3.4528, 1.7264, 0};
        //double par1eta[3] = {5.49805, 2.74902, 0}; 
        //double par2eta[3] = {0, 0.9, 1.18}; 
        //double par0pi0[3] = {7.30661, 3.6533, 0}; 
        //double par1pi0[3] = {30.05331, 15.2665, 0}; 
        //double par2pi0[3] = {0, 0.400002, 0.800003}; 

	//double peakWidtheta[2] = { 0.546184, 0.0221398};
	//double par0eta[3] = { 5.38593, 2.69296, 0};
	//double par1eta[3] = { 1.11649, 0.558246, 0};
	//double par2eta[3] = { 0, 0.9, 1.8};
	//
	// BW for the eta
	// double peakWidtheta[2] =  {0.548531, 0.0466235};
	// double par0eta[3] = { -0.449614, -0.224807, 0};
	// double par1eta[3] = { 1.90839, 0.954194, 0};
	// double par2eta[3] = { 0, 0.0901173, 0.180235};
	// // These are for a guassian pi0
	// double peakWidthpi0[2] = { 0.135198, 0.00789097};
	// double par0pi0[3] = { 6.69442, 3.34721, 0};
	// double par1pi0[3] = { 35.0678, 17.5339, 0};
	// double par2pi0[3] = { 0, 0.400004, 0.800007};

	// ********************** THIS IS FOR THE MEAS VALUES *********************
		// kDim = 200
	// old
	//double par0eta[3] = { 0.0110229, 0.00551144, 0};
	//double par1eta[3] = { 0.0345038, 0.0172519, 0};
	//double par2eta[3] = { 0, 0.202743, 0.405486};
	// Newly calculated 11/14/19 5000 Bins
	//double par0eta[3] = { 0.0288, 0.01438, 0};
	//double par1eta[3] = { 0.09, 0.045, 0};
	//double par2eta[3] = { 0, 0.005722, 0.01144};
	// Newly calculated 11/16/19 5000 Bins
	//double peakWidtheta[2] = {0.540358, 0.0234706};
	//double peakWidthpi0[2] = { 0.134273, 0.00781488};
	//double par0eta[3] = { 3, 1.5, 0};
	//double par1eta[3] = { 9, 4.5, 0 };
	//double par2eta[3] = { 0, 0.574, 1.14 };
	double par0pi0[3] = { 8.58665, 4.29333, 0};
	double par1pi0[3] = { 21.0513, 10.5256, 0};
	double par2pi0[3] = { 0, 54.8886, 109.777};

	// all - Newly calculated 11/21/19 for double gaussian
	//double peakWidtheta[2] = {0.54625, 0.00994868};
	//double ampRatio = 463.361/131.364;
	//double widthRatio = 0.023968/0.00994868;
	//double par0eta[3] = { 0.71, 0.3547 , 0};
	//double par1eta[3] = { 0.946, 0.472, 0 };
	//double par2eta[3] = { 0, 0.022, 0.04414 };
	// (WORKS OK) bcal - calculated 12/9/18 - 400 neighbors 
	//double peakWidtheta[2] = {0.545679, 0.0242801};
	//double ampRatio = 1.0715;
	//double widthRatio = 0.48342;
	//double par0eta[3] = { 5.9, 2.95 , 0};
	//double par1eta[3] = { -2.92 , -1.462, 0 };
	//double par2eta[3] = { 0, 0.6, 1.2 };
	//
	// () fcal - calculated 12/9/18 - 400 neighbors
	double peakWidtheta[2] = {0.54655, 0.0262712};
	double ampRatio = 0.474188;
	double widthRatio = 0.398465;
	//double par0eta[3] = { 5.18 , 2.59, 0 };
	//double par1eta[3] = { 1.017 , 0.5085, 0 };
	//double par2eta[3] = { 0, 0.65, 1.3 };
	double par0eta[3] = { 5.18 , 0, 0 };
	double par1eta[3] = { 1.017 , 0, 0 };
	double par2eta[3] = { 1.3, 0, 0 };

        // with 5000 bins
        //double peakWidtheta[2] = {0.545949, 0.0194813};
        //double peakWidthpi0[2] = {0.13541, 0.00747303};
        //double par0eta[3] = {0.0348218, 0.0174109, 0};
        //double par1eta[3] = {0.0544581, 0.0272291, 0}; 
        //double par2eta[3] = {0, 0.009, 0.018}; 
        //double par0pi0[3] = {0.0785842, 0.0392921, 0}; 
        //double par1pi0[3] = {0.264456, 0.132228, 0}; 
        //double par2pi0[3] = {0, 0.00400001, 0.00800002}; 

        showInit[0]->SetParameters(par0eta[0],par1eta[0],par2eta[0],peakWidtheta[0],peakWidtheta[1], ampRatio, widthRatio);
        showInit[1]->SetParameters(par0eta[1],par1eta[1],par2eta[1],peakWidtheta[0],peakWidtheta[1], ampRatio, widthRatio);
        showInit[2]->SetParameters(par0eta[2],par1eta[2],par2eta[2],peakWidtheta[0],peakWidtheta[1], ampRatio, widthRatio);

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
	int saveN_badEvents=1;
	int savedN_badEvents=0;
        for (int ientry=lowest_nentry; ientry<largest_nentry; ientry++){ 
		legend_init->Clear();
		legend_conv->Clear();
		flatEntryNumber=ientry;
		cumulativeStd stdCalc(kDim);
		cumulativeStd stdCalc2(kDim);
		
		if(verbose) { cout << "Getting next event!\n--------------------------------\n" << endl;  }
		auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
		auto duration_beginEvent = std::chrono::high_resolution_clock::now();
		if(verbose2){logFile << "Starting event " << ientry << "/" << largest_nentry << " ---- Time: " << duration2 << "ms" << endl; }
		
		// clear data from previous events
		mapDistToJ.clear();
		distances.clear();
		allCanvases->Clear();
		allCanvases_badFit->Clear();
		allCanvases->Divide(2,1);
		discriminatorHist = new TH1F("","",binRange[0],binRange[1],binRange[2]);
		discriminatorHist2 = new TH1F("","",binRange2[0],binRange2[1],binRange2[2]);
		
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
		if(verbose2){logFile << "	Found neighbors: " << duration2 << "ms" << endl; }
		if (distKNN.kNN.size() != kDim){ cout << "size of distKNN is not equal to kDim! size,kDim="<< distKNN.kNN.size() << "," << kDim 
		    << "\n    -- if size is 1 less than kDim it is probably because kDim=nentries and event i cannot be a neighbor to itself" << 
		    "\n    -- if size != kDim it could also mean that the number of spectroscopically unique neighbors reduces the number of poential neighbors below kDim" << endl;}
		while ( distKNN.kNN.empty() == false ){
		        newPair = distKNN.kNN.top();
		        distKNN.kNN.pop();
		        discriminatorHist->Fill(Metas[newPair.second]);//,AccWeights[newPair.second]);
		        discriminatorHist2->Fill(Mpi0s[newPair.second]);//,AccWeights[newPair.second]);
		        stdCalc.insertValue(Metas[newPair.second]);
		        stdCalc2.insertValue(Mpi0s[newPair.second]);
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
		
		bestChiSq_eta=DBL_MAX;
		int best_iFit=0;
		bestChiSq_pi0=DBL_MAX;
		int best_iFit2=0;
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
			fit = new TF1("fit",fitFunc,fitRange[0],fitRange[1],numDOFbkg+numDOFsig);
			bkgFit = new TF1("bkgFit",background,fitRange[0],fitRange[1],numDOFbkg);
			sigFit = new TF1("sigFit",signalDG,fitRange[0],fitRange[1],numDOFsig);
			//fit2 = new TF1("fit2",fitFunc,fitRange2[0],fitRange2[1],numDOFbkg+numDOFsig);
			//bkgFit2 = new TF1("bkgFit2",background,fitRange2[0],fitRange2[1],numDOFbkg);
			//sigFit2 = new TF1("sigFit2",signal,fitRange2[0],fitRange2[1],numDOFsig);
			
			//dataFit.fitToData(iFit);
			
			// Should use getInitParams.C whenever we get a new dataset to initialize the peak and width of the pi0 and eta
			fit->SetParameters(par0eta[iFit],par1eta[iFit],par2eta[iFit],peakWidtheta[0],peakWidtheta[1],ampRatio,widthRatio);
			fit->SetParLimits(3,peakWidtheta[0]*0.95, peakWidtheta[0]*1.05);
			fit->SetParLimits(4,peakWidtheta[1]*0.95, peakWidtheta[1]*1.05); 
			fit->SetParLimits( 5,ampRatio*0.95,ampRatio*1.05 );
			fit->SetParLimits( 6,widthRatio*0.95, widthRatio*1.05 );

			// we have to calculate the q-value for the eta distribution to check if it is between 0 and 1 
			discriminatorHist->Fit("fit","RQBNL"); // B will enforce the bounds, N will be no draw
			fit->GetParameters(par);
			bkgFit->SetParameters(par);
			sigFit->SetParameters(&par[numDOFbkg]);
			qvalue_eta=sigFit->Eval(Metas[ientry])/fit->Eval(Metas[ientry]);
			chiSq_eta = fit->GetChisquare()/(fit->GetNDF());
			
			// save some fit and the chiSq_eta for all the fits
			showConv[iFit]->SetParameters(par); // save all the fits
			chiSqs[iFit]=chiSq_eta;
			
			// The pi0 fit
			//discriminatorHist2->Fit("fit2","RQBNL"); // B will enforce the bounds, N will be no draw
			//chiSq_pi0 = fit2->GetChisquare()/(fit2->GetNDF());
			if (verbose2) { logFile << "current ChiSq, best ChiSq: " << chiSq_eta << ", " << bestChiSq_eta << endl; }
			//if (verbose2) { logFile << "current ChiSq2, best ChiSq2: " << chiSq_pi0 << ", " << bestChiSq_pi0 << endl; }
			
			if (chiSq_eta < bestChiSq_eta){
				best_qvalue_eta = qvalue_eta;
				bestChiSq_eta=chiSq_eta;
				best_iFit=iFit;
				for (int i=0; i < sizeof(par)/sizeof(Double_t); ++i){
					parBest[i]=par[i];
				}
				if(qvalue_eta>1 || qvalue_eta<0) {
				      cout << "qvalue_eta out of bounds!\n-------------" << endl;
				      for (double parVal : par){
				          cout << parVal << endl;
				      } 
				      cout << "---- BREAKING ----" << endl;
				      //exit(0);
				}
			} 
		}


		// /////////////////////////////////////////
		// need to calcuclate new q-value since it is out of bounds
		// /////////////////////////////////////////
		if (best_qvalue_eta>1 || best_qvalue_eta<0){
			cout << "Using flat fit instead of linear on event:" << ientry << endl;
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
			
		        fit->SetParameters(par0eta[0],0,par2eta[0],peakWidtheta[0],peakWidtheta[1],ampRatio,widthRatio);
			fit->FixParameter(1,0); 
			discriminatorHist->Fit("fit","RQBNL"); // B will enforce the bounds, N will be no draw
			fit->GetParameters(par);
			bkgFit->SetParameters(par);
			sigFit->SetParameters(&par[numDOFbkg]);

			qvalue_eta=sigFit->Eval(Metas[ientry])/fit->Eval(Metas[ientry]);
			if (qvalue_eta>1 || qvalue_eta<0){
			    cout << "Not sure why qvalue_eta is still >1 or <0. Need to fix this!" << endl;
			    cout << "These are the parameters:"<<endl;
			    for ( double parVal : par ){
					cout << " " << parVal << endl;
			    }
			    cout << " **************** BREAKING ****************** " << endl;
			    //exit(0);
			}

			best_qvalue_eta = qvalue_eta;
			chiSq_eta = fit->GetChisquare()/(fit->GetNDF());
			bestChiSq_eta=chiSq_eta;
			best_iFit=4;
			for (int i=0; i < sizeof(par)/sizeof(Double_t); ++i){
				parBest[i]=par[i];
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
		if ( best_iFit2 == 0){ ++nBest100Bkg2; }
		else if ( best_iFit2 == 1){ ++nBest100Sig2; }
		else if ( best_iFit2 == 2){ ++nBest50Bkg50Sig2; }
		
		// Here we draw the histograms that were randomly selected
		if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end()) {
			if(verbose) { cout << "Keeping this event" << endl; }
			discriminatorHist->SetTitle(("BEST VALS:  QValue="+to_string(best_qvalue_eta)+"  ChiSq="+to_string(bestChiSq_eta)+"     Std="+to_string(comboStd)+"    iFit="+to_string(best_iFit)).c_str() );
			//discriminatorHist2->SetTitle(("BEST VALS:  QValue="+to_string(best_qvalue_pi0)+"  ChiSq="+to_string(bestChiSq_pi0)+"     Std="+to_string(comboStd2)+"    iFit="+to_string(best_iFit2)).c_str());
			etaLine = new TLine(Metas[ientry],0,Metas[ientry],kDim);
			pi0Line = new TLine(Mpi0s[ientry],0,Mpi0s[ientry],kDim);
			etaLine->SetLineColor(kOrange);
			pi0Line->SetLineColor(kOrange);
			
                  	fit->SetParameters(parBest);
                  	//fit2->GetParameters(par2);
	          	bkgFit->SetParameters(parBest);
	          	//bkgFit2->SetParameters(par2);
	          	sigFit->SetParameters(&parBest[numDOFbkg]);
	          	//sigFit2->SetParameters(&par2[numDOFbkg]);
                  	fit->SetLineColor(kRed);
                  	//fit2->SetLineColor(kRed);
  	          	bkgFit->SetFillColor(kMagenta);
                  	bkgFit->SetLineColor(kMagenta);
  	          	bkgFit->SetFillStyle(3004);
  	          	//bkgFit2->SetFillColor(kMagenta);
                  	//bkgFit2->SetLineColor(kMagenta);
  	          	//bkgFit2->SetFillStyle(3004);
  	          	sigFit->SetFillColor(kBlue);
                  	sigFit->SetLineColor(kBlue);
  	          	sigFit->SetFillStyle(3005);
  	          	//sigFit2->SetFillColor(kBlue);
                  	//sigFit2->SetLineColor(kBlue);
  	          	//sigFit2->SetFillStyle(3005);

			qSigLine = new TLine(0,sigFit->Eval(Metas[ientry]),binRange[2],sigFit->Eval(Metas[ientry]));
			qBkgLine = new TLine(0,bkgFit->Eval(Metas[ientry]),binRange[2],bkgFit->Eval(Metas[ientry]));
                  	qSigLine->SetLineColor(kBlue);
			qSigLine->SetLineStyle(9);
                  	qBkgLine->SetLineColor(kMagenta);
			qBkgLine->SetLineStyle(9);


                  	allCanvases->cd(1);
	          	discriminatorHist->Draw();
                  	drawText(parBest,numDOFbkg+numDOFsig,"par",sigFit->Eval(Metas[ientry]),bkgFit->Eval(Metas[ientry]),fit->Eval(Metas[ientry]));
	          	etaLine->Draw("same");
                  	fit->Draw("SAME");
  	          	bkgFit->Draw("SAME FC");
  	          	sigFit->Draw("SAME FC");
			qBkgLine->Draw("SAME");
			qSigLine->Draw("SAME");
                  	allCanvases->cd(2);
	          	discriminatorHist2->Draw();
                  	//drawText(par2,numDOFbkg+numDOFsig,"par");
	          	pi0Line->Draw("same");
                  	//fit2->Draw("SAME");
  	          	//bkgFit2->Draw("SAME FC");
  	          	//sigFit2->Draw("SAME FC");

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
			drawText(chiSqs,3,"chiSq_eta",0,0,0);
        		legend_conv->Draw();
			//allCanvases->SaveAs(("histograms/fitCheck-event"+std::to_string(ientry)+".pdf").c_str());

			//if(verbose2){
			//     duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
			//     logFile << "	Best ChiSq (" << bestChiSq_eta << ") with Q-value: " << qvalue_eta << " -- time: " << duration2 <<  "ms" << endl;
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
	auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
	//logFile << "Total time: " << duration2 << " ms" << endl;
	//logFile << "Time Per Event: " << ( std::clock() - start ) / (double) CLOCKS_PER_SEC / nentries << "ns" << endl;
	logFile.close();
        return 0;
}

