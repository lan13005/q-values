#include "main.h"

//int test(int kDim, string varString, int numberEventsToSavePerProcess, int nProcess, int seedShift, Long64_t nentries, bool override_nentries, bool verbose){
int main( int argc, char* argv[] ){
	ROOT::EnableThreadSafety();
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
	TH1::AddDirectory(kFALSE);

	// This suppresses all the "info" messages from roofit. Like saying that it will use numeric integrator for the generic pdf we define
	// https://root-forum.cern.ch/t/suppressing-info-messages/14642/6
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
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
        int numberEventsToSavePerProcess= atoi(argv[3]);
	int nProcess=atoi(argv[4]);
        int seedShift= atoi(argv[5]);
        Long64_t nentries= atoi(argv[6]);
        bool override_nentries;
        cout << "Num Vars="<< dim << endl;
        if ( atoi(argv[7]) ==1 ){ 
            override_nentries=true;
        }
        else{ 
            override_nentries=false;
        }
        bool verbose;
        if ( atoi(argv[8])==1 ){ verbose=true;}
        else{ verbose=false;}
        cout << "----------------------------" << endl;
        cout << "kDim: " << kDim << endl;
        cout << "numberEventsToSavePerProcess: " << numberEventsToSavePerProcess << endl;
        cout << "seedShift: " << seedShift << endl;
        cout << "nentries: " << nentries << endl; 
        cout << "override_nentries: " << override_nentries << endl;
        cout << "verbose: " << verbose  << endl; 
	cout << "varString: " << varString << endl; 
        cout << "nProcess: " << nProcess << endl;
        cout << "----------------------------" << endl;
    
	QFactorAnalysis analysisControl(kDim, varString, numberEventsToSavePerProcess, nProcess, seedShift, nentries, override_nentries, verbose);
	string rootFileLocation = "degALL_"+detector+"_treeFlat_DSelector.root";
	string rootTreeName = "degALL_"+detector+"_tree_flat";
	string fitLocation = "fitResults/etaFitNoAccSub_"+detector+".txt";
	analysisControl.loadTree(rootFileLocation, rootTreeName);
	analysisControl.loadFitParameters(fitLocation);
	analysisControl.loadData();
	analysisControl.runQFactorThreaded();
        return 0;
}


