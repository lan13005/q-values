//The goal of this program is to read in all the q-value results and organizes them to be used when creating any histogram
#include "helperFuncs.h"

void mergeQresults(){
        cout << "\n========================" << endl;
        cout << "Merging QFactors results with the original root file" << endl;
        cout << "1. Load qfactors result tree" << endl;
        cout << "2. Load input Data tree -- Clone -- Add new branches -- Fill" << endl;
        cout << "========================" << endl;



	// ---------------------------------------------------------------------------
	// ------------------------------Settting branch addresses for Q-Value Results
	// ---------------------------------------------------------------------------
	// Read in the qvalue data
	string preQFileName = "logs"+runTag+"/"+fileTag+"/resultsMerged_"+fileTag+".root";
	cout << "\n------------------\nOpening " << preQFileName << endl;
	TFile* qResultsFile = new TFile((preQFileName).c_str());
        TTree* qResultsTree;
	qResultsFile->GetObject("resultsTree",qResultsTree);
        cout << "Loaded the Q-factor results file..." << endl;

	double qvalue;
	double conjugate_qvalue;
	double bestNLL;
	double worstNLL;
        double worst_qvalue;
        double qvalueBS_std;
        double eff_nentries;
        int kDim1; // kDim1 is an int whereas kDim is a const int so we can make an array out of it
        int neighbors[kDim];

	qResultsTree->SetBranchAddress("qvalue",&qvalue);
	qResultsTree->SetBranchAddress("bestNLL",&bestNLL);
	qResultsTree->SetBranchAddress("worstNLL",&worstNLL);
        qResultsTree->SetBranchAddress("worst_qvalue",&worst_qvalue);
        qResultsTree->SetBranchAddress("qvalueBS_std",&qvalueBS_std);
        qResultsTree->SetBranchAddress("eff_nentries",&eff_nentries);
        qResultsTree->SetBranchAddress("kDim",&kDim1);
        qResultsTree->SetBranchAddress("neighbors",neighbors);
        cout << "Setting branch address..." << endl;

        int nentries;
	nentries=qResultsTree->GetEntries();
        cout << "Total number of events in Q-factor results: " << nentries << endl;

	// ---------------------------------------------------------------------------
	// ----------------------------------------------------- Loading the datafile data
	// and clone the datafile and add new branches to track the qvalue, NLL, etc
        // Now we have a root tree that has all the data from the DSelector and the q-factors analysis
	// ---------------------------------------------------------------------------
	string inputFileLoc = rootFileLoc;
	cout << "\n------------------\nOpening " << preQFileName << endl;
	TFile* dataFile=new TFile((inputFileLoc).c_str());
	TTree *dataTree;
	dataFile->GetObject((rootTreeName).c_str(),dataTree);
        cout << "Entries in input root tree: " << dataTree->GetEntries()  << endl;
	

	string postQFileName = "logs"+runTag+"/"+fileTag+"/postQVal_flatTree_"+fileTag+".root";	
	TFile *qd_dataFile = TFile::Open((postQFileName).c_str(),"RECREATE"); 
	TTree *outputTree = dataTree->CloneTree(-1,"fast"); 
	cout << "Cloned input to output file" << postQFileName << endl;

	TBranch* b_qvalue;
	TBranch* b_NLLBest;
	TBranch* b_NLLWorst;
        TBranch* b_worst_qvalue;
        TBranch* b_qvalueBS_std;
        TBranch* b_eff_nentries;
        TBranch* b_kDim;
        TBranch* b_neighbors;
	b_qvalue = outputTree->Branch("qvalue",&qvalue,"qvalue/D");
	b_NLLBest = outputTree->Branch("NLLBest",&bestNLL,"qvalue_NLLBest/D");
	b_NLLWorst = outputTree->Branch("NLLWorst",&worstNLL,"qvalue_NLLWorst/D");
	b_worst_qvalue = outputTree->Branch("worst_qvalue",&worst_qvalue,"worst_qvalue/D");
	b_qvalueBS_std = outputTree->Branch("qvalueBS_std",&qvalueBS_std,"qvalueBS_std/D");
        b_eff_nentries = outputTree->Branch("eff_nentries",&eff_nentries,"eff_nentries/D");
        b_kDim = outputTree->Branch("kDim",&kDim1,"kDim/I");
        b_neighbors = outputTree->Branch("neighbors",&neighbors,"neighbors[kDim]/I");
        cout << "Added new branches to output file..." <<endl;

	for (int ientry=0; ientry<nentries; ientry++)
	{
            qResultsTree->GetEntry(ientry);
	    b_qvalue->Fill();
	    b_NLLBest->Fill();
	    b_NLLWorst->Fill();
	    b_worst_qvalue->Fill();
	    b_qvalueBS_std->Fill();
            b_eff_nentries->Fill();
            b_kDim->Fill();
            b_neighbors->Fill();
	}
        cout << "Filled new output branches with qfactor results" << endl;
	qResultsFile->Close();
	qd_dataFile->cd();
	outputTree->Write(); 
        cout << "FINIHSED MERGING Q RESULTS!"<<endl;
}

