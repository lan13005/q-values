// "main" outputs results files that contain the q-values. The goal of tihs program is to read in all the q-value results and organizes them to be used when creating any histogram


#include <ctime>
#include <math.h> 
#include "makeDiagnosticHists.h"
#include "helperFuncs.h"


bool verbose = true;
string fileTag="bcal";
string rootFileLoc="/home/lawrence/Desktop/gluex/q-values/degALL_bcal_treeFlat_DSelector_UTweights.root";
string rootTreeName="degALL_bcal_tree_flat";
string weightingScheme="as"; // "" or "as*bs"
string s_accWeight="AccWeight";
string s_discrimVar="Meta";
string s_sideBandVar="Mpi0";
string s_uniquenessTracking="weighted"; // "default" or "weighted". Any neighbor is possible in the weighted setting

void makeDiagnosticHists(){
	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
        gStyle->SetOptStat(0);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	
	TH1F* dHist_bestChiSq= new TH1F( "dHist_bestChiSq", "#Chi^2; Events/0.01", 50, 0, 50);
	TH1F* dHist_deltaChiSq= new TH1F( "dHist_deltaChiSq", "Best #Chi^2 - Worst #Chi^2; Events/0.5", 50, -1, 0);

	TH1F* dHist_qvalues = new TH1F( "dHist_qvalues", "Q-Values; Events/0.01", 100, 0, 1);

	// ---------------------------------------------------------------------------
	// ------------------------------Settting branch addresses for Q-Value Results
        // AFTER THIS SECTION, ALL THE Q-FACTORS SHOULD BE ORDERED, WE JUST HAVE TO 
        // READ IN THE DATA AND AND FILL THEM. YOU PROABLY DONT HAVE TO EDIT THIS SECTION
        // JUST DEFINE YOUR HISTOGRAMS AND DRAW THEM
	// ---------------------------------------------------------------------------
        

	// Read in the qvalue data
	string preQFileName = "diagnosticPlots/"+fileTag+"/qvalResults_"+fileTag+".root";
	cout << "Opening " << preQFileName << endl;
	TFile* qResultsFile = new TFile((preQFileName).c_str());
        TTree* qResultsTree;
	qResultsFile->GetObject("resultsTree",qResultsTree);
        cout << "Loaded the Q-factor results file" << endl;

	Double_t qvalue;
	Double_t conjugate_qvalue;
	Double_t bestChiSq;
	Double_t worstChiSq;
        Bool_t bool_MetaFlat;
	ULong64_t flatEntryNumber;

	qResultsTree->SetBranchAddress("qvalue",&qvalue);
	qResultsTree->SetBranchAddress("bestChiSq",&bestChiSq);
	qResultsTree->SetBranchAddress("worstChiSq",&worstChiSq);
	qResultsTree->SetBranchAddress("flatEntryNumber",&flatEntryNumber);


	ULong64_t nentries;
	nentries=qResultsTree->GetEntries();
	const int c_nentriesResults = (const int)nentries;
        cout << "Total number of events in Q-factor results: " << c_nentriesResults << endl;


	// make a vector of 0s such that we can fill the q-values in order to unscramble to multiprocessing effects
	std::vector< double > qvalues(c_nentriesResults,0);
	std::vector< double > bestChiSqs(c_nentriesResults,0);
	std::vector< double > worstChiSqs(c_nentriesResults,0);
        std::vector<ULong64_t> flatEntryNumbers(c_nentriesResults,0);
	for (int ientry=0; ientry<c_nentriesResults; ientry++)
	{
        	qResultsTree->GetEntry(ientry);
        	//cout << "ientry, flatEntryNumber: " << ientry << "," << flatEntryNumber << endl;            
        	
        	qvalues[flatEntryNumber]=qvalue;
        	bestChiSqs[flatEntryNumber]=bestChiSq;
        	worstChiSqs[flatEntryNumber]=worstChiSq;
        	flatEntryNumbers[flatEntryNumber] = flatEntryNumber;
        	
        	dHist_qvalues->Fill(qvalue);
        	dHist_bestChiSq->Fill(bestChiSq);
        	dHist_deltaChiSq->Fill(bestChiSq-worstChiSq);
        }
	qResultsFile->Close();
        cout << "Loaded all the q-factor results from the results file" << endl;

        allCanvases->Clear();
	dHist_qvalues->Draw();
        allCanvases->SaveAs(("diagnosticPlots/"+fileTag+"/qvalues.png").c_str());
        allCanvases->Clear();
	dHist_bestChiSq->Draw();
        allCanvases->SaveAs(("diagnosticPlots/"+fileTag+"/bestChiSq.png").c_str());
        allCanvases->Clear();
        gPad->SetLogy(1);
	dHist_deltaChiSq->Draw();
        allCanvases->SaveAs(("diagnosticPlots/"+fileTag+"/deltaChiSq.png").c_str());
        gPad->SetLogy(0);
        cout << "Making some diagnostic distributions of the Q-factors results" << endl;

	// ---------------------------------------------------------------------------
	// ------------------------------------------------MANUALLY DEFINE HISTOGRAMS
	// ---------------------------------------------------------------------------
	TH1F* Meta_tot[2];
	TH1F* Meta_sig[2];
	TH1F* Meta_bkg[2];
	TH1F* Mpi0_tot[2];
	TH1F* Mpi0_sig[2];
	TH1F* Mpi0_bkg[2];
	TH1F* Mpi0eta_tot[2];
	TH1F* Mpi0eta_sig[2];
	TH1F* Mpi0eta_bkg[2];
	TH1F* cosThetaEta_GJ_tot[2];
	TH1F* cosThetaEta_GJ_sig[2];
	TH1F* cosThetaEta_GJ_bkg[2];
	TH1F* cosThetaX_CM_tot[2];
	TH1F* cosThetaX_CM_sig[2];
	TH1F* cosThetaX_CM_bkg[2];
	TH1F* phiEta_GJ_tot[2];
	TH1F* phiEta_GJ_sig[2];
	TH1F* phiEta_GJ_bkg[2];

	TH1F *Mpi0g_tot = new TH1F( "Mpi0g_tot", "M(#pi^{0}g); M(#pi^{0}g) GeV; Events/0.001 GeV", 350, 0, 3.5 );
	TH1F *Mpi0g_sig = new TH1F( "Mpi0g_sig", "M(#pi^{0}g); M(#pi^{0}g) GeV; Events/0.001 GeV", 350, 0, 3.5 );
	TH1F *Mpi0g_bkg = new TH1F( "Mpi0g_bkg", "M(#pi^{0}g); M(#pi^{0}g) GeV; Events/0.001 GeV", 350, 0, 3.5 );
	for (int i=0; i<2; i++){
		string tag="";
		if (i%2==0){ tag="meas"; }
		else { tag="kin"; }
		// WE SHOULD TRY TO MATCH THE BINS USED TO FIT, USED IN THE Q-VALUE CALCULATION, AND HERE. USING IT IN THE FIT AND DURING GRAPHS WILL MAKE FOR GOOD COMPARISION
		// AND USING FOR IN THE FIT AND THE Q-VALUE WOULD ALSO BE GOOD SINCE THE FIT COULD CHANGE DEPENDING ON THE BINNING
		Meta_tot[i] = new TH1F( ("Meta_tot_"+tag).c_str(), "M(#eta); M(#eta) GeV; Events/0.003 GeV", 200, 0.25, 0.85);//Events/0.002 GeV", 300, 0.25, 0.85 );
		Meta_sig[i] = new TH1F( ("Meta_sig_"+tag).c_str(), "M(#eta); M(#eta) GeV; Events/0.003 GeV", 200, 0.25, 0.85);//Events/0.002 GeV", 300, 0.25, 0.85 );
		Meta_bkg[i] = new TH1F( ("Meta_bkg_"+tag).c_str(), "M(#eta); M(#eta) GeV; Events/0.003 GeV", 200, 0.25, 0.85);//Events/0.002 GeV", 300, 0.25, 0.85 );
		Mpi0_tot[i] = new TH1F( ("Mpi0_tot_"+tag).c_str(), "M(#pi^{0}); M(#pi^{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
		Mpi0_sig[i] = new TH1F( ("Mpi0_sig_"+tag).c_str(), "M(#pi^{0}); M(#pi^{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
		Mpi0_bkg[i] = new TH1F( ("Mpi0_bkg_"+tag).c_str(), "M(#pi^{0}); M(#pi^{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
		Mpi0eta_tot[i] = new TH1F( ("Mpi0eta_tot_"+tag).c_str(), "M(#pi^{0}#eta); M(#pi^{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
		Mpi0eta_sig[i] = new TH1F( ("Mpi0eta_sig_"+tag).c_str(), "M(#pi^{0}#eta); M(#pi^{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
		Mpi0eta_bkg[i] = new TH1F( ("Mpi0eta_bkg_"+tag).c_str(), "M(#pi^{0}#eta); M(#pi^{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
		cosThetaEta_GJ_tot[i] = new TH1F( ("cosThetaEta_GJ_tot_"+tag).c_str(), "cos(#theta) GJ of #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, -1, 1 );
		cosThetaEta_GJ_sig[i] = new TH1F( ("cosThetaEta_GJ_sig_"+tag).c_str(), "cos(#theta) GJ of #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, -1, 1 );
		cosThetaEta_GJ_bkg[i] = new TH1F( ("cosThetaEta_GJ_bkg_"+tag).c_str(), "cos(#theta) GJ of #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, -1, 1 );
		cosThetaX_CM_tot[i] = new TH1F( ("cosThetaX_CM_tot_"+tag).c_str(), "cos(#theta) of CM #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, 0.9, 1 );
		cosThetaX_CM_sig[i] = new TH1F( ("cosThetaX_CM_sig_"+tag).c_str(), "cos(#theta) of CM #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, 0.9, 1 );
		cosThetaX_CM_bkg[i] = new TH1F( ("cosThetaX_CM_bkg_"+tag).c_str(), "cos(#theta) of CM #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, 0.9, 1 );
		phiEta_GJ_tot[i] = new TH1F( ("phiEta_GJ_tot_"+tag).c_str(), "#phi of #eta;#phi of #eta ;Events/0.02 GeV", 90, -180, 180 );
		phiEta_GJ_sig[i] = new TH1F( ("phiEta_GJ_sig_"+tag).c_str(), "#phi of #eta;#phi of #eta ;Events/0.02 GeV", 90, -180, 180 );
		phiEta_GJ_bkg[i] = new TH1F( ("phiEta_GJ_bkg_"+tag).c_str(), "#phi of #eta;#phi of #eta ;Events/0.02 GeV", 90, -180, 180 );
	}
        cout << "Defined all histograms" << endl;



	// ---------------------------------------------------------------------------
	// ----------------------------------------------------- Set up branch address of the data file 
	// But fill it after we load the q-value data so we can write the data to the tree
	// I have a felling it might be bad to fill the q-values as we calculate it during the main program
	// since we could have 36+ processes all calculating and writing to the same branch
	// ---------------------------------------------------------------------------
	//
	// setting up some basic root stuff and getting the file and tree
	string inputFileLoc = rootFileLoc;
	cout << "Loading the entire data root file" << inputFileLoc << endl;
	TFile* dataFile=new TFile((inputFileLoc).c_str());
	TTree *dataTree;
	dataFile->GetObject((rootTreeName).c_str(),dataTree);
	
        TLine* etaLine;
        double AccWeight;
	double Meta;
	double Mpi0;
	double Mpi0g1;
	double Mpi0g2;
	double Mpi0eta;
	double cosTheta_X_cm;
	double cosTheta_eta_gj;
	double phi_eta_gj;
	double Meta_meas;
	double Mpi0_meas;
	double Mpi0eta_meas;
	double cosTheta_X_cm_meas;
	double cosTheta_eta_gj_meas;
	double phi_eta_gj_meas;
        double utWeight;

        dataTree->SetBranchAddress("AccWeight",&AccWeight);
	dataTree->SetBranchAddress("Meta_meas",&Meta_meas);
	dataTree->SetBranchAddress("Mpi0_meas",&Mpi0_meas);
	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("Mpi0",&Mpi0);

	dataTree->SetBranchAddress("Mpi0eta_meas",&Mpi0eta_meas);
        dataTree->SetBranchAddress("cosTheta_X_cm_meas", &cosTheta_X_cm_meas); 
        dataTree->SetBranchAddress("cosTheta_eta_gj_meas",&cosTheta_eta_gj_meas);
        dataTree->SetBranchAddress("phi_eta_gj_meas",&phi_eta_gj_meas); 
	dataTree->SetBranchAddress("Mpi0g1",&Mpi0g1);
	dataTree->SetBranchAddress("Mpi0g2",&Mpi0g2);
	dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
        dataTree->SetBranchAddress("cosTheta_X_cm", &cosTheta_X_cm); 
        dataTree->SetBranchAddress("cosTheta_eta_gj",&cosTheta_eta_gj);
        dataTree->SetBranchAddress("phi_eta_gj",&phi_eta_gj); 

        // continue to do uniqueness tracking for the histograms we care about
        bool isUniqueEtaB;
        bool isUniquePi0B;
        bool isUniquePi0g1B;
        bool isUniquePi0g2B;
        bool isUniquePi0EtaB;
        dataTree->SetBranchAddress("isNotRepeated_eta",&isUniqueEtaB);
        dataTree->SetBranchAddress("isNotRepeated_pi0",&isUniquePi0B);
        dataTree->SetBranchAddress("isNotRepeated_pi0g1",&isUniquePi0g1B);
        dataTree->SetBranchAddress("isNotRepeated_pi0g2",&isUniquePi0g2B);
        dataTree->SetBranchAddress("isNotRepeated_pi0eta",&isUniquePi0EtaB);
        dataTree->SetBranchAddress("uniqunessTrackingWeights",&utWeight);

        cout << "Finished setting up the tree" << endl;
	
	nentries=dataTree->GetEntries();
	const int c_nentriesData = (const int)nentries;
        cout << "c_nentriesData: " << c_nentriesData << endl;
	if ( c_nentriesData != c_nentriesResults ){
		cout << "Entries are not equal between dataFile and resultsFile. Ran over a subset" << endl;
                cout << "Only making histograms for the subset" << endl;
	}
	else{
		cout << "Entries are equal. Ran over the full dataset" << endl;
	}

        std::vector<double> AccWeights; AccWeights.reserve(c_nentriesResults); 
        std::vector<bool> isUniqueEtaBs; isUniqueEtaBs.reserve(c_nentriesResults); 
        std::vector<bool> isUniquePi0Bs; isUniquePi0Bs.reserve(c_nentriesResults); 
        std::vector<bool> isUniquePi0g1Bs; isUniquePi0g1Bs.reserve(c_nentriesResults); 
        std::vector<bool> isUniquePi0g2Bs; isUniquePi0g2Bs.reserve(c_nentriesResults); 
        std::vector<bool> isUniquePi0EtaBs; isUniquePi0EtaBs.reserve(c_nentriesResults); 
        std::vector<double> utWeights; utWeights.reserve(c_nentriesResults);
        std::vector<double> Metas_meas; Metas_meas.reserve(c_nentriesResults); 
        std::vector<double> Mpi0s_meas; Mpi0s_meas.reserve(c_nentriesResults); 
        std::vector<double> Mpi0etas_meas; Mpi0etas_meas.reserve(c_nentriesResults); 
	std::vector<double> cosTheta_X_cms_meas; cosTheta_X_cms_meas.reserve(c_nentriesResults);
	std::vector<double> cosTheta_eta_gjs_meas; cosTheta_eta_gjs_meas.reserve(c_nentriesResults);
	std::vector<double> phi_eta_gjs_meas; phi_eta_gjs_meas.reserve(c_nentriesResults);
        std::vector<double> Metas; Metas.reserve(c_nentriesResults); 
        std::vector<double> Mpi0s; Mpi0s.reserve(c_nentriesResults); 
        std::vector<double> Mpi0g1s; Mpi0g1s.reserve(c_nentriesResults); 
        std::vector<double> Mpi0g2s; Mpi0g2s.reserve(c_nentriesResults); 
        std::vector<double> Mpi0etas; Mpi0etas.reserve(c_nentriesResults); 
	std::vector<double> cosTheta_X_cms; cosTheta_X_cms.reserve(c_nentriesResults);
	std::vector<double> cosTheta_eta_gjs; cosTheta_eta_gjs.reserve(c_nentriesResults);
	std::vector<double> phi_eta_gjs; phi_eta_gjs.reserve(c_nentriesResults);

	std::vector< double > sbWeights; sbWeights.reserve(c_nentriesResults);
	double sbWeight;

        // ***************** CHECK 1 ********************
        for (int iEntry=0; iEntry<c_nentriesResults; ++iEntry){
            if ( flatEntryNumbers[iEntry] != iEntry){
                cout << "flatEntryNumbers[iEntry] != iEntry: flatEntryNumbers[iEntry], iEntry: " << flatEntryNumbers[iEntry] << ", " << iEntry << endl;
            }
        }

        // **************** CHECK 2 ********************
        //std::set<int> difference;
        //std::set_difference(all_ientries.begin(), all_ientries.end(), all_flatEntryNumber.begin(), all_flatEntryNumber.end(),
        //    std::inserter(difference, difference.end()));
        //if ( difference.size() != 0 ) {
        //    cout << "These elements are not common in both sets! Fix me!\n\t\t";
        //    for (auto it=difference.begin(); it != difference.end(); ++it){ 
        //        cout << ' ' << *it; 
        //    }
        //    cout << endl;
        //}
        //else { cout << "The flatEntryNumber set is complete!" << endl; }

        //cout << "Finished consistency check on results file and data" << endl;

	// ---------------------------------------------------------------------------
	// ----------------------------------------------------- Loading the datafile data
	// and clone the datafile and add 3 new branches to track the qvalue, chiSq, flatEntryNumber
        // Now we have a root tree that has all the data from the DSelector and the q-factors analysis
	// ---------------------------------------------------------------------------

	string postQFileName = "diagnosticPlots/"+fileTag+"/postQ_"+fileTag+"_flatTree.root";	
	cout << "Remaking " << postQFileName << endl;
	TFile *qd_dataFile = TFile::Open((postQFileName).c_str(),"RECREATE"); 
	TTree *outputTree = dataTree->CloneTree(-1,"fast"); 
	TBranch* b_qvalue;
	TBranch* b_qvalue_chisq;
	TBranch* b_flatEntryNumber;
	b_qvalue = outputTree->Branch("qvalue",&qvalue,"qvalue/D");
	b_qvalue_chisq = outputTree->Branch("qvalue_chisqBest",&bestChiSq,"qvalue_chisqBest/D");
	b_qvalue_chisq = outputTree->Branch("qvalue_chisqWorst",&worstChiSq,"qvalue_chisqWorst/D");
	b_flatEntryNumber = outputTree->Branch("flatEntryNumber",&flatEntryNumber,"flatEntryNumber/l");

	for (int ientry=0; ientry<c_nentriesResults; ientry++)
	{
		outputTree->GetEntry(ientry);
                AccWeights.push_back(AccWeight);

                getSBWeight(Mpi0,&sbWeight,weightingScheme);
		sbWeights.push_back(sbWeight);

		Metas_meas.push_back(Meta_meas);
		Mpi0s_meas.push_back(Mpi0_meas);
		Mpi0etas_meas.push_back(Mpi0eta_meas);
		cosTheta_X_cms_meas.push_back(cosTheta_X_cm_meas);
		cosTheta_eta_gjs_meas.push_back(cosTheta_eta_gj_meas);
		phi_eta_gjs_meas.push_back(phi_eta_gj_meas);
		Metas.push_back(Meta);
		Mpi0s.push_back(Mpi0);
		Mpi0g1s.push_back(Mpi0g1);
		Mpi0g2s.push_back(Mpi0g2);
		Mpi0etas.push_back(Mpi0eta);
		cosTheta_X_cms.push_back(cosTheta_X_cm);
		cosTheta_eta_gjs.push_back(cosTheta_eta_gj);
		phi_eta_gjs.push_back(phi_eta_gj);

                isUniqueEtaBs.push_back(isUniqueEtaB);
                isUniquePi0Bs.push_back(isUniquePi0B);
                isUniquePi0g1Bs.push_back(isUniquePi0g1B);
                isUniquePi0g2Bs.push_back(isUniquePi0g2B);
                isUniquePi0EtaBs.push_back(isUniquePi0EtaB);
                utWeights.push_back(utWeight);
                
		// will take the opportunity to save the now-ordered qvalues into the root file as we load the rest of the data
		qvalue = qvalues[ientry];
		bestChiSq = bestChiSqs[ientry];
		worstChiSq = worstChiSqs[ientry];
		flatEntryNumber = flatEntryNumbers[ientry];
		b_qvalue->Fill();
		b_qvalue_chisq->Fill();
		b_flatEntryNumber->Fill();
	}

        cout << "Imported all the tree data into arrays" << endl;
        cout << "Size of array: " << c_nentriesData << endl;

	qd_dataFile->cd();
	outputTree->Write(); 

    
        // We can now fill the histograms, properly weighted
	double sigWeight;
	double totWeight;
	double bkgWeight;
	int numNan=0;
	for (int ientry=0; ientry<c_nentriesResults; ientry++){
		qvalue = qvalues[ientry];
		conjugate_qvalue = 1-qvalue;
                AccWeight=AccWeights[ientry];
                sbWeight=sbWeights[ientry];
                if(qvalue != qvalue){
                    //cout << "ientry is nan: " << ientry << endl;
		    ++numNan;
                }
                double weight; 
                if (weightingScheme==""){ weight=1; }
                if (weightingScheme=="as"){ weight=AccWeights[ientry]; }
                if (weightingScheme=="as*bs"){ weight=AccWeights[ientry]*sbWeights[ientry]; }
		sigWeight = qvalue*weight;
		totWeight = weight;
		bkgWeight = conjugate_qvalue*weight;
                //cout << "ientry,qVal,conj_qVal,AccWeight,sbWeight: " << ientry << "," << qvalue << "," << conjugate_qvalue << ", " << AccWeight << ", " << sbWeight << endl;
                if (s_uniquenessTracking=="default"){
                    if ( isUniqueEtaBs[ientry] ) {
		    	cosThetaEta_GJ_sig[0]->Fill(cosTheta_eta_gjs_meas[ientry], sigWeight);
		    	cosThetaEta_GJ_tot[0]->Fill(cosTheta_eta_gjs_meas[ientry], totWeight);
		    	cosThetaEta_GJ_bkg[0]->Fill(cosTheta_eta_gjs_meas[ientry], bkgWeight);
		    	cosThetaX_CM_sig[0]->Fill(cosTheta_X_cms_meas[ientry], sigWeight);
		    	cosThetaX_CM_tot[0]->Fill(cosTheta_X_cms_meas[ientry], totWeight);
		    	cosThetaX_CM_bkg[0]->Fill(cosTheta_X_cms_meas[ientry], bkgWeight);
		    	Meta_sig[0]->Fill(Metas_meas[ientry],sigWeight);
		    	Meta_tot[0]->Fill(Metas_meas[ientry],totWeight);
		    	Meta_bkg[0]->Fill(Metas_meas[ientry],bkgWeight);

		    	cosThetaEta_GJ_sig[1]->Fill(cosTheta_eta_gjs[ientry], sigWeight);
		    	cosThetaEta_GJ_tot[1]->Fill(cosTheta_eta_gjs[ientry], totWeight);
		    	cosThetaEta_GJ_bkg[1]->Fill(cosTheta_eta_gjs[ientry], bkgWeight);
		    	cosThetaX_CM_sig[1]->Fill(cosTheta_X_cms[ientry], sigWeight);
		    	cosThetaX_CM_tot[1]->Fill(cosTheta_X_cms[ientry], totWeight);
		    	cosThetaX_CM_bkg[1]->Fill(cosTheta_X_cms[ientry], bkgWeight);
		    	Meta_sig[1]->Fill(Metas[ientry],sigWeight);
		    	Meta_tot[1]->Fill(Metas[ientry],totWeight);
		    	Meta_bkg[1]->Fill(Metas[ientry],bkgWeight);
                    }
                    if ( isUniquePi0g1Bs[ientry] ) { 
		    	Mpi0g_sig->Fill(Mpi0g1s[ientry],sigWeight);
		    	Mpi0g_tot->Fill(Mpi0g1s[ientry],totWeight);
		    	Mpi0g_bkg->Fill(Mpi0g1s[ientry],bkgWeight);
		    }
                    if ( isUniquePi0g2Bs[ientry] ) { 
		    	Mpi0g_sig->Fill(Mpi0g2s[ientry],sigWeight);
		    	Mpi0g_tot->Fill(Mpi0g2s[ientry],totWeight);
		    	Mpi0g_bkg->Fill(Mpi0g2s[ientry],bkgWeight);
		    }
                    if ( isUniquePi0Bs[ientry] ) { 
		    	Mpi0_sig[0]->Fill(Mpi0s_meas[ientry],sigWeight);
		    	Mpi0_tot[0]->Fill(Mpi0s_meas[ientry],totWeight);
		    	Mpi0_bkg[0]->Fill(Mpi0s_meas[ientry], bkgWeight);

		    	Mpi0_sig[1]->Fill(Mpi0s[ientry],sigWeight);
		    	Mpi0_tot[1]->Fill(Mpi0s[ientry],totWeight);
		    	Mpi0_bkg[1]->Fill(Mpi0s[ientry], bkgWeight);
                    }
                    if ( isUniquePi0EtaBs[ientry] ) { 
		    	phiEta_GJ_sig[0]->Fill(phi_eta_gjs_meas[ientry], sigWeight);
		    	phiEta_GJ_tot[0]->Fill(phi_eta_gjs_meas[ientry], totWeight);
		    	phiEta_GJ_bkg[0]->Fill(phi_eta_gjs_meas[ientry], bkgWeight);
		    	Mpi0eta_sig[0]->Fill(Mpi0etas_meas[ientry],sigWeight);
		    	Mpi0eta_tot[0]->Fill(Mpi0etas_meas[ientry],totWeight);
		    	Mpi0eta_bkg[0]->Fill(Mpi0etas_meas[ientry], bkgWeight);

		    	phiEta_GJ_sig[1]->Fill(phi_eta_gjs[ientry], sigWeight);
		    	phiEta_GJ_tot[1]->Fill(phi_eta_gjs[ientry], totWeight);
		    	phiEta_GJ_bkg[1]->Fill(phi_eta_gjs[ientry], bkgWeight);
		    	Mpi0eta_sig[1]->Fill(Mpi0etas[ientry],sigWeight);
		    	Mpi0eta_tot[1]->Fill(Mpi0etas[ientry],totWeight);
		    	Mpi0eta_bkg[1]->Fill(Mpi0etas[ientry], bkgWeight);
                    }
	        }
                else if (s_uniquenessTracking=="weighted"){
		        sigWeight = sigWeight*utWeights[ientry];
		        totWeight = totWeight*utWeights[ientry]; 
		        bkgWeight = bkgWeight*utWeights[ientry];
		    	cosThetaEta_GJ_sig[0]->Fill(cosTheta_eta_gjs_meas[ientry], sigWeight);
		    	cosThetaEta_GJ_tot[0]->Fill(cosTheta_eta_gjs_meas[ientry], totWeight);
		    	cosThetaEta_GJ_bkg[0]->Fill(cosTheta_eta_gjs_meas[ientry], bkgWeight);
		    	cosThetaX_CM_sig[0]->Fill(cosTheta_X_cms_meas[ientry], sigWeight);
		    	cosThetaX_CM_tot[0]->Fill(cosTheta_X_cms_meas[ientry], totWeight);
		    	cosThetaX_CM_bkg[0]->Fill(cosTheta_X_cms_meas[ientry], bkgWeight);
		    	Meta_sig[0]->Fill(Metas_meas[ientry],sigWeight);
		    	Meta_tot[0]->Fill(Metas_meas[ientry],totWeight);
		    	Meta_bkg[0]->Fill(Metas_meas[ientry],bkgWeight);

		    	cosThetaEta_GJ_sig[1]->Fill(cosTheta_eta_gjs[ientry], sigWeight);
		    	cosThetaEta_GJ_tot[1]->Fill(cosTheta_eta_gjs[ientry], totWeight);
		    	cosThetaEta_GJ_bkg[1]->Fill(cosTheta_eta_gjs[ientry], bkgWeight);
		    	cosThetaX_CM_sig[1]->Fill(cosTheta_X_cms[ientry], sigWeight);
		    	cosThetaX_CM_tot[1]->Fill(cosTheta_X_cms[ientry], totWeight);
		    	cosThetaX_CM_bkg[1]->Fill(cosTheta_X_cms[ientry], bkgWeight);
		    	Meta_sig[1]->Fill(Metas[ientry],sigWeight);
		    	Meta_tot[1]->Fill(Metas[ientry],totWeight);
		    	Meta_bkg[1]->Fill(Metas[ientry],bkgWeight);

		    	Mpi0g_sig->Fill(Mpi0g1s[ientry],sigWeight);
		    	Mpi0g_tot->Fill(Mpi0g1s[ientry],totWeight);
		    	Mpi0g_bkg->Fill(Mpi0g1s[ientry],bkgWeight);

		    	Mpi0g_sig->Fill(Mpi0g2s[ientry],sigWeight);
		    	Mpi0g_tot->Fill(Mpi0g2s[ientry],totWeight);
		    	Mpi0g_bkg->Fill(Mpi0g2s[ientry],bkgWeight);

		    	Mpi0_sig[0]->Fill(Mpi0s_meas[ientry],sigWeight);
		    	Mpi0_tot[0]->Fill(Mpi0s_meas[ientry],totWeight);
		    	Mpi0_bkg[0]->Fill(Mpi0s_meas[ientry], bkgWeight);

		    	Mpi0_sig[1]->Fill(Mpi0s[ientry],sigWeight);
		    	Mpi0_tot[1]->Fill(Mpi0s[ientry],totWeight);
		    	Mpi0_bkg[1]->Fill(Mpi0s[ientry], bkgWeight);

		    	phiEta_GJ_sig[0]->Fill(phi_eta_gjs_meas[ientry], sigWeight);
		    	phiEta_GJ_tot[0]->Fill(phi_eta_gjs_meas[ientry], totWeight);
		    	phiEta_GJ_bkg[0]->Fill(phi_eta_gjs_meas[ientry], bkgWeight);
		    	Mpi0eta_sig[0]->Fill(Mpi0etas_meas[ientry],sigWeight);
		    	Mpi0eta_tot[0]->Fill(Mpi0etas_meas[ientry],totWeight);
		    	Mpi0eta_bkg[0]->Fill(Mpi0etas_meas[ientry], bkgWeight);

		    	phiEta_GJ_sig[1]->Fill(phi_eta_gjs[ientry], sigWeight);
		    	phiEta_GJ_tot[1]->Fill(phi_eta_gjs[ientry], totWeight);
		    	phiEta_GJ_bkg[1]->Fill(phi_eta_gjs[ientry], bkgWeight);
		    	Mpi0eta_sig[1]->Fill(Mpi0etas[ientry],sigWeight);
		    	Mpi0eta_tot[1]->Fill(Mpi0etas[ientry],totWeight);
		    	Mpi0eta_bkg[1]->Fill(Mpi0etas[ientry], bkgWeight);
                }
                else {
                    cout << "uniqueness tracking setting not valid!" << endl;
                    exit(0);
                }
        }
        cout << "Made the histograms" << endl;

	// HERE WE WILL JUST DRAW SOME OF THE HISTOGRAMS WITH THE BKG FILLED IN TO SEE THEIR CONTRIBUTION
	makeStackedHist(Mpi0g_tot,Mpi0g_sig,Mpi0g_bkg,"Mpi0gkin", "diagnosticPlots/"+fileTag);
	for (int i=0; i<2; i++){
		string tag="";
		if (i%2==0){ tag="meas"; }
		else { tag="kin"; } 
		makeStackedHist(Meta_tot[i],Meta_sig[i],Meta_bkg[i],"Meta"+tag, "diagnosticPlots/"+fileTag);
		makeStackedHist(Mpi0_tot[i],Mpi0_sig[i],Mpi0_bkg[i],"Mpi0"+tag, "diagnosticPlots/"+fileTag);
		makeStackedHist(Mpi0eta_tot[i],Mpi0eta_sig[i],Mpi0eta_bkg[i],"Mpi0eta"+tag, "diagnosticPlots/"+fileTag);
		makeStackedHist(cosThetaEta_GJ_tot[i],cosThetaEta_GJ_sig[i],cosThetaEta_GJ_bkg[i],"cosThetaEta_GJ"+tag, "diagnosticPlots/"+fileTag);	
		makeStackedHist(cosThetaX_CM_tot[i],cosThetaX_CM_sig[i],cosThetaX_CM_bkg[i],"cosThetaX_CM"+tag, "diagnosticPlots/"+fileTag);	
		makeStackedHist(phiEta_GJ_tot[i],phiEta_GJ_sig[i],phiEta_GJ_bkg[i],"phiEta_GJ"+tag, "diagnosticPlots/"+fileTag);	
	}

        cout << "FINIHSED!"<<endl;
	cout << "There were " << numNan << " nan values where q-value was not calculated. Fix me if nonzero!" << endl;

        // we can also directly save all the histograms to a root file
	TFile* dataFile3 = new TFile(("diagnosticPlots/"+fileTag+"/postQValHists_"+fileTag+".root").c_str(),"RECREATE");
	for (int i=0; i<2; i++){
        	Meta_bkg[i]->Write();
        	Meta_sig[i]->Write();
        	Meta_tot[i]->Write();
        	Mpi0_sig[i]->Write();
        	Mpi0_bkg[i]->Write();
        	Mpi0_tot[i]->Write();
        	Mpi0eta_sig[i]->Write();
        	Mpi0eta_bkg[i]->Write();
        	Mpi0eta_tot[i]->Write();
	}
	Mpi0g_sig->Write();
	Mpi0g_bkg->Write();
	Mpi0g_tot->Write();
        
}

