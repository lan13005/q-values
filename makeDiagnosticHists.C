#include <ctime>
#include <math.h> 
#include "makeDiagnosticHists.h"
bool verbose = true;
string detector="split";

void makeDiagnosticHists(){
	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	
	TH1F* dHist_chisq= new TH1F( "dHist_chisq", "#Chi^2; Events/0.1", 50, 0, 5);
	TH1F* dHist_chisq01= new TH1F( "dHist_chisq01", "#Chi^2; Events/0.1", 50, 0, 5);
	TH1F* dHist_chisq02= new TH1F( "dHist_chisq02", "#Chi^2; Events/0.1", 50, 0, 5);

        TH1F* dHist_std = new TH1F("dHist_std", "#sigma; Event/0.0033", 60, 0, 0.2);
	TH1F* dHist_qvalues = new TH1F( "dHist_qvalues", "Q-Values; Events/0.01", 100, -1, 1);

	// ---------------------------------------------------------------------------
	// --------------------------------------------- Settting branch addresses for Q-Value Results
	// ---------------------------------------------------------------------------

	TFile* dataFile2 = new TFile(("qvalResults_"+detector+".root").c_str());
        TTree* dataTree2;
	dataFile2->GetObject("resultsTree",dataTree2);
        cout << "Loaded the results file" << endl;

	Double_t qvalue;
	Double_t conjugate_qvalue;
	Double_t chisq;
        Double_t combostd;
        Bool_t bool_MetaFlat;
	ULong64_t flatEntryNumber;
	Double_t chiSq_eta_01;
	Double_t chiSq_eta_02;

	dataTree2->SetBranchAddress("qvalue",&qvalue);
	dataTree2->SetBranchAddress("chisq_eta",&chisq);
	dataTree2->SetBranchAddress("chisq_eta_01", &chiSq_eta_01);
	dataTree2->SetBranchAddress("chisq_eta_02", &chiSq_eta_02);
	dataTree2->SetBranchAddress("combostd",&combostd);
	dataTree2->SetBranchAddress("flatEntryNumber",&flatEntryNumber);


	ULong64_t nentries;
	nentries=dataTree2->GetEntries();
	const int c_nentries2 = (const int)nentries;
        cout << "c_nentries2: " << c_nentries2 << endl;


	// make a comparison set to see if we have all the combos
        // set<ULong64_t> all_ientries;
        // set<ULong64_t> all_flatEntryNumber;
        // cout << "built comparision sets..." << endl;

	// make a vector of 0s such that we can fill the q-values in order to unscramble to multiprocessing effects
	std::vector< double > qvalues(c_nentries2,0);
	std::vector< double > chisqs(c_nentries2,0);
	std::vector< double > combostds(c_nentries2,0);
        std::vector<ULong64_t> flatEntryNumbers(c_nentries2,0);
	for (int ientry=0; ientry<c_nentries2; ientry++)
	{
        	dataTree2->GetEntry(ientry);
        	//cout << "ientry, flatEntryNumber: " << ientry << "," << flatEntryNumber << endl;            
        	
        	//all_ientries.insert(ientry);
        	//all_flatEntryNumber.insert(flatEntryNumber);
        	qvalues[flatEntryNumber]=qvalue;
        	chisqs[flatEntryNumber]=chisq;
        	combostds[flatEntryNumber]=combostd;
        	flatEntryNumbers[flatEntryNumber] = flatEntryNumber;
        	
        	dHist_qvalues->Fill(qvalue);
        	dHist_chisq->Fill(chisq);
		dHist_chisq01->Fill(chiSq_eta_01);
		dHist_chisq02->Fill(chiSq_eta_02);
        	dHist_std->Fill(combostd);
        }
	dataFile2->Close();
        cout << "Loaded all the results file" << endl;


	// ---------------------------------------------------------------------------
	// ----------------------------------------------------- Define the histograms 
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

	// ---------------------------------------------------------------------------
	// ----------------------------------------------------- Set up branch address of the data file 
	// But fill it after we load the q-value data so we can write the data to the tree
	// I have a felling it might be bad to fill the q-values as we calculate it during the main program
	// since we could have 36+ processes all calculating and writing to the same branch
	// ---------------------------------------------------------------------------
	//
	// setting up some basic root stuff and getting the file and tree
	TFile* dataFile=new TFile(("pi0eta_"+detector+"_treeFlat_DSelector.root").c_str());
	//TFile* dataFile=new TFile("pi0eta_a0_recotreeFlat_DSelector.root");
	TTree *dataTree;
	dataFile->GetObject(("pi0eta_"+detector+"_tree_flat").c_str(),dataTree);
	
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
	double Meta2;
	double Mpi02;
	double Mpi0eta2;
	double cosTheta_X_cm2;
	double cosTheta_eta_gj2;
	double phi_eta_gj2;

        dataTree->SetBranchAddress("AccWeight",&AccWeight);
	dataTree->SetBranchAddress("Meta_meas",&Meta2);
	dataTree->SetBranchAddress("Mpi0_meas",&Mpi02);
	dataTree->SetBranchAddress("Mpi0eta_meas",&Mpi0eta2);
        dataTree->SetBranchAddress("cosTheta_X_cm_meas", &cosTheta_X_cm2); 
        dataTree->SetBranchAddress("cosTheta_eta_gj_meas",&cosTheta_eta_gj2);
        dataTree->SetBranchAddress("phi_eta_gj_meas",&phi_eta_gj2); 
	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("Mpi0",&Mpi0);
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

        cout << "Finished setting up" << endl;
	
	nentries=dataTree->GetEntries();
	const int c_nentries = (const int)nentries;
        cout << "c_nentries: " << c_nentries << endl;
	if ( c_nentries != c_nentries2 ){
		cout << "Entries are not equal between dataFile and resultsFile. Ran over a subset" << endl;
	}
	else{
		cout << "Entries are equal. Ran over the full dataset" << endl;
	}

        std::vector<double> AccWeights; AccWeights.reserve(c_nentries2); 
        std::vector<bool> isUniqueEtaBs; isUniqueEtaBs.reserve(c_nentries2); 
        std::vector<bool> isUniquePi0Bs; isUniquePi0Bs.reserve(c_nentries2); 
        std::vector<bool> isUniquePi0g1Bs; isUniquePi0g1Bs.reserve(c_nentries2); 
        std::vector<bool> isUniquePi0g2Bs; isUniquePi0g2Bs.reserve(c_nentries2); 
        std::vector<bool> isUniquePi0EtaBs; isUniquePi0EtaBs.reserve(c_nentries2); 
        std::vector<double> Metas2; Metas2.reserve(c_nentries2); 
        std::vector<double> Mpi0s2; Mpi0s2.reserve(c_nentries2); 
        std::vector<double> Mpi0etas2; Mpi0etas2.reserve(c_nentries2); 
	std::vector<double> cosTheta_X_cms2; cosTheta_X_cms2.reserve(c_nentries2);
	std::vector<double> cosTheta_eta_gjs2; cosTheta_eta_gjs2.reserve(c_nentries2);
	std::vector<double> phi_eta_gjs2; phi_eta_gjs2.reserve(c_nentries2);
        std::vector<double> Metas; Metas.reserve(c_nentries2); 
        std::vector<double> Mpi0s; Mpi0s.reserve(c_nentries2); 
        std::vector<double> Mpi0g1s; Mpi0g1s.reserve(c_nentries2); 
        std::vector<double> Mpi0g2s; Mpi0g2s.reserve(c_nentries2); 
        std::vector<double> Mpi0etas; Mpi0etas.reserve(c_nentries2); 
	std::vector<double> cosTheta_X_cms; cosTheta_X_cms.reserve(c_nentries2);
	std::vector<double> cosTheta_eta_gjs; cosTheta_eta_gjs.reserve(c_nentries2);
	std::vector<double> phi_eta_gjs; phi_eta_gjs.reserve(c_nentries2);

	std::vector< double > sbWeights; sbWeights.reserve(c_nentries2);
	double sbWeight;
	double sbRL = 0.09; // Right sideband left line
	double sbRR = 0.105; // Right sideband right line
	double sigL = 0.12;
	double sigR = 0.15;
	double sbLL = 0.165;
	double sbLR = 0.18;

        // ***************** CHECK 1 ********************
        for (int iEntry=0; iEntry<c_nentries2; ++iEntry){
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
	// ---------------------------------------------------------------------------

	TFile *qd_dataFile = TFile::Open(("postQ_"+detector+"_flatTree.root").c_str(),"RECREATE"); 
	TTree *outputTree = dataTree->CloneTree(-1,"fast"); 
	TBranch* b_qvalue;
	TBranch* b_qvalue_chisq;
	TBranch* b_flatEntryNumber;
	b_qvalue = outputTree->Branch("qvalue",&qvalue,"qvalue/D");
	b_qvalue_chisq = outputTree->Branch("qvalue_chisq",&chisq,"qvalue_chisq/D");
	b_flatEntryNumber = outputTree->Branch("flatEntryNumber",&flatEntryNumber,"flatEntryNumber/l");

	for (int ientry=0; ientry<c_nentries2; ientry++)
	{
		outputTree->GetEntry(ientry);
                AccWeights.push_back(AccWeight);

		if ( Mpi0 > sbRL && Mpi0 < sbRR ) { sbWeight = -1; } 
		else if ( Mpi0 > sbLL && Mpi0 < sbLR ) { sbWeight = -1; } 
		else if ( Mpi0 > sigL && Mpi0 < sigR ) { sbWeight = 1; } 
		else { sbWeight = 0; }
		sbWeights.push_back(sbWeight);

		Metas2.push_back(Meta2);
		Mpi0s2.push_back(Mpi02);
		Mpi0etas2.push_back(Mpi0eta2);
		cosTheta_X_cms2.push_back(cosTheta_X_cm2);
		cosTheta_eta_gjs2.push_back(cosTheta_eta_gj2);
		phi_eta_gjs2.push_back(phi_eta_gj2);
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
                
		// will take the opportunity to save the now-ordered qvalues into the root file as we load the rest of the data
		qvalue = qvalues[ientry];
		chisq = chisqs[ientry];
		flatEntryNumber = flatEntryNumbers[ientry];
		b_qvalue->Fill();
		b_qvalue_chisq->Fill();
		b_flatEntryNumber->Fill();
	}

        cout << "Imported all the tree data into arrays" << endl;
        cout << "Size of array: " << c_nentries << endl;

	qd_dataFile->cd();
	outputTree->Write(); 

	
        allCanvases->Clear();
	dHist_chisq01->Draw();
        allCanvases->SaveAs(("diagnosticPlots/"+detector+"/chisq01.png").c_str());
        allCanvases->Clear();
	dHist_chisq02->Draw();
        allCanvases->SaveAs(("diagnosticPlots/"+detector+"/chisq02.png").c_str());

        allCanvases->Clear();
        dHist_std->Draw();
        allCanvases->SaveAs(("diagnosticPlots/"+detector+"/combostds.png").c_str());
	THStack* stackedHists = new THStack("stackedHists","");
	dHist_qvalues->Draw();
        allCanvases->SaveAs(("diagnosticPlots/"+detector+"/qvalues.png").c_str());
        allCanvases->Clear();
	dHist_chisq->Draw();
        allCanvases->SaveAs(("diagnosticPlots/"+detector+"/chisq.png").c_str());
    
        cout << "Drew the imported data" << endl;

	double sigWeightNoSB;
	double totWeightNoSB;
	double bkgWeightNoSB;
	double sigWeight;
	double totWeight;
	double bkgWeight;
	int numNan=0;
	for (int ientry=0; ientry<c_nentries2; ientry++){
		qvalue = qvalues[ientry];
		conjugate_qvalue = 1-qvalue;
                if(qvalue != qvalue){
                    //cout << "ientry is nan: " << ientry << endl;
		    ++numNan;
                }
		sigWeight = qvalue*AccWeights[ientry];//*sbWeights[ientry];
		totWeight = AccWeights[ientry];//*sbWeights[ientry];
		bkgWeight = conjugate_qvalue*AccWeights[ientry];//*sbWeights[ientry];
                //cout << "ientry,qVal,conj_qVal: " << ientry << "," << qvalue << "," << conjugate_qvalue << endl;
                if ( isUniqueEtaBs[ientry] ) {
			cosThetaEta_GJ_sig[0]->Fill(cosTheta_eta_gjs2[ientry], sigWeight);
			cosThetaEta_GJ_tot[0]->Fill(cosTheta_eta_gjs2[ientry], totWeight);
			cosThetaEta_GJ_bkg[0]->Fill(cosTheta_eta_gjs2[ientry], bkgWeight);
			cosThetaX_CM_sig[0]->Fill(cosTheta_X_cms2[ientry], sigWeight);
			cosThetaX_CM_tot[0]->Fill(cosTheta_X_cms2[ientry], totWeight);
			cosThetaX_CM_bkg[0]->Fill(cosTheta_X_cms2[ientry], bkgWeight);
			Meta_sig[0]->Fill(Metas2[ientry],sigWeight);
			Meta_tot[0]->Fill(Metas2[ientry],totWeight);
			Meta_bkg[0]->Fill(Metas2[ientry],bkgWeight);

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
			Mpi0_sig[0]->Fill(Mpi0s2[ientry],sigWeight);
			Mpi0_tot[0]->Fill(Mpi0s2[ientry],totWeight);
			Mpi0_bkg[0]->Fill(Mpi0s2[ientry], bkgWeight);

			Mpi0_sig[1]->Fill(Mpi0s[ientry],sigWeight);
			Mpi0_tot[1]->Fill(Mpi0s[ientry],totWeight);
			Mpi0_bkg[1]->Fill(Mpi0s[ientry], bkgWeight);
                }
                if ( isUniquePi0EtaBs[ientry] ) { 
			phiEta_GJ_sig[0]->Fill(phi_eta_gjs2[ientry], sigWeight);
			phiEta_GJ_tot[0]->Fill(phi_eta_gjs2[ientry], totWeight);
			phiEta_GJ_bkg[0]->Fill(phi_eta_gjs2[ientry], bkgWeight);
			Mpi0eta_sig[0]->Fill(Mpi0etas2[ientry],sigWeight);
			Mpi0eta_tot[0]->Fill(Mpi0etas2[ientry],totWeight);
			Mpi0eta_bkg[0]->Fill(Mpi0etas2[ientry], bkgWeight);

			phiEta_GJ_sig[1]->Fill(phi_eta_gjs[ientry], sigWeight);
			phiEta_GJ_tot[1]->Fill(phi_eta_gjs[ientry], totWeight);
			phiEta_GJ_bkg[1]->Fill(phi_eta_gjs[ientry], bkgWeight);
			Mpi0eta_sig[1]->Fill(Mpi0etas[ientry],sigWeight);
			Mpi0eta_tot[1]->Fill(Mpi0etas[ientry],totWeight);
			Mpi0eta_bkg[1]->Fill(Mpi0etas[ientry], bkgWeight);
                }
	}
        cout << "Made the histograms" << endl;

	// HERE WE WILL JUST DRAW SOME OF THE HISTOGRAMS WITH THE BKG FILLED IN TO SEE THEIR CONTRIBUTION
	makeStackedHist(Mpi0g_tot,Mpi0g_sig,Mpi0g_bkg,"Mpi0gkin", "diagnosticPlots/"+detector);
	for (int i=0; i<2; i++){
		string tag="";
		if (i%2==0){ tag="meas"; }
		else { tag="kin"; } 
		makeStackedHist(Meta_tot[i],Meta_sig[i],Meta_bkg[i],"Meta"+tag, "diagnosticPlots/"+detector);
		makeStackedHist(Mpi0_tot[i],Mpi0_sig[i],Mpi0_bkg[i],"Mpi0"+tag, "diagnosticPlots/"+detector);
		makeStackedHist(Mpi0eta_tot[i],Mpi0eta_sig[i],Mpi0eta_bkg[i],"Mpi0eta"+tag, "diagnosticPlots/"+detector);
		makeStackedHist(cosThetaEta_GJ_tot[i],cosThetaEta_GJ_sig[i],cosThetaEta_GJ_bkg[i],"cosThetaEta_GJ"+tag, "diagnosticPlots/"+detector);	
		makeStackedHist(cosThetaX_CM_tot[i],cosThetaX_CM_sig[i],cosThetaX_CM_bkg[i],"cosThetaX_CM"+tag, "diagnosticPlots/"+detector);	
		makeStackedHist(phiEta_GJ_tot[i],phiEta_GJ_sig[i],phiEta_GJ_bkg[i],"phiEta_GJ"+tag, "diagnosticPlots/"+detector);	
	}

        cout << "FINIHSED!"<<endl;
	cout << "There were " << numNan << " nan values meaning q-value was not calculated. Fix me if nonzero!" << endl;

	TFile* dataFile3 = new TFile(("postQValHists_"+detector+".root").c_str(),"RECREATE");
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

