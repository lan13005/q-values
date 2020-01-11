#include <ctime>
#include <math.h> 
#include "makeDiagnosticHists.h"
bool verbose = true;

void makeDiagnosticHists(Long64_t nentries, string tag){

	// setting up some basic root stuff and getting the file and tree
	TFile* dataFile=new TFile("pi0eta_fcal_treeFlat_DSelector.root");
	//TFile* dataFile=new TFile("pi0eta_a0_recotreeFlat_DSelector.root");
	TTree *dataTree;
	dataFile->GetObject("pi0eta_fcaltree_flat",dataTree);
	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	TH1F* dHist_chisq= new TH1F( "dHist_chisq", "#Chi^2; Events/0.1", 50, 0, 5);
	TH1F* dHist_chisq01= new TH1F( "dHist_chisq01", "#Chi^2; Events/0.1", 50, 0, 5);
	TH1F* dHist_chisq02= new TH1F( "dHist_chisq02", "#Chi^2; Events/0.1", 50, 0, 5);

        TH1F* dHist_std = new TH1F("dHist_std", "#sigma; Event/0.0033", 60, 0, 0.2);
	TH1F* dHist_qvalues = new TH1F( "dHist_qvalues", "Q-Values; Events/0.01", 100, -1, 1);


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

	for (int i=0; i<2; i++){
		string tag="";
		if (i==0){ tag="meas"; }
		else { tag="kin"; } 
		Meta_tot[i] = new TH1F( ("Meta_tot_"+tag).c_str(), "M(#eta) no Q-Value Weighting; M(#eta) GeV; Events/0.002 GeV", 300, 0.25, 0.85 );
		Meta_sig[i] = new TH1F( ("Meta_sig_"+tag).c_str(), "M(#eta) after Q-Value Weighting; M(#eta) GeV; Events/0.002 GeV", 300, 0.25, 0.85 );
		Meta_bkg[i] = new TH1F( ("Meta_bkg_"+tag).c_str(), "M(#eta) after Q-Value Weighting; M(#eta) GeV; Events/0.002 GeV", 300, 0.25, 0.85 );
		Mpi0_tot[i] = new TH1F( ("Mpi0_tot_"+tag).c_str(), "M(#pi_{0}) no Q-Value Weighting; M(#pi_{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
		Mpi0_sig[i] = new TH1F( ("Mpi0_sig_"+tag).c_str(), "M(#pi_{0}) after Q-Value Weighting; M(#pi_{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
		Mpi0_bkg[i] = new TH1F( ("Mpi0_bkg_"+tag).c_str(), "M(#pi_{0}) after Q-Value Weighting; M(#pi_{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
		Mpi0eta_tot[i] = new TH1F( ("Mpi0eta_tot_"+tag).c_str(), "M(#pi_{0}#eta) no Q-Value Weighting; M(#pi_{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
		Mpi0eta_sig[i] = new TH1F( ("Mpi0eta_sig_"+tag).c_str(), "M(#pi_{0}#eta) after Q-Value Weighting; M(#pi_{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
		Mpi0eta_bkg[i] = new TH1F( ("Mpi0eta_bkg_"+tag).c_str(), "M(#pi_{0}#eta) after Q-Value Weighting; M(#pi_{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
		cosThetaEta_GJ_tot[i] = new TH1F( ("cosThetaEta_GJ_tot_"+tag).c_str(), "cos(#theta) GJ of #eta no Q-Value Weighting; cos(#theta) of #eta; Events/0.02 GeV", 100, -1, 1 );
		cosThetaEta_GJ_sig[i] = new TH1F( ("cosThetaEta_GJ_sig_"+tag).c_str(), "cos(#theta) GJ of #eta after Q-Value Weighting; cos(#theta) of #eta; Events/0.02 GeV", 100, -1, 1 );
		cosThetaEta_GJ_bkg[i] = new TH1F( ("cosThetaEta_GJ_bkg_"+tag).c_str(), "cos(#theta) GJ of #eta after Q-Value Weighting; cos(#theta) of #eta; Events/0.02 GeV", 100, -1, 1 );
		cosThetaX_CM_tot[i] = new TH1F( ("cosThetaX_CM_tot_"+tag).c_str(), "cos(#theta) of CM #eta no Q-Value Weighting; cos(#theta) of #eta; Events/0.02 GeV", 100, 0.9, 1 );
		cosThetaX_CM_sig[i] = new TH1F( ("cosThetaX_CM_sig_"+tag).c_str(), "cos(#theta) of CM #eta after Q-Value Weighting; cos(#theta) of #eta; Events/0.02 GeV", 100, 0.9, 1 );
		cosThetaX_CM_bkg[i] = new TH1F( ("cosThetaX_CM_bkg_"+tag).c_str(), "cos(#theta) of CM #eta after Q-Value Weighting; cos(#theta) of #eta; Events/0.02 GeV", 100, 0.9, 1 );
		phiEta_GJ_tot[i] = new TH1F( ("phiEta_GJ_tot_"+tag).c_str(), "#phi of #eta no Q-Value Weighting;#phi of #eta ;Events/0.02 GeV", 90, -180, 180 );
		phiEta_GJ_sig[i] = new TH1F( ("phiEta_GJ_sig_"+tag).c_str(), "#phi of #eta after Q-Value Weighting;#phi of #eta ;Events/0.02 GeV", 90, -180, 180 );
		phiEta_GJ_bkg[i] = new TH1F( ("phiEta_GJ_bkg_"+tag).c_str(), "#phi of #eta after Q-Value Weighting;#phi of #eta ;Events/0.02 GeV", 90, -180, 180 );
	}

        TLine* etaLine;
        double AccWeight;
	double Meta;
	double Mpi0;
	double Mpi0eta;
	double cosTheta_X_cm;
	double cosTheta_eta_gj;
	double phi_eta_gj;
	double Meta_kin;
	double Mpi0_kin;
	double Mpi0eta_kin;
	double cosTheta_X_cm_kin;
	double cosTheta_eta_gj_kin;
	double phi_eta_gj_kin;

        dataTree->SetBranchAddress("AccWeight",&AccWeight);
	dataTree->SetBranchAddress("Meta_meas",&Meta);
	dataTree->SetBranchAddress("Mpi0_meas",&Mpi0);
	dataTree->SetBranchAddress("Mpi0eta_meas",&Mpi0eta);
        dataTree->SetBranchAddress("cosTheta_X_cm_meas", &cosTheta_X_cm); 
        dataTree->SetBranchAddress("cosTheta_eta_gj_meas",&cosTheta_eta_gj);
        dataTree->SetBranchAddress("phi_eta_gj_meas",&phi_eta_gj); 
	dataTree->SetBranchAddress("Meta",&Meta_kin);
	dataTree->SetBranchAddress("Mpi0",&Mpi0_kin);
	dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta_kin);
        dataTree->SetBranchAddress("cosTheta_X_cm", &cosTheta_X_cm_kin); 
        dataTree->SetBranchAddress("cosTheta_eta_gj",&cosTheta_eta_gj_kin);
        dataTree->SetBranchAddress("phi_eta_gj",&phi_eta_gj_kin); 

        // continue to do uniqueness tracking for the histograms we care about
        bool isUniqueEtaB;
        bool isUniquePi0B;
        bool isUniquePi0EtaB;
        dataTree->SetBranchAddress("isNotRepeated_eta",&isUniqueEtaB);
        dataTree->SetBranchAddress("isNotRepeated_pi0",&isUniquePi0B);
        dataTree->SetBranchAddress("isNotRepeated_pi0eta",&isUniquePi0EtaB);

        cout << "Finished setting up" << endl;
	
	const int c_nentries = (const int)nentries;
        cout << "c_nentries: " << c_nentries << endl;

        std::vector<double> AccWeights; AccWeights.reserve(c_nentries); 
        std::vector<bool> isUniqueEtaBs; isUniqueEtaBs.reserve(c_nentries); 
        std::vector<bool> isUniquePi0Bs; isUniquePi0Bs.reserve(c_nentries); 
        std::vector<bool> isUniquePi0EtaBs; isUniquePi0EtaBs.reserve(c_nentries); 
        std::vector<double> Metas; Metas.reserve(c_nentries); 
        std::vector<double> Mpi0s; Mpi0s.reserve(c_nentries); 
        std::vector<double> Mpi0etas; Mpi0etas.reserve(c_nentries); 
	std::vector<double> cosTheta_X_cms; cosTheta_X_cms.reserve(c_nentries);
	std::vector<double> cosTheta_eta_gjs; cosTheta_eta_gjs.reserve(c_nentries);
	std::vector<double> phi_eta_gjs; phi_eta_gjs.reserve(c_nentries);
        std::vector<double> Metas_kin; Metas_kin.reserve(c_nentries); 
        std::vector<double> Mpi0s_kin; Mpi0s_kin.reserve(c_nentries); 
        std::vector<double> Mpi0etas_kin; Mpi0etas_kin.reserve(c_nentries); 
	std::vector<double> cosTheta_X_cms_kin; cosTheta_X_cms_kin.reserve(c_nentries);
	std::vector<double> cosTheta_eta_gjs_kin; cosTheta_eta_gjs_kin.reserve(c_nentries);
	std::vector<double> phi_eta_gjs_kin; phi_eta_gjs_kin.reserve(c_nentries);

	for (int ientry=0; ientry<nentries; ientry++)
	{
		dataTree->GetEntry(ientry);
                AccWeights.push_back(AccWeight);
		Metas.push_back(Meta);
		Mpi0s.push_back(Mpi0);
		Mpi0etas.push_back(Mpi0eta);
		cosTheta_X_cms.push_back(cosTheta_X_cm);
		cosTheta_eta_gjs.push_back(cosTheta_eta_gj);
		phi_eta_gjs.push_back(phi_eta_gj);
		Metas_kin.push_back(Meta_kin);
		Mpi0s_kin.push_back(Mpi0_kin);
		Mpi0etas_kin.push_back(Mpi0eta_kin);
		cosTheta_X_cms_kin.push_back(cosTheta_X_cm_kin);
		cosTheta_eta_gjs_kin.push_back(cosTheta_eta_gj_kin);
		phi_eta_gjs_kin.push_back(phi_eta_gj_kin);

                isUniqueEtaBs.push_back(isUniqueEtaB);
                isUniquePi0Bs.push_back(isUniquePi0B);
                isUniquePi0EtaBs.push_back(isUniquePi0EtaB);
                
	}
        cout << "Imported all the tree data into arrays" << endl;
        nentries = AccWeights.size();
        cout << "Size of array: " << nentries << endl;
	//cout << "\033[1;31m THE ANALYSIS DEPENDS ON THE CONSISTENCY OF READING A TTREE FROM DIFFERENT FILES. IN ONE FILE WE READ IN THE DATA AND SAVED THE ENTRY NUMBER RAN THIS PROGRAM IN PARALLEL SO THE WAY WE MERGE THE ROOT FILES MIGHT HAVE SCRAMBLED SOME BLOCKS OF DATA. IF WE USE THE ENTRY NUMBER IT SHOULD GIVE US CONSISTENT DATA IF WE COMPARE TO THE DATA WE JUST IMPORTED IN THE PREVIOUS STEP. \033[0m\n";
	TFile* dataFile2 = new TFile("qvalResults.root");
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


	nentries=dataTree2->GetEntries();
	const int c_nentries2 = (const int)nentries;
        cout << "c_nentries2: " << c_nentries2 << endl;

        if (nentries!=c_nentries2){
            cout << "\033[1;31m THE IMPORTED SIZES ARE NOT THE SAME .... EXITING... \033[0m\n"; 
            cout << "size of dataFile, resultsFile: " << nentries << "," << c_nentries2 << endl;
            exit(0);
        }

        set<ULong64_t> all_ientries;
        set<ULong64_t> all_flatEntryNumber;
        cout << "built comparision sets..." << endl;

	std::vector< double > qvalues(c_nentries,0);
	std::vector< double > chisqs(c_nentries,0);
	std::vector< double > combostds(c_nentries,0);
        std::vector<ULong64_t> flatEntryNumbers(c_nentries,0);
	for (int ientry=0; ientry<c_nentries; ientry++)
	{
        	dataTree2->GetEntry(ientry);
        	//cout << "ientry, flatEntryNumber: " << ientry << "," << flatEntryNumber << endl;            
        	
        	all_ientries.insert(ientry);
        	all_flatEntryNumber.insert(flatEntryNumber);
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

        cout << "Loaded all the results file" << endl;


        // ***************** CHECK 1 ********************
        for (int iEntry=0; iEntry<nentries; ++iEntry){
            if ( flatEntryNumbers[iEntry] != iEntry){
                cout << "flatEntryNumbers[iEntry] != iEntry: flatEntryNumbers[iEntry], iEntry: " << flatEntryNumbers[iEntry] << ", " << iEntry << endl;
            }
        }

        // **************** CHECK 2 ********************
        std::set<int> difference;
        std::set_difference(all_ientries.begin(), all_ientries.end(), all_flatEntryNumber.begin(), all_flatEntryNumber.end(),
            std::inserter(difference, difference.end()));
        if ( difference.size() != 0 ) {
            cout << "These elements are not common in both sets! Fix me!\n\t\t";
            for (auto it=difference.begin(); it != difference.end(); ++it){ 
                cout << ' ' << *it; 
            }
            cout << endl;
        }
        else { cout << "The flatEntryNumber set is complete!" << endl; }

        cout << "Finished consistency check on results file and data" << endl;

	
        allCanvases->Clear();
	dHist_chisq01->Draw();
        allCanvases->SaveAs("diagnosticPlots/chisq01.png");
        allCanvases->Clear();
	dHist_chisq02->Draw();
        allCanvases->SaveAs("diagnosticPlots/chisq02.png");

        allCanvases->Clear();
        dHist_std->Draw();
        allCanvases->SaveAs("diagnosticPlots/combostds.png");
	THStack* stackedHists = new THStack("stackedHists","");
	dHist_qvalues->Draw();
        allCanvases->SaveAs("diagnosticPlots/qvalues.png");
        allCanvases->Clear();
	dHist_chisq->Draw();
        allCanvases->SaveAs("diagnosticPlots/chisq.png");
    
        cout << "Drew the imported data" << endl;

	for (int ientry=0; ientry<nentries; ientry++){
		qvalue = qvalues[ientry];
		conjugate_qvalue = 1-qvalue;
                if(qvalue != qvalue){
                    cout << "ientry is nan: " << ientry << endl;
                }
                //cout << "ientry,qVal,conj_qVal: " << ientry << "," << qvalue << "," << conjugate_qvalue << endl;
                if ( isUniqueEtaBs[ientry] ) {
			cosThetaEta_GJ_sig[0]->Fill(cosTheta_eta_gjs[ientry], qvalue*AccWeights[ientry]);
			cosThetaEta_GJ_tot[0]->Fill(cosTheta_eta_gjs[ientry], AccWeights[ientry]);
			cosThetaEta_GJ_bkg[0]->Fill(cosTheta_eta_gjs[ientry], conjugate_qvalue*AccWeights[ientry]);
			cosThetaX_CM_sig[0]->Fill(cosTheta_X_cms[ientry], qvalue*AccWeights[ientry]);
			cosThetaX_CM_tot[0]->Fill(cosTheta_X_cms[ientry], AccWeights[ientry]);
			cosThetaX_CM_bkg[0]->Fill(cosTheta_X_cms[ientry], conjugate_qvalue*AccWeights[ientry]);
			Meta_sig[0]->Fill(Metas[ientry],qvalue*AccWeights[ientry]);
			Meta_tot[0]->Fill(Metas[ientry],AccWeights[ientry]);
			Meta_bkg[0]->Fill(Metas[ientry],conjugate_qvalue*AccWeights[ientry]);

			cosThetaEta_GJ_sig[1]->Fill(cosTheta_eta_gjs_kin[ientry], qvalue*AccWeights[ientry]);
			cosThetaEta_GJ_tot[1]->Fill(cosTheta_eta_gjs_kin[ientry], AccWeights[ientry]);
			cosThetaEta_GJ_bkg[1]->Fill(cosTheta_eta_gjs_kin[ientry], conjugate_qvalue*AccWeights[ientry]);
			cosThetaX_CM_sig[1]->Fill(cosTheta_X_cms_kin[ientry], qvalue*AccWeights[ientry]);
			cosThetaX_CM_tot[1]->Fill(cosTheta_X_cms_kin[ientry], AccWeights[ientry]);
			cosThetaX_CM_bkg[1]->Fill(cosTheta_X_cms_kin[ientry], conjugate_qvalue*AccWeights[ientry]);
			Meta_sig[1]->Fill(Metas_kin[ientry],qvalue*AccWeights[ientry]);
			Meta_tot[1]->Fill(Metas_kin[ientry], AccWeights[ientry]);
			Meta_bkg[1]->Fill(Metas_kin[ientry],conjugate_qvalue*AccWeights[ientry]);
                }
                if ( isUniquePi0Bs[ientry] ) { 
			Mpi0_sig[0]->Fill(Mpi0s[ientry],qvalue*AccWeights[ientry]);
			Mpi0_tot[0]->Fill(Mpi0s[ientry],AccWeights[ientry]);
			Mpi0_bkg[0]->Fill(Mpi0s[ientry], conjugate_qvalue*AccWeights[ientry]);

			Mpi0_sig[1]->Fill(Mpi0s_kin[ientry],qvalue*AccWeights[ientry]);
			Mpi0_tot[1]->Fill(Mpi0s_kin[ientry],AccWeights[ientry]);
			Mpi0_bkg[1]->Fill(Mpi0s_kin[ientry], conjugate_qvalue*AccWeights[ientry]);
                }
                if ( isUniquePi0EtaBs[ientry] ) { 
			phiEta_GJ_sig[0]->Fill(phi_eta_gjs[ientry], qvalue*AccWeights[ientry]);
			phiEta_GJ_tot[0]->Fill(phi_eta_gjs[ientry], AccWeights[ientry]);
			phiEta_GJ_bkg[0]->Fill(phi_eta_gjs[ientry], conjugate_qvalue*AccWeights[ientry]);
			Mpi0eta_sig[0]->Fill(Mpi0etas[ientry],qvalue*AccWeights[ientry]);
			Mpi0eta_tot[0]->Fill(Mpi0etas[ientry],AccWeights[ientry]);
			Mpi0eta_bkg[0]->Fill(Mpi0etas[ientry], conjugate_qvalue*AccWeights[ientry]);

			phiEta_GJ_sig[1]->Fill(phi_eta_gjs_kin[ientry], qvalue*AccWeights[ientry]);
			phiEta_GJ_tot[1]->Fill(phi_eta_gjs_kin[ientry], AccWeights[ientry]);
			phiEta_GJ_bkg[1]->Fill(phi_eta_gjs_kin[ientry], conjugate_qvalue*AccWeights[ientry]);
			Mpi0eta_sig[1]->Fill(Mpi0etas_kin[ientry],qvalue*AccWeights[ientry]);
			Mpi0eta_tot[1]->Fill(Mpi0etas_kin[ientry],AccWeights[ientry]);
			Mpi0eta_bkg[1]->Fill(Mpi0etas_kin[ientry], conjugate_qvalue*AccWeights[ientry]);
                }
	}
        cout << "Made the histograms" << endl;

	// HERE WE WILL JUST DRAW SOME OF THE HISTOGRAMS WITH THE BKG FILLED IN TO SEE THEIR CONTRIBUTION
	for (int i=0; i<2; i++){
		string tag="";
		if (i==0){ tag="meas"; }
		else { tag="kin"; } 
		makeStackedHist(Meta_tot[i],Meta_sig[i],Meta_bkg[i],"Meta"+tag, "diagnosticPlots");
		makeStackedHist(Mpi0_tot[i],Mpi0_sig[i],Mpi0_bkg[i],"Mpi0"+tag, "diagnosticPlots");
		makeStackedHist(Mpi0eta_tot[i],Mpi0eta_sig[i],Mpi0eta_bkg[i],"Mpi0eta"+tag, "diagnosticPlots");
		makeStackedHist(cosThetaEta_GJ_tot[i],cosThetaEta_GJ_sig[i],cosThetaEta_GJ_bkg[i],"cosThetaEta_GJ"+tag, "diagnosticPlots");	
		makeStackedHist(cosThetaX_CM_tot[i],cosThetaX_CM_sig[i],cosThetaX_CM_bkg[i],"cosThetaX_CM"+tag, "diagnosticPlots");	
		makeStackedHist(phiEta_GJ_tot[i],phiEta_GJ_sig[i],phiEta_GJ_bkg[i],"phiEta_GJ"+tag, "diagnosticPlots");	
	}

        cout << "FINIHSED!"<<endl;

	TFile* dataFile3 = new TFile("postQValHists.root","RECREATE");
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
        
}

