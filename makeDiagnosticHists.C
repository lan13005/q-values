#include <ctime>
bool verbose = true;

void makeDiagnosticHists(){

	// setting up some basic root stuff and getting the file and tree
	TFile* dataFile=new TFile("pi0eta_a0_recotreeFlat_DSelector.root");
	TTree *dataTree;
	dataFile->GetObject("pi0eta_a0_recotree_flat",dataTree);
	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	TH1F* dHist_chisq= new TH1F( "dHist_chisq", "#Chi^2; Events/1", 100, 0, 100);
        TH1F* dHist_std = new TH1F("dHist_std", "#sigma; Event/0.002", 60, 0, 0.12);
	TH1F* dHist_qvalues = new TH1F( "dHist_qvalues_", "Q-Values; Events/0.01", 100, -1, 1);

	TH1F* Meta_Qd = new TH1F( "Meta_Qd", "M(#eta) after Q-Value Weighting; M(#eta) GeV; Events/0.002 GeV", 300, 0.25, 0.85 );
	TH1F* Meta_bkg_Qd = new TH1F( "Meta_bkg_Qd", "M(#eta) after Q-Value Weighting; M(#eta) GeV; Events/0.002 GeV", 300, 0.25, 0.85 );
	TH1F* Mpi0_Qd = new TH1F( "Mpi0_Qd", "M(#pi_{0}) after Q-Value Weighting; M(#pi_{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
	TH1F* Mpi0_bkg_Qd = new TH1F( "Mpi0_bkg_Qd", "M(#pi_{0}) after Q-Value Weighting; M(#pi_{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
	TH1F* Mpi0eta_Qd = new TH1F( "Mpi0eta_Qd", "M(#pi_{0}#eta) after Q-Value Weighting; M(#pi_{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
	TH1F* Mpi0eta_bkg_Qd = new TH1F( "Mpi0eta_bkg_Qd", "M(#pi_{0}#eta) after Q-Value Weighting; M(#pi_{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
        
        TLine* etaLine;
	double Meta;
	double Mpi0;
	double Mpi0eta;
        double AccWeight;
	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("Mpi0",&Mpi0);
	dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
        dataTree->SetBranchAddress("AccWeight",&AccWeight);

        // continue to do uniqueness tracking for the histograms we care about
        bool isUniqueEtaB;
        bool isUniquePi0B;
        bool isUniquePi0EtaB;
        dataTree->SetBranchAddress("isNotRepeated_eta",&isUniqueEtaB);
        dataTree->SetBranchAddress("isNotRepeated_pi0",&isUniquePi0B);
        dataTree->SetBranchAddress("isNotRepeated_pi0eta",&isUniquePi0EtaB);
	Long64_t nentries=dataTree->GetEntries();

        cout << "Finished setting up" << endl;
	
	const int c_nentries = (const int)nentries;
        cout << "c_nentries: " << c_nentries << endl;
	std::vector<double> Metas(c_nentries,0);
	std::vector<double> Mpi0s(c_nentries,0);
	std::vector<double> Mpi0etas(c_nentries,0);
        std::vector<double> AccWeights(c_nentries,0);

        std::vector<bool> isUniqueEtaBs(c_nentries,0);
        std::vector<bool> isUniquePi0Bs(c_nentries,0);
        std::vector<bool> isUniquePi0EtaBs(c_nentries,0);

	for (int ientry=0; ientry<nentries; ientry++)
	{
		dataTree->GetEntry(ientry);
		Metas[ientry]=Meta;
		Mpi0s[ientry]=Mpi0;
		Mpi0etas[ientry]=Mpi0eta;
                AccWeights[ientry]=AccWeight;

                isUniqueEtaBs[ientry]=isUniqueEtaB;
                isUniquePi0Bs[ientry]=isUniquePi0B;
                isUniquePi0EtaBs[ientry]=isUniquePi0EtaB;
	}
        
        cout << "Imported all the tree data into arrays" << endl;
	//cout << "\033[1;31m THE ANALYSIS DEPENDS ON THE CONSISTENCY OF READING A TTREE FROM DIFFERENT FILES. IN ONE FILE WE READ IN THE DATA AND SAVED THE ENTRY NUMBER RAN THIS PROGRAM IN PARALLEL SO THE WAY WE MERGE THE ROOT FILES MIGHT HAVE SCRAMBLED SOME BLOCKS OF DATA. IF WE USE THE ENTRY NUMBER IT SHOULD GIVE US CONSISTENT DATA IF WE COMPARE TO THE DATA WE JUST IMPORTED IN THE PREVIOUS STEP. \033[0m\n";

	TFile* dataFile2 = new TFile("qvalResults.root");
        TTree* dataTree2;
	dataFile2->GetObject("resultsTree",dataTree2);
        cout << "Loaded the results file" << endl;

	Double_t qvalue;
	Double_t conjugate_qvalue;
	Double_t chisq;
        Double_t combostd;
	ULong64_t flatEntryNumber;

	dataTree2->SetBranchAddress("qvalue",&qvalue);
	dataTree2->SetBranchAddress("chisq",&chisq);
	dataTree2->SetBranchAddress("combostd",&combostd);
	dataTree2->SetBranchAddress("flatEntryNumber",&flatEntryNumber);

	nentries=dataTree2->GetEntries();

	const int c_nentries2 = (const int)nentries;
        cout << "c_nentries2: " << c_nentries2 << endl;

        if (c_nentries!=c_nentries2){
            cout << "\033[1;31m THE IMPORTED SIZES ARE NOT THE SAME .... EXITING... \033[0m\n"; 
            cout << "size of dataFile, resultsFile: " << c_nentries << "," << c_nentries2 << endl;
            exit(0);
        }

        set<ULong64_t> all_ientries;
        set<ULong64_t> all_flatEntryNumber;
        cout << "built comparision sets..." << endl;

	std::vector< double > qvalues(c_nentries,0);
	std::vector< double > chisqs(c_nentries,0);
	std::vector< double > combostds(c_nentries,0);
        //ULong64_t flatEntryNumbers[c_nentries];
        std::vector<ULong64_t> flatEntryNumbers(c_nentries,0);
	for (int ientry=0; ientry<nentries; ientry++)
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

        cout << "Finished consistency check on results file and data" << endl;

        allCanvases->Clear();
        dHist_std->Draw();
        allCanvases->SaveAs("combostds.png");
	THStack* stackedHists = new THStack("stackedHists","");
        //allCanvases->Clear();
        //dHist_qvalues_bkgRegion->SetFillColorAlpha(kBlue,0.3);
        //dHist_qvalues_sigRegion->SetFillColorAlpha(kMagenta,0.3);
        //stackedHists->Add(dHist_qvalues_bkgRegion,"HIST");
        //stackedHists->Add(dHist_qvalues_sigRegion,"HIST");
        //stackedHists->Draw();
	//stackedHists->GetXaxis()->SetTitle(dHist_qvalues_bkgRegion->GetXaxis()->GetTitle());
	//stackedHists->SetTitle("Q-Values");
	dHist_qvalues->Draw();
        allCanvases->SaveAs("qvalues.png");
        allCanvases->Clear();
	//stackedHists = new THStack("stackedHists","");
        //dHist_chisq_bkgRegion->SetFillColorAlpha(kBlue,0.3);
        //dHist_chisq_sigRegion->SetFillColorAlpha(kMagenta,0.3);
        //stackedHists->Add(dHist_chisq_bkgRegion,"HIST");
        //stackedHists->Add(dHist_chisq_sigRegion,"HIST");
        //stackedHists->Draw();
	//stackedHists->GetXaxis()->SetTitle(dHist_chisq_bkgRegion->GetXaxis()->GetTitle());
	//stackedHists->SetTitle("#Chi^2");
	dHist_chisq->Draw();
        allCanvases->SaveAs("chisq.png");
    
        cout << "Drew the imported data" << endl;

	for (int ientry=0; ientry<nentries; ientry++){
		qvalue = qvalues[ientry];
		conjugate_qvalue = 1-qvalue;
                //cout << "ientry,qVal,conj_qVal: " << ientry << "," << qvalue << "," << conjugate_qvalue << endl;
                if ( isUniqueEtaBs[ientry] ) {
		    Meta_Qd->Fill(Metas[ientry]);//,AccWeights[ientry]);
		    Meta_bkg_Qd->Fill(Metas[ientry], conjugate_qvalue);//*AccWeights[ientry]);
                }
                if ( isUniquePi0Bs[ientry] ) { 
		    Mpi0_Qd->Fill(Mpi0s[ientry]);//,AccWeights[ientry]);
		    Mpi0_bkg_Qd->Fill(Mpi0s[ientry], conjugate_qvalue);//*AccWeights[ientry]);
                }
                if ( isUniquePi0EtaBs[ientry] ) { 
		    Mpi0eta_Qd->Fill(Mpi0etas[ientry]);//,AccWeights[ientry]);
		    Mpi0eta_bkg_Qd->Fill(Mpi0etas[ientry], conjugate_qvalue);//*AccWeights[ientry]);
                }
		//Meta_Qd->Fill(Metas[ientry],AccWeights[ientry]);
		//Meta_bkg_Qd->Fill(Metas[ientry], conjugate_qvalue*AccWeights[ientry]);
		//Mpi0_Qd->Fill(Mpi0s[ientry],AccWeights[ientry]);
		//Mpi0_bkg_Qd->Fill(Mpi0s[ientry], conjugate_qvalue*AccWeights[ientry]);
		//Mpi0eta_Qd->Fill(Mpi0etas[ientry],AccWeights[ientry]);
		//Mpi0eta_bkg_Qd->Fill(Mpi0etas[ientry], conjugate_qvalue*AccWeights[ientry]);
	}
        cout << "Made the histograms" << endl;

	// HERE WE WILL JUST DRAW SOME OF THE HISTOGRAMS WITH THE BKG FILLED IN TO SEE THEIR CONTRIBUTION
	// ----------------- Meta 
	stackedHists = new THStack("stackedHists","");
	allCanvases->Clear();
	Meta_bkg_Qd->SetFillColorAlpha(kMagenta,0.5);
	Meta_bkg_Qd->SetLineColorAlpha(kMagenta,0);
	stackedHists->Add(Meta_Qd,"HIST");
	stackedHists->Add(Meta_bkg_Qd,"HIST");
	stackedHists->Draw("nostack");
	stackedHists->GetXaxis()->SetTitle(Meta_Qd->GetXaxis()->GetTitle());
	stackedHists->GetYaxis()->SetTitle(Meta_Qd->GetYaxis()->GetTitle());
	stackedHists->SetTitle(Meta_Qd->GetTitle());
	allCanvases->SaveAs("Meta_Qd.png");
	// ----------------- Mpi0 
	stackedHists = new THStack("stackedHists","");
	allCanvases->Clear();
	Mpi0_bkg_Qd->SetFillColorAlpha(kMagenta,0.5);
	Mpi0_bkg_Qd->SetLineColorAlpha(kMagenta,0);
	stackedHists->Add(Mpi0_Qd,"HIST");
	stackedHists->Add(Mpi0_bkg_Qd,"HIST");
	stackedHists->Draw("nostack");
	stackedHists->GetXaxis()->SetTitle(Mpi0_Qd->GetXaxis()->GetTitle());
	stackedHists->GetYaxis()->SetTitle(Mpi0_Qd->GetYaxis()->GetTitle());
	stackedHists->SetTitle(Mpi0_Qd->GetTitle());
	allCanvases->SaveAs("Mpi0_Qd.png");
	// ----------------- Mpi0eta 
	stackedHists = new THStack("stackedHists","");
	allCanvases->Clear();
	Mpi0eta_bkg_Qd->SetFillColorAlpha(kMagenta,0.5);
	Mpi0eta_bkg_Qd->SetLineColorAlpha(kMagenta,0);
	stackedHists->Add(Mpi0eta_Qd,"HIST");
	stackedHists->Add(Mpi0eta_bkg_Qd,"HIST");
	stackedHists->Draw("nostack");
	stackedHists->GetXaxis()->SetTitle(Mpi0eta_Qd->GetXaxis()->GetTitle());
	stackedHists->GetYaxis()->SetTitle(Mpi0eta_Qd->GetYaxis()->GetTitle());
	stackedHists->SetTitle(Mpi0eta_Qd->GetTitle());
	allCanvases->SaveAs("Mpi0eta_Qd.png");

        cout << "FINIHSED!"<<endl;


	TFile* dataFile3 = new TFile("postQValHists.root","RECREATE");
        Meta_bkg_Qd->Write();
        Meta_Qd->Write();
        Mpi0_Qd->Write();
        Mpi0_bkg_Qd->Write();
        Mpi0eta_Qd->Write();
        Mpi0eta_bkg_Qd->Write();
        







}

