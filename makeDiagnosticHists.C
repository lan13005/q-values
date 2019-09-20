#include <ctime>
bool verbose = false;

void makeDiagnosticHists(){

	// setting up some basic root stuff and getting the file and tree
	TFile* dataFile=new TFile("/d/grid15/ln16/pi0eta/121818/z_pi0eta_a0_reco/pi0eta_a0_recotreeFlat_DSelector.root");
	TTree *dataTree;
	dataFile->GetObject("pi0eta_a0_recotree_flat",dataTree);
	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	TH1F* dHist_qvalues_bkgRegion = new TH1F( "dHist_qvalues_bkgRegion", "Q-Values; Events/0.01", 100, 0, 1);
	TH1F* dHist_qvalues_sigRegion = new TH1F( "dHist_qvalues_signRegion", "Q-Values; Events/0.01", 100, 0, 1);
	TH1F* dHist_chisq_bkgRegion = new TH1F( "dHist_chisq_bkgRegion", "#Chi^2; Events/1", 100, 0, 100);
	TH1F* dHist_chisq_sigRegion = new TH1F( "dHist_chisq_sigRegion", "#Chi^2; Events/1", 100, 0, 100);
        TH1F* dHist_std = new TH1F("dHist_std", "#sigma; Event/0.002", 60, 0, 0.12);

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
	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("Mpi0",&Mpi0);
	dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
	Long64_t nentries=dataTree->GetEntries();
	
	// ***** MAKE SURE THE NENTRIES MATCH THE ONES IN RUN_PARALLEL.C ******
	//string segment;
        //ifstream file ("run.sh");
        //while (getline(file, segment, '='))
        //{
        //    if (segment=="nentries"){
        //        getline(file, segment, '=');
        //        nentries=stoi(segment);
        //    }
        //}
        //cout << nentries << endl;
	//nentries=10000;

	const int c_nentries = (const int)nentries;
	double Metas[c_nentries];
	double Mpi0s[c_nentries];
	double Mpi0etas[c_nentries];
	for (int ientry=0; ientry<nentries; ientry++)
	{
		dataTree->GetEntry(ientry);
		Metas[ientry]=Meta;
		Mpi0s[ientry]=Mpi0;
		Mpi0etas[ientry]=Mpi0eta;
	}

	double qvalues[c_nentries];
	double chiSqs[c_nentries];
	double stds[c_nentries];
	double qvalue;
	double conjugate_qvalue;
	double chiSq;
        double std;
	int event;
	// we will not save Meta since we already have it
	// the order of the array[c_entries] are ordered by event which starts at 0 to nentries-1. If we just fill qvalues and chiSqs with "value" at array index "event" 
	// 	then it should be properly ordered.
	ifstream diagnosticlogFile ("diagnostic_logs.txt");
	if ( diagnosticlogFile.is_open() ) {
		while ( diagnosticlogFile >> event >> qvalue >> chiSq >> std) {
			if (verbose) { cout << event << " " << qvalue << " " << chiSq << " " << std << endl; }
                        if ( qvalue < 0 || qvalue > 1){
                            cout << event << " " << qvalue << " " << chiSq << " " << std << endl;
                        }
			chiSqs[event]=chiSq;
                        stds[event]=std;      
                        //if (std>0.09){ qvalue = 0; }
			qvalues[event]=qvalue;
                        if ( (0.115 < Mpi0s[event]) && (Mpi0s[event] < .16) ) { 
                            dHist_qvalues_sigRegion->Fill(qvalue);
                            dHist_chisq_sigRegion->Fill(chiSq);
                        }
                        else {
                            dHist_qvalues_bkgRegion->Fill(qvalue);
                            dHist_chisq_bkgRegion->Fill(chiSq);
                        }
                        dHist_std->Fill(std);
		}
	}
        allCanvases->Clear();
        dHist_std->Draw();
        allCanvases->SaveAs("stds.png");
	THStack* stackedHists = new THStack("stackedHists","");
        allCanvases->Clear();
        dHist_qvalues_bkgRegion->SetFillColorAlpha(kBlue,0.3);
        dHist_qvalues_sigRegion->SetFillColorAlpha(kMagenta,0.3);
        stackedHists->Add(dHist_qvalues_bkgRegion,"HIST");
        stackedHists->Add(dHist_qvalues_sigRegion,"HIST");
        stackedHists->Draw();
	stackedHists->GetXaxis()->SetTitle(dHist_qvalues_bkgRegion->GetXaxis()->GetTitle());
	stackedHists->SetTitle("Q-Values");
        allCanvases->SaveAs("qvalues.png");
        allCanvases->Clear();
	stackedHists = new THStack("stackedHists","");
        dHist_chisq_bkgRegion->SetFillColorAlpha(kBlue,0.3);
        dHist_chisq_sigRegion->SetFillColorAlpha(kMagenta,0.3);
        stackedHists->Add(dHist_chisq_bkgRegion,"HIST");
        stackedHists->Add(dHist_chisq_sigRegion,"HIST");
        stackedHists->Draw();
	stackedHists->GetXaxis()->SetTitle(dHist_chisq_bkgRegion->GetXaxis()->GetTitle());
	stackedHists->SetTitle("#Chi^2");
        allCanvases->SaveAs("chisq.png");
    

	for (int ientry=0; ientry<nentries; ientry++){
		qvalue = qvalues[ientry];
		conjugate_qvalue = 1-qvalue;
		Meta_Qd->Fill(Metas[ientry]);
		Meta_bkg_Qd->Fill(Metas[ientry], conjugate_qvalue);
		Mpi0_Qd->Fill(Mpi0s[ientry]);
		Mpi0_bkg_Qd->Fill(Mpi0s[ientry], conjugate_qvalue);
		Mpi0eta_Qd->Fill(Mpi0etas[ientry]);
		Mpi0eta_bkg_Qd->Fill(Mpi0etas[ientry], conjugate_qvalue);
	}

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
}

