#include "main.h"
#include "/d/grid15/ln16/pi0eta/092419/makeGraphs.h"

bool isEta2g=true;

void getInitParams_step1(){
	//TFile* dataFile=new TFile("pi0eta_a0_recotreeFlat_DSelector.root");

	gStyle->SetOptFit(111);
	gStyle->SetOptStat(0);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
	TTree *dataTree;

	string detectorNames[3] = {"fcal","bcal","split"};
	for (int i=0; i<3; ++i){
		if (isEta2g) {
			TFile* dataFile=new TFile(("pi0eta_"+detectorNames[i]+"_treeFlat_DSelector.root").c_str());
			dataFile->GetObject(("pi0eta_"+detectorNames[i]+"_tree_flat").c_str(),dataTree);
		}
		else {
			TFile* dataFile=new TFile("pi0eta_reco_3pi0treeFlat_DSelector.root");
			dataFile->GetObject("pi0eta_reco_3pi0tree_flat",dataTree);
		}
    		ofstream logFile_eta;
    		ofstream logFile_pi0;
    		logFile_eta.open(("fitResults/etaFitNoAccSub_"+detectorNames[i]+".txt").c_str());
    		logFile_pi0.open(("fitResults/pi0FitNoAccSub_"+detectorNames[i]+".txt").c_str());

    		TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
        	bool isUniqueEtaB;
        	bool isUniquePi0B;
        	bool isUniquePi0EtaB;
		double Meta;
		double Mpi0;
		double Mpi0eta;
		double AccWeight;
		dataTree->SetBranchAddress("Meta",&Meta);
		dataTree->SetBranchAddress("Mpi0",&Mpi0);
		dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
		dataTree->SetBranchAddress("AccWeight",&AccWeight);
        	dataTree->SetBranchAddress("isNotRepeated_eta",&isUniqueEtaB);
        	dataTree->SetBranchAddress("isNotRepeated_pi0",&isUniquePi0B);
        	dataTree->SetBranchAddress("isNotRepeated_pi0eta",&isUniquePi0EtaB);
		long long nentries=dataTree->GetEntries();

        	TH1F *massHistEta; 
        	TH1F *massHistPi0; 
        	TH1F *massHistPi0Eta; 
		TF1* fit;
        	TF1* bkgFit;
        	TF1* sigFit;
        	double par[7];
		
        	cout <<"Initialized" << endl;
        	string namePar[8] = {"nentries","const","linear","amp1","mass","sigma1","ampRatio","sigmaRatio"};

		allCanvases->Clear();
		std::vector<double> binRangeEta;
		std::vector<double> fitRangeEta;
		double binWidthEta;
		std::vector<double> binRangePi0;
		std::vector<double> fitRangePi0;
		double binWidthPi0;

		double integralBKG;
		double integralSIG;
		double weightedSignma;

		// ///////////////////////////////////////////////////
		// START ETA FIT
		// ///////////////////////////////////////////////////
		//// eta -> 2g
		//binRangeEta={100,0.25,0.8};
		//fitRangeEta={0.38,0.65};
		//binRangePi0={100,0.05,0.25};
		//fitRangePi0={0.1,0.17};
		// eta->3pi0
		binRangeEta={200,0.25,0.85};
		fitRangeEta={0.44,0.65};
		binRangePi0={200,0.05,0.25};
		fitRangePi0={0.1,0.17};

		binWidthEta=(binRangeEta[2]-binRangeEta[1])/binRangeEta[0];
		binWidthPi0=(binRangePi0[2]-binRangePi0[1])/binRangePi0[0];



		fit = new TF1("fit",fitFunc,fitRangeEta[0],fitRangeEta[1],numDOFbkg+numDOFsig);
		bkgFit = new TF1("bkgFit",background,fitRangeEta[0],fitRangeEta[1],numDOFbkg);
		sigFit = new TF1("sigFit",signalDG,fitRangeEta[0],fitRangeEta[1],numDOFsig);

		// eta->gg
		if (isEta2g) {
			fit->SetParameters(11500,250,1750,0.547,0.003,3,5);
		}
		// eta->3pi0
		else {
			fit->SetParameters(100,250,1750,0.547,0.003,1,1);
		}

		massHistEta = new TH1F("","",binRangeEta[0],binRangeEta[1],binRangeEta[2]);
		massHistPi0 = new TH1F("","",binRangePi0[0],binRangePi0[1],binRangePi0[2]);
		massHistPi0Eta = new TH1F("","", 350, 0, 3.5);
		cout << "Initialized for a specific mass (eta/pi0) fit" << endl;
		
		double sbRL = 0.09; // Right sideband left line
		double sbRR = 0.105; // Right sideband right line
		double sigL = 0.12;
		double sigR = 0.15;
		double sbLL = 0.165;
		double sbLR = 0.18;
		double sbWeight;
		for (int ientry=0; ientry<nentries; ientry++)
		{
			dataTree->GetEntry(ientry);
			if ( Mpi0 > sbRL && Mpi0 < sbRR ) { sbWeight = -1; } 
			else if ( Mpi0 > sbLL && Mpi0 < sbLR ) { sbWeight = -1; } 
			else if ( Mpi0 > sigL && Mpi0 < sigR ) { sbWeight = 1; } 
			else { sbWeight = 0; }
                	if ( isUniqueEtaB ) {
		        	massHistEta->Fill(Meta,AccWeight);//*sbWeight);
			}
                	if ( isUniquePi0B ) {
		        	massHistPi0->Fill(Mpi0,AccWeight); /////////////////////////////////////////// NOT WEIGHTED SINCE WE WONT BE ABLE TO FIT IT PROPERLY 
			}
                	if ( isUniquePi0EtaB ) {
				massHistPi0Eta->Fill(Mpi0eta,AccWeight);//*sbWeight);
			}
		}
		cout << "Filled all entries into histogram for a specific fit" << endl;

		massHistPi0Eta->Draw("HIST");
		allCanvases->SaveAs(("fitResults/Mpi0eta_fit_"+detectorNames[i]+".png").c_str());
		allCanvases->Clear();
		
		massHistEta->Fit("fit","RQB"); // B will enforce the bounds
		fit->GetParameters(par);
		bkgFit->SetParameters(par);
		sigFit->SetParameters(&par[numDOFbkg]);
		massHistEta->SetTitle(";M(#eta) (GeV)");
		logFile_eta << namePar[0] << " " << nentries << endl;
		
		for (int iPar=0; iPar<numDOFbkg+numDOFsig; ++iPar){
		        logFile_eta << namePar[iPar+1] << " " << par[iPar] << endl;
		}
		
		integralBKG = bkgFit->Integral(fitRangeEta[0],fitRangeEta[1]);
		integralSIG = sigFit->Integral(fitRangeEta[0],fitRangeEta[1]);
		cout << "IntegralBKG before scaling: " << integralBKG << endl;
		cout << "IntegralSIG before scaling: " << integralSIG << endl;
		integralBKG *= 1/binWidthEta;
		integralSIG *= 1/binWidthEta;

		cout << "IntegralBKG after scaling: " << integralBKG << endl;
		cout << "IntegralSIG after scaling: " << integralSIG << endl;
		cout << "nentries: " << nentries << endl;

		weightedSigma = 1.0/(1+par[5]/par[6])*par[4]+1.0/(1+par[6]/par[5])*par[6]*par[4];
		logFile_eta << "weightedSigma: " << weightedSigma << endl;
		
		massHistEta->Draw();
		massHistEta->SetTitle(("Peak: "+to_string(par[4])+"    width: "+to_string(par[5])).c_str());
		allCanvases->SaveAs(("fitResults/Meta_fit_"+detectorNames[i]+".png").c_str());
		cout << "Saved for a specific fit!" << endl;
        	

		// ///////////////////////////////////////////////////
		// START Pi0 FIT
		// ///////////////////////////////////////////////////
		fit = new TF1("fit",fitFunc,fitRangePi0[0],fitRangePi0[1],numDOFbkg+numDOFsig);
		bkgFit = new TF1("bkgFit",background,fitRangePi0[0],fitRangePi0[1],numDOFbkg);
		sigFit = new TF1("sigFit",signalDG,fitRangePi0[0],fitRangePi0[1],numDOFsig);

		// eta->gg
		if (isEta2g) {
			fit->SetParameters(11500,250,1750,0.134,0.001,5,2);
		}
		// eta->3pi0
		else { 
			fit->SetParameters(700,0,300,0.134,0.01,1,1);
		}

		massHistPi0->Fit("fit","RQB"); // B will enforce the bounds
		fit->GetParameters(par);
		bkgFit->SetParameters(par);
		sigFit->SetParameters(&par[numDOFbkg]);
		massHistPi0->SetTitle(";M(#pi^{0}) (GeV)");
		logFile_pi0 << namePar[0] << " " << nentries << endl;
		
		for (int iPar=0; iPar<numDOFbkg+numDOFsig; ++iPar){
		        logFile_pi0 << namePar[iPar+1] << " " << par[iPar] << endl;
		}
		
		integralBKG = bkgFit->Integral(fitRangePi0[0],fitRangePi0[1]);
		integralSIG = sigFit->Integral(fitRangePi0[0],fitRangePi0[1]);
		cout << "IntegralBKG before scaling: " << integralBKG << endl;
		cout << "IntegralSIG before scaling: " << integralSIG << endl;
		integralBKG *= 1/binWidthPi0;
		integralSIG *= 1/binWidthPi0;

		cout << "IntegralBKG after scaling: " << integralBKG << endl;
		cout << "IntegralSIG after scaling: " << integralSIG << endl;
		cout << "nentries: " << nentries << endl;

		weightedSigma = 1.0/(1+par[5]/par[6])*par[4]+1.0/(1+par[6]/par[5])*par[6]*par[4];
		logFile_pi0 << "weightedSigma: " << weightedSigma << endl;
		
		massHistPi0->Draw();
		massHistPi0->SetTitle(("Peak: "+to_string(par[4])+"    width: "+to_string(par[5])).c_str());
		TLine *line = new TLine( sbRL, 0, sbRL, massHistPi0->GetMaximum());
		line->SetLineColor(kRed);
		line->Draw("SAME");
		line->DrawLine( sbRR, 0, sbRR, massHistPi0->GetMaximum());
		line->DrawLine( sigL, 0, sigL, massHistPi0->GetMaximum());
		line->DrawLine( sigR, 0, sigR, massHistPi0->GetMaximum());
		line->DrawLine( sbLL, 0, sbLL, massHistPi0->GetMaximum());
		line->DrawLine( sbLR, 0, sbLR, massHistPi0->GetMaximum());

		allCanvases->SaveAs(("fitResults/Mpi0_fit_"+detectorNames[i]+".png").c_str());
		cout << "Saved for a specific fit!" << endl;
	}
}




































