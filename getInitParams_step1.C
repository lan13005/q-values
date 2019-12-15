#include "main.h"
#include "/d/grid15/ln16/pi0eta/092419/makeGraphs.h"

void getInitParams_step1(){
	//TFile* dataFile=new TFile("pi0eta_a0_recotreeFlat_DSelector.root");
    	ofstream logFile_eta;
    	ofstream logFile_pi0;
    	logFile_eta.open("fitResults/etaFitNoAccSub.txt");
    	logFile_pi0.open("fitResults/pi0FitNoAccSub.txt");

	gStyle->SetOptFit(111);
	gStyle->SetOptStat(0);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
	TFile* dataFile=new TFile("pi0eta_fcal_tLT1_treeFlat_DSelector.root");
	TTree *dataTree;
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	dataFile->GetObject("pi0eta_fcaltree_flat",dataTree);
	double Meta;
	double Mpi0;
	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("Mpi0",&Mpi0);
	nentries=dataTree->GetEntries();

        TH1F *massHist; 
	TF1* fit;
        TF1* bkgFit;
        TF1* sigFit;
        std::vector<double> binRange;
	double binWidth;
        double par[5];
	
        cout <<"Initialized" << endl;
        string namePar[8] = {"nentries","const","linear","amp1","mass","sigma1","ampRatio","sigmaRatio"};

        for (int i=0; i<1; ++i){
		allCanvases->Clear();

		std::vector<double> binRange;
		std::vector<double> fitRange;
		binRange={100,0.25,0.8};
		//binRange={300,0.25,0.8};
		fitRange={0.38,0.65};

		binWidth=(binRange[2]-binRange[1])/binRange[0];
		fit = new TF1("fit",fitFunc,fitRange[0],fitRange[1],numDOFbkg+numDOFsig);
		bkgFit = new TF1("bkgFit",background,fitRange[0],fitRange[1],numDOFbkg);
		sigFit = new TF1("sigFit",signalDG,fitRange[0],fitRange[1],numDOFsig);

		//fit->SetParameters(49965.9,7668.6,2000,0.547,0.02,0.2,5);
		fit->SetParameters(11500,250,1750,0.547,0.003,10,5);
		//fit->SetParameters(229,5,34,0.547,0.007,10,5);

		massHist = new TH1F("","",binRange[0],binRange[1],binRange[2]);
		cout << "Initialized for a specific mass (eta/pi0) fit" << endl;
		
		for (int ientry=0; ientry<nentries; ientry++)
		{
			dataTree->GetEntry(ientry);
		        massHist->Fill(Meta);
		}
		cout << "Filled all entries into histogram for a specific fit" << endl;
		
		massHist->Fit("fit","RQB"); // B will enforce the bounds
		fit->GetParameters(par);
		bkgFit->SetParameters(par);
		sigFit->SetParameters(&par[numDOFbkg]);
		massHist->SetTitle(";M(#eta) (GeV)");
		logFile_eta << namePar[0] << " " << nentries << endl;
		
		for (int iPar=0; iPar<numDOFbkg+numDOFsig; ++iPar){
		        logFile_eta << namePar[iPar+1] << " " << par[iPar] << endl;
		}
		
		double integralBKG = bkgFit->Integral(fitRange[0],fitRange[1]);
		double integralSIG = sigFit->Integral(fitRange[0],fitRange[1]);
		cout << "IntegralBKG before scaling: " << integralBKG << endl;
		cout << "IntegralSIG before scaling: " << integralSIG << endl;
		integralBKG *= 1/binWidth;
		integralSIG *= 1/binWidth;

		cout << "IntegralBKG after scaling: " << integralBKG << endl;
		cout << "IntegralSIG after scaling: " << integralSIG << endl;
		cout << "nentries: " << nentries << endl;

		
		massHist->Draw();
		massHist->SetTitle(("Peak: "+to_string(par[2])+"    width: "+to_string(par[3])).c_str());
		allCanvases->SaveAs("fitResults/Meta_fit.png");
		cout << "Saved for a specific fit!" << endl;
        }


        
}
