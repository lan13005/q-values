#include "main.h"

void getInitParams_step1(){
	//TFile* dataFile=new TFile("pi0eta_a0_recotreeFlat_DSelector.root");
    	ofstream logFile_eta;
    	ofstream logFile_pi0;
    	logFile_eta.open("fitResults/etaFitNoAccSub.txt");
    	logFile_pi0.open("fitResults/pi0FitNoAccSub.txt");

	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
	TFile* dataFile=new TFile("pi0eta_datatreeFlat_DSelector.root");
	TTree *dataTree;
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	dataFile->GetObject("pi0eta_datatree_flat",dataTree);
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
        bool useEta;
        double par[5];

        cout <<"Initialized" << endl;
        string namePar[6] = {"nentries","par0","par1","par2","par3","par4"};

        for (int i=0; i<2; ++i){
            allCanvases->Clear();
            if (i==0) {useEta=true;}
            else {useEta=false; }

            std::vector<double> binRange;
            std::vector<double> fitRange;
            if (useEta){
                binRange={5000,0.35,0.8};
                fitRange={0.425,0.7};
            } 
            else{ 
                binRange={5000,0.05,0.25};
                fitRange={0.1,0.17};
            }
	    fit = new TF1("fit",fitFunc,fitRange[0],fitRange[1],numDOFbkg+numDOFsig);
            if (useEta){
	        fit->SetParameters(5,100,60,0.55,0.025);
            }
            else {
	        fit->SetParameters(10,0,140,0.136,0.01);
            }
	    massHist = new TH1F("","",binRange[0],binRange[1],binRange[2]);
            cout << "Initialized for a specific mass (eta/pi0) fit" << endl;


	    for (int ientry=0; ientry<nentries; ientry++)
	    {
	    	dataTree->GetEntry(ientry);
                if ( useEta){
                    massHist->Fill(Meta);
                }
                else{
                    massHist->Fill(Mpi0);
                }
            }
            cout << "Filled all entries into histogram for a specific fit" << endl;

	    massHist->Fit("fit","RQB"); // B will enforce the bounds
	    fit->GetParameters(par);
            if (useEta){
                logFile_eta << namePar[0] << " " << nentries << endl;
            } else {
                logFile_pi0 << namePar[0] << " " << nentries << endl;
            }
            for (int iPar=0; iPar<numDOFbkg+numDOFsig; ++iPar){
                if ( useEta ) {
                    logFile_eta << namePar[iPar+1] << " " << par[iPar] << endl;
                }
                else { 
                    logFile_pi0 << namePar[iPar+1] << " " << par[iPar] << endl;
                }
            }
            massHist->Draw();
            massHist->SetTitle(("Peak: "+to_string(par[2])+"    width: "+to_string(par[3])).c_str());
            //fit->Draw();
            if (useEta){
                allCanvases->SaveAs("fitResults/Meta_fit.png");
            }
            else{
                allCanvases->SaveAs("fitResults/Mpi0_fit.png");
            }
            cout << "Saved for a specific fit!" << endl;
        }


        
}
