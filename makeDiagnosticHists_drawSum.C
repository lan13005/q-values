#include "makeDiagnosticHists.h"

void makeDiagnosticHists_drawSum(){
	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	
	TFile* dataFile2 = new TFile(("diagnosticPlots"+runTag+"/postQVal_hists.root").c_str());
	TH1F* totHist;
	TH1F* bkgHist;
	TH1F* sigHist;
	TH1F* bkgHist_sb;
	TH1F* sigHist_sb;

	// HERE WE WILL JUST DRAW SOME OF THE HISTOGRAMS WITH THE BKG FILLED IN TO SEE THEIR CONTRIBUTION
	string particles[6] = {"Meta", "Mpi0", "Mpi0eta","cosThetaEta_GJ","cosThetaX_CM","phiEta_GJ"};
	string sigTypes[3] = {"bkg","sig","tot"};
	string measurementTypes[2] = {"meas","kin"};
	for (auto particle : particles){ 
		for (auto measurementType : measurementTypes){
			for (auto sigType : sigTypes){
				cout << "Grabbing "+particle+"_"+sigType+"_"+measurementType << endl;
				if (sigType == "bkg"){ 
                                    dataFile2->GetObject((particle+"_"+sigType+"_"+measurementType).c_str(),bkgHist); 
                                    dataFile2->GetObject((particle+"_"+sigType+"_"+measurementType+"_sb").c_str(),bkgHist_sb); 
                                }  
				if (sigType == "sig"){ 
                                    dataFile2->GetObject((particle+"_"+sigType+"_"+measurementType).c_str(),sigHist); 
                                    dataFile2->GetObject((particle+"_"+sigType+"_"+measurementType+"_sb").c_str(),sigHist_sb); 
                                }  
				if (sigType == "tot"){ dataFile2->GetObject((particle+"_"+sigType+"_"+measurementType).c_str(),totHist); }  
			}
			cout << "Grabbed tot/bkg/sig... \n Making stacked histograms" << endl;
			makeStackedHist(totHist,sigHist,bkgHist,sigHist_sb,bkgHist_sb,particle+measurementType, ("diagnosticPlots"+runTag).c_str());
			cout << "----------" << endl;
		}
	}
	dataFile2->GetObject("Mpi0g_bkg", bkgHist);
	dataFile2->GetObject("Mpi0g_sig", sigHist);
	dataFile2->GetObject("Mpi0g_bkg_sb", bkgHist_sb);
	dataFile2->GetObject("Mpi0g_sig_sb", sigHist_sb);
	dataFile2->GetObject("Mpi0g_tot", totHist);
	cout << "Grabbed tot/bkg/sig... \n Making stacked histograms" << endl;
	makeStackedHist(totHist,sigHist,bkgHist,sigHist_sb,bkgHist_sb,"Mpi0gkin", ("diagnosticPlots"+runTag).c_str());
	cout << "----------" << endl;


}

