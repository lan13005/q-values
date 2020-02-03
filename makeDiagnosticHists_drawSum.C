#include <ctime>
#include <math.h> 
#include "makeDiagnosticHists.h"
bool verbose = true;

void makeDiagnosticHists_drawSum(){
	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	
	TFile* dataFile2 = new TFile("postQVal.root");
	TH1F* totHist;
	TH1F* bkgHist;
	TH1F* sigHist;

	// HERE WE WILL JUST DRAW SOME OF THE HISTOGRAMS WITH THE BKG FILLED IN TO SEE THEIR CONTRIBUTION
	string particles[3] = {"Meta", "Mpi0", "Mpi0eta"};
	string sigTypes[3] = {"bkg","sig","tot"};
	string measurementTypes[2] = {"meas","kin"};
	for (auto particle : particles){ 
		for (auto measurementType : measurementTypes){
			for (auto sigType : sigTypes){
				cout << "Grabbing "+particle+"_"+sigType+"_"+measurementType << endl;
				if (sigType == "bkg"){ dataFile2->GetObject((particle+"_"+sigType+"_"+measurementType).c_str(),bkgHist); }  
				if (sigType == "sig"){ dataFile2->GetObject((particle+"_"+sigType+"_"+measurementType).c_str(),sigHist); }  
				if (sigType == "tot"){ dataFile2->GetObject((particle+"_"+sigType+"_"+measurementType).c_str(),totHist); }  
			}
			cout << "Grabbed tot/bkg/sig... \n Making stacked histograms" << endl;
			makeStackedHist(totHist,sigHist,bkgHist,particle+measurementType, "diagnosticPlots");
			cout << "----------" << endl;
		}
	}
	dataFile2->GetObject("Mpi0g_bkg", bkgHist);
	dataFile2->GetObject("Mpi0g_sig", sigHist);
	dataFile2->GetObject("Mpi0g_tot", totHist);
	cout << "Grabbed tot/bkg/sig... \n Making stacked histograms" << endl;
	makeStackedHist(totHist,sigHist,bkgHist,"Mpi0gkin", "diagnosticPlots");
	cout << "----------" << endl;


}

