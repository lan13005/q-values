#include "run.h"
#include <ctime>
bool verbose = true;

void run(){
	clock_t start;
	double duration;
	start = clock();
	
	TFile* dataFile=new TFile("/d/grid15/ln16/pi0eta/121818/z_pi0eta_a0_reco/pi0eta_a0_recotreeFlat_DSelector.root");
	TTree *dataTree;
	dataFile->GetObject("pi0eta_a0_recotree_flat",dataTree);

    	ofstream logFile;
    	logFile.open("logEventChiSqQValue.txt");
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	TH1F* discriminatorHist = new TH1F("","",50,0.35,0.8);
	
	double Meta;
	double Mpi0;
	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("Mpi0",&Mpi0);

	Long64_t nentries=dataTree->GetEntries();
	nentries=20;
	cout << "Total Entries: " << nentries << endl;
	const int c_nentries = (const int)nentries;
	double Metas[c_nentries];
	double Mpi0s[c_nentries];
	for (int ientry=0; ientry<nentries; ientry++)
	{
		dataTree->GetEntry(ientry);
		Metas[ientry]=Meta;
		Mpi0s[ientry]=Mpi0;
	}
	cout << "LOADED ALL THE DATA INTO MEMORY" << endl;

	map<double, int> mapDistToJ;
	set<double> distances;
	double phasePoint1[dim];
	double phasePoint2[dim];
	double distance;

	for (int ientry=0; ientry<nentries; ientry++){
		if(verbose) { cout << "Getting next event!\n--------------------------------\n" << endl;  }
		duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << "Starting event " << ientry << "/" << nentries << " ---- Time: " << duration << "s" << endl;
		mapDistToJ.clear();
		distances.clear();

		phasePoint1[0] = Mpi0s[ientry];
		for (int jentry=0; jentry<nentries; jentry++){
			phasePoint2[0] = Mpi0s[jentry];
			if (jentry != ientry){
				distance = calc_distance(phasePoint1,phasePoint2);
				if ( distances.size() >= kDim ) {
					if ( distance < *distances.rbegin() ) {
						mapDistToJ.erase(*distances.rbegin());
						mapDistToJ[distance] = jentry;
						distances.erase(*distances.rbegin());
						distances.insert(distance);
					}
				}
				else{
					mapDistToJ[distance] = jentry;
					distances.insert(distance);
				}
			}
			if ( verbose) { 
				cout << "CURRENT SET: " << endl;
				for(auto elem : mapDistToJ){
					std::cout << elem.first << " " << elem.second << "\n";
				}
			}
		}
	}


	//TF1* fit = new 	TF1("fit",fitFunc,0.52,0.58,numDOFbkg+numDOFsig);
	//TF1* bkgFit = new TF1("bkgFit",background,0.52,0.58,numDOFbkg);
	//TF1* sigFit = new TF1("sigFit",signal,0.52,0.58,numDOFsig);
	//bkgFit->SetLineColor(kMagenta);
	//sigFit->SetLineColor(kBlue);
	//Double_t par[numDOFbkg+numDOFsig];
	//fit->SetParameters(5,20,0.545,0.02);
	//discriminatorHist->Fit("fit","RQ");
	//fit->GetParameters(par);
	//bkgFit->SetParameters(par);
	//bkgFit->Draw("same");
	//sigFit->SetParameters(&par[numDOFbkg]);
	//sigFit->Draw("same");
	
	logFile.close();
}
	//allCanvases->SaveAs(("/d/grid15/ln16/pi0eta/121818/readRootTrees/fitResults/pi0eta_datamc_overlaid_"+names_vertex[i_vertex]+".png").c_str());

