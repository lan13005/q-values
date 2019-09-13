#include "run.h"
#include <ctime>
bool verbose = false;

void run(){
	// Starting timing
	clock_t start;
	double duration;
	start = clock();
	
	// setting up some basic root stuff and getting the file and tree
	TFile* dataFile=new TFile("/d/grid15/ln16/pi0eta/121818/z_pi0eta_a0_reco/pi0eta_a0_recotreeFlat_DSelector.root");
	TTree *dataTree;
	dataFile->GetObject("pi0eta_a0_recotree_flat",dataTree);
	const int dim=5;
	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	TH1F* discriminatorHist;
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
	double cosTheta_X_cm;
	double phi_X_cm;
	double cosTheta_eta_gj;
	double phi_eta_gj;
	double cosThetaHighestEphotonIneta_gj;
	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("Mpi0",&Mpi0);
	dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
        dataTree->SetBranchAddress("cosTheta_X_cm", &cosTheta_X_cm); 
        dataTree->SetBranchAddress("phi_X_cm",&phi_X_cm); 
        dataTree->SetBranchAddress("cosTheta_eta_gj",&cosTheta_eta_gj);
        dataTree->SetBranchAddress("phi_eta_gj",&phi_eta_gj); 
        dataTree->SetBranchAddress("cosThetaHighestEphotonIneta_gj",&cosThetaHighestEphotonIneta_gj);
	Long64_t nentries=dataTree->GetEntries();
	nentries=1000;
	cout << "Total Entries: " << nentries << endl;

	// opening a file to write my log data to
    	ofstream logFile;
	ofstream entireLogFile;
    	logFile.open(("logEventChiSqQValue_nentries"+to_string(nentries)+".txt").c_str());
    	entireLogFile.open(("logEventChiSqQValue_nentries"+to_string(nentries)+"-all.txt").c_str());
	logFile << "Event\tQ-Value\tChiSq\tMpi0" << endl;
	entireLogFile << "Event\tQ-Value\tChiSq\tMpi0" << endl;
	
	// randomly select some events to write histograms for 
	int numberEventsToSave=100;
	const int kDim=500;
	set<int> selectRandomIdxToSave;
	const int c_nentries = (const int)nentries;

	// importing all the data to RAM instead of reading from root file
	double Metas[c_nentries];
	double Mpi0s[c_nentries];
	double Mpi0etas[c_nentries];
	double cosTheta_X_cms[c_nentries];
	double phi_X_cms[c_nentries];
	double cosTheta_eta_gjs[c_nentries];
	double phi_eta_gjs[c_nentries];
	double cosThetaHighestEphotonIneta_gjs[c_nentries];
	for (int ientry=0; ientry<nentries; ientry++)
	{
		dataTree->GetEntry(ientry);
		Metas[ientry]=Meta;
		Mpi0s[ientry]=Mpi0;
		Mpi0etas[ientry]=Mpi0eta;
		cosTheta_X_cms[ientry]=cosTheta_X_cm;
		phi_X_cms[ientry]=phi_X_cm;
		cosTheta_eta_gjs[ientry]=cosTheta_eta_gj;
		phi_eta_gjs[ientry]=phi_eta_gj;
		cosThetaHighestEphotonIneta_gjs[ientry]=cosThetaHighestEphotonIneta_gj;	 
		if ( selectRandomIdxToSave.size() < numberEventsToSave){
			selectRandomIdxToSave.insert(rand() % nentries);
		}
		
	}
	// outputting the results before and after standardizeArray will show that it works
	// for(auto& cosTheta_X_cm1 : cosTheta_X_cms){ cout << cosTheta_X_cm1 << endl; }
	standardizeArray(cosTheta_X_cms, nentries, "cosTheta_X_cms"); 
	standardizeArray(phi_X_cms, nentries, "phi_X_cms"); 
	standardizeArray(cosTheta_eta_gjs,nentries,"cosTheta_eta_gjs");
	standardizeArray(phi_eta_gjs,nentries,"phi_eta_gjs");
	standardizeArray(cosThetaHighestEphotonIneta_gjs,nentries,"cosThetaHighestEphotonIneta_gjs");

	// cout << "LOADED ALL THE DATA INTO MEMORY" << endl;

	// defining some variables we will use in the main loop to get the distances and then the q-values
	map<double, int> mapDistToJ;
	set<double> distances;
	double phasePoint1[dim];
	double phasePoint2[dim];
	double distance;
	double qvalue;
	double conjugate_qvalue;

	// the main loop where we loop through all events in a double for loop to calculate dij. Iterating through all j we can find the k nearest neighbors to event i.
	// Then we can plot the k nearest neighbors in the discriminating distribution which happens to be a double gaussian and flat bkg. Calculate Q-Value from event
	// i's discriminating variable's value. This value will be plugged into the signal PDF and the total PDF. The ratio of these two values are taken which is the Q-Value.
	for (int ientry=0; ientry<nentries; ientry++){
		if(verbose) { cout << "Getting next event!\n--------------------------------\n" << endl;  }
		duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << "Starting event " << ientry << "/" << nentries << " ---- Time: " << duration << "s" << endl;
		mapDistToJ.clear();
		distances.clear();
		allCanvases->Clear();
		discriminatorHist = new TH1F("","",75,0.35,0.8);

		phasePoint1[0] = cosTheta_X_cms[ientry];
		phasePoint1[1] = phi_X_cms[ientry];
		phasePoint1[2] = cosTheta_eta_gjs[ientry];
		phasePoint1[3] = phi_eta_gjs[ientry];
		phasePoint1[4] = cosThetaHighestEphotonIneta_gjs[ientry];
		for (int jentry=0; jentry<nentries; jentry++){
			phasePoint2[0] = cosTheta_X_cms[jentry];
			phasePoint2[1] = phi_X_cms[jentry];
			phasePoint2[2] = cosTheta_eta_gjs[jentry];
			phasePoint2[3] = phi_eta_gjs[jentry];
			phasePoint2[4] = cosThetaHighestEphotonIneta_gjs[jentry];
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
		// Filling the discriminatorHist with all the nearest neighbors 
		for(auto elem : mapDistToJ){
			discriminatorHist->Fill(Metas[elem.second]);
		}
	
		// Building the fit functions. We have to set some parameter limits to try to guide Minuit. We choose a double gaussian since for whatever reason the Meta has a asymmetry
		// such that there are more events to the left side of the peak. The second gaussian will have lower amplitude and a large sigma with a left shifted mean. We also want to 
		// fix the amplitudes and constnat background to be strictly positive and choose kDim as the max value ( which is reached only if all values are filled into a single bin)
		TF1* fit = new 	TF1("fit",fitFunc,0.425,0.7,numDOFbkg+numDOFsig);
		TF1* bkgFit = new TF1("bkgFit",background,0.425,0.7,numDOFbkg);
		TF1* sigFit = new TF1("sigFit",signal,0.425,0.7,numDOFsig);
		bkgFit->SetLineColor(kMagenta);
		sigFit->SetLineColor(kBlue);
		Double_t par[numDOFbkg+numDOFsig];
		fit->SetParName(0,"const");
		fit->SetParName(1,"Amp_Gaus1");
		fit->SetParName(2,"Mean_Gaus1");
		fit->SetParName(3,"Sigma_Gaus1");
		fit->SetParName(4,"Amp_Gaus2");
		fit->SetParName(5,"Mean_Gaus2");
		fit->SetParName(6,"Sigma_Gaus2");
		fit->SetParameters(5,70,0.545,0.02,10,0.51,0.04);
		fit->SetParLimits(0, 0, kDim);
		fit->SetParLimits(1, 0, kDim);
		fit->SetParLimits(4, 0, kDim);
		fit->SetParLimits(6, 0.02, 0.06);
		discriminatorHist->Fit("fit","RQB"); // B will enforce the bounds
		fit->GetParameters(par);
		bkgFit->SetParameters(par);
		sigFit->SetParameters(&par[numDOFbkg]);
		qvalue=sigFit->Eval(Metas[ientry])/fit->Eval(Metas[ientry]);
		conjugate_qvalue = 1-qvalue;
		Meta_Qd->Fill(Metas[ientry]);
		Meta_bkg_Qd->Fill(Metas[ientry], conjugate_qvalue);
		Mpi0_Qd->Fill(Mpi0s[ientry]);
		Mpi0_bkg_Qd->Fill(Mpi0s[ientry], conjugate_qvalue);
		Mpi0eta_Qd->Fill(Mpi0etas[ientry]);
		Mpi0eta_bkg_Qd->Fill(Mpi0etas[ientry], conjugate_qvalue);

		// Here we draw the histograms that were randomly selected
		if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end()) {
			logFile << ientry << "\t" << qvalue << "\t" << fit->GetChisquare()<< "\t" << Mpi0s[ientry] << endl;
			discriminatorHist->SetTitle(("QValue="+to_string(qvalue)+"  ChiSq="+to_string(fit->GetChisquare())).c_str());
			discriminatorHist->Draw();
        		etaLine = new TLine(Metas[ientry],0,Metas[ientry],discriminatorHist->GetMaximum());
			etaLine->SetLineColor(kOrange);
			//bkgFit->Draw("same");
			//sigFit->Draw("same");
  			bkgFit->SetFillColor(kMagenta);
  			bkgFit->SetFillStyle(3004);
  			bkgFit->Draw("SAME FC");
  			sigFit->SetFillColor(kBlue);
  			sigFit->SetFillStyle(3005);
  			sigFit->Draw("SAME FC");
			etaLine->Draw("same");
			allCanvases->SaveAs(("Meta-event"+std::to_string(ientry)+".png").c_str());
		}
		entireLogFile << ientry << "\t" << qvalue << "\t" << fit->GetChisquare() << "\t" << Mpi0s[ientry] << endl;
	}

	// HERE WE WILL JUST DRAW SOME OF THE HISTOGRAMS WITH THE BKG FILLED IN TO SEE THEIR CONTRIBUTION
	// ----------------- Meta 
	THStack* stackedHists = new THStack("stackedHists","");
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


	// Finish the log files by including an elapsed time and finally closing the file
	logFile << "Total time: " << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << "s" << endl;
	logFile << "Time Per Event: " << ( std::clock() - start ) / (double) CLOCKS_PER_SEC / nentries << "s" << endl;
	logFile.close();
	entireLogFile << "Total time: " << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << "s" << endl;
	entireLogFile << "Time Per Event: " << ( std::clock() - start ) / (double) CLOCKS_PER_SEC / nentries << "s" << endl;
	entireLogFile.close();
}

