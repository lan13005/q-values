#include "main.h"

//bool verbose = false;
//int kDim=200;
//int numberEventsToSave=5;
//bool override_nentries=true;
//Long64_t nentries=100;
//int nProcess=5;
//int seedShift=123125;

bool useEta=false;

//void main(int iProcess, int kDim, int numberEventsToSave, int nProcess, int seedShift, Long64_t nentries, bool override_nentries, bool verbose){
int main( int argc, char* argv[] ){
        //int iProcess = std::stoi(argv[0]);
        //int kDim= std::stoi(argv[1]);
        //int numberEventsToSave= std::stoi(argv[2]);
        //int nProcess= std::stoi(argv[3]);
        //int seedShift= std::stoi(argv[4]);
        //Long64_t nentries= std::stoll(argv[5]);
        int iProcess = atoi(argv[1]);
        int kDim= atoi(argv[2]);
        int numberEventsToSave= atoi(argv[3]);
        int nProcess= atoi(argv[4]);
        int seedShift= atoi(argv[5]);
        Long64_t nentries= atoi(argv[6]);
        bool override_nentries;
        if ( atoi(argv[7]) ==1 ){ 
            override_nentries=true;
        }
        else{ 
            override_nentries=false;
        }
        bool verbose;
        if ( atoi(argv[8])==1 ){ verbose=true;}
        else{ verbose=false;}
        cout << "----------------------------" << endl;
        cout << "iProcess: " << iProcess << endl;
        cout << "kDim: " << kDim << endl;
        cout << "numberEventsToSave: " << numberEventsToSave << endl;
        cout << "nProcess: " << nProcess << endl;
        cout << "seedShift: " << seedShift << endl;
        cout << "nentries: " << nentries << endl; 
        cout << "override_nentries: " << override_nentries << endl;
        cout << "verbose: " << verbose  << endl; 
        cout << "----------------------------" << endl;
    
	// Starting timing
	//clock_t start;
	//double duration;
	auto start2 = std::chrono::high_resolution_clock::now();
	//start = clock();
	
	// setting up some basic root stuff and getting the file and tree
	TFile* dataFile=new TFile("/d/grid15/ln16/pi0eta/121818/z_pi0eta_a0_reco/pi0eta_a0_recotreeFlat_DSelector.root");
	TTree *dataTree;
	dataFile->GetObject("pi0eta_a0_recotree_flat",dataTree);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
        TLine* etaLine;
	TH1F* discriminatorHist;
	double Meta;
	double Mpi0;
	double Mpi0eta;
	double cosTheta_X_cm;
	double phi_X_cm;
	double cosTheta_eta_gj;
	double phi_eta_gj;
	double cosThetaHighestEphotonIneta_gj;
	double cosThetaHighestEphotonInpi0_cm;
        ULong64_t eventNumber;
        double uniqueComboID;
	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("Mpi0",&Mpi0);
	dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
        dataTree->SetBranchAddress("cosTheta_X_cm", &cosTheta_X_cm); 
        dataTree->SetBranchAddress("phi_X_cm",&phi_X_cm); 
        dataTree->SetBranchAddress("cosTheta_eta_gj",&cosTheta_eta_gj);
        dataTree->SetBranchAddress("phi_eta_gj",&phi_eta_gj); 
        dataTree->SetBranchAddress("cosThetaHighestEphotonIneta_gj",&cosThetaHighestEphotonIneta_gj);
        dataTree->SetBranchAddress("cosThetaHighestEphotonInpi0_cm",&cosThetaHighestEphotonInpi0_cm);
        dataTree->SetBranchAddress("uniqueComboID",&uniqueComboID);
        dataTree->SetBranchAddress("event",&eventNumber);
	if (!override_nentries){
		nentries=dataTree->GetEntries();
	}
	if(verbose){cout << "Chosen Total Entries: " << nentries << endl;}

	double batchEntries = nentries/nProcess;
	double lowest_nentry = iProcess*batchEntries;
	double largest_nentry = (iProcess+1)*batchEntries;
	if(verbose){cout << "nentries we will use for this process: " << lowest_nentry << ", " << largest_nentry << endl;}

	// opening a file to write my log data to
    	ofstream logFile;
    	logFile.open(("logs/logEventChiSqQValue_process"+to_string(iProcess)+".txt").c_str());
	//logFile << "Event\tQ-Value\tChiSq\tMpi0" << endl;
	
	// randomly select some events to write histograms for 
	set<int> selectRandomIdxToSave;
	int randomEvent;
	srand(iProcess+seedShift);
	for (int i=0; i<numberEventsToSave/nProcess; i++){
		randomEvent = rand() % (int)batchEntries;
		randomEvent += lowest_nentry;
		selectRandomIdxToSave.insert( randomEvent );
	}

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
	double cosThetaHighestEphotonInpi0_cms[c_nentries];

	// We will use a ientry to keep track of which entries we will get from the tree. We will simply use ientry when filling the arrays.  
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
		cosThetaHighestEphotonInpi0_cms[ientry]=cosThetaHighestEphotonInpi0_cm;	 
	}
	dataFile->Close();
	
        if ( verbose_outputDistCalc ) {
            cout << "Before standarization" << endl;
            for ( int ientry=0 ; ientry < nentries; ientry++){
                cout << phi_X_cms[ientry] << "," << cosTheta_eta_gjs[ientry] << endl;
            }
        }
	//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
	//logFile << "Time marker to load data: " << duration << "ms" << endl;

	// outputting the results before and after standardizeArray will show that it works
	// for(auto& cosTheta_X_cm1 : cosTheta_X_cms){ cout << cosTheta_X_cm1 << endl; }
	standardizeArray(cosTheta_X_cms, nentries, "cosTheta_X_cms"); 
	standardizeArray(phi_X_cms, nentries, "phi_X_cms"); 
	standardizeArray(cosTheta_eta_gjs,nentries,"cosTheta_eta_gjs");
	standardizeArray(phi_eta_gjs,nentries,"phi_eta_gjs");
	standardizeArray(cosThetaHighestEphotonIneta_gjs,nentries,"cosThetaHighestEphotonIneta_gjs");
	standardizeArray(cosThetaHighestEphotonInpi0_cms,nentries,"cosThetaHighestEphotonInpi0_cms");
	//
	//
        if ( verbose_outputDistCalc ) {
	    cout << "After standardization" << endl;
            for ( int ientry=0 ; ientry < nentries; ientry++){
                cout << cosTheta_X_cms[ientry] << "," << phi_X_cms[ientry] << endl;
            }
        }


	// defining some variables we will use in the main loop to get the distances and then the q-values
	map<double, int> mapDistToJ;
	set<double> distances;
	double phasePoint1[dim];
	double phasePoint2[dim];
	double distance;
	double qvalue;
	double conjugate_qvalue;

        distSort_kNN distKNN(kDim);
        pair<double,int> newPair;

	// the main loop where we loop through all events in a double for loop to calculate dij. Iterating through all j we can find the k nearest neighbors to event i.
	// Then we can plot the k nearest neighbors in the discriminating distribution which happens to be a double gaussian and flat bkg. Calculate Q-Value from event
	// i's discriminating variable's value. This value will be plugged into the signal PDF and the total PDF. The ratio of these two values are taken which is the Q-Value.
	//logFile << std::fixed << std::setprecision(6);
	for (int ientry=lowest_nentry; ientry<largest_nentry; ientry++){ 
		if(verbose) { cout << "Getting next event!\n--------------------------------\n" << endl;  }
		//duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
		auto duration_beginEvent = std::chrono::high_resolution_clock::now();
		if(verbose){cout << "Starting event " << ientry << "/" << largest_nentry << " ---- Time: " << duration2 << "ms" << endl; }
		if(verbose2){logFile << "Starting event " << ientry << "/" << largest_nentry << " ---- Time: " << duration2 << "ms" << endl; }
		mapDistToJ.clear();
		distances.clear();
		allCanvases->Clear();
                if ( useEta) { 
		    discriminatorHist = new TH1F("","",50,0.35,0.8);
		}
                else {
		    discriminatorHist = new TH1F("","",50,0.05,0.25);
                }

		phasePoint1[0] = cosTheta_X_cms[ientry];
		phasePoint1[1] = phi_X_cms[ientry];
		phasePoint1[2] = cosTheta_eta_gjs[ientry];
		phasePoint1[3] = phi_eta_gjs[ientry];
		phasePoint1[4] = cosThetaHighestEphotonIneta_gjs[ientry];
		phasePoint1[5] = cosThetaHighestEphotonInpi0_cms[ientry];
		for (int jentry=0; jentry<nentries; jentry++){
                        if ( verbose_outputDistCalc ) { cout << "event i,j = " << ientry << "," << jentry << endl;} 
			phasePoint2[0] = cosTheta_X_cms[jentry];
			phasePoint2[1] = phi_X_cms[jentry];
			phasePoint2[2] = cosTheta_eta_gjs[jentry];
			phasePoint2[3] = phi_eta_gjs[jentry];
			phasePoint2[4] = cosThetaHighestEphotonIneta_gjs[jentry];
			phasePoint2[5] = cosThetaHighestEphotonInpi0_cms[jentry];
			if (jentry != ientry){
		    	        distance = calc_distance(phasePoint1,phasePoint2);
                                //distance = rgen.Uniform(nentries);
                                distKNN.insertPair(make_pair(distance,jentry));
			}
			//if ( verbose) { 
			//	cout << "CURRENT SET: " << endl;
			//	for(auto elem : mapDistToJ){
			//		std::cout << elem.first << " " << elem.second << "\n";
			//	}
			//}
		}
		duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
		if(verbose2){logFile << "	Found neighbors: " << duration2 << "ms" << endl; }
		//// Filling the discriminatorHist with all the nearest neighbors 
		//for(auto elem : mapDistToJ){
		//	discriminatorHist->Fill(Metas[elem.second]);
		//}
		//for (int i=0; i<200; i++){
		//	discriminatorHist->Fill(rgen.Gaus(0.55,0.15));
		//}
		if (distKNN.kNN.size() != kDim){ cout << "size of distKNN is not equal to kDim!" << endl; exit(0); }
                //cout << "New Event\n" ;
                while ( distKNN.kNN.empty() == false ){
                        newPair = distKNN.kNN.top();
                        distKNN.kNN.pop();
                        if ( useEta ){
                            discriminatorHist->Fill(Metas[newPair.second]);
                        }
                        else {
                            discriminatorHist->Fill(Mpi0s[newPair.second]);
                        }
                        //cout << "(" << newPair.first << ", " << newPair.second << ")"; 
                        //cout << endl; 
                }
		//duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
		//if(verbose2){logFile << "	Filled neighbors: " << duration2 << "ms" << endl;}
	
		// Building the fit functions. We have to set some parameter limits to try to guide Minuit. We choose a double gaussian since for whatever reason the Meta has a asymmetry
		// such that there are more events to the left side of the peak. The second gaussian will have lower amplitude and a large sigma with a left shifted mean. We also want to 
		// fix the amplitudes and constnat background to be strictly positive and choose kDim as the max value ( which is reached only if all values are filled into a single bin)
		
    
	        TF1* fit;
                TF1* bkgFit;
                TF1* sigFit;
                if (useEta){
		    fit = new TF1("fit",fitFunc,0.425,0.7,numDOFbkg+numDOFsig);
		    bkgFit = new TF1("bkgFit",background,0.425,0.7,numDOFbkg);
		    sigFit = new TF1("sigFit",signal,0.425,0.7,numDOFsig);
		}
                else { 
		    fit = new TF1("fit",fitFunc,0.1,0.17,numDOFbkg+numDOFsig);
		    bkgFit = new TF1("bkgFit",background,0.1,0.17,numDOFbkg);
		    sigFit = new TF1("sigFit",signal,0.1,0.17,numDOFsig);
                }
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
                if (useEta) { 
		    fit->SetParameters(5,70,0.545,0.02,10,0.51,0.04);
		    fit->SetParLimits(0, 0, kDim);
		    fit->SetParLimits(1, 0, kDim);
		    fit->SetParLimits(2,0.53,0.56);
		    fit->SetParLimits(4, 0, kDim);
		    fit->SetParLimits(5,0.49,0.55);
		    fit->SetParLimits(6, 0.02, 0.06);
                }
                else {
		    fit->SetParameters(5,70,0.135,0.01,10,0,0,0);
                    fit->SetParLimits(0, 0, kDim);
                    fit->SetParLimits(1, 0, kDim);
                    fit->SetParLimits(2, 0.12, 0.1425);
                    fit->SetParLimits(3, 0.001, 0.01);
                    fit->SetParLimits(4,0,0);
                    fit->SetParLimits(5,0,0);
                    fit->SetParLimits(6,0,0);
		}
		discriminatorHist->Fit("fit","RQB"); // B will enforce the bounds
		fit->GetParameters(par);
		bkgFit->SetParameters(par);
		sigFit->SetParameters(&par[numDOFbkg]);
                if ( useEta) { 
		    qvalue=sigFit->Eval(Metas[ientry])/fit->Eval(Metas[ientry]);
                }
                else {
		    qvalue=sigFit->Eval(Mpi0s[ientry])/fit->Eval(Mpi0s[ientry]);
                }
		conjugate_qvalue = 1-qvalue;
		if(verbose2){
			duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
			logFile << "	Got Q-value: " << duration2 <<  "ms" << endl;
		}

		// Here we draw the histograms that were randomly selected
		logFile << ientry << "\t" << qvalue << "\t" << fit->GetChisquare()<< "\t" << Mpi0s[ientry] << endl;
		if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end()) {
			discriminatorHist->SetTitle(("QValue="+to_string(qvalue)+"  ChiSq="+to_string(fit->GetChisquare())).c_str());
			discriminatorHist->Draw();
                        if ( useEta) { 
        		    etaLine = new TLine(Metas[ientry],0,Metas[ientry],discriminatorHist->GetMaximum());
                        }
                        else { 
        		    etaLine = new TLine(Mpi0s[ientry],0,Mpi0s[ientry],discriminatorHist->GetMaximum());
                        }
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
			allCanvases->SaveAs(("histograms/Mass-event"+std::to_string(ientry)+".png").c_str());
		}
		if(verbose2){
			duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
			logFile << "	Delta T to finish event: " << duration2 <<  "ms" << endl;
		}
            
	}

	// Finish the log files by including an elapsed time and finally closing the file
	auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
	//logFile << "Total time: " << duration2 << " ms" << endl;
	//logFile << "Time Per Event: " << ( std::clock() - start ) / (double) CLOCKS_PER_SEC / nentries << "ns" << endl;
	logFile.close();
        return 0;
}

