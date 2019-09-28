#include "main.h"

//bool verbose = false;
//int kDim=200;
//int numberEventsToSavePerProcess=5;
//bool override_nentries=true;
//Long64_t nentries=100;
//int nProcess=5;
//int seedShift=123125;

bool useEta=true;

//void main(int iProcess, int kDim, int numberEventsToSavePerProcess, int nProcess, int seedShift, Long64_t nentries, bool override_nentries, bool verbose){
int main( int argc, char* argv[] ){
	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
        //int iProcess = std::stoi(argv[0]);
        //int kDim= std::stoi(argv[1]);
        //int numberEventsToSavePerProcess= std::stoi(argv[2]);
        //int nProcess= std::stoi(argv[3]);
        //int seedShift= std::stoi(argv[4]);
        //Long64_t nentries= std::stoll(argv[5]);
        int iProcess = atoi(argv[1]);
        int kDim= atoi(argv[2]);
        int numberEventsToSavePerProcess= atoi(argv[3]);
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
        std::string varString=argv[9];
        cout << "----------------------------" << endl;
        cout << "iProcess: " << iProcess << endl;
        cout << "kDim: " << kDim << endl;
        cout << "numberEventsToSavePerProcess: " << numberEventsToSavePerProcess << endl;
        cout << "nProcess: " << nProcess << endl;
        cout << "seedShift: " << seedShift << endl;
        cout << "nentries: " << nentries << endl; 
        cout << "override_nentries: " << override_nentries << endl;
        cout << "verbose: " << verbose  << endl; 
    
        parseVarString parse(varString);
        parse.parseString();
        for (int iVar=0; iVar<parse.varStringSet.size(); ++iVar){
            cout << "var" << std::to_string(iVar) << ": " << parse.varStringSet[iVar] << endl;
        }
        cout << "----------------------------" << endl;


	// Starting timing
	//clock_t start;
	//double duration;
	auto start2 = std::chrono::high_resolution_clock::now();
	//start = clock();
	
	// setting up some basic root stuff and getting the file and tree
	TFile* dataFile=new TFile("/d/grid15/ln16/pi0eta/092419/pi0eta_a0_recotreeFlat_DSelector.root");
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
        double vanHove_x;
        double vanHove_y;
        double vanHove_omega;
        double pi0_energy;
        double mandelstam_tp;
        ULong64_t eventNumber;
        double uniqueComboID;
        double AccWeight;

	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("Mpi0",&Mpi0);
	dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
        dataTree->SetBranchAddress("cosTheta_X_cm", &cosTheta_X_cm); 
        dataTree->SetBranchAddress("phi_X_cm",&phi_X_cm); 
        dataTree->SetBranchAddress("cosTheta_eta_gj",&cosTheta_eta_gj);
        dataTree->SetBranchAddress("phi_eta_gj",&phi_eta_gj); 
        dataTree->SetBranchAddress("cosThetaHighestEphotonIneta_gj",&cosThetaHighestEphotonIneta_gj);
        dataTree->SetBranchAddress("cosThetaHighestEphotonInpi0_cm",&cosThetaHighestEphotonInpi0_cm);
        dataTree->SetBranchAddress("vanHove_x",&vanHove_x);
        dataTree->SetBranchAddress("vanHove_y",&vanHove_y);
        dataTree->SetBranchAddress("vanHove_omega",&vanHove_omega);
        dataTree->SetBranchAddress("pi0_energy", &pi0_energy);
        dataTree->SetBranchAddress("mandelstam_tp", &mandelstam_tp);
        dataTree->SetBranchAddress("uniqueComboID",&uniqueComboID);
        dataTree->SetBranchAddress("event",&eventNumber);
        dataTree->SetBranchAddress("AccWeight",&AccWeight);

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
	for (int i=0; i<numberEventsToSavePerProcess; i++){
		randomEvent = rand() % (int)batchEntries;
		randomEvent += lowest_nentry;
		selectRandomIdxToSave.insert( randomEvent );
	}

	const int c_nentries = (const int)nentries;

	// importing all the data to RAM instead of reading from root file
	std::vector<double> Metas; Metas.reserve(c_nentries);
        std::vector<double> Mpi0s; Mpi0s.reserve(c_nentries);
        std::vector<double> Mpi0etas; Mpi0etas.reserve(c_nentries);
        std::vector<double> cosTheta_X_cms; cosTheta_X_cms.reserve(c_nentries);
        std::vector<double> phi_X_cms; phi_X_cms.reserve(c_nentries);
        std::vector<double> cosTheta_eta_gjs; cosTheta_eta_gjs.reserve(c_nentries);
        std::vector<double> phi_eta_gjs; phi_eta_gjs.reserve(c_nentries);
        std::vector<double> cosThetaHighestEphotonIneta_gjs; cosThetaHighestEphotonIneta_gjs.reserve(c_nentries);
        std::vector<double> cosThetaHighestEphotonInpi0_cms; cosThetaHighestEphotonInpi0_cms.reserve(c_nentries);
        std::vector<double> pi0_energies; pi0_energies.reserve(c_nentries);
        std::vector<double> mandelstam_tps; mandelstam_tps.reserve(c_nentries);
	std::vector<double> vanHove_xs; vanHove_xs.reserve(c_nentries);
	std::vector<double> vanHove_ys; vanHove_ys.reserve(c_nentries);
	std::vector<double> vanHove_omegas; vanHove_omegas.reserve(c_nentries);
        std::vector<double> AccWeights; AccWeights.reserve(c_nentries);

	// We will use a ientry to keep track of which entries we will get from the tree. We will simply use ientry when filling the arrays.  
	for (int ientry=0; ientry<nentries; ientry++)
	{
		dataTree->GetEntry(ientry);
		Metas.push_back(Meta);
		Mpi0s.push_back(Mpi0);
		Mpi0etas.push_back(Mpi0eta);
		cosTheta_X_cms.push_back(cosTheta_X_cm);
		phi_X_cms.push_back(phi_X_cm);
		cosTheta_eta_gjs.push_back(cosTheta_eta_gj);
		phi_eta_gjs.push_back(phi_eta_gj);
		cosThetaHighestEphotonIneta_gjs.push_back(cosThetaHighestEphotonIneta_gj);	 
		cosThetaHighestEphotonInpi0_cms.push_back(cosThetaHighestEphotonInpi0_cm);	 
                vanHove_xs.push_back(vanHove_x);
                vanHove_ys.push_back(vanHove_y);
                vanHove_omegas.push_back(vanHove_omega);
                pi0_energies.push_back(pi0_energy);
                mandelstam_tps.push_back(mandelstam_tp);
                AccWeights.push_back(AccWeight);
	}
	dataFile->Close();
	
        if ( verbose_outputDistCalc ) {
            cout << "Before standarization" << endl;
            for ( int ientry=0 ; ientry < nentries; ientry++){
                cout << cosTheta_X_cms[ientry] << endl;//"," << phi_X_cms[ientry] << endl;
                //cout << phi_X_cms[ientry] <<endl;// "," << cosTheta_eta_gjs[ientry] << endl;
            }
        }
	//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
	//logFile << "Time marker to load data: " << duration << "ms" << endl;

	// outputting the results before and after standardizeArray will show that it works
	// for(auto& cosTheta_X_cm1 : cosTheta_X_cms){ cout << cosTheta_X_cm1 << endl; }
	
        standardizeArray class_cosTheta_X_cms(cosTheta_X_cms,nentries);
        standardizeArray class_phi_X_cms(phi_X_cms,nentries);
        standardizeArray class_cosTheta_eta_gjs(cosTheta_eta_gjs,nentries);
        standardizeArray class_phi_eta_gjs(phi_eta_gjs,nentries);
        standardizeArray class_Mpi0s(Mpi0s,nentries);
        standardizeArray class_cosThetaHighestEphotonIneta_gjs(cosThetaHighestEphotonIneta_gjs,nentries);
        standardizeArray class_cosThetaHighestEphotonInpi0_cms(cosThetaHighestEphotonInpi0_cms,nentries);
        standardizeArray class_vanHove_xs(vanHove_xs,nentries);
        standardizeArray class_vanHove_ys(vanHove_ys,nentries);
        standardizeArray class_vanHove_omegas(vanHove_omegas,nentries);
        standardizeArray class_pi0_energies(pi0_energies, nentries);
        standardizeArray class_mandelstam_tps(mandelstam_tps, nentries);

        class_cosTheta_X_cms.stdevStandardization();
        class_phi_X_cms.stdevStandardization();
        class_cosTheta_eta_gjs.stdevStandardization();
        class_phi_eta_gjs.stdevStandardization();
        class_Mpi0s.stdevStandardization();
        class_cosThetaHighestEphotonIneta_gjs.stdevStandardization();
        class_cosThetaHighestEphotonInpi0_cms.stdevStandardization();
        class_vanHove_xs.stdevStandardization();
        class_vanHove_ys.stdevStandardization();
        class_vanHove_omegas.stdevStandardization();
        class_pi0_energies.stdevStandardization();
        class_mandelstam_tps.stdevStandardization();

        cosTheta_X_cms = class_cosTheta_X_cms.getVector();
        phi_X_cms=class_phi_X_cms.getVector();
        cosTheta_eta_gjs=class_cosTheta_eta_gjs.getVector();
        phi_eta_gjs=class_phi_eta_gjs.getVector();
        Mpi0s=class_Mpi0s.getVector();
        cosThetaHighestEphotonIneta_gjs=class_cosThetaHighestEphotonIneta_gjs.getVector();
        cosThetaHighestEphotonInpi0_cms=class_cosThetaHighestEphotonInpi0_cms.getVector();
        vanHove_xs=class_vanHove_xs.getVector();
        vanHove_ys=class_vanHove_xs.getVector();
        vanHove_omegas=class_vanHove_omegas.getVector();
        pi0_energies=class_pi0_energies.getVector();
        mandelstam_tps=class_mandelstam_tps.getVector();

        map<std::string, std::vector<double>> nameToVec;
        nameToVec["cosTheta_X_cms"] = cosTheta_X_cms;
        nameToVec["phi_X_cms"] = phi_X_cms; 
        nameToVec["cosTheta_eta_gjs"] = cosTheta_eta_gjs; 
        nameToVec["phi_eta_gjs"] = phi_eta_gjs; 
        nameToVec["Mpi0s"] = Mpi0s; 
        nameToVec["cosThetaHighestEphotonIneta_gjs"] = cosThetaHighestEphotonIneta_gjs; 
        nameToVec["cosThetaHighestEphotonInpi0_cms"] = cosThetaHighestEphotonInpi0_cms; 
        nameToVec["vanHove_omegas"] = vanHove_omegas; 
        nameToVec["vanHove_xs"] = vanHove_xs; 
        nameToVec["vanHove_ys"] = vanHove_ys; 
        nameToVec["pi0_energies"] = pi0_energies; 
        nameToVec["mandelstam_tps"] = mandelstam_tps; 


        //int nameMapCounter=0;
	//for(auto elem : nameToVec){
        //    phasePoint1[nameMapCounter] = elem.second[ientry];
        //    ++nameMapCounter;
	//}

        if ( verbose_outputDistCalc ) {
	    cout << "After standardization" << endl;
            for ( int ientry=0 ; ientry < nentries; ientry++){
                cout << cosTheta_X_cms[ientry] << endl;//"," << phi_X_cms[ientry] << endl;
                //cout << phi_X_cms[ientry] << endl;//"," << cosTheta_eta_gjs[ientry] << endl;
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


        // It is much slower to constantly read a map to get the vector rather than just importing it all into a vector first.
        std::vector< std::vector< double > > varVector;
        int numVars=parse.varStringSet.size();
        for (int iVar=0; iVar<numVars; ++iVar){
            varVector.push_back(nameToVec[parse.varStringSet[iVar]]);
        }

        double comboStd; 
	// the main loop where we loop through all events in a double for loop to calculate dij. Iterating through all j we can find the k nearest neighbors to event i.
	// Then we can plot the k nearest neighbors in the discriminating distribution which happens to be a double gaussian and flat bkg. Calculate Q-Value from event
	// i's discriminating variable's value. This value will be plugged into the signal PDF and the total PDF. The ratio of these two values are taken which is the Q-Value.
	//logFile << std::fixed << std::setprecision(6);
	for (int ientry=lowest_nentry; ientry<largest_nentry; ientry++){ 
                cumulativeStd std(kDim);
                //std = new cumulativeStd(kDim);
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

                for ( int iVar=0; iVar<numVars; ++iVar ){
                    phasePoint1[iVar] = varVector[iVar][ientry];
                }
                
		//phasePoint1[0] = cosTheta_X_cms[ientry];
		//phasePoint1[1] = phi_X_cms[ientry];
		//phasePoint1[2] = cosTheta_eta_gjs[ientry];
		//phasePoint1[3] = phi_eta_gjs[ientry];
		//phasePoint1[4] = cosThetaHighestEphotonIneta_gjs[ientry];
		//phasePoint1[5] = cosThetaHighestEphotonInpi0_cms[ientry];
                //phasePoint1[6] = mandelstam_tps[ientry];
                //phasePoint1[6] = pi0_energies[ientry];
                //phasePoint1[6] = Mpi0s[ientry];
                //phasePoint1[6] = vanHove_omegas[ientry];
                //phasePoint1[6] = vanHove_xs[ientry];
                //phasePoint1[7] = vanHove_ys[ientry];
		for (int jentry=0; jentry<nentries; jentry++){
                        if ( verbose_outputDistCalc ) { cout << "event i,j = " << ientry << "," << jentry << endl;} 
        
                        for ( int iVar=0; iVar<numVars; ++iVar ){
                            phasePoint2[iVar] = varVector[iVar][jentry];
                        }
                        
		        //phasePoint2[0] = cosTheta_X_cms[jentry];
			//phasePoint2[1] = phi_X_cms[jentry];
			//phasePoint2[2] = cosTheta_eta_gjs[jentry];
			//phasePoint2[3] = phi_eta_gjs[jentry];
			//phasePoint2[4] = cosThetaHighestEphotonIneta_gjs[jentry];
			//phasePoint2[5] = cosThetaHighestEphotonInpi0_cms[jentry];
                        //phasePoint1[6] = mandelstam_tps[ientry];
                        //phasePoint1[6] = pi0_energies[ientry];
			//phasePoint2[6] = Mpi0s[jentry];
                        //phasePoint2[6] = vanHove_omegas[jentry];
                        //phasePoint2[6] = vanHove_xs[jentry];
                        //phasePoint2[7] = vanHove_ys[jentry];
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
                        std.insertValue(Metas[newPair.second]);
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
                comboStd = std.calcStd();
                
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
		//fit->SetParName(5,"Mean_Gaus2");
		fit->SetParName(5,"Sigma_Gaus2");
                if (useEta) { 
		    fit->SetParameters(5,50,0.545,0.01,20,0.03);
		    fit->SetParLimits(0, 0, kDim);
		    fit->SetParLimits(1, 0, kDim);
		    //fit->SetParLimits(2,0.53,0.56);
		    //fit->SetParLimits(3,0.009,0.1);
		    fit->SetParLimits(2,0.50,0.56);
		    fit->SetParLimits(3,0.005,0.1);
		    //fit->SetParLimits(4, 0, 0);// kDim);
		    fit->SetParLimits(4, 0, kDim);//0.49,0.55);
		    fit->SetParLimits(5, 0.02, 0.04);//0.02, 0.06);
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
		logFile << ientry << " " << qvalue << " " << fit->GetChisquare()<< " " << comboStd << endl;//"\t" << Metas[ientry] << "\t" << Mpi0s[ientry] << endl;
		if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end()) {
			discriminatorHist->SetTitle(("QValue="+to_string(qvalue)+"  ChiSq="+to_string(fit->GetChisquare())+"     Std="+to_string(comboStd)).c_str());
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

