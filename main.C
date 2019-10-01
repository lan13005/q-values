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
	//TFile* dataFile=new TFile("pi0eta_a0_recotreeFlat_DSelector.root");
	TFile* dataFile=new TFile("pi0eta_datatreeFlat_DSelector.root");
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
	double largest_nentry;
        if (iProcess!=(nProcess-1)) {
            largest_nentry  = (iProcess+1)*batchEntries;
        }
        else {
            largest_nentry = nentries; 
        }

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
                //if ( ientry != eventNumber){ cout << "ientry != eventNumber. The events are in the root tree are out of order or missing!" <<
                //        "\n ientry,eventNumber: " << ientry << ", " << eventNumber << endl; break; }
                //                ***** THE COMBINATIONS COME IN AT RANDOM ORDERS! ***** DOESNT MATTER I THINK
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
        cout << "Size of each array: " << AccWeights.size() << endl;
	
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
	// **** WE FIRST IMPORT THE DATA INTO A CLASS, SAVES AN INTERNAL COPY
        standardizeArray class_cosTheta_X_cms(cosTheta_X_cms,nentries);
        standardizeArray class_phi_X_cms(phi_X_cms,nentries);
        standardizeArray class_cosTheta_eta_gjs(cosTheta_eta_gjs,nentries);
        standardizeArray class_phi_eta_gjs(phi_eta_gjs,nentries);
        standardizeArray class_cosThetaHighestEphotonIneta_gjs(cosThetaHighestEphotonIneta_gjs,nentries);
        standardizeArray class_cosThetaHighestEphotonInpi0_cms(cosThetaHighestEphotonInpi0_cms,nentries);
        standardizeArray class_vanHove_xs(vanHove_xs,nentries);
        standardizeArray class_vanHove_ys(vanHove_ys,nentries);
        standardizeArray class_vanHove_omegas(vanHove_omegas,nentries);
        standardizeArray class_pi0_energies(pi0_energies, nentries);
        standardizeArray class_mandelstam_tps(mandelstam_tps, nentries);
        //standardizeArray class_Mpi0s(Mpi0s,nentries);
        //standardizeArray class_Metas(Metas,nentries);

        // ***** STANDARDIZES THE INTERAL COPY 
        class_cosTheta_X_cms.stdevStandardization();
        class_phi_X_cms.stdevStandardization();
        class_cosTheta_eta_gjs.stdevStandardization();
        class_phi_eta_gjs.stdevStandardization();
        class_cosThetaHighestEphotonIneta_gjs.stdevStandardization();
        class_cosThetaHighestEphotonInpi0_cms.stdevStandardization();
        class_vanHove_xs.stdevStandardization();
        class_vanHove_ys.stdevStandardization();
        class_vanHove_omegas.stdevStandardization();
        class_pi0_energies.stdevStandardization();
        class_mandelstam_tps.stdevStandardization();
        //class_Mpi0s.stdevStandardization();
        //class_Metas.stdevStandardization();
    
        // ******* REPLACES THE ORIGNAL WITH THE STANDARDIZED COPY INSIDE THE ARRAY
        cosTheta_X_cms = class_cosTheta_X_cms.getVector();
        phi_X_cms=class_phi_X_cms.getVector();
        cosTheta_eta_gjs=class_cosTheta_eta_gjs.getVector();
        phi_eta_gjs=class_phi_eta_gjs.getVector();
        cosThetaHighestEphotonIneta_gjs=class_cosThetaHighestEphotonIneta_gjs.getVector();
        cosThetaHighestEphotonInpi0_cms=class_cosThetaHighestEphotonInpi0_cms.getVector();
        vanHove_xs=class_vanHove_xs.getVector();
        vanHove_ys=class_vanHove_xs.getVector();
        vanHove_omegas=class_vanHove_omegas.getVector();
        pi0_energies=class_pi0_energies.getVector();
        mandelstam_tps=class_mandelstam_tps.getVector();
        //Mpi0s=class_Mpi0s.getVector();
        //Metas=class_Metas.getVector();

        map<std::string, std::vector<double>> nameToVec;
        nameToVec["cosTheta_X_cms"] = cosTheta_X_cms;
        nameToVec["phi_X_cms"] = phi_X_cms; 
        nameToVec["cosTheta_eta_gjs"] = cosTheta_eta_gjs; 
        nameToVec["phi_eta_gjs"] = phi_eta_gjs; 
        nameToVec["cosThetaHighestEphotonIneta_gjs"] = cosThetaHighestEphotonIneta_gjs; 
        nameToVec["cosThetaHighestEphotonInpi0_cms"] = cosThetaHighestEphotonInpi0_cms; 
        nameToVec["vanHove_omegas"] = vanHove_omegas; 
        nameToVec["vanHove_xs"] = vanHove_xs; 
        nameToVec["vanHove_ys"] = vanHove_ys; 
        nameToVec["pi0_energies"] = pi0_energies; 
        nameToVec["mandelstam_tps"] = mandelstam_tps; 
        nameToVec["Mpi0s"] = Mpi0s; 
        nameToVec["Metas"] = Metas; 

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

        distSort_kNN distKNN(kDim);
        pair<double,int> newPair;


        // It is much slower to constantly read a map to get the vector rather than just importing it all into a vector first.
        std::vector< std::vector< double > > varVector;
        int numVars=parse.varStringSet.size();
        for (int iVar=0; iVar<numVars; ++iVar){
            varVector.push_back(nameToVec[parse.varStringSet[iVar]]);
        }

        double comboStd; 
        double chiSq;
        ULong64_t flatEntryNumber;
        double bestQ;
        double bestChiSq;

        TFile *resultsFile = new TFile(("logs/results"+to_string(iProcess)+".root").c_str(),"RECREATE");
        TTree* resultsTree = new TTree("resultsTree","results");
        resultsTree->Branch("flatEntryNumber",&flatEntryNumber,"flatEntryNumber/l");
        resultsTree->Branch("qvalue",&bestQ,"qvalue/D");
        resultsTree->Branch("chisq",&bestChiSq,"chisq/D");
        resultsTree->Branch("combostd",&comboStd,"combostd/D");

        std::vector<double> binRange;
        std::vector<double> fitRange;
	TF1* fit;
        TF1* bkgFit;
        TF1* sigFit;
        if (useEta){
            binRange={50,0.35,0.8};
            fitRange={0.425,0.7};
        } 
        else{ 
            binRange={50,0.05,0.25};
            fitRange={0.1,0.17};
        }

        int nBest100Bkg=0;
        int nBest100Sig=0;
        int nBest50Bkg50Sig=0;
        // Getting the scaling of the flat bkg is a bit tricky. Have to count how much bins we have in our fit range. kDim divided by the bins in fit range is the height of the flat function.
        double numBinsInFit = (fitRange[1]-fitRange[0])/((binRange[2]-binRange[1])/binRange[0]);
        double binSize=((binRange[2]-binRange[1])/binRange[0]);
        // this will be a unit gaussian with peak location and width given in peakWidth array. 
        double ampPeakWidth[3]; 
        if (useEta){
            ampPeakWidth[0]=1; ampPeakWidth[1]=peakWidth_eta[0]; ampPeakWidth[2]=peakWidth_eta[1];
        }
        else {
            ampPeakWidth[0]=1; ampPeakWidth[1]=peakWidth_pi0[0]; ampPeakWidth[2]=peakWidth_pi0[1];
        }
        double constAmp = (double)kDim/numBinsInFit;
        double flatAmpInit[3]={constAmp,0,constAmp/2};
        double gausAmp=gausAmplitude(fitRange[0],binSize,(int)numBinsInFit,kDim,&ampPeakWidth[0]); // see header file for info and motivation
        double gausAmpInit[3]={0,gausAmp,gausAmp/2};
        cout << "\n";
        cout << "Num Bins In Fit: " << (int)numBinsInFit << endl;
        cout << "BinSize: " << binSize << endl;
        cout << "GausAmp that would have k total entries in the bins: " << gausAmp << endl;
        cout << "   ^- Corresponds to a total amplitude of: " << gausAmp*1/TMath::Sqrt(2*TMath::Pi())/ampPeakWidth[2] << endl;
        cout << "Flat: Amp1, Amp2, Amp3: " << flatAmpInit[0] << ", " << flatAmpInit[1] << ", " << flatAmpInit[2]<< endl;
        cout << "Gaus: Amp1, Amp2, Amp3: " << gausAmpInit[0] << ", " << gausAmpInit[1] << ", " << gausAmpInit[2]<< endl;
        cout << "\n" << endl;


	// the main loop where we loop through all events in a double for loop to calculate dij. Iterating through all j we can find the k nearest neighbors to event i.
	// Then we can plot the k nearest neighbors in the discriminating distribution which happens to be a double gaussian and flat bkg. Calculate Q-Value from event
	// i's discriminating variable's value. This value will be plugged into the signal PDF and the total PDF. The ratio of these two values are taken which is the Q-Value.
	//logFile << std::fixed << std::setprecision(6);
        for (int ientry=lowest_nentry; ientry<largest_nentry; ientry++){ 
                flatEntryNumber=ientry;
                cumulativeStd stdCalc(kDim);

	        if(verbose) { cout << "Getting next event!\n--------------------------------\n" << endl;  }
	        auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
	        auto duration_beginEvent = std::chrono::high_resolution_clock::now();
	        if(verbose2){logFile << "Starting event " << ientry << "/" << largest_nentry << " ---- Time: " << duration2 << "ms" << endl; }

                // clear data from previous events
	        mapDistToJ.clear();
	        distances.clear();
	        allCanvases->Clear();
	        discriminatorHist = new TH1F("","",binRange[0],binRange[1],binRange[2]);

                for ( int iVar=0; iVar<numVars; ++iVar ){
                    phasePoint1[iVar] = varVector[iVar][ientry];
                }
	        for (int jentry=0; jentry<nentries; jentry++){
                     if ( verbose_outputDistCalc ) { cout << "event i,j = " << ientry << "," << jentry << endl;} 
        
                     for ( int iVar=0; iVar<numVars; ++iVar ){
                         phasePoint2[iVar] = varVector[iVar][jentry];
                     }
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
	        if (distKNN.kNN.size() != kDim){ cout << "size of distKNN is not equal to kDim! size,kDim="<< distKNN.kNN.size() << "," << kDim 
                    << " -- if size is 1 less than kDim it is probably because kDim=nentries and event i cannot be a neighbor to itself" << endl;}
                while ( distKNN.kNN.empty() == false ){
                        newPair = distKNN.kNN.top();
                        distKNN.kNN.pop();
                        if ( useEta ){
                            discriminatorHist->Fill(Metas[newPair.second]);//,AccWeights[newPair.second]);
                            stdCalc.insertValue(Metas[newPair.second]);
                        }
                        else {
                            discriminatorHist->Fill(Mpi0s[newPair.second]);//,AccWeights[newPair.second]);
                            stdCalc.insertValue(Mpi0s[newPair.second]);
                        }
                        //cout << "(" << newPair.first << ", " << newPair.second << ")"; 
                        //cout << endl; 
                }
                comboStd = stdCalc.calcStd();
                
	        //duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
	        //if(verbose2){logFile << "	Filled neighbors: " << duration2 << "ms" << endl;}
	
	        // Building the fit functions. We have to set some parameter limits to try to guide Minuit. We choose a double gaussian since for whatever reason the Meta has a asymmetry
	        // such that there are more events to the left side of the peak. The second gaussian will have lower amplitude and a large sigma with a left shifted mean. We also want to 
	        // fix the amplitudes and constnat background to be strictly positive and choose kDim as the max value ( which is reached only if all values are filled into a single bin)
	        
                
                // We will always do 3 fits per event. In the original paper on Q-factors they use 100% bkg, %100 signal, 50% bkg and 50% signal.
     
                //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
	        Double_t par[numDOFbkg+numDOFsig]; // needed to calculate the qvalue
	        Double_t parBest[numDOFbkg+numDOFsig]; // needed to save the best initialization

                bestQ=0;
                bestChiSq=DBL_MAX;
                int best_iFit=0;

                for (UInt_t iFit=0; iFit<3; ++iFit){
                    // We use a normalized gaussian and a flat function. 
	            fit = new TF1("fit",fitFunc,fitRange[0],fitRange[1],numDOFbkg+numDOFsig);
	            bkgFit = new TF1("bkgFit",background,fitRange[0],fitRange[1],numDOFbkg);
	            sigFit = new TF1("sigFit",signal,fitRange[0],fitRange[1],numDOFsig);

                    // Should use getInitParams.C whenever we get a new dataset to initialize the peak and width of the pi0 and eta
                    if (useEta) { 
	                fit->SetParameters(flatAmpInit[iFit],gausAmpInit[iFit],peakWidth_eta[0],peakWidth_eta[1]);//,20,0.03);
                    }
                    else {
	                fit->SetParameters(flatAmpInit[iFit],gausAmpInit[iFit],peakWidth_pi0[0],peakWidth_pi0[1]);
	            }
                    // we have to enforce the functions to be positive. Easiest way is to make min=0 and max=kDim, the number of neighbors
                    fit->SetParLimits(0,0,kDim); 
                    fit->SetParLimits(1,0,kDim); 

	            //discriminatorHist->Fit("fit","RQBWL"); // WL for weighted histogram fitting
	            discriminatorHist->Fit("fit","RQBN"); // B will enforce the bounds, N will be no draw
                    chiSq = fit->GetChisquare()/(fit->GetNDF());
                    if (verbose2) { logFile << "current ChiSq, best ChiSq: " << chiSq << ", " << bestChiSq << endl; }

                    if (chiSq < bestChiSq){
	                fit->GetParameters(par);
	                bkgFit->SetParameters(par);
	                sigFit->SetParameters(&par[numDOFbkg]);
                        if ( useEta) { 
	                    qvalue=sigFit->Eval(Metas[ientry])/fit->Eval(Metas[ientry]);
                        }
                        else {
	                    qvalue=sigFit->Eval(Mpi0s[ientry])/fit->Eval(Mpi0s[ientry]);
                        }
                        bestChiSq=chiSq;
                        bestQ=qvalue;
                        best_iFit=iFit;
                        for( int iPar=0; iPar<numDOFbkg+numDOFsig; ++iPar){
                            parBest[iPar]=par[iPar];
                        }
                    } 
	        }
                 
                if ( best_iFit == 0){ ++nBest100Bkg; }
                if ( best_iFit == 1){ ++nBest100Sig; }
                if ( best_iFit == 2){ ++nBest50Bkg50Sig; }

	       // Here we draw the histograms that were randomly selected
	       logFile << ientry << " " << bestQ << " " << bestChiSq << " " << comboStd << endl;//"\t" << Metas[ientry] << "\t" << Mpi0s[ientry] << endl;
	       if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end()) {
	           discriminatorHist->SetTitle(("BEST VALS:  QValue="+to_string(bestQ)+"  ChiSq="+to_string(bestChiSq)+"     Std="+to_string(comboStd)+"    iFit="+to_string(best_iFit)).c_str());
	           discriminatorHist->Draw();
                   if ( useEta) { 
               	        etaLine = new TLine(Metas[ientry],0,Metas[ientry],discriminatorHist->GetMaximum());
                   }
                   else { 
               	        etaLine = new TLine(Mpi0s[ientry],0,Mpi0s[ientry],discriminatorHist->GetMaximum());
                   }
	           etaLine->SetLineColor(kOrange);

                   // bkgFit and sigFit has parameters that are set by SetParameters in the for loop above. The SetParameters function is only called when the new chiSq is better than
                   // the old. the "fit" function has parameters that are overwritten every time we do a fit so if we wanted to plot the best fit function we have to save the 
                   // best paramters but for bkgFit and sigFit I think that the best params are inherently saved by the above setup.
	           fit->GetParameters(parBest);
	           bkgFit->SetParameters(parBest);
	           sigFit->SetParameters(&parBest[numDOFbkg]);
	           fit->SetParName(0,"const");
	           fit->SetParName(1,"Amp_Gaus1");
	           fit->SetParName(2,"Mean_Gaus1");
	           fit->SetParName(3,"Sigma_Gaus1");

                   fit->SetLineColor(kRed);
                   fit->Draw("SAME");
  	           bkgFit->SetFillColor(kMagenta);
                   bkgFit->SetLineColor(kMagenta);
  	           bkgFit->SetFillStyle(3004);
  	           bkgFit->Draw("SAME FC");
  	           sigFit->SetFillColor(kBlue);
                   sigFit->SetLineColor(kBlue);
  	           sigFit->SetFillStyle(3005);
  	           sigFit->Draw("SAME FC");

	           etaLine->Draw("same");
	           allCanvases->SaveAs(("histograms/Mass-event"+std::to_string(ientry)+".png").c_str());
	       }

	       if(verbose2){
		    duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
		    logFile << "	Best ChiSq (" << bestChiSq << ") with Q-value: " << bestQ << " -- time: " << duration2 <<  "ms" << endl;
		    logFile << "	Delta T to finish event: " << duration2 <<  "ms" << endl;
	       }
               resultsTree->Fill();
        } 
        resultsFile->cd();
        resultsTree->Write();
        //resultsFile->Write();
        cout << "nentries: " << nentries << endl;
        cout << "Number of times 100% bkg initialization was the best: " << nBest100Bkg << endl;
        cout << "Number of times 100% sig initialization was the best: " << nBest100Sig << endl;
        cout << "Number of times 50/50 bkg/sig initialization was the best: " << nBest50Bkg50Sig << endl;
        
        
	// Finish the log files by including an elapsed time and finally closing the file
	auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start2).count();
	//logFile << "Total time: " << duration2 << " ms" << endl;
	//logFile << "Time Per Event: " << ( std::clock() - start ) / (double) CLOCKS_PER_SEC / nentries << "ns" << endl;
	logFile.close();
        return 0;
}

