// It is not trivial to get enter get the number of combos that pass the event selections into the flat trees. Since the DSelector runs combo by combo we have to save the values and fill them all.
// This would not be hard but just kind of annoying to implement. It is best to just read in the data and count the number of combos that pass the event. This is what we do here.
// This is probably a more straight forward way to implement uniqueness tracking into q-factors and takes away the ambiguity how how to do tracking. i.e. should we track all particles in a combo
// when using combo based variables? i.e. should we only track photons 1,2 and beam when we fill Meta histograms?
// This also makes it easier to transport to other analyses if needed

bool verbose=false;

void addUTWeightsBranch(string rootFileLoc, string rootFileName, string treeName){
    TFile* dataFile = new TFile((rootFileLoc+rootFileName+".root").c_str());
    TTree* dataTree;
    dataFile->GetObject((treeName).c_str(),dataTree);

    Long64_t nentries = (Long64_t)dataTree->GetEntries();
    ULong64_t event;
    dataTree->SetBranchAddress("event",&event);

    TFile *ut_dataFile = TFile::Open((rootFileLoc+rootFileName+"_UTweights.root").c_str(),"RECREATE"); 
    TTree *ut_dataTree = dataTree->CloneTree(-1,"fast"); 
    double uniquenessTrackingWeight;
    TBranch* uniquenessTrackingWeights =ut_dataTree->Branch("uniqunessTrackingWeights",&uniquenessTrackingWeight,"uniquenessTrackingWeights/D");
    
    vector<int> countEvents;

    int count=1;
    //nentries=20;
    ULong64_t previousEvent=-1;
    for(Long64_t ientry=0; ientry<nentries; ientry++)
    {
    	dataTree->GetEntry(ientry);
        if (previousEvent==event){
            ++count;
            if (verbose){
                cout << "event: " << event << " at count " << count << endl;
            }
        }
        else {
            previousEvent=event;
            if(ientry>0){
                countEvents.push_back(count);
            }
            count=1;
            if (verbose) {
                cout << "event: " << event << " at count " << count << endl;
            }
        }
    }
    countEvents.push_back(count); // need to fill it one more time or we dont get all the events

    int totalCounts=0;
    for (auto numCount : countEvents){
        totalCounts += numCount;
        if(verbose){
            cout << numCount << " ";
        }
        for(int iCount=0; iCount<numCount; ++iCount){
            uniquenessTrackingWeight=1.0/numCount;
            uniquenessTrackingWeights->Fill();
        }
    }
    if (verbose){
        cout << endl;
    }
    cout << "Total Counts: " << totalCounts << ", asked for: " << nentries << endl;
    
    ut_dataFile->cd();
    ut_dataTree->Write(); 
}


void getUniquenessWeights(){
    //string rootFileLoc = "/d/grid15/ln16/pi0eta/q-values/";
    string rootFileLoc="/home/lawrence/Desktop/gluex/q-values/";
    string rootFileName = "degALL_bcal_treeFlat_DSelector";
    string treeName = "degALL_bcal_tree_flat";
    addUTWeightsBranch(rootFileLoc, rootFileName, treeName);
}
