// The goal of this code is to convert q-values flat trees into a format that is amptools ready
// Normally one would use tree_to_amptools but q-values code requires flat tree input and thus outputs in flat tree format
//   whereas tree_to_amptools requires tree format
// For Zlm amplitudes we need to split the data into various polarizations which this code can do in you have the polarization angle saved in the flat tree
// Amptools can also accept a background file so you can import the weights appropriately and include them: i.e. AccWeight*qvalue

void flat_to_amptools_pol(string pol){
    // FIRST DEFINE IF YOU WANT TO USE THROWN VALUES OR KIN FIT OR MEASURED VALUES
    // TYPICALL IT SHOULD BE THROWN VALUES FOR ACC AND KIN FOR DATA
    string sourceTag="kin"; // {true, kin, meas}

    map<string,int> mapPolToInt;
    mapPolToInt["000"]=0;
    mapPolToInt["045"]=45;
    mapPolToInt["090"]=90;
    mapPolToInt["135"]=135;
    mapPolToInt["AMO"]=-1;

    // Define reaction and storage variables for trees
    string beam="beam";
    const int cNumFinalState=3;
    int NumFinalState=(int)cNumFinalState;
    int finalStates_PID[cNumFinalState]={14,7,17};
    vector<string> finalState1={"p"};
    vector<string> finalState2={"g1","g2"};
    vector<string> finalState3={"g3","g4"};
    vector<TLorentzVector*> finalState1_p4={new TLorentzVector()};
    vector<TLorentzVector*> finalState2_p4={new TLorentzVector(), new TLorentzVector()};
    vector<TLorentzVector*> finalState3_p4={new TLorentzVector(), new TLorentzVector()};
    vector<vector<TLorentzVector*>> finalStates_p4={finalState1_p4,finalState2_p4,finalState3_p4};
    string accidentalBranch="AccWeight";
    string asbsBranch="weightASBS";
    string qvalueBranch="qvalue";
    string polAngleBranch="BeamAngle";
    string rfTimeBranch="rfTime";
    string tBranch="mandelstam_tp";
    double accidental;
    double qvalue;
    double weightASBS;
    int polAngle;
    double rfTime;
    double mandelstam_tp;
    float beam_e;
    float beam_px;
    float beam_py;
    float beam_pz;
    float targetMass=0.9382719;
    float finalStates_e[cNumFinalState];
    float finalStates_px[cNumFinalState];
    float finalStates_py[cNumFinalState];
    float finalStates_pz[cNumFinalState];
    float weight;
    TLorentzVector* beam_p4=new TLorentzVector();

    // Load input root file/tree
    TTree* tree;
    TFile* file;
    string fileLoc;
    string baseFolder="./results_flat_2017_allPols_40/";
    if (pol=="all")
        fileLoc=baseFolder+"postQVal_flatTree.root";
    else
        fileLoc=baseFolder+"logs/"+pol+"/postQVal_flatTree"+"_"+pol+".root";
    string treeName="degALL_flat_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat";
    //string treeName="degALL_2017_8288_chi13_tree_flat";
    
    //fileLoc="./degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat_pol"+pol+".root"; 
    //treeName="degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat";

    file=TFile::Open(fileLoc.c_str());
    file->GetObject(treeName.c_str(),tree);
    cout << "Opening file: " << fileLoc << endl;
    cout << "Tree name: " << treeName << endl;

    // Turn on branches we care about in the input file
    tree->SetBranchStatus("*",0);
    for (auto s: finalState1)
        tree->SetBranchStatus((s+"_p4_"+sourceTag).c_str(),1);
    for (auto s: finalState2)
        tree->SetBranchStatus((s+"_p4_"+sourceTag).c_str(),1);
    for (auto s: finalState3)
        tree->SetBranchStatus((s+"_p4_"+sourceTag).c_str(),1);
    tree->SetBranchStatus((beam+"_p4_"+sourceTag).c_str(),1);
    tree->SetBranchStatus(polAngleBranch.c_str(),1);
    tree->SetBranchStatus(accidentalBranch.c_str(),1);
    tree->SetBranchStatus(asbsBranch.c_str(),1);
    tree->SetBranchStatus(qvalueBranch.c_str(),1);
    tree->SetBranchStatus(rfTimeBranch.c_str(),1);

    // Set branch addresses of storage variables for the input tree
    for (auto i=0; i<finalState1.size(); ++i)
        tree->SetBranchAddress((finalState1[i]+"_p4_"+sourceTag).c_str(),&finalState1_p4[i]);
    for (auto i=0; i<finalState2.size(); ++i)
        tree->SetBranchAddress((finalState2[i]+"_p4_"+sourceTag).c_str(),&finalState2_p4[i]);
    for (auto i=0; i<finalState3.size(); ++i)
        tree->SetBranchAddress((finalState3[i]+"_p4_"+sourceTag).c_str(),&finalState3_p4[i]);
    tree->SetBranchAddress(accidentalBranch.c_str(),&accidental);
    tree->SetBranchAddress(asbsBranch.c_str(),&weightASBS);
    tree->SetBranchAddress(qvalueBranch.c_str(),&qvalue);
    tree->SetBranchAddress((beam+"_p4_"+sourceTag).c_str(),&beam_p4);
    tree->SetBranchAddress(polAngleBranch.c_str(),&polAngle);
    tree->SetBranchAddress(rfTimeBranch.c_str(),&rfTime);
    tree->SetBranchAddress(tBranch.c_str(),&mandelstam_tp);
    
    // Create output root file/tree
    //string outputFileTag="degALL_flat_8288_chi13_2_tree_flat_tot";
    string outputFileTag=treeName+"_"+"sig";
    //string outputFileTag="degALL_malte_mc_pol035_tree_flat_sb";
    TFile* outFile; 
    TTree* outTree;
    outFile = new TFile(("AmpToolsInput_"+pol+"_"+outputFileTag+".root").c_str(), "RECREATE");
    outTree = new TTree("kin", "kin");
    outTree->Branch("Weight", new float, "Weight/F");
    outTree->Branch("E_Beam", new float, "E_Beam/F");
    outTree->Branch("Px_Beam", new float, "Px_Beam/F");
    outTree->Branch("Py_Beam", new float, "Py_Beam/F");
    outTree->Branch("Pz_Beam", new float, "Pz_Beam/F");
    outTree->Branch("Target_Mass", new float, "Target_Mass/F");
    outTree->Branch("NumFinalState", new int, "NumFinalState/I");
    outTree->Branch("PID_FinalState", new int[cNumFinalState], "PID_FinalState[NumFinalState]/I");
    outTree->Branch("E_FinalState", new float[cNumFinalState], "E_FinalState[NumFinalState]/F");
    outTree->Branch("Px_FinalState", new float[cNumFinalState], "Px_FinalState[NumFinalState]/F");
    outTree->Branch("Py_FinalState", new float[cNumFinalState], "Py_FinalState[NumFinalState]/F");
    outTree->Branch("Pz_FinalState", new float[cNumFinalState], "Pz_FinalState[NumFinalState]/F");

    // set branch addresses for output tree 
    outTree->SetBranchAddress("NumFinalState", &NumFinalState);
    outTree->SetBranchAddress("Target_Mass", &targetMass);
    outTree->SetBranchAddress("PID_FinalState", finalStates_PID);
    outTree->SetBranchAddress("E_FinalState", finalStates_e);
    outTree->SetBranchAddress("Px_FinalState", finalStates_px);
    outTree->SetBranchAddress("Py_FinalState", finalStates_py);
    outTree->SetBranchAddress("Pz_FinalState", finalStates_pz);
    outTree->SetBranchAddress("E_Beam", &beam_e);
    outTree->SetBranchAddress("Px_Beam", &beam_px);
    outTree->SetBranchAddress("Py_Beam", &beam_py);
    outTree->SetBranchAddress("Pz_Beam", &beam_pz);
    outTree->SetBranchAddress("Weight", &weight);
    
    
    // Read input tree and dump into amptools output format
    Long64_t nentries=tree->GetEntries();
    //for (Long64_t ientry=0; ientry<400; ++ientry){
    for (Long64_t ientry=0; ientry<nentries; ++ientry){
        tree->GetEntry(ientry);
        //if (ientry%10000==0)
        //    cout << "current entry: " << ientry << endl;
        cout << ientry << " ";
        int ith_finalState=0;
        for (auto v: finalStates_p4){ 
            TLorentzVector p4_sum;
            for (auto p4: v)
                p4_sum += *p4;
            finalStates_e[ith_finalState]=p4_sum.E();
            finalStates_px[ith_finalState]=p4_sum.Px();
            finalStates_py[ith_finalState]=p4_sum.Py();
            finalStates_pz[ith_finalState]=p4_sum.Pz();
            ++ith_finalState; 
            cout << p4_sum.M() << " ";
        }
        beam_e=beam_p4->E();
        beam_px=beam_p4->Px();
        beam_py=beam_p4->Py();
        beam_pz=beam_p4->Pz();

        //if (weightASBS!=0){
        // For FLAT MC
        weight=accidental*qvalue;
        //weight=weightASBS;
        // FOR DATA TOTAL
        //weight=1;
        // FOR DATA SB
        //weight=1-accidental*qvalue;
        //weight=1-weightASBS;
        //if (abs(weight)<0.0000001) 
        // recall that we have weights that can be negative. We just want weights that are close to zero to be set to some small number
        // IT IS PROBABLY SAFE TO NOT FILL WEIGHT=0 EVEN FOR DATA SB WHICH WOULD BE HIGH CONFIDENCE SIGNAL. THE OUTPUT DATASET
        //      FOR DATA SB WOULD HOPEFULLY BE SMALLER THAN DATA TOTAL SINCE THERE SHOULD BE MORE SIGNAL EVENTS WHICH WOULD HAVE BKG WEIGHT = 0
        //    weight=0.0000001;
        if (abs(weight)>0.0001)
            outTree->Fill();
        //}
        cout << weight << " " << endl;
    }

    // Completed! Just write the files and output the number of entries in each orientation 
    cout << "Total number of entries: " << nentries << endl;
    Long64_t nentries_postSelection=0;
    outFile->Write("kin");
    cout << "Entries with polAngle " << pol << ": " << outTree->GetEntries() << endl;
    nentries_postSelection+=outTree->GetEntries();
    cout << "Post-Selection number of entries: " << nentries_postSelection << endl;
}

void flat_to_amptools(){
    vector<string> pols={"000","045","090","135","AMO"};
    for ( auto pol: pols){
        flat_to_amptools_pol(pol);
    }
}


