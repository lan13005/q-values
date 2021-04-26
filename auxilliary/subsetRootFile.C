void subsetRootFile() {
   //Get old file, old tree and set top branch address
   TFile *oldfile = new TFile("/d/grid13/ln16/q-values-2/degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_treeFlat_DSelector.root");
   TTree *oldtree = (TTree*)oldfile->Get("degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat");
   Long64_t nentries = oldtree->GetEntries();
   //Create a new file + a clone of old tree in new file
   TFile *newfile[5];
   TTree *newtree[5]; 

   vector<string> pols={"000","045","090","135","AMO"};
   for (int idata=0; idata<(int)pols.size(); ++idata){
       newfile[idata] = new TFile(("/d/grid13/ln16/q-values-2/degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat_pol"+pols[idata]+".root").c_str(),"recreate");
       newtree[idata] = oldtree->CloneTree(0);
   }
   map<int,int> mapPolToIdx;
   mapPolToIdx[0]=0;
   mapPolToIdx[45]=1;
   mapPolToIdx[90]=2;
   mapPolToIdx[135]=3;
   mapPolToIdx[-1]=4;

   int beamAngle;
   double mandelstam_tp;
   oldtree->SetBranchAddress("BeamAngle",&beamAngle);
   oldtree->SetBranchAddress("mandelstam_tp",&mandelstam_tp);

   // THIS IS A BIT TRICKY IF YOU ARE DOING THIS FOR ACCEPTANCE CORRECTION
   // YOU MUST SOMEHOW DO THIS CONSISTENTLY WITH THE GEN TREES ALSO
   double percToKeep=100;
   for (Long64_t i=0;i<nentries; i++) {
        if (i%1000==0)
            cout << "current entry: " << i << endl;
        if (rand()%100<percToKeep){
            oldtree->GetEntry(i);
            if ( (mandelstam_tp<0.5) )
                newtree[mapPolToIdx[beamAngle]]->Fill();
                //newtree[0]->Fill();
        }
   }
   for (int idata=0; idata<(int)pols.size(); ++idata){
       cout << pols[idata] << " has " << newtree[idata]->GetEntries() << " entries" << endl;
       newtree[idata]->AutoSave();
//       delete oldfile[idata];
//       delete newfile[idata];
   }
}
