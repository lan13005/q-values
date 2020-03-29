void makeProjections(){
        TFile* dataFile=new TFile("degALL_BCAL_treeFlat_DSelector.root");
        TTree *dataTree;
        dataFile->GetObject("degALL_BCALtree_flat",dataTree);
        TCanvas *allCanvases = new TCanvas("anyCanvas","",1440,900);
        TCanvas *allCanvasesX = new TCanvas("anyCanvasX","",1440,900);
        TCanvas *allCanvasesY = new TCanvas("anyCanvasY","",1440,900);

        std::vector<double> binRange;
        std::vector<double> binRange2;
        binRange={100,0.35,0.8};
        binRange2={100,0.08,0.205};

	TH2D* discriminatorHist;
	TH1D* discriminatorHistX;
	TH1D* discriminatorHistY;
        discriminatorHist = new TH2D("2DHist","",binRange2[0],binRange2[1],binRange2[2],binRange[0],binRange[1],binRange[2]);
        discriminatorHistX = new TH1D("1DHistX","",binRange2[0],binRange2[1],binRange2[2]);
        discriminatorHistY = new TH1D("1DHistY","",binRange[0],binRange[1],binRange[2]);
	TH1D* projection;

	double Meta;
        double Mpi0;

        dataTree->SetBranchAddress("Meta",&Meta);
        dataTree->SetBranchAddress("Mpi0",&Mpi0);

	nentries=dataTree->GetEntries();
        for (Long64_t ientry=0; ientry<nentries; ientry++)
        {
                dataTree->GetEntry(ientry);
        	discriminatorHist->Fill(Mpi0,Meta);
        	discriminatorHistX->Fill(Mpi0);
        	discriminatorHistY->Fill(Meta);
        }

	allCanvases->cd();
	discriminatorHist->Draw("COLZ");
	allCanvases->SaveAs("projections/2DMassHist.png");
	allCanvases->Clear();
	discriminatorHistX->Draw("COLZ");
	allCanvases->SaveAs("projections/1DMassHistX.png");
	allCanvases->Clear();
	discriminatorHistY->Draw("COLZ");
	allCanvases->SaveAs("projections/1DMassHistY.png");

	allCanvasesX->Divide(3,3);
	allCanvasesY->Divide(3,3);

	for (int i=0; i<9; ++i){
		cout << "Drawing on pad: " << i << endl;
		cout << "Using bin: " << 11*i << " to " << 11*(i+1) << endl;
		allCanvasesX->cd(i+1);
		if ( i < 9 ) {
			discriminatorHist->ProjectionX(("_px"+to_string(i)).c_str(),11*i,11*(i+1),"")->Draw();
		}
		else {
			discriminatorHist->ProjectionX(("_px"+to_string(i)).c_str(),11*i,-1,"")->Draw();
		}
		allCanvasesY->cd(i+1);
		if ( i < 9 ) {
			discriminatorHist->ProjectionY(("_py"+to_string(i)).c_str(),11*i,11*(i+1),"")->Draw();
		}
		else {
			discriminatorHist->ProjectionY(("_py"+to_string(i)).c_str(),11*i,-1,"")->Draw();
		}
	}
	allCanvasesX->SaveAs("projections/xproj.png");
	allCanvasesY->SaveAs("projections/yproj.png");
}



















