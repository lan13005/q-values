// This program loads data from postQ root file and draws a bunch of histograms

#include "makeDiagnosticHists.h"


void makeDiagnosticHists(){
	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
        gStyle->SetOptStat(0);

        ///////////////////////////////////////
        // LOAD THE DATA FROM THE POST-Q ROOT FILE
        ///////////////////////////////////////

	string postQFileName = "logs"+runTag+"/"+fileTag+"/postQVal_flatTree_"+fileTag+".root";	
	TFile* postQFile = TFile::Open((postQFileName).c_str()); 
        TTree* dataTree;
	postQFile->GetObject(rootTreeName.c_str(),dataTree);

        // SETUP VARIALBES TO LOAD QFACTOR RELATED VARIABLES
        int nentries=dataTree->GetEntries();
	double qvalue;
        double conjugate_qvalue;
	double bestNLL;
	double worstNLL;
        double worst_qvalue;
        double qvalueBS_std;
        double eff_nentries;
        int kDim1; // kDim1 is a regular int whereas kDim needs to be a const int to make an array out of it
        int neighbors[kDim];
	dataTree->SetBranchAddress("qvalue",&qvalue);
	dataTree->SetBranchAddress("NLLBest",&bestNLL);
	dataTree->SetBranchAddress("NLLWorst",&worstNLL);
	dataTree->SetBranchAddress("worst_qvalue",&worst_qvalue);
	dataTree->SetBranchAddress("qvalueBS_std",&qvalueBS_std);
        dataTree->SetBranchAddress("eff_nentries",&eff_nentries);
        dataTree->SetBranchAddress("kDim",&kDim1);
        dataTree->SetBranchAddress("neighbors",&neighbors);
	std::vector< double > qvalues;
	std::vector< double > bestNLLs;
	std::vector< double > worstNLLs;
	std::vector< double > worst_qvalues;
	std::vector< double > qvalueBS_stds;
	std::vector< double > eff_nentrieses;
        std::vector< int > kDims;
        std::vector< std::array<int,kDim> > neighborses;
        std::array<int,kDim> copyableNeighbors;


        // SETUP VARIALBES TO LOAD INPUT DATA RELATED VARIABLES
        double AccWeight;
        double sbWeight;
        double sbWeightPi0;
        double sbWeightEta;
	double Meta;
	double Mpi0;
	double Mpi0g1;
	double Mpi0g2;
	double Mpi0eta;
	double cosTheta_X_cm;
	double cosTheta_eta_gj;
	double phi_eta_gj;
	double Meta_meas;
	double Mpi0_meas;
	double Mpi0eta_meas;
	double cosTheta_X_cm_meas;
	double cosTheta_eta_gj_meas;
	double phi_eta_gj_meas;
        double utWeight;

        dataTree->SetBranchAddress(s_accWeight.c_str(),&AccWeight);
        if(!s_sbWeight.empty()){
            dataTree->SetBranchAddress(s_sbWeight.c_str(),&sbWeight);
            dataTree->SetBranchAddress("weightBSpi0",&sbWeightPi0);
            dataTree->SetBranchAddress("weightBSeta",&sbWeightEta);
        }
        else{
            sbWeight=1;
            sbWeightPi0=1;
            sbWeightEta=1;
        }
	dataTree->SetBranchAddress("Meta_meas",&Meta_meas);
	dataTree->SetBranchAddress("Mpi0_meas",&Mpi0_meas);
	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("Mpi0",&Mpi0);

	dataTree->SetBranchAddress("Mpi0eta_meas",&Mpi0eta_meas);
        dataTree->SetBranchAddress("cosTheta_X_cm_meas", &cosTheta_X_cm_meas); 
        dataTree->SetBranchAddress("cosTheta_eta_gj_meas",&cosTheta_eta_gj_meas);
        dataTree->SetBranchAddress("phi_eta_gj_meas",&phi_eta_gj_meas); 
	dataTree->SetBranchAddress("Mpi0g1",&Mpi0g1);
	dataTree->SetBranchAddress("Mpi0g2",&Mpi0g2);
	dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
        dataTree->SetBranchAddress("cosTheta_X_cm", &cosTheta_X_cm); 
        dataTree->SetBranchAddress("cosTheta_eta_gj",&cosTheta_eta_gj);
        dataTree->SetBranchAddress("phi_eta_gj",&phi_eta_gj); 
        if(!s_utBranch.empty()){
            dataTree->SetBranchAddress(s_utBranch.c_str(),&utWeight);
        }
        else {
            utWeight=1;
        }
        std::vector<double> AccWeights; AccWeights.reserve(nentries); 
        std::vector<double> utWeights; utWeights.reserve(nentries);
        std::vector<double> Metas_meas; Metas_meas.reserve(nentries); 
        std::vector<double> Mpi0s_meas; Mpi0s_meas.reserve(nentries); 
        std::vector<double> Mpi0etas_meas; Mpi0etas_meas.reserve(nentries); 
	std::vector<double> cosTheta_X_cms_meas; cosTheta_X_cms_meas.reserve(nentries);
	std::vector<double> cosTheta_eta_gjs_meas; cosTheta_eta_gjs_meas.reserve(nentries);
	std::vector<double> phi_eta_gjs_meas; phi_eta_gjs_meas.reserve(nentries);
        std::vector<double> Metas; Metas.reserve(nentries); 
        std::vector<double> Mpi0s; Mpi0s.reserve(nentries); 
        std::vector<double> Mpi0g1s; Mpi0g1s.reserve(nentries); 
        std::vector<double> Mpi0g2s; Mpi0g2s.reserve(nentries); 
        std::vector<double> Mpi0etas; Mpi0etas.reserve(nentries); 
	std::vector<double> cosTheta_X_cms; cosTheta_X_cms.reserve(nentries);
	std::vector<double> cosTheta_eta_gjs; cosTheta_eta_gjs.reserve(nentries);
	std::vector<double> phi_eta_gjs; phi_eta_gjs.reserve(nentries);
	std::vector< double > sbWeights; sbWeights.reserve(nentries);
	std::vector< double > sbWeightsPi0; sbWeightsPi0.reserve(nentries);
	std::vector< double > sbWeightsEta; sbWeightsEta.reserve(nentries);
	

	for (int ientry=0; ientry<nentries; ientry++)
	{
		dataTree->GetEntry(ientry);

                AccWeights.push_back(AccWeight);
		sbWeights.push_back(sbWeight);
		sbWeightsPi0.push_back(sbWeightPi0);
		sbWeightsEta.push_back(sbWeightEta);
                utWeights.push_back(utWeight);

		qvalues.push_back(qvalue);
		bestNLLs.push_back(bestNLL);
		worstNLLs.push_back(worstNLL);
                worst_qvalues.push_back(worst_qvalue);
                qvalueBS_stds.push_back(qvalueBS_std);
                eff_nentrieses.push_back(eff_nentries);
                kDims.push_back(kDim1);
                std::copy(std::begin(copyableNeighbors),std::end(copyableNeighbors),std::begin(neighbors));
                neighborses.push_back(copyableNeighbors);

		Metas_meas.push_back(Meta_meas);
		Mpi0s_meas.push_back(Mpi0_meas);
		Mpi0etas_meas.push_back(Mpi0eta_meas);
		cosTheta_X_cms_meas.push_back(cosTheta_X_cm_meas);
		cosTheta_eta_gjs_meas.push_back(cosTheta_eta_gj_meas);
		phi_eta_gjs_meas.push_back(phi_eta_gj_meas);
		Metas.push_back(Meta);
		Mpi0s.push_back(Mpi0);
		Mpi0g1s.push_back(Mpi0g1);
		Mpi0g2s.push_back(Mpi0g2);
		Mpi0etas.push_back(Mpi0eta);
		cosTheta_X_cms.push_back(cosTheta_X_cm);
		cosTheta_eta_gjs.push_back(cosTheta_eta_gj);
		phi_eta_gjs.push_back(phi_eta_gj);
	}
        cout << "Imported all the tree data into arrays" << endl;


        // -----------------------------------
        // Defining Histograms
        // -----------------------------------
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);

        // The array holds {kin,measured} values. We plot them both even though phase space for the kNN is dependent on only one of them. Might be insightful
	TH1F* Meta_tot[2];
	TH1F* Meta_sig[2];
	TH1F* Meta_bkg[2];
	TH1F* Mpi0_tot[2];
	TH1F* Mpi0_sig[2];
	TH1F* Mpi0_bkg[2];
	TH1F* Mpi0eta_tot[2];
	TH1F* Mpi0eta_sig[2];
	TH1F* Mpi0eta_bkg[2];
	TH1F* cosThetaEta_GJ_tot[2];
	TH1F* cosThetaEta_GJ_sig[2];
	TH1F* cosThetaEta_GJ_bkg[2];
	TH1F* cosThetaX_CM_tot[2];
	TH1F* cosThetaX_CM_sig[2];
	TH1F* cosThetaX_CM_bkg[2];
	TH1F* phiEta_GJ_tot[2];
	TH1F* phiEta_GJ_sig[2];
	TH1F* phiEta_GJ_bkg[2];
	TH1F* Meta_sig_sb[2];
	TH1F* Meta_bkg_sb[2];
	TH1F* Mpi0_sig_sb[2];
	TH1F* Mpi0_bkg_sb[2];
	TH1F* Mpi0eta_sig_sb[2];
	TH1F* Mpi0eta_bkg_sb[2];
	TH1F* cosThetaEta_GJ_sig_sb[2];
	TH1F* cosThetaEta_GJ_bkg_sb[2];
	TH1F* cosThetaX_CM_sig_sb[2];
	TH1F* cosThetaX_CM_bkg_sb[2];
	TH1F* phiEta_GJ_sig_sb[2];
	TH1F* phiEta_GJ_bkg_sb[2];

	TH1F *Mpi0g_tot = new TH1F( "Mpi0g_tot", "M(#pi^{0}g); M(#pi^{0}g) GeV; Events/0.001 GeV", 350, 0, 3.5 );
	TH1F *Mpi0g_sig = new TH1F( "Mpi0g_sig", "M(#pi^{0}g); M(#pi^{0}g) GeV; Events/0.001 GeV", 350, 0, 3.5 );
	TH1F *Mpi0g_bkg = new TH1F( "Mpi0g_bkg", "M(#pi^{0}g); M(#pi^{0}g) GeV; Events/0.001 GeV", 350, 0, 3.5 );
	TH1F *Mpi0g_sig_sb = new TH1F( "Mpi0g_sig_sb", "M(#pi^{0}g); M(#pi^{0}g) GeV; Events/0.001 GeV", 350, 0, 3.5 );
	TH1F *Mpi0g_bkg_sb = new TH1F( "Mpi0g_bkg_sb", "M(#pi^{0}g); M(#pi^{0}g) GeV; Events/0.001 GeV", 350, 0, 3.5 );

        double gj_lowMass=0.7;
        double gj_uppMass=2.0;
        const int nBins=50;
        double gj_width=(gj_uppMass-gj_lowMass)/nBins;
	TH1F* cosThetaEta_GJ_tot_Mbinned[nBins];
	TH1F* cosThetaEta_GJ_sig_Mbinned[nBins];
	TH1F* cosThetaEta_GJ_bkg_Mbinned[nBins];
	TH1F* cosThetaEta_GJ_sig_sb_Mbinned[nBins];
	TH1F* cosThetaEta_GJ_bkg_sb_Mbinned[nBins];
        for (int imass=0; imass<nBins; ++imass){
            double lmass=imass*gj_width+gj_lowMass;
            double umass=(imass+1)*gj_width+gj_lowMass;
            char clowMass[20];
            char cuppMass[20];
            sprintf(clowMass,"%.2lf",lmass);
            sprintf(cuppMass,"%.2lf",umass);
            string slowMass=clowMass;
            string suppMass=cuppMass;
	    cosThetaEta_GJ_tot_Mbinned[imass] = new TH1F( ("cosThetaEta_GJ_tot-"+to_string(imass)).c_str(), 
                    ("cos(#theta) GJ of #eta; cos(#theta) of #eta "+slowMass+"<Mpi0eta<"+suppMass+"; Events/0.02 GeV").c_str(), 100, -1, 1 );
	    cosThetaEta_GJ_sig_Mbinned[imass] = new TH1F( ("cosThetaEta_GJ_sig-"+to_string(imass)).c_str(), 
                    ("cos(#theta) GJ of #eta; cos(#theta) of #eta "+slowMass+"<Mpi0eta<"+suppMass+"; Events/0.02 GeV").c_str(), 100, -1, 1 );
	    cosThetaEta_GJ_bkg_Mbinned[imass] = new TH1F( ("cosThetaEta_GJ_bkg-"+to_string(imass)).c_str(), 
                    ("cos(#theta) GJ of #eta; cos(#theta) of #eta "+slowMass+"<Mpi0eta<"+suppMass+"; Events/0.02 GeV").c_str(), 100, -1, 1 );
	    cosThetaEta_GJ_sig_sb_Mbinned[imass] = new TH1F( ("cosThetaEta_GJ_sig_sb-"+to_string(imass)).c_str(),
                    ("cos(#theta) GJ of #eta; cos(#theta) of #eta "+slowMass+"<Mpi0eta<"+suppMass+"; Events/0.02 GeV").c_str(), 100, -1, 1 );
	    cosThetaEta_GJ_bkg_sb_Mbinned[imass] = new TH1F( ("cosThetaEta_GJ_bkg_sb-"+to_string(imass)).c_str(),
                    ("cos(#theta) GJ of #eta; cos(#theta) of #eta "+slowMass+"<Mpi0eta<"+suppMass+"; Events/0.02 GeV").c_str(), 100, -1, 1 );
        }

        cout << "Defined all histograms" << endl;

	for (int i=0; i<2; i++){
		string tag="";
		if (i%2==0){ tag="meas"; }
		else { tag="kin"; }
		// WE SHOULD TRY TO MATCH THE BINS USED TO FIT, USED IN THE Q-VALUE CALCULATION, AND HERE. USING IT IN THE FIT AND DURING GRAPHS WILL MAKE FOR GOOD COMPARISION
		// AND USING FOR IN THE FIT AND THE Q-VALUE WOULD ALSO BE GOOD SINCE THE FIT COULD CHANGE DEPENDING ON THE BINNING
		Meta_tot[i] = new TH1F( ("Meta_tot_"+tag).c_str(), "M(#eta); M(#eta) GeV; Events/0.003 GeV", 200, 0.25, 0.85);//Events/0.002 GeV", 300, 0.25, 0.85 );
		Meta_sig[i] = new TH1F( ("Meta_sig_"+tag).c_str(), "M(#eta); M(#eta) GeV; Events/0.003 GeV", 200, 0.25, 0.85);//Events/0.002 GeV", 300, 0.25, 0.85 );
		Meta_bkg[i] = new TH1F( ("Meta_bkg_"+tag).c_str(), "M(#eta); M(#eta) GeV; Events/0.003 GeV", 200, 0.25, 0.85);//Events/0.002 GeV", 300, 0.25, 0.85 );
		Mpi0_tot[i] = new TH1F( ("Mpi0_tot_"+tag).c_str(), "M(#pi^{0}); M(#pi^{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
		Mpi0_sig[i] = new TH1F( ("Mpi0_sig_"+tag).c_str(), "M(#pi^{0}); M(#pi^{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
		Mpi0_bkg[i] = new TH1F( ("Mpi0_bkg_"+tag).c_str(), "M(#pi^{0}); M(#pi^{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
		Mpi0eta_tot[i] = new TH1F( ("Mpi0eta_tot_"+tag).c_str(), "M(#pi^{0}#eta); M(#pi^{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
		Mpi0eta_sig[i] = new TH1F( ("Mpi0eta_sig_"+tag).c_str(), "M(#pi^{0}#eta); M(#pi^{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
		Mpi0eta_bkg[i] = new TH1F( ("Mpi0eta_bkg_"+tag).c_str(), "M(#pi^{0}#eta); M(#pi^{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
		cosThetaEta_GJ_tot[i] = new TH1F( ("cosThetaEta_GJ_tot_"+tag).c_str(), "cos(#theta) GJ of #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, -1, 1 );
		cosThetaEta_GJ_sig[i] = new TH1F( ("cosThetaEta_GJ_sig_"+tag).c_str(), "cos(#theta) GJ of #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, -1, 1 );
		cosThetaEta_GJ_bkg[i] = new TH1F( ("cosThetaEta_GJ_bkg_"+tag).c_str(), "cos(#theta) GJ of #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, -1, 1 );
		cosThetaX_CM_tot[i] = new TH1F( ("cosThetaX_CM_tot_"+tag).c_str(), "cos(#theta) of CM #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, 0.9, 1 );
		cosThetaX_CM_sig[i] = new TH1F( ("cosThetaX_CM_sig_"+tag).c_str(), "cos(#theta) of CM #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, 0.9, 1 );
		cosThetaX_CM_bkg[i] = new TH1F( ("cosThetaX_CM_bkg_"+tag).c_str(), "cos(#theta) of CM #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, 0.9, 1 );
		phiEta_GJ_tot[i] = new TH1F( ("phiEta_GJ_tot_"+tag).c_str(), "#phi of #eta;#phi of #eta ;Events/0.02 GeV", 90, -180, 180 );
		phiEta_GJ_sig[i] = new TH1F( ("phiEta_GJ_sig_"+tag).c_str(), "#phi of #eta;#phi of #eta ;Events/0.02 GeV", 90, -180, 180 );
		phiEta_GJ_bkg[i] = new TH1F( ("phiEta_GJ_bkg_"+tag).c_str(), "#phi of #eta;#phi of #eta ;Events/0.02 GeV", 90, -180, 180 );

		Meta_sig_sb[i] = new TH1F( ("Meta_sig_"+tag+"_sb").c_str(), "M(#eta); M(#eta) GeV; Events/0.003 GeV", 200, 0.25, 0.85);//Events/0.002 GeV", 300, 0.25, 0.85 );
		Meta_bkg_sb[i] = new TH1F( ("Meta_bkg_"+tag+"_sb").c_str(), "M(#eta); M(#eta) GeV; Events/0.003 GeV", 200, 0.25, 0.85);//Events/0.002 GeV", 300, 0.25, 0.85 );
		Mpi0_sig_sb[i] = new TH1F( ("Mpi0_sig_"+tag+"_sb").c_str(), "M(#pi^{0}); M(#pi^{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
		Mpi0_bkg_sb[i] = new TH1F( ("Mpi0_bkg_"+tag+"_sb").c_str(), "M(#pi^{0}); M(#pi^{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
		Mpi0eta_sig_sb[i] = new TH1F( ("Mpi0eta_sig_"+tag+"_sb").c_str(), "M(#pi^{0}#eta); M(#pi^{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
		Mpi0eta_bkg_sb[i] = new TH1F( ("Mpi0eta_bkg_"+tag+"_sb").c_str(), "M(#pi^{0}#eta); M(#pi^{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
		cosThetaEta_GJ_sig_sb[i] = new TH1F( ("cosThetaEta_GJ_sig_"+tag+"_sb").c_str(), "cos(#theta) GJ of #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, -1, 1 );
		cosThetaEta_GJ_bkg_sb[i] = new TH1F( ("cosThetaEta_GJ_bkg_"+tag+"_sb").c_str(), "cos(#theta) GJ of #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, -1, 1 );
		cosThetaX_CM_sig_sb[i] = new TH1F( ("cosThetaX_CM_sig_"+tag+"_sb").c_str(), "cos(#theta) of CM #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, 0.9, 1 );
		cosThetaX_CM_bkg_sb[i] = new TH1F( ("cosThetaX_CM_bkg_"+tag+"_sb").c_str(), "cos(#theta) of CM #eta; cos(#theta) of #eta; Events/0.02 GeV", 100, 0.9, 1 );
		phiEta_GJ_sig_sb[i] = new TH1F( ("phiEta_GJ_sig_"+tag+"_sb").c_str(), "#phi of #eta;#phi of #eta ;Events/0.02 GeV", 90, -180, 180 );
		phiEta_GJ_bkg_sb[i] = new TH1F( ("phiEta_GJ_bkg_"+tag+"_sb").c_str(), "#phi of #eta;#phi of #eta ;Events/0.02 GeV", 90, -180, 180 );
	}
        cout << "Initialized all histograms" << endl;


        // We can now fill the histograms, properly weighted
	double sigWeight;
	double totWeight;
	double bkgWeight;
	double sigWeight_sb;
	double bkgWeight_sb;
	double sigWeight_sbPi0;
	double bkgWeight_sbPi0;
	double sigWeight_sbEta;
	double bkgWeight_sbEta;
        double baseWeight; 
	for (int ientry=0; ientry<nentries; ientry++){
		qvalue = qvalues[ientry];
		conjugate_qvalue = 1-qvalue;
                AccWeight=AccWeights[ientry];
                sbWeight=sbWeights[ientry];
                sbWeightPi0=sbWeightsPi0[ientry];
                sbWeightEta=sbWeightsEta[ientry];
                if (weightingScheme==""){ baseWeight=1; }
                if (weightingScheme=="as"){ baseWeight=AccWeight; }
                ////////////////////////////////////
                // Multiply q-factor and accidetal weights if requested
                ////////////////////////////////////
		sigWeight = qvalue*baseWeight;
		totWeight = baseWeight;
		bkgWeight = conjugate_qvalue*baseWeight;

                //////////////////////////////
                // Multply sideband and accidental weights if requested 
                //////////////////////////////
                sigWeight_sb = baseWeight*sbWeight;
                bkgWeight_sb = baseWeight*(1-sbWeight);
                sigWeight_sbPi0 = baseWeight*sbWeightPi0;
                bkgWeight_sbPi0 = baseWeight*(1-sbWeightPi0);
                sigWeight_sbEta = baseWeight*sbWeightEta;
                bkgWeight_sbEta = baseWeight*(1-sbWeightEta);
            
                //////////////////////////////
                // Include tracking weights
                //////////////////////////////
		sigWeight = sigWeight*utWeights[ientry];
		totWeight = totWeight*utWeights[ientry]; 
		bkgWeight = bkgWeight*utWeights[ientry];

                sigWeight_sb = sigWeight_sb*utWeights[ientry];
                bkgWeight_sb = bkgWeight_sb*utWeights[ientry];
                sigWeight_sbPi0 = sigWeight_sbPi0*utWeights[ientry];
                bkgWeight_sbPi0 = bkgWeight_sbPi0*utWeights[ientry];
                sigWeight_sbEta = sigWeight_sbEta*utWeights[ientry];
                bkgWeight_sbEta = bkgWeight_sbEta*utWeights[ientry];

                ////////////////////////////////////
                // Fill histograms
                ////////////////////////////////////
		cosThetaEta_GJ_sig[0]->Fill(cosTheta_eta_gjs_meas[ientry], sigWeight);
		cosThetaEta_GJ_tot[0]->Fill(cosTheta_eta_gjs_meas[ientry], totWeight);
		cosThetaEta_GJ_bkg[0]->Fill(cosTheta_eta_gjs_meas[ientry], bkgWeight);
		cosThetaX_CM_sig[0]->Fill(cosTheta_X_cms_meas[ientry], sigWeight);
		cosThetaX_CM_tot[0]->Fill(cosTheta_X_cms_meas[ientry], totWeight);
		cosThetaX_CM_bkg[0]->Fill(cosTheta_X_cms_meas[ientry], bkgWeight);
		Meta_sig[0]->Fill(Metas_meas[ientry],sigWeight);
		Meta_tot[0]->Fill(Metas_meas[ientry],totWeight);
		Meta_bkg[0]->Fill(Metas_meas[ientry],bkgWeight);
		cosThetaEta_GJ_sig_sb[0]->Fill(cosTheta_eta_gjs_meas[ientry], sigWeight_sb);
		cosThetaEta_GJ_bkg_sb[0]->Fill(cosTheta_eta_gjs_meas[ientry], bkgWeight_sb);
		cosThetaX_CM_sig_sb[0]->Fill(cosTheta_X_cms_meas[ientry], sigWeight_sb);
		cosThetaX_CM_bkg_sb[0]->Fill(cosTheta_X_cms_meas[ientry], bkgWeight_sb);
		Meta_sig_sb[0]->Fill(Metas_meas[ientry],sigWeight_sbPi0);
		Meta_bkg_sb[0]->Fill(Metas_meas[ientry],bkgWeight_sbPi0);

		cosThetaEta_GJ_sig[1]->Fill(cosTheta_eta_gjs[ientry], sigWeight);
		cosThetaEta_GJ_tot[1]->Fill(cosTheta_eta_gjs[ientry], totWeight);
		cosThetaEta_GJ_bkg[1]->Fill(cosTheta_eta_gjs[ientry], bkgWeight);
		cosThetaX_CM_sig[1]->Fill(cosTheta_X_cms[ientry], sigWeight);
		cosThetaX_CM_tot[1]->Fill(cosTheta_X_cms[ientry], totWeight);
		cosThetaX_CM_bkg[1]->Fill(cosTheta_X_cms[ientry], bkgWeight);
		Meta_sig[1]->Fill(Metas[ientry],sigWeight);
		Meta_tot[1]->Fill(Metas[ientry],totWeight);
		Meta_bkg[1]->Fill(Metas[ientry],bkgWeight);
		cosThetaEta_GJ_sig_sb[1]->Fill(cosTheta_eta_gjs[ientry], sigWeight_sb);
		cosThetaEta_GJ_bkg_sb[1]->Fill(cosTheta_eta_gjs[ientry], bkgWeight_sb);
		cosThetaX_CM_sig_sb[1]->Fill(cosTheta_X_cms[ientry], sigWeight_sb);
		cosThetaX_CM_bkg_sb[1]->Fill(cosTheta_X_cms[ientry], bkgWeight_sb);
		Meta_sig_sb[1]->Fill(Metas[ientry],sigWeight_sbPi0);
		Meta_bkg_sb[1]->Fill(Metas[ientry],bkgWeight_sbPi0);

		Mpi0g_sig->Fill(Mpi0g1s[ientry],sigWeight);
		Mpi0g_tot->Fill(Mpi0g1s[ientry],totWeight);
		Mpi0g_bkg->Fill(Mpi0g1s[ientry],bkgWeight);
		Mpi0g_sig_sb->Fill(Mpi0g1s[ientry],sigWeight_sb);
		Mpi0g_bkg_sb->Fill(Mpi0g1s[ientry],bkgWeight_sb);

		Mpi0g_sig->Fill(Mpi0g2s[ientry],sigWeight);
		Mpi0g_tot->Fill(Mpi0g2s[ientry],totWeight);
		Mpi0g_bkg->Fill(Mpi0g2s[ientry],bkgWeight);
		Mpi0g_sig_sb->Fill(Mpi0g2s[ientry],sigWeight_sb);
		Mpi0g_bkg_sb->Fill(Mpi0g2s[ientry],bkgWeight_sb);

		Mpi0_sig[0]->Fill(Mpi0s_meas[ientry],sigWeight);
		Mpi0_tot[0]->Fill(Mpi0s_meas[ientry],totWeight);
		Mpi0_bkg[0]->Fill(Mpi0s_meas[ientry], bkgWeight);
		Mpi0_sig_sb[0]->Fill(Mpi0s_meas[ientry],sigWeight_sbEta);
		Mpi0_bkg_sb[0]->Fill(Mpi0s_meas[ientry],bkgWeight_sbEta);

		Mpi0_sig[1]->Fill(Mpi0s[ientry],sigWeight);
		Mpi0_tot[1]->Fill(Mpi0s[ientry],totWeight);
		Mpi0_bkg[1]->Fill(Mpi0s[ientry], bkgWeight);
		Mpi0_sig_sb[1]->Fill(Mpi0s[ientry],sigWeight_sbEta);
		Mpi0_bkg_sb[1]->Fill(Mpi0s[ientry],bkgWeight_sbEta);

		phiEta_GJ_sig[0]->Fill(phi_eta_gjs_meas[ientry], sigWeight);
		phiEta_GJ_tot[0]->Fill(phi_eta_gjs_meas[ientry], totWeight);
		phiEta_GJ_bkg[0]->Fill(phi_eta_gjs_meas[ientry], bkgWeight);
		Mpi0eta_sig[0]->Fill(Mpi0etas_meas[ientry],sigWeight);
		Mpi0eta_tot[0]->Fill(Mpi0etas_meas[ientry],totWeight);
		Mpi0eta_bkg[0]->Fill(Mpi0etas_meas[ientry], bkgWeight);
		phiEta_GJ_sig_sb[0]->Fill(phi_eta_gjs_meas[ientry], sigWeight_sb);
		phiEta_GJ_bkg_sb[0]->Fill(phi_eta_gjs_meas[ientry], bkgWeight_sb);
		Mpi0eta_sig_sb[0]->Fill(Mpi0etas_meas[ientry],sigWeight_sb);
		Mpi0eta_bkg_sb[0]->Fill(Mpi0etas_meas[ientry],bkgWeight_sb);

		phiEta_GJ_sig[1]->Fill(phi_eta_gjs[ientry], sigWeight);
		phiEta_GJ_tot[1]->Fill(phi_eta_gjs[ientry], totWeight);
		phiEta_GJ_bkg[1]->Fill(phi_eta_gjs[ientry], bkgWeight);
		Mpi0eta_sig[1]->Fill(Mpi0etas[ientry],sigWeight);
		Mpi0eta_tot[1]->Fill(Mpi0etas[ientry],totWeight);
		Mpi0eta_bkg[1]->Fill(Mpi0etas[ientry], bkgWeight);
		phiEta_GJ_sig_sb[1]->Fill(phi_eta_gjs[ientry], sigWeight_sb);
		phiEta_GJ_bkg_sb[1]->Fill(phi_eta_gjs[ientry], bkgWeight_sb);
		Mpi0eta_sig_sb[1]->Fill(Mpi0etas[ientry],sigWeight_sb);
		Mpi0eta_bkg_sb[1]->Fill(Mpi0etas[ientry],bkgWeight_sb);

                int massBin=floor( (Mpi0etas[ientry]-gj_lowMass)/gj_width );
                if ((massBin>=0)*(massBin<nBins)){
		    cosThetaEta_GJ_tot_Mbinned[massBin]->Fill(cosTheta_eta_gjs[ientry], totWeight);
		    cosThetaEta_GJ_sig_Mbinned[massBin]->Fill(cosTheta_eta_gjs[ientry], sigWeight);
		    cosThetaEta_GJ_bkg_Mbinned[massBin]->Fill(cosTheta_eta_gjs[ientry], bkgWeight);
		    cosThetaEta_GJ_sig_sb_Mbinned[massBin]->Fill(cosTheta_eta_gjs[ientry], sigWeight_sb);
		    cosThetaEta_GJ_bkg_sb_Mbinned[massBin]->Fill(cosTheta_eta_gjs[ientry], bkgWeight_sb);
                }
        }
        cout << "Made the histograms" << endl;

	// HERE WE WILL JUST DRAW SOME OF THE HISTOGRAMS WITH THE BKG FILLED IN TO SEE THEIR CONTRIBUTION
	makeStackedHist(Mpi0g_tot,Mpi0g_sig,Mpi0g_bkg,Mpi0g_sig_sb,Mpi0g_bkg_sb,"Mpi0gkin", "diagnosticPlots"+runTag+"/"+fileTag);
	for (int i=0; i<2; i++){
		string tag="";
		if (i%2==0){ tag="meas"; }
		else { tag="kin"; } 
		makeStackedHist(Meta_tot[i],Meta_sig[i],Meta_bkg[i],Meta_sig_sb[i], Meta_bkg_sb[i],"Meta"+tag, "diagnosticPlots"+runTag+"/"+fileTag);
		makeStackedHist(Mpi0_tot[i],Mpi0_sig[i],Mpi0_bkg[i],Mpi0_sig_sb[i], Mpi0_bkg_sb[i],"Mpi0"+tag, "diagnosticPlots"+runTag+"/"+fileTag);
		makeStackedHist(Mpi0eta_tot[i],Mpi0eta_sig[i],Mpi0eta_bkg[i],Mpi0eta_sig_sb[i],Mpi0eta_bkg_sb[i],"Mpi0eta"+tag, "diagnosticPlots"+runTag+"/"+fileTag);
		makeStackedHist(cosThetaEta_GJ_tot[i],cosThetaEta_GJ_sig[i],cosThetaEta_GJ_bkg[i],cosThetaEta_GJ_sig_sb[i],cosThetaEta_GJ_bkg_sb[i],"cosThetaEta_GJ"+tag, "diagnosticPlots"+runTag+"/"+fileTag);	
		makeStackedHist(cosThetaX_CM_tot[i],cosThetaX_CM_sig[i],cosThetaX_CM_bkg[i],cosThetaX_CM_sig_sb[i],cosThetaX_CM_bkg_sb[i],"cosThetaX_CM"+tag, "diagnosticPlots"+runTag+"/"+fileTag);	
		makeStackedHist(phiEta_GJ_tot[i],phiEta_GJ_sig[i],phiEta_GJ_bkg[i],phiEta_GJ_sig_sb[i],phiEta_GJ_bkg_sb[i],"phiEta_GJ"+tag, "diagnosticPlots"+runTag+"/"+fileTag);	
	}
        for (int imass=0; imass<nBins; ++imass){
	    makeStackedHist(cosThetaEta_GJ_tot_Mbinned[imass],cosThetaEta_GJ_sig_Mbinned[imass],cosThetaEta_GJ_bkg_Mbinned[imass],cosThetaEta_GJ_sig_sb_Mbinned[imass],cosThetaEta_GJ_bkg_sb_Mbinned[imass],"cosThetaEta_GJ-"+to_string(imass), "diagnosticPlots"+runTag+"/"+fileTag);	
        }

        cout << "FINIHSED!"<<endl;

        // we can also directly save all the histograms to a root file
	TFile* dataFile3 = new TFile(("diagnosticPlots"+runTag+"/"+fileTag+"/postQVal_hists_"+fileTag+".root").c_str(),"RECREATE");
	for (int i=0; i<2; i++){
        	Meta_bkg[i]->Write();
        	Meta_sig[i]->Write();
        	Meta_tot[i]->Write();
        	Mpi0_sig[i]->Write();
        	Mpi0_bkg[i]->Write();
        	Mpi0_tot[i]->Write();
        	Mpi0eta_sig[i]->Write();
        	Mpi0eta_bkg[i]->Write();
        	Mpi0eta_tot[i]->Write();
                cosThetaEta_GJ_tot[i]->Write();
                cosThetaEta_GJ_sig[i]->Write();
                cosThetaEta_GJ_bkg[i]->Write();
                cosThetaX_CM_tot[i]->Write();
                cosThetaX_CM_bkg[i]->Write();
                cosThetaX_CM_sig[i]->Write();
                phiEta_GJ_tot[i]->Write();
                phiEta_GJ_sig[i]->Write();
                phiEta_GJ_bkg[i]->Write();

        	Meta_bkg_sb[i]->Write();
        	Meta_sig_sb[i]->Write();
        	Mpi0_sig_sb[i]->Write();
        	Mpi0_bkg_sb[i]->Write();
        	Mpi0eta_sig_sb[i]->Write();
        	Mpi0eta_bkg_sb[i]->Write();
                cosThetaEta_GJ_sig_sb[i]->Write();
                cosThetaEta_GJ_bkg_sb[i]->Write();
                cosThetaX_CM_bkg_sb[i]->Write();
                cosThetaX_CM_sig_sb[i]->Write();
                phiEta_GJ_sig_sb[i]->Write();
                phiEta_GJ_bkg_sb[i]->Write();

	}
	Mpi0g_sig->Write();
	Mpi0g_bkg->Write();
	Mpi0g_tot->Write();
	Mpi0g_sig_sb->Write();
	Mpi0g_bkg_sb->Write();
        for (int imass=0; imass<nBins; ++imass){
	    cosThetaEta_GJ_tot_Mbinned[imass]->Write();
	    cosThetaEta_GJ_sig_Mbinned[imass]->Write();
	    cosThetaEta_GJ_bkg_Mbinned[imass]->Write();
	    cosThetaEta_GJ_sig_sb_Mbinned[imass]->Write();
	    cosThetaEta_GJ_bkg_sb_Mbinned[imass]->Write();
        }
}

