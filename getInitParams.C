// README
// The purpose of this code is to fit the distribution of the discriminating variables and extract the fit parameters
// We can then pass these fitted parameters, scale them, and use them as initializations for the q-values

#include <string>
#include "helperFuncs.h"


string rootFileLoc="/d/grid15/ln16/pi0eta/q-values/degALL_fcal_treeFlat_DSelector.root";
string rootTreeName="degALL_fcal_tree_flat";
string fileTag="fcal";
string weightingScheme="as"; // "" or "as*bs"
string s_accWeight="AccWeight";
string s_discrimVar="Meta";
string s_sideBandVar="Mpi0";

void getInitParams(){
	gStyle->SetOptFit(111);
	gStyle->SetOptStat(0);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);
	TTree *dataTree;

	static const int nDetectors = 1;
	for (int i=0; i<nDetectors; ++i){
		cout << "File Name: " << rootFileLoc << endl;
		cout << "Tree Name: " << rootTreeName << endl;
		TFile* dataFile=new TFile(rootFileLoc.c_str());
		dataFile->GetObject(rootTreeName.c_str(),dataTree);

    		ofstream logFile_discrimVar1;
    		ofstream logFile_discrimVar2;
    		logFile_discrimVar1.open(("fitResults/"+fileTag+"/discrimVarFit_toMain_"+fileTag+".txt").c_str());
    		logFile_discrimVar2.open(("fitResults/"+fileTag+"/pi0Fit_toMain_"+fileTag+".txt").c_str());

    		TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
        	bool isUniqueEtaB;
        	bool isUniquePi0B;
        	bool isUniquePi0EtaB;
		double discrimVar;
		double sideBandVar;
		double Mpi0eta;
		double AccWeight;
		dataTree->SetBranchAddress(s_discrimVar.c_str(),&discrimVar);
		dataTree->SetBranchAddress(s_sideBandVar.c_str(),&sideBandVar);
		dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
		dataTree->SetBranchAddress(s_accWeight.c_str(),&AccWeight);
        	dataTree->SetBranchAddress("isNotRepeated_eta",&isUniqueEtaB);
        	dataTree->SetBranchAddress("isNotRepeated_pi0",&isUniquePi0B);
        	dataTree->SetBranchAddress("isNotRepeated_pi0eta",&isUniquePi0EtaB);
		long long nentries=dataTree->GetEntries();

        	TH1F *massHistEta; 
        	TH1F *massHistPi0; 
        	TH1F *massHistPi0Eta; 
		TF1* fit;
        	TF1* bkgFit;
        	TF1* sigFit;
        	double par[numDOFbkg+numDOFsig];
		
        	cout <<"Initialized" << endl;
        	string namePar[numDOFbkg+numDOFsig] = {"bern01","bern11","amp","mass","sigma"};//,"ampRatio","sigmaRatio"};
                int meanIndex=3;
                int widthIndex=4;


		allCanvases->Clear();
		std::vector<double> binRangeEta;
		std::vector<double> fitRangeEta;
		double binWidthEta;
		std::vector<double> binRangePi0;
		std::vector<double> fitRangePi0;
		double binWidthPi0;

		double integralBKG;
		double integralSIG;

		// ///////////////////////////////////////////////////
		// START ETA FIT
		// ///////////////////////////////////////////////////
                double _SET_DiscrimBinLower = 0.25;
                double _SET_DiscrimBinUpper = 0.85;
                double _SET_DiscrimBinNum = 200;
                double _SET_DiscrimFitLower = 0.4;
                double _SET_DiscrimFitUpper = 0.7;

		binRangeEta={_SET_DiscrimBinNum,_SET_DiscrimBinLower,_SET_DiscrimBinUpper};
		fitRangeEta={_SET_DiscrimFitLower,_SET_DiscrimFitUpper};
		binRangePi0={200,0.005,0.25};
		fitRangePi0={0.09,0.18};

		binWidthEta=(binRangeEta[2]-binRangeEta[1])/binRangeEta[0];
		binWidthPi0=(binRangePi0[2]-binRangePi0[1])/binRangePi0[0];


		fit = new TF1("fit",fitFunc,fitRangeEta[0],fitRangeEta[1],numDOFbkg+numDOFsig);
		bkgFit = new TF1("bkgFit",background,fitRangeEta[0],fitRangeEta[1],numDOFbkg);
		sigFit = new TF1("sigFit",signal,fitRangeEta[0],fitRangeEta[1],numDOFsig);

		fit->SetParameters(3000,2000,15000,0.547,0.003);//,3,5);
                fit->SetParLimits(0,0,4000);
                fit->SetParLimits(1,0,3000);
                fit->SetParLimits(2,0,20000);
                fit->SetParLimits(3,0.5,0.6);
                fit->SetParLimits(4,0,0.05);

		massHistEta = new TH1F("","",binRangeEta[0],binRangeEta[1],binRangeEta[2]);
		massHistPi0 = new TH1F("","",binRangePi0[0],binRangePi0[1],binRangePi0[2]);
		massHistPi0Eta = new TH1F("","", 350, 0, 3.5);
		cout << "Initialized for a specific mass (eta/pi0) fit" << endl;
		
		double sbWeight;
		for (int ientry=0; ientry<nentries; ientry++)
		{
			dataTree->GetEntry(ientry);
                        getSBWeight(sideBandVar,&sbWeight,weightingScheme);
                        double weight;
                        if (weightingScheme==""){ weight=1; }
                        if (weightingScheme=="as"){ weight=AccWeight; }
                        if (weightingScheme=="as*bs"){ weight=AccWeight*sbWeight; }
                	if ( isUniqueEtaB ) {
		        	massHistEta->Fill(discrimVar,weight);//*sbWeight);
			}
                	if ( isUniquePi0B ) {
		        	massHistPi0->Fill(sideBandVar,weight); /////////////////////////////////////////// NOT WEIGHTED SINCE WE WONT BE ABLE TO FIT IT PROPERLY 
			}
                	if ( isUniquePi0EtaB ) {
				massHistPi0Eta->Fill(Mpi0eta,weight);//*sbWeight);
			}
		}
		cout << "Filled all entries into histogram for a specific fit" << endl;

		massHistPi0Eta->Draw("HIST");
		massHistPi0Eta->GetXaxis()->SetTitleSize(0.04);
		massHistPi0Eta->GetYaxis()->SetTitleSize(0.04);
		allCanvases->SaveAs(("fitResults/"+fileTag+"/Mpi0eta_fit_"+fileTag+".png").c_str());
		allCanvases->Clear();
		
		massHistEta->Fit("fit","RQB"); // B will enforce the bounds
		fit->GetParameters(par);
                cout << "Setting par[" << widthIndex << "] = width to be positive" << endl;
                if(par[widthIndex] < 0){ par[widthIndex] = abs(par[widthIndex]); }
                fit->SetParameters(par);
		bkgFit->SetParameters(par);
		sigFit->SetParameters(&par[numDOFbkg]);
		massHistEta->SetTitle(";M(#eta) (GeV)");
		logFile_discrimVar1 << "#nentries " << nentries << endl;
		
		for (int iPar=0; iPar<numDOFbkg+numDOFsig; ++iPar){
		        logFile_discrimVar1 << namePar[iPar] << " " << par[iPar] << endl;
		}
		
                cout << "IF SIDEBAND SUBTRACTING ON THIS, THE INTEGRALS WILL OBVIOULSY BE WEIRD" << endl;
                cout << "----------------------------------------------------------------------" << endl;
		integralBKG = bkgFit->Integral(fitRangeEta[0],fitRangeEta[1]);
		integralSIG = sigFit->Integral(fitRangeEta[0],fitRangeEta[1]);
		cout << "IntegralBKG before scaling: " << integralBKG << endl;
		cout << "IntegralSIG before scaling: " << integralSIG << endl;
		integralBKG *= 1/binWidthEta;
		integralSIG *= 1/binWidthEta;

		cout << "IntegralBKG after scaling: " << integralBKG << endl;
		cout << "IntegralSIG after scaling: " << integralSIG << endl;
		cout << "nentries: " << nentries << endl;

		logFile_discrimVar1 << "#integralBKG " << integralBKG << endl;
		logFile_discrimVar1 << "#integralSIG " << integralSIG << endl;

		double eventRatioSigToBkg = integralSIG/integralBKG;
		logFile_discrimVar1 << "#eventRatioSigToBkg " << eventRatioSigToBkg << endl;

		int nSig=3;
		//double weightedSigma = par[2]/(par[2]+par[2]*par[5])*par[4]+(par[2]*par[5])/(par[2]+par[2]*par[5])*par[6]*par[4];
                double weightedSigma = par[widthIndex];
                cout << "Checking if our calculated integrals for the bkg and signal add up to the total" << endl;
		logFile_discrimVar1 << "#weightedSigma " << weightedSigma << endl;
		TLine *line = new TLine(par[meanIndex]-nSig*weightedSigma,0,par[meanIndex]-nSig*weightedSigma,massHistEta->GetMaximum());
		line->SetLineColor(kMagenta);
		double integralBKG_nsig = bkgFit->Integral(par[meanIndex]-nSig*weightedSigma,par[meanIndex]+nSig*weightedSigma)/binWidthEta;
		double integralSIG_nsig = sigFit->Integral(par[meanIndex]-nSig*weightedSigma,par[meanIndex]+nSig*weightedSigma)/binWidthEta;
		logFile_discrimVar1 << "#integralBKG_nSig " << integralBKG_nsig << endl;
		logFile_discrimVar1 << "#integralSIG_nSig " << integralSIG_nsig << endl;
		Int_t binx1 = massHistEta->GetXaxis()->FindBin(par[meanIndex]-nSig*weightedSigma);
		Int_t binx2 = massHistEta->GetXaxis()->FindBin(par[meanIndex]+nSig*weightedSigma);
		logFile_discrimVar1 << "#-- 3sigma BinLower, BinUpper = " << binx1 << ", " << binx2 << endl;
		logFile_discrimVar1 << "#-- Actual counts in the 3sigma range " << massHistEta->Integral(binx1,binx2) << endl;
		logFile_discrimVar1 << "#-- integralBKG_nSig+integralSIG_nSig " << integralBKG_nsig+integralSIG_nsig << endl;
		if ( abs(1-massHistEta->Integral(binx1,binx2)/(integralBKG_nsig+integralSIG_nsig)) < 0.05 ) {
			logFile_discrimVar1 << "#--There is agreement within 5%" << endl;
		}
		else {
			logFile_discrimVar1 << "#--Percent off " << abs(1-massHistEta->Integral(binx1,binx2)/(integralBKG_nsig+integralSIG_nsig)) << endl;
		}
		double purity = integralSIG_nsig/(integralBKG_nsig+integralSIG_nsig);
		logFile_discrimVar1 << "#purity " << purity << endl;

		
		massHistEta->Draw();
		line->DrawLine(par[meanIndex]-nSig*weightedSigma,0,par[meanIndex]-nSig*weightedSigma,massHistEta->GetMaximum());
		line->DrawLine(par[meanIndex]+nSig*weightedSigma,0,par[meanIndex]+nSig*weightedSigma,massHistEta->GetMaximum());
		massHistEta->GetXaxis()->SetTitleSize(0.04);
		massHistEta->GetYaxis()->SetTitleSize(0.04);
		massHistEta->SetTitle(("Peak: "+to_string(par[meanIndex])+"    width: "+to_string(weightedSigma)).c_str());
		allCanvases->SaveAs(("fitResults/"+fileTag+"/"+s_discrimVar+"_fit_"+fileTag+".png").c_str());
		cout << "Saved for a specific fit!" << endl;
        	

		// ///////////////////////////////////////////////////
		// START Pi0 FIT
		// ///////////////////////////////////////////////////
		fit = new TF1("fit",fitFunc,fitRangePi0[0],fitRangePi0[1],numDOFbkg+numDOFsig);
		bkgFit = new TF1("bkgFit",background,fitRangePi0[0],fitRangePi0[1],numDOFbkg);
		sigFit = new TF1("sigFit",signal,fitRangePi0[0],fitRangePi0[1],numDOFsig);

		fit->SetParameters(4000,3000,14000,0.134,0.015);//,5,2);
                fit->SetParLimits(0,0,7000);
                fit->SetParLimits(1,0,7000);
                fit->SetParLimits(2,0,20000);
                fit->SetParLimits(3,0.13,0.14);
                fit->SetParLimits(4,0.001,0.01);


		massHistPi0->Fit("fit","RQB"); // B will enforce the bounds
		fit->GetParameters(par);
                cout << "Setting par[" << widthIndex << "] = width to be positive" << endl;
                if(par[widthIndex] < 0){ par[widthIndex] = abs(par[widthIndex]); }
                fit->SetParameters(par);
		bkgFit->SetParameters(par);
		sigFit->SetParameters(&par[numDOFbkg]);
		massHistPi0->SetTitle(";M(#pi^{0}) (GeV)");
		logFile_discrimVar2 << "#nentries " << nentries << endl;
		
		for (int iPar=0; iPar<numDOFbkg+numDOFsig; ++iPar){
		        logFile_discrimVar2 << namePar[iPar] << " " << par[iPar] << endl;
		}
		
                cout << "IF SIDEBAND SUBTRACTING ON THIS, THE INTEGRALS WILL OBVIOULSY BE WEIRD" << endl;
                cout << "----------------------------------------------------------------------" << endl;
		integralBKG = bkgFit->Integral(fitRangePi0[0],fitRangePi0[1]);
		integralSIG = sigFit->Integral(fitRangePi0[0],fitRangePi0[1]);
		cout << "IntegralBKG before scaling: " << integralBKG << endl;
		cout << "IntegralSIG before scaling: " << integralSIG << endl;
		integralBKG *= 1/binWidthPi0;
		integralSIG *= 1/binWidthPi0;
		logFile_discrimVar2 << "#integralBKG " << integralBKG << endl;
		logFile_discrimVar2 << "#integralSIG " << integralSIG << endl;

		cout << "IntegralBKG after scaling: " << integralBKG << endl;
		cout << "IntegralSIG after scaling: " << integralSIG << endl;
		cout << "nentries: " << nentries << endl;

		//weightedSigma = par[2]/(par[2]+par[2]*par[5])*par[4]+(par[2]*par[5])/(par[2]+par[2]*par[5])*par[6]*par[4];
                weightedSigma=par[widthIndex];
		logFile_discrimVar2 << "#weightedSigma: " << weightedSigma << endl;

                cout << "Checking if our calculated integrals for the bkg and signal add up to the total" << endl;
		integralBKG_nsig = bkgFit->Integral(par[meanIndex]-nSig*weightedSigma,par[meanIndex]+nSig*weightedSigma)/binWidthPi0;
		integralSIG_nsig = sigFit->Integral(par[meanIndex]-nSig*weightedSigma,par[meanIndex]+nSig*weightedSigma)/binWidthPi0;
		logFile_discrimVar2 << "#integralBKG_nSig " << integralBKG_nsig << endl;
		logFile_discrimVar2 << "#integralSIG_nSig " << integralSIG_nsig << endl;
		binx1 = massHistPi0->GetXaxis()->FindBin(par[meanIndex]-nSig*weightedSigma);
		binx2 = massHistPi0->GetXaxis()->FindBin(par[meanIndex]+nSig*weightedSigma);
		logFile_discrimVar2 << "#-- 3sigma BinLower, BinUpper = " << binx1 << ", " << binx2 << endl;
		logFile_discrimVar2 << "#-- Actual counts in the 3sigma range " << massHistPi0->Integral(binx1,binx2) << endl;
		logFile_discrimVar2 << "#-- integralBKG_nSig+integralSIG_nSig " << integralBKG_nsig+integralSIG_nsig << endl;
		if ( abs(1-massHistPi0->Integral(binx1,binx2)/(integralBKG_nsig+integralSIG_nsig)) < 0.05 ) {
			logFile_discrimVar2 << "#--There is agreement within 5%" << endl;
		}
		else {
			logFile_discrimVar2 << "#--Percent off " << abs(1-massHistPi0->Integral(binx1,binx2)/(integralBKG_nsig+integralSIG_nsig)) << endl;
		}
		purity = integralSIG_nsig/(integralBKG_nsig+integralSIG_nsig);
		logFile_discrimVar2 << "#purity " << purity << endl;
		
		massHistPi0->Draw();
		line->DrawLine(par[meanIndex]-nSig*weightedSigma,0,par[meanIndex]-nSig*weightedSigma,massHistEta->GetMaximum());
		line->DrawLine(par[meanIndex]+nSig*weightedSigma,0,par[meanIndex]+nSig*weightedSigma,massHistEta->GetMaximum());
		massHistPi0->GetXaxis()->SetTitleSize(0.04);
		massHistPi0->GetYaxis()->SetTitleSize(0.04);
		massHistPi0->SetTitle(("Peak: "+to_string(par[meanIndex])+"    width: "+to_string(weightedSigma)).c_str());
		//TLine *line = new TLine( sbRL, 0, sbRL, massHistPi0->GetMaximum());
		//line->SetLineColor(kRed);
		//line->Draw("SAME");
		//line->DrawLine( sbRR, 0, sbRR, massHistPi0->GetMaximum());
		//line->DrawLine( sigL, 0, sigL, massHistPi0->GetMaximum());
		//line->DrawLine( sigR, 0, sigR, massHistPi0->GetMaximum());
		//line->DrawLine( sbLL, 0, sbLL, massHistPi0->GetMaximum());
		//line->DrawLine( sbLR, 0, sbLR, massHistPi0->GetMaximum());

		allCanvases->SaveAs(("fitResults/"+fileTag+"/"+s_sideBandVar+"_fit_"+fileTag+".png").c_str());
		cout << "Saved for a specific fit!" << endl;
	}
}




































