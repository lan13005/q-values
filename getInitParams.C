// README
// The purpose of this code is to fit the distribution of the discriminating variables and extract the fit parameters
// We can then pass these fitted parameters, scale them, and use them as initializations for the q-values

#include <string>
#include "helperFuncs.h"


// This takes in a bunch of arguments and outputs the things we need to a file which can then be read into the "main" program 
void fitAndOutput(TH1F* inputHist1, TF1* fit1, TF1* bkgFit1, TF1* sigFit1, double lowerFitRange1, double upperFitRange1, double binWidth1, 
                    int nentries1, double* par1, ofstream* outputFile1, TCanvas* allCanvases1, string varTag1){ 
    // Fit the histogram and extract parameters
    cout << "\n\nStarting new fit output" << endl;
    inputHist1->Fit("fit","RQB"); // B will enforce the bounds
    fit1->GetParameters(par1);
    cout << "Setting par[" << widthIndex << "] = width to be positive" << endl;
    if(par1[widthIndex] < 0){ par1[widthIndex] = abs(par1[widthIndex]); }
    fit1->SetParameters(par1);
    bkgFit1->SetParameters(par1);
    sigFit1->SetParameters(&par1[numDOFbkg]);
    *outputFile1 << "#nentries " << nentries1 << endl;
    
    for (int iPar=0; iPar<numDOFbkg+numDOFsig; ++iPar){
            *outputFile1 << namePar[iPar] << " " << par1[iPar] << endl;
    }
    
    // Want to to extract the yields of signal and bkg
    cout << "IF SIDEBAND SUBTRACTING ON THIS, THE INTEGRALS WILL OBVIOULSY BE WEIRD" << endl;
    cout << "----------------------------------------------------------------------" << endl;
    double integralBKG = bkgFit1->Integral(lowerFitRange1, upperFitRange1);
    double integralSIG = sigFit1->Integral(lowerFitRange1, upperFitRange1);
    cout << "IntegralBKG before scaling: " << integralBKG << endl;
    cout << "IntegralSIG before scaling: " << integralSIG << endl;
    integralBKG *= 1/binWidth1;
    integralSIG *= 1/binWidth1;
    cout << "IntegralBKG after scaling: " << integralBKG << endl;
    cout << "IntegralSIG after scaling: " << integralSIG << endl;
    cout << "nentries: " << nentries1 << endl;
    *outputFile1 << "#integralBKG " << integralBKG << endl;
    *outputFile1 << "#integralSIG " << integralSIG << endl;

    // Checking to see if the integral of the fits is the same/similar as the integral of the histogram
    int nSig=3;
    //double weightedSigma = par1[2]/(par1[2]+par1[2]*par1[5])*par1[4]+(par1[2]*par1[5])/(par1[2]+par1[2]*par1[5])*par1[6]*par1[4]; // for a double gaussian
    double weightedSigma = par1[widthIndex]; // not actually weighted in this case
    cout << "Checking if our calculated integrals for the bkg and signal add up to the total" << endl;
    *outputFile1 << "#weightedSigma " << weightedSigma << endl;
    TLine *line = new TLine(par1[meanIndex]-nSig*weightedSigma,0,par1[meanIndex]-nSig*weightedSigma,inputHist1->GetMaximum());
    line->SetLineColor(kMagenta);
    double integralBKG_nsig = bkgFit1->Integral(par1[meanIndex]-nSig*weightedSigma,par1[meanIndex]+nSig*weightedSigma);
    double integralSIG_nsig = sigFit1->Integral(par1[meanIndex]-nSig*weightedSigma,par1[meanIndex]+nSig*weightedSigma);
    *outputFile1 << "#integralBKG_nSig " << integralBKG_nsig << endl;
    *outputFile1 << "#integralSIG_nSig " << integralSIG_nsig << endl;
    // Calculate eventRatioSigToBkg.
    double eventRatioSigToBkg = integralSIG_nsig/integralBKG_nsig;
    *outputFile1 << "#eventRatioSigToBkg " << eventRatioSigToBkg << endl;

    Int_t binx1 = inputHist1->GetXaxis()->FindBin(par1[meanIndex]-nSig*weightedSigma);
    Int_t binx2 = inputHist1->GetXaxis()->FindBin(par1[meanIndex]+nSig*weightedSigma);
    *outputFile1 << "#-- 3sigma BinLower, BinUpper = " << binx1 << ", " << binx2 << endl;
    *outputFile1 << "#-- Actual counts in the 3sigma range " << inputHist1->Integral(binx1,binx2) << endl;
    *outputFile1 << "#-- integralBKG_nSig+integralSIG_nSig " << integralBKG_nsig+integralSIG_nsig << endl;
    if ( abs(1-inputHist1->Integral(binx1,binx2)/(integralBKG_nsig+integralSIG_nsig)) < 0.05 ) {
    	*outputFile1 << "#--There is agreement within 5%" << endl;
    }
    else {
    	*outputFile1 << "#--Percent off " << abs(1-inputHist1->Integral(binx1,binx2)/(integralBKG_nsig+integralSIG_nsig)) << endl;
    }
    double purity = integralSIG_nsig/(integralBKG_nsig+integralSIG_nsig);
    *outputFile1 << "#purity " << purity << endl;

    // Finally, we draw some histograms to see what the fit looks like
    inputHist1->Draw();
    line->DrawLine(par1[meanIndex]-nSig*weightedSigma,0,par1[meanIndex]-nSig*weightedSigma,inputHist1->GetMaximum());
    line->DrawLine(par1[meanIndex]+nSig*weightedSigma,0,par1[meanIndex]+nSig*weightedSigma,inputHist1->GetMaximum());
    inputHist1->GetXaxis()->SetTitleSize(0.04);
    inputHist1->GetYaxis()->SetTitleSize(0.04);
    inputHist1->SetTitle(("Peak: "+to_string(par1[meanIndex])+"    width: "+to_string(weightedSigma)).c_str());
    allCanvases1->SaveAs(("fitResults/"+fileTag+"/"+varTag1+"_fit_"+fileTag+".png").c_str());
    cout << "Saved for a specific fit!" << endl;
}
        	
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

    		TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
		double discrimVar;
		double sideBandVar;
		double Mpi0eta;
		double AccWeight;
                double sbWeight;
		dataTree->SetBranchAddress(s_discrimVar.c_str(),&discrimVar);
		dataTree->SetBranchAddress(s_sideBandVar.c_str(),&sideBandVar);
		dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
		dataTree->SetBranchAddress(s_accWeight.c_str(),&AccWeight);
                dataTree->SetBranchAddress(s_sbWeight.c_str(),&sbWeight);

                // Uniqueness tracking, the default way or using weighted by number of combos in an event 
        	bool isUniqueEtaB;
        	bool isUniquePi0B;
        	bool isUniquePi0EtaB;
        	dataTree->SetBranchAddress("isNotRepeated_eta",&isUniqueEtaB);
        	dataTree->SetBranchAddress("isNotRepeated_pi0",&isUniquePi0B);
        	dataTree->SetBranchAddress("isNotRepeated_pi0eta",&isUniquePi0EtaB);
                double utWeight;
                dataTree->SetBranchAddress("uniqunessTrackingWeights",&utWeight);

		long long nentries=dataTree->GetEntries();

		// ///////////////////////////////////////////////////
		// FILL YOUR HISTOGRAMS
		// ///////////////////////////////////////////////////
		std::vector<double> binRangeEta;
		std::vector<double> binRangePi0;
		double binWidthEta;
		double binWidthPi0;
		binRangeEta={200,0.25,0.85};
		binWidthEta=(binRangeEta[2]-binRangeEta[1])/binRangeEta[0];
		binRangePi0={200,0.005,0.25};
		binWidthPi0=(binRangePi0[2]-binRangePi0[1])/binRangePi0[0];

		TH1F* massHistEta = new TH1F("",";M(#eta) (GeV)",binRangeEta[0],binRangeEta[1],binRangeEta[2]);
		TH1F* massHistPi0 = new TH1F("",";M(#pi) (GeV)",binRangePi0[0],binRangePi0[1],binRangePi0[2]);
		TH1F* massHistPi0Eta = new TH1F("","", 350, 0, 3.5);
		cout << "Initialized for a specific mass (eta/pi0) fit" << endl;
		
		for (int ientry=0; ientry<nentries; ientry++)
		{
			dataTree->GetEntry(ientry);
                        //getSBWeight(sideBandVar,&sbWeight,weightingScheme);
                        double weight;
                        if (weightingScheme==""){ weight=1; }
                        if (weightingScheme=="as"){ weight=AccWeight; }
                        if (weightingScheme=="as*bs"){ weight=AccWeight*sbWeight; }
                        if (s_uniquenessTracking=="default"){
                	    if ( isUniqueEtaB ) {
		            	massHistEta->Fill(discrimVar,weight);
			    }
                	    if ( isUniquePi0B ) {
		            	massHistPi0->Fill(sideBandVar,weight); /////////////////////////////////////////// NOT WEIGHTED SINCE WE WONT BE ABLE TO FIT IT PROPERLY 
			    }
                	    if ( isUniquePi0EtaB ) {
			    	massHistPi0Eta->Fill(Mpi0eta,weight);//
			    }
                        }
                        else if (s_uniquenessTracking=="weighted"){
                                weight=weight*utWeight;
		            	massHistEta->Fill(discrimVar,weight);//
		            	massHistPi0->Fill(sideBandVar,weight); /////////////////////////////////////////// NOT WEIGHTED SINCE WE WONT BE ABLE TO FIT IT PROPERLY 
			    	massHistPi0Eta->Fill(Mpi0eta,weight);//
                        }
                        else {
                            cout << "uniqueness tracking setting not valid!" << endl;
                            exit(0);
                        }
		}
		cout << "Filled all entries into histogram for a specific fit" << endl;

		massHistPi0Eta->Draw("HIST");
		massHistPi0Eta->GetXaxis()->SetTitleSize(0.04);
		massHistPi0Eta->GetYaxis()->SetTitleSize(0.04);
		allCanvases->SaveAs(("fitResults/"+fileTag+"/Mpi0eta_fit_"+fileTag+".png").c_str());
		allCanvases->Clear();
        	cout <<"Initialized" << endl;

		TF1* fit;
        	TF1* bkgFit;
        	TF1* sigFit;
        	double par[numDOFbkg+numDOFsig];
		



		// ///////////////////////////////////////////////////
		// FIT YOUR HISTOGRAMS. fitAndOutput WILL OUTPUT THE NEEDED STUFF
		// ///////////////////////////////////////////////////
                //
    		ofstream logFile;
                // -------------------
		allCanvases->Clear();
		std::vector<double> fitRange;
		fitRange={0.4,0.7};
		fit = new TF1("fit",fitFunc,fitRange[0],fitRange[1],numDOFbkg+numDOFsig);
		bkgFit = new TF1("bkgFit",background,fitRange[0],fitRange[1],numDOFbkg);
		sigFit = new TF1("sigFit",signal,fitRange[0],fitRange[1],numDOFsig);

		fit->SetParameters(3000,2000,15000,0.547,0.003);//,3,5);
                fit->SetParLimits(0,0,4000);
                fit->SetParLimits(1,0,3000);
                fit->SetParLimits(2,0,20000);
                fit->SetParLimits(3,0.5,0.6);
                fit->SetParLimits(4,0,0.05);

    		logFile.open(("fitResults/"+fileTag+"/discrimVarFit_toMain_"+fileTag+".txt").c_str());
                string varTag="Meta";
                fitAndOutput(massHistEta, fit, bkgFit, sigFit, fitRange[0], fitRange[1], binWidthEta, nentries, &par[0], &logFile, allCanvases, varTag);
                logFile.close();

                // -------------------
		allCanvases->Clear();
		fitRange={0.09,0.18};
		fit = new TF1("fit",fitFunc,fitRange[0],fitRange[1],numDOFbkg+numDOFsig);
		bkgFit = new TF1("bkgFit",background,fitRange[0],fitRange[1],numDOFbkg);
		sigFit = new TF1("sigFit",signal,fitRange[0],fitRange[1],numDOFsig);

		fit->SetParameters(4000,3000,14000,0.134,0.015);//,5,2);
                fit->SetParLimits(0,0,7000);
                fit->SetParLimits(1,0,7000);
                fit->SetParLimits(2,0,20000);
                fit->SetParLimits(3,0.13,0.14);
                fit->SetParLimits(4,0.001,0.03);

    		logFile.open(("fitResults/"+fileTag+"/pi0Fit_toMain_"+fileTag+".txt").c_str());
                varTag="Mpi0";
                fitAndOutput(massHistPi0, fit, bkgFit, sigFit, fitRange[0], fitRange[1], binWidthPi0, nentries, &par[0], &logFile, allCanvases, varTag);
                logFile.close();

	}
}




































