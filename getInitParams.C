// README
// The purpose of this code is to fit the distribution of the discriminating variables and extract the fit parameters
// We can then pass these fitted parameters, scale them, and use them as initializations for the q-values

#include <string>
#include "getInitParams.h"
bool do2Dfit=true;

// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// This takes in a bunch of arguments and outputs the things we need to a file which can then be read into the "main" program 
//   also makes some pretty plots
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
void fitAndOutput(TH1F* inputHist1, TF1* fit1, TF1* bkgFit1, TF1* sigFit1, double lowerFitRange1, double upperFitRange1, double binWidth1, 
                    int nentries1, double* par1, ofstream* outputFile1, TCanvas* allCanvases1, string varTag1){ 
    // Fit the histogram and extract parameters
    cout << "\n\nStarting new fit output" << endl;
    inputHist1->Fit(fit1->GetName(),"RQB"); // B will enforce the bounds
    fit1->GetParameters(par1);

    int meanIndex=meanIndices[0];
    int widthIndex=widthIndices[0];
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
    cout << "nentries" << nentries1 << endl;
    cout << "ChiSqPerDOF: " << fit1->GetChisquare()/(fit1->GetNDF()) << endl;
    *outputFile1 << "#integralBKG " << integralBKG << endl;
    *outputFile1 << "#integralSIG " << integralSIG << endl;

    // Checking to see if the integral of the fits is the same/similar as the integral of the histogram
    int nSig=3;
    //double weightedSigma = par1[2]/(par1[2]+par1[2]*par1[5])*par1[4]+(par1[2]*par1[5])/(par1[2]+par1[2]*par1[5])*par1[6]*par1[4]; // for a double gaussian
    double weightedSigma = par1[widthIndex]; // not actually weighted in this case
    cout << "Checking if our calculated integrals for the bkg and signal add up to the total" << endl;
    *outputFile1 << "#weightedSigma " << weightedSigma << endl;
    TLine* line = new TLine(par1[meanIndex]-nSig*weightedSigma,0,par1[meanIndex]-nSig*weightedSigma,inputHist1->GetMaximum());
    line->SetLineColor(kMagenta);
    double integralBKG_nsig = bkgFit1->Integral(par1[meanIndex]-nSig*weightedSigma,par1[meanIndex]+nSig*weightedSigma);
    double integralSIG_nsig = sigFit1->Integral(par1[meanIndex]-nSig*weightedSigma,par1[meanIndex]+nSig*weightedSigma);
    *outputFile1 << "#integralBKG_nSig " << integralBKG_nsig << endl;
    *outputFile1 << "#integralSIG_nSig " << integralSIG_nsig << endl;
    integralBKG_nsig *= 1/binWidth1;
    integralSIG_nsig *= 1/binWidth1;
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
    allCanvases1->SaveAs(("fitResults"+runTag+"/"+fileTag+"/"+varTag1+"_fit_"+fileTag+".png").c_str());
    cout << "Saved for a specific fit!" << endl;
}

// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// 2D version of the above. Incase the sideband is 2D which the background in pi0eta is
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
//
void fitAndOutput2D(TH2F* inputHist1, TF2* fit1, TF2* bkgFit1, TF2* sigFit1, double lowerFitRange1, double upperFitRange1, double lowerFitRange2, double upperFitRange2, vector<double> binRangePi0, vector<double> binRangeEta,
                    TF1* proj_var1, TF1* proj_var2, int nentries1, double* par1, ofstream* outputFile1, TCanvas* allCanvases1, string varTag1){ 
    // Fit the histogram and extract parameters
    cout << "\n\nStarting new fit output" << endl;
    inputHist1->Fit(fit1->GetName(),"RQB"); // B will enforce the bounds
    fit1->GetParameters(par1);

    double binWidthPi01=(binRangePi0[2]-binRangePi0[1])/binRangePi0[0];
    double binWidthEta1=(binRangeEta[2]-binRangeEta[1])/binRangeEta[0];

    int meanIndex1=meanIndices2[0]; // pi0 is 0th index, eta is 1st
    int meanIndex2=meanIndices2[1];
    int widthIndex1=widthIndices2[0];
    int widthIndex2=widthIndices2[1];
    cout << "Setting par[" << widthIndex1 << "] = width to be positive" << endl;
    cout << "Setting par[" << widthIndex2 << "] = width to be positive" << endl;
    par1[widthIndex1] = abs(par1[widthIndex1]);
    par1[widthIndex2] = abs(par1[widthIndex2]);
    fit1->SetParameters(par1);
    bkgFit1->SetParameters(par1);
    sigFit1->SetParameters(&par1[numDOFbkg2]);
    *outputFile1 << "#nentries " << nentries1 << endl;
    
    for (int iPar=0; iPar<numDOFbkg2+numDOFsig2; ++iPar){
            *outputFile1 << namePar2[iPar] << " " << par1[iPar] << endl;
    }
    
    // Want to to extract the yields of signal and bkg
    cout << "IF SIDEBAND SUBTRACTING ON THIS, THE INTEGRALS WILL OBVIOULSY BE WEIRD" << endl;
    cout << "----------------------------------------------------------------------" << endl;
    double integralBKG = bkgFit1->Integral(binRangePi0[1], binRangePi0[2], binRangeEta[1], binRangeEta[2]);
    double integralSIG = sigFit1->Integral(binRangePi0[1], binRangePi0[2], binRangeEta[1], binRangeEta[2]);
    cout << "IntegralBKG before scaling: " << integralBKG << endl;
    cout << "IntegralSIG before scaling: " << integralSIG << endl;
    integralBKG *= 1/binWidthPi01/binWidthEta1;
    integralSIG *= 1/binWidthPi01/binWidthEta1;
    cout << "IntegralBKG after scaling: " << integralBKG << endl;
    cout << "IntegralSIG after scaling: " << integralSIG << endl;
    cout << "Summed Integrals: " << integralBKG+integralSIG << endl;
    cout << "Entries " << nentries1 << endl;
    cout << "ChiSqPerDOF: " << fit1->GetChisquare()/(fit1->GetNDF()) << endl;
    *outputFile1 << "#integralBKG " << integralBKG << endl;
    *outputFile1 << "#integralSIG " << integralSIG << endl;
    
    // Checking to see if the integral of the fits is the same/similar as the integral of the histogram
    int nSig=3;
    //double weightedSigma = par1[2]/(par1[2]+par1[2]*par1[5])*par1[4]+(par1[2]*par1[5])/(par1[2]+par1[2]*par1[5])*par1[6]*par1[4]; // for a double gaussian
    double weightedSigma1 = par1[widthIndex1]; // not actually weighted in this case
    double weightedSigma2 = par1[widthIndex2]; // not actually weighted in this case
    cout << "Checking if our calculated integrals for the bkg and signal add up to the total" << endl;
    *outputFile1 << "#weightedSigma1 " << weightedSigma1 << endl;
    cout << "3Sig pi0 range: " << par1[meanIndex1]-nSig*weightedSigma1 << ", " << par1[meanIndex1]+nSig*weightedSigma1 << endl;
    cout << "3Sig eta range: " << par1[meanIndex2]-nSig*weightedSigma2 << ", " << par1[meanIndex2]+nSig*weightedSigma2 << endl;
    double integralBKG_nsig = bkgFit1->Integral(par1[meanIndex1]-nSig*weightedSigma1,par1[meanIndex1]+nSig*weightedSigma1, par1[meanIndex2]-nSig*weightedSigma2,par1[meanIndex2]+nSig*weightedSigma2);
    double integralSIG_nsig = sigFit1->Integral(par1[meanIndex1]-nSig*weightedSigma1,par1[meanIndex1]+nSig*weightedSigma1, par1[meanIndex2]-nSig*weightedSigma2,par1[meanIndex2]+nSig*weightedSigma2);
    integralBKG_nsig *= 1/binWidthPi01/binWidthEta1;
    integralSIG_nsig *= 1/binWidthPi01/binWidthEta1;
    *outputFile1 << "#integralBKG_nSig " << integralBKG_nsig << endl;
    *outputFile1 << "#integralSIG_nSig " << integralSIG_nsig << endl;
    
    // Calculate eventRatioSigToBkg.
    double eventRatioSigToBkg = integralSIG_nsig/integralBKG_nsig;
    *outputFile1 << "#eventRatioSigToBkg " << eventRatioSigToBkg << endl;

    Int_t binx1 = inputHist1->GetXaxis()->FindBin(par1[meanIndex1]-nSig*weightedSigma1);
    Int_t binx2 = inputHist1->GetXaxis()->FindBin(par1[meanIndex1]+nSig*weightedSigma1);
    Int_t biny1 = inputHist1->GetYaxis()->FindBin(par1[meanIndex2]-nSig*weightedSigma2);
    Int_t biny2 = inputHist1->GetYaxis()->FindBin(par1[meanIndex2]+nSig*weightedSigma2);
    *outputFile1 << "#-- 3sigma BinXLower, BinXUpper = " << binx1 << ", " << binx2 << endl;
    *outputFile1 << "#-- 3sigma BinYLower, BinYUpper = " << biny1 << ", " << biny2 << endl;
    *outputFile1 << "#-- Actual counts in the 3sigma range " << inputHist1->Integral(binx1,binx2,biny1,biny2) << endl;
    *outputFile1 << "#-- Actual total counts " << inputHist1->Integral(0,binRangePi0[0],0,binRangeEta[0]) << endl;
    *outputFile1 << "#-- integralBKG_nSig+integralSIG_nSig " << integralBKG_nsig+integralSIG_nsig << endl;
    if ( abs(1-inputHist1->Integral(binx1,binx1,biny1,biny2)/(integralBKG_nsig+integralSIG_nsig)) < 0.05 ) {
    	*outputFile1 << "#--There is agreement within 5%" << endl;
    }
    else {
    	*outputFile1 << "#--Percent off " << abs(1-inputHist1->Integral(binx1,binx2,biny1,biny2)/(integralBKG_nsig+integralSIG_nsig)) << endl;
    }
    double purity = integralSIG_nsig/(integralBKG_nsig+integralSIG_nsig);
    *outputFile1 << "#purity " << purity << endl;

    // Finally, we draw some histograms to see what the fit looks like
    allCanvases->Clear();
    allCanvases->Divide(3,1);
    allCanvases->cd(2);
    inputHist1->Draw("SURF1");

    // Setting up the projections then drawing
    allCanvases->cd(1);
    for (int iPar=0; iPar<numDOFsig2+numDOFbkg2; ++iPar){
        par2D_proj1[iPar] = par1[iPar];
    }
    proj_var1->SetParameters(par2D_proj1);
    TH1F *dHist_var1 = (TH1F *)inputHist1->ProjectionX("proj_pi0", 0, -1);
    dHist_var1->GetXaxis()->SetTitleSize(0.04);
    dHist_var1->GetYaxis()->SetTitleSize(0.04);
    dHist_var1->SetTitle(("Peak: "+to_string(par1[meanIndex1])+"    width: "+to_string(weightedSigma1)).c_str());
    dHist_var1->Draw("HIST");
    proj_var1->Draw("SAME");
    TLine *line1 = new TLine(par1[meanIndex1]-nSig*weightedSigma1,0,par1[meanIndex1]-nSig*weightedSigma1,dHist_var1->GetMaximum());
    line1->SetLineColor(kMagenta);
    line1->DrawLine(par1[meanIndex1]-nSig*weightedSigma1,0,par1[meanIndex1]-nSig*weightedSigma1,dHist_var1->GetMaximum());
    line1->DrawLine(par1[meanIndex1]+nSig*weightedSigma1,0,par1[meanIndex1]+nSig*weightedSigma1,dHist_var1->GetMaximum());
    
    // Setting up the projections then drawing
    allCanvases->cd(3);
    for (int iPar=0; iPar<numDOFsig2+numDOFbkg2; ++iPar){
        par2D_proj2[iPar] = par1[iPar];
    }
    proj_var2->SetParameters(par2D_proj2);
    TH1F *dHist_var2 = (TH1F *)inputHist1->ProjectionY("proj_eta", 0, -1);
    dHist_var2->GetXaxis()->SetTitleSize(0.04);
    dHist_var2->GetYaxis()->SetTitleSize(0.04);
    dHist_var2->SetTitle(("Peak: "+to_string(par1[meanIndex2])+"    width: "+to_string(weightedSigma2)).c_str());
    TLine *line2 = new TLine(par1[meanIndex2]-nSig*weightedSigma2,0,par1[meanIndex2]-nSig*weightedSigma2,dHist_var2->GetMaximum());
    line2->SetLineColor(kMagenta);
    dHist_var2->Draw("HIST");
    proj_var2->Draw("SAME");
    line2->DrawLine(par1[meanIndex2]-nSig*weightedSigma2,0,par1[meanIndex2]-nSig*weightedSigma2,dHist_var2->GetMaximum());
    line2->DrawLine(par1[meanIndex2]+nSig*weightedSigma2,0,par1[meanIndex2]+nSig*weightedSigma2,dHist_var2->GetMaximum());
    allCanvases1->SaveAs(("fitResults"+runTag+"/"+fileTag+"/"+varTag1+"_fit_"+fileTag+".png").c_str());
    cout << "Saved for a specific fit!" << endl;
}
        	
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// Load the data, fill the histograms, fitAndOuput the results
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
void getInitParams(){
	gStyle->SetOptFit(111);
	gStyle->SetOptStat(0);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);

	static const int nDetectors = 1;
	for (int i=0; i<nDetectors; ++i){
		cout << "File Name: " << rootFileLoc << endl;
		cout << "Tree Name: " << rootTreeName << endl;
		TFile* dataFile=new TFile(rootFileLoc.c_str());
		dataFile->GetObject(rootTreeName.c_str(),dataTree);

                allCanvases = new TCanvas("anyHists","",1440,900);
                
                parseVarString parseDiscrimVars;
                parseDiscrimVars.updateString(s_discrimVar); 
	        parseDiscrimVars.parseString();
                std::vector<double> discrimVar;
                // we have to push_back first for w.e. reason
	        for (int iVar=0; iVar<parseDiscrimVars.varStringSet.size(); ++iVar){
                    discrimVar.push_back(0);
                    discrimVar.push_back(0);
                }
	        for (int iVar=0; iVar<parseDiscrimVars.varStringSet.size(); ++iVar){
                    cout << "Loading discrimVar: " << parseDiscrimVars.varStringSet[iVar] << endl;
	            dataTree->SetBranchAddress(parseDiscrimVars.varStringSet[iVar].c_str(), &discrimVar[iVar]);
                }
		dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
		dataTree->SetBranchAddress(s_accWeight.c_str(),&AccWeight);

                if(!s_sbWeight.empty()){
                    dataTree->SetBranchAddress(s_sbWeight.c_str(),&sbWeight);
                }
                else{
                    sbWeight=1;
                }

                // Uniqueness tracking, the default way or using weighted by number of combos in an event 
                if (!s_utBranch.empty()){
                    dataTree->SetBranchAddress(s_utBranch.c_str(),&utWeight);
                }
                else{
                    utWeight=1;
                }

		long long nentries=dataTree->GetEntries();

		// ///////////////////////////////////////////////////
		// FILL YOUR HISTOGRAMS
		// ///////////////////////////////////////////////////

                massHistEta = new TH1F("",";M(#eta) (GeV)",binRangeEta[0],binRangeEta[1],binRangeEta[2]);
                massHistPi0 = new TH1F("",";M(#pi) (GeV)",binRangePi0[0],binRangePi0[1],binRangePi0[2]);
                massHistPi0VsEta = new TH2F("",";M(#pi) (GeV);M(#eta) (GeV)", binRangePi0[0],binRangePi0[1],binRangePi0[2],binRangeEta[0],binRangeEta[1],binRangeEta[2]);
                massHistPi0Eta = new TH1F("","", 350, 0, 3.5);
		cout << "Initialized for a specific mass (eta/pi0) fit" << endl;
		
		for (int ientry=0; ientry<nentries; ientry++)
		{
			dataTree->GetEntry(ientry);
                        //getSBWeight(discrimVar[1],&sbWeight,weightingScheme);
                        if (weightingScheme==""){ weight=1; }
                        if (weightingScheme=="as"){ weight=AccWeight; }
                        if (weightingScheme=="as*bs"){ weight=AccWeight*sbWeight; }
                        weight=weight*utWeight;

			massHistPi0Eta->Fill(Mpi0eta,weight);//
		        massHistEta->Fill(discrimVar[1],weight);//
		        massHistPi0->Fill(discrimVar[0],weight); /////////////////////////////////////////// NOT WEIGHTED SINCE WE WONT BE ABLE TO FIT IT PROPERLY 
                        massHistPi0VsEta->Fill(discrimVar[0],discrimVar[1],weight);
		}
		cout << "Filled all entries into histogram for a specific fit" << endl;

		massHistPi0Eta->Draw("HIST");
		massHistPi0Eta->GetXaxis()->SetTitleSize(0.04);
		massHistPi0Eta->GetYaxis()->SetTitleSize(0.04);
		allCanvases->SaveAs(("fitResults"+runTag+"/"+fileTag+"/Mpi0eta_fit_"+fileTag+".png").c_str());
		allCanvases->Clear();
		massHistPi0->Draw("HIST");
		massHistPi0->GetXaxis()->SetTitleSize(0.04);
		massHistPi0->GetYaxis()->SetTitleSize(0.04);
		allCanvases->SaveAs(("fitResults"+runTag+"/"+fileTag+"/Mpi0_fit_"+fileTag+".png").c_str());
		allCanvases->Clear();
		massHistEta->Draw("HIST");
		massHistEta->GetXaxis()->SetTitleSize(0.04);
		massHistEta->GetYaxis()->SetTitleSize(0.04);
		allCanvases->SaveAs(("fitResults"+runTag+"/"+fileTag+"/Meta_fit_"+fileTag+".png").c_str());
		allCanvases->Clear();
        	cout <<"Initialized" << endl;

		// ///////////////////////////////////////////////////
		// FIT YOUR HISTOGRAMS. fitAndOutput WILL OUTPUT THE NEEDED STUFF
		// ///////////////////////////////////////////////////
    		ofstream logFile;

                // 1D Fits
                if (!do2Dfit){
                    // -------------------
		    allCanvases->Clear();
		    fit = new TF1("fit",fitFunc,fitRange2[0],fitRange2[1],numDOFbkg+numDOFsig);
		    bkgFit = new TF1("bkgFit",background,fitRange2[0],fitRange2[1],numDOFbkg);
		    sigFit = new TF1("sigFit",signal,fitRange2[0],fitRange2[1],numDOFsig);

		    fit->SetParameters(3000,2000,15000,0.547,0.05);//,3,5);
                    fit->SetParLimits(0,0,4000);
                    fit->SetParLimits(1,0,3000);
                    fit->SetParLimits(2,0,20000);
                    fit->SetParLimits(3,0.5,0.6);
                    fit->SetParLimits(4,0,0.05);

    		    logFile.open(("fitResults"+runTag+"/"+fileTag+"/var2Fit_toMain_"+fileTag+".txt").c_str());
                    string varTag="Meta";
                    fitAndOutput(massHistEta, fit, bkgFit, sigFit, fitRange2[0], fitRange2[1], binWidthEta, nentries, &par[0], &logFile, allCanvases, varTag);
                    logFile.close();

                    // -------------------
		    allCanvases->Clear();
		    fit = new TF1("fit",fitFunc,fitRange1[0],fitRange1[1],numDOFbkg+numDOFsig);
		    bkgFit = new TF1("bkgFit",background,fitRange1[0],fitRange1[1],numDOFbkg);
		    sigFit = new TF1("sigFit",signal,fitRange1[0],fitRange1[1],numDOFsig);

		    fit->SetParameters(4000,3000,14000,0.134,0.015);//,5,2);
                    fit->SetParLimits(0,0,7000);
                    fit->SetParLimits(1,0,7000);
                    fit->SetParLimits(2,0,20000);
                    fit->SetParLimits(3,0.13,0.14);
                    fit->SetParLimits(4,0.001,0.03);

    		    logFile.open(("fitResults"+runTag+"/"+fileTag+"/var1Fit_toMain_"+fileTag+".txt").c_str());
                    varTag="Mpi0";
                    fitAndOutput(massHistPi0, fit, bkgFit, sigFit, fitRange1[0], fitRange1[1], binWidthPi0, nentries, &par[0], &logFile, allCanvases, varTag);
                    logFile.close();
                }
                else {
                    // -------------------
		    allCanvases->Clear();
		    fit2D = new TF2("fit2D",fitFunc2,fitRange1[0],fitRange1[1], fitRange2[0], fitRange2[1], numDOFbkg2+numDOFsig2);
		    bkgFit2D = new TF2("bkgFit2D",background2,fitRange1[0],fitRange1[1], fitRange2[0], fitRange2[1], numDOFbkg2);
		    sigFit2D = new TF2("sigFit2D",signal2,fitRange1[0],fitRange1[1], fitRange2[0], fitRange2[1], numDOFsig2);
                    proj_var1 = new TF1("proj_var1",fitFunc2_projectVar1,fitRange1[0],fitRange1[1],numDOFsig2+numDOFbkg2);
                    proj_var2 = new TF1("proj_var2",fitFunc2_projectVar2,fitRange2[0],fitRange2[1],numDOFsig2+numDOFbkg2);

		    fit2D->SetParameters(50,50,50,50,1,0.135,0.01,0.547,0.02);
                    fit2D->SetParLimits(0,0,1000);
                    fit2D->SetParLimits(1,0,1000);
                    fit2D->SetParLimits(2,0,1000);
                    fit2D->SetParLimits(3,0,1000);
                    fit2D->SetParLimits(4,0,1000);
                    fit2D->SetParLimits(5,0.13,0.14);
                    fit2D->SetParLimits(6,0.001,0.03);
                    fit2D->SetParLimits(7,0.5,0.6);
                    fit2D->SetParLimits(8,0,0.05);

    		    logFile.open(("fitResults"+runTag+"/"+fileTag+"/var1Vsvar2Fit_toMain_"+fileTag+".txt").c_str());
                    string varTag="Mpi0VsMeta";
                    fitAndOutput2D(massHistPi0VsEta, fit2D, bkgFit2D, sigFit2D, fitRange1[0], fitRange1[1], fitRange2[0], fitRange2[1], binRangePi0, binRangeEta, proj_var1,proj_var2,
                            nentries, &par2D[0], &logFile, allCanvases, varTag);
                    logFile.close();
                }

	}
}




































