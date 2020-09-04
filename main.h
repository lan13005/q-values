#ifndef MAIN_H
#define MAIN_H

#include "helperFuncs.h"
#include <thread>
#include "TThread.h"
#include "Math/MinimizerOptions.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <TCanvas.h>
#include <TRandom.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPolyLine3D.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGenericPdf.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooConstVar.h>
#include <RooArgList.h>
#include <RooPlot.h>
#include <RooMinuit.h>
#include <RooFitResult.h>
#include <RooChi2Var.h>
#include <RooWorkspace.h>
//#include "RooAddPdf.h"
//#include "RooFormulaVar.h"
//
//#include "Math/MinimizerOptions.h"


const int dim=5; // will get replaced by run.py
bool verbose_outputDistCalc=false;
TRandom rgen;

using namespace std;

// OUT OF DATED CODE THAT USES ROOFIT TO DO UNBINNED MAX LIKELIHOOD FIT. WILL PROBABLY NEED TO REIMPLEMENT THIS
//class rooFitML{
//	private:
//        	double peakWidtheta[2] = {0.545928, 0.0196892};
//        	double peakWidthpi0[2] = {0.135399, 0.00760648};
//        	double par0eta[3] = {3.4528, 1.7264, 0};
//        	double par1eta[3] = {5.49805, 2.74902, 0}; 
//        	double par2eta[3] = {0, 0.9, 1.18}; 
//        	double par0pi0[3] = {7.30661, 3.6533, 0}; 
//        	double par1pi0[3] = {30.05331, 15.2665, 0}; 
//        	double par2pi0[3] = {0, 0.400002, 0.800003}; 
//
//		/// we had to make everything deal with pointers since we couldn't define and initialize at the same time here. It can only be done with static const primitive types
//		RooRealVar par0_eta("par0_eta","par0_eta", 0, 1, 5.25);
//		RooRealVar par1_eta("par1_eta","par1_eta", 0, -8.25, 8.25);
//		RooRealVar par2_eta("par2_eta","par2_eta", 0, 0, 1.77);
//		RooRealVar peak_eta("peak_eta","peak_eta",0.545928,0.52,0.58);
//		RooRealVar x("x","x",0, 1);
//		RooRealVar width_eta("width_eta","width_eta",0.0196892,0.017,0.027);
//		RooDataSet data("data","data",RooArgSet(x));
//
//		RooPlot xframe = x->frame();
//	public:
//		//RooGaussian gauss("gauss","gauss(x,m,s)", x, peak_eta, width_eta) ;
//		//RooFormulaVar poly("poly","polynomial",x,RooArgList(par0_eta,par1_eta)) ;
//		//RooAddPdf sigWithBkg("sigWithBkg","gaus+polyBkg",RooArgList(gauss,poly),RooArgList(scale
//		RooGenericPdf sigWithBkg("sigWithBkg","sigWithBkg","par0_eta+par1_eta*x+par2_eta*exp(-0.5*((x-peak_eta)/width_eta)**2)", RooArgSet(x,par0_eta,par1_eta,par2_eta,peak_eta,width_eta));
//
//
//		void fillValues( double mass ){
//			data.add(RooArgSet(mass));
//		}
//
//		void fitToData( UInt_t iFit ){
//			par0_eta.setVal(par0eta[iFit]);
//			par1_eta.setVal(par1eta[iFit]);
//			par2_eta.setVal(par2eta[iFit]); 
//			// not sure what these flags do
//			sigWithBkg.fitTo(data,"mhvr");
//		}
//
//		void drawFitData(){
//			data.plotOn(xframe);
//			sigWithBkg.plotOn(xframe);
//			xframe.Draw();
//		}
//
//		void clear(){
//			// think this deletes the pool of data
//			data.Clear();
//			// we might also have to reinitialize the RooRealVars if the fitTo function decides to change things, idky it would though
//		}
//
//};


class QFactorAnalysis{
	private:
		// variables to help load the data
		TTree *dataTree;
		Long64_t total_nentries;
	
		// we define the variables that we need in the thread functions
		int kDim;
		int numberEventsToSavePerProcess;
		int nProcess;
		int seedShift;
		Long64_t nentries;
                int nRndRepSubset;
                bool saveBShistsAlso;
                int nBS;
		bool override_nentries;
		bool verbose;
		string varString;
                string standardizationType;
                bool redistributeBkgSigFits;
                bool doKRandomNeighbors;
		std::chrono::time_point<std::chrono::high_resolution_clock> start2;
                
                // These block of variables will be used to hold the initialization parameters. In the Q-factor paper they use 3 different initializations which
                // correspond to 100% bkg, 50/50, and 100% sig. If we want to do this here, the yields in the bkg and signal need to be modified. These vectors
                // will hold that information
                map<string, double> initializationParMap;
                std::vector<double> redistributeFactorBkg;
                std::vector<double> redistributeFactorSig;
                parameterLimits parLimits;
                parameterLimits2 parLimits2;
                std::vector<double> sigFracs;

                // initialize vectors to hold the discriminating and phase space variables
                parseVarString parseDiscrimVars;
	        parseVarString parsePhaseSpace;
                std::vector<std::vector<double>> discrimVars; 
                std::vector<double> discrimVar; // will contain the discriminating variables for just one entry whereas the above will hold all the data

                std::vector<std::vector<double>> phaseSpaceVars;
		std::vector<double> AccWeights; 
                std::vector<double> utWeights;
		std::vector<ULong64_t> spectroscopicComboIDs; 
                // Not all combinations will be a valid pairing. Suppose we only care about spectroscopically unique pairs, then we can fill phasePoint2PotentialNeighbor with
                // only unique combos.
		std::vector<int> phasePoint2PotentialNeighbor; 
	
	public:
		QFactorAnalysis(int kDim1, string varString1, string standardizationType1, bool redistributeBkgSigFits1, bool doKRandomNeighbors1, int numberEventsToSavePerProcess1, int nProcess1,
                                int seedShift1, Long64_t nentries1, int nRndRepSubset1, int nBS1, bool saveBShistsAlso1, bool override_nentries1, bool verbose1){ 
			cout << "Constructed QFactorAnalysis class..." << endl;
			kDim=kDim1;
			numberEventsToSavePerProcess=numberEventsToSavePerProcess1;
			nProcess = nProcess1;
			seedShift=seedShift1;
			nentries=nentries1;
                        nRndRepSubset=nRndRepSubset1;
                        nBS=nBS1;
                        saveBShistsAlso=saveBShistsAlso1;
			override_nentries=override_nentries1;
			verbose=verbose1;
			start2 = std::chrono::high_resolution_clock::now();
			varString=varString1;
                        standardizationType=standardizationType1;
                        redistributeBkgSigFits=redistributeBkgSigFits1;
                        doKRandomNeighbors=doKRandomNeighbors1;
		}
		void loadTree(string rootFileLoc, string rootTreeName);
		void loadFitParameters(string fitLocation);
		void loadData();
		void runQFactorThreaded(int iProcess);

};


#endif
