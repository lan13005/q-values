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
#include <TSystem.h>
#include <TPolyLine3D.h>
#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
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
#include <RooMinimizer.h>
#include <RooTrace.h>
//#include "RooAddPdf.h"
//#include "RooFormulaVar.h"
//
//#include "Math/MinimizerOptions.h"


const int dim=7; // will get replaced by run.py
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
                string alwaysSaveTheseEvents;
                bool saveBranchOfNeighbors;
                string standardizationType;
                string cwd;
                bool redistributeBkgSigFits;
                bool doKRandomNeighbors;
		std::chrono::time_point<std::chrono::high_resolution_clock> start2;
                
                // These block of variables will be used to hold the initialization parameters. In the Q-factor paper they use 3 different initializations which
                // correspond to 100% bkg, 50/50, and 100% sig. If we want to do this here, the yields in the bkg and signal need to be modified. These vectors
                // will hold that information
                map<string, double> initializationParMap;
                std::vector<double> sigFracs;

                // initialize vectors to hold the discriminating and phase space variables
                parseVarString parseDiscrimVars;
	        parseVarString parsePhaseSpace;
                parseVarString parseEventsToSave;
                std::vector<std::vector<double>> discrimVars; 
                std::vector<double> discrimVar; // will contain the discriminating variables for just one entry whereas the above will hold all the data

                std::vector<std::vector<double>> phaseSpaceVars;
		std::vector<double> AccWeights; 
		std::vector<int> mcprocesses; 
                std::vector<double> utWeights;
                // Not all combinations will be a valid pairing. Suppose we only care about spectroscopically unique pairs, then we can fill phasePoint2PotentialNeighbor with
                // only unique combos.
		std::vector<int> phasePoint2PotentialNeighbor; 
	
	public:
		QFactorAnalysis(int kDim1, string varString1, string cwd1, string standardizationType1, bool redistributeBkgSigFits1, bool doKRandomNeighbors1,
                                int numberEventsToSavePerProcess1, int nProcess1, int seedShift1, Long64_t nentries1, int nRndRepSubset1, 
                                int nBS1, bool saveBShistsAlso1, bool override_nentries1, bool verbose1, string alwaysSaveTheseEvents1, 
                                bool saveBranchOfNeighbors1){ 
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
                        alwaysSaveTheseEvents=alwaysSaveTheseEvents1;
                        saveBranchOfNeighbors=saveBranchOfNeighbors1;
                        cwd=cwd1;
                        standardizationType=standardizationType1;
                        redistributeBkgSigFits=redistributeBkgSigFits1;
                        doKRandomNeighbors=doKRandomNeighbors1;
		}
		void loadTree(string rootFileLoc, string rootTreeName);
		void loadFitParameters(string fitLocation,string cwd);
		void loadData();
		void runQFactorThreaded(int iProcess);

};

double calculateStd(int nentries, double* input){
    double mean=0;
    for(int ientry=0; ientry<nentries; ++ientry){ 
        mean += input[ientry];
    }
    mean /= nentries;
    double diff;
    double std=0;
    for(int ientry=0; ientry<nentries; ++ientry){ 
        diff = input[ientry]-mean;
        std += diff*diff;
    }
    std /= nentries-1;
    return sqrt(std);
}


// Not used anymore. I thought it would interseting to check the stdev of the k nearest neighbors. 
// Used to calculate the current std as we stream/insert more data
class cumulativeStd{
    public:
        std::vector<double> inputVector;
        // kDim would be the typical size of the calculation. But sometimes we will choose nentries < kDim when testing quickly (this is never the case in real example though)
        cumulativeStd ( int kDim ){ _kDim=kDim; inputVector.reserve(_kDim); }
        void insertValue ( double value ) {
            inputVector[_timesCalled] = value;
            ++_timesCalled;
        }
        double calcStd(){
            for (int i=0; i<_timesCalled; ++i){
                _sum+=inputVector[i];
            }
            _sum /= _timesCalled;
            _mean=_sum; _sum=0;
            //std::cout << "mean: " << _mean << std::endl;
            for (int i=0; i<_timesCalled; ++i){
                _diff = (inputVector[i]-_mean);
                _sum += _diff*_diff;
                //std::cout << "sum: " << _sum << std::endl;
            }
            _sum /= _timesCalled;
            return sqrt(_sum);
        }


    private:
        int _timesCalled=0;
        double _sum=0;
        double _mean=0;
        double _diff;
        UInt_t _kDim;

};


// The following two classes will be used to keep track of pairs of distances and index, sorted by distance. priority_queue in stl requires three arguments
// which are type, container for type, and a comparator. we will use pair as the type which is held in a vector container. compareDist is the comparator
// which compares the first elements of two pairs. The first element is the distance, the second is the index j (in dij when calculating the distances). 
// distSort_kNN will setup our priority_queue that keeps a maximum of kDim elements simply by popping and pushing data. 
// This type of sorting should have k*log(k) sorting, I think. If we do this N times then the complexit ~ N*k*log(k). Probably have to check me on this
class compareDist
{
    public:
        bool operator() (const std::pair<double,int>& p1, const std::pair<double,int>& p2)
        {
            // < for min heap
            // maybe > for max heap?
            return p1.first < p2.first;
        }
};
class distSort_kNN
{
    public:
        // constructor with basic output and setting kDim
        distSort_kNN ( UInt_t kDim ) {
            //std::cout << "Priority Queue is set up to sort the {distance,index} pairs keeping " << kDim << " neighbors" << std::endl;
            _kDim = kDim;
        }

        std::priority_queue <std::pair<double,int>, vector<std::pair<double,int>>, compareDist > kNN;
        
        // depending on the size and the distanceis we will push, or pop+push
        void insertPair ( std::pair<double,int> newPair ){
            if ( kNN.size() >= _kDim ){
                // > for min heap and maybe < for max heap 
                if (kNN.top().first > newPair.first) {
                    //std::cout << "\nPOPPING OUT " << kNN.top().first << std::endl;
                    kNN.pop();
                    kNN.push( newPair );
                }
            } 
            else {      
                //std::cout << "\nPUSHING " << newPair.first << std::endl;
                kNN.push( newPair );
            }
        }

    private:
        UInt_t _kDim;
        std::pair<double,int> _pair;
};

// This class will be used to standardize the phase space variables. 
// Can either do stddev or range standardization
class standardizeArray{
	public:
    		double max_inputVector;
         	double min_inputVector;
		
		void rangeStandardization(std::vector<double> &inputVector, long long nentries){
    			max_inputVector = DBL_MIN;
         		min_inputVector = DBL_MAX;
                        std::cout << "\nStarting range std" << std::endl;
			for (int ientry=0; ientry<nentries; ++ientry){
				if (inputVector[ientry] > max_inputVector){
					max_inputVector = inputVector[ientry];
				}
				if (inputVector[ientry] < min_inputVector){
					min_inputVector = inputVector[ientry];
				}
			}
			for (int ientry=0; ientry<nentries; ++ientry){
				inputVector[ientry] = (inputVector[ientry]-min_inputVector)/(max_inputVector-min_inputVector);
			}
                        std::cout << "Max,min: " << max_inputVector << "," << min_inputVector << std::endl;
                        std::cout << "--Finished Range standardizing " << std::endl;
		}
		
		double calcStd(std::vector<double> &inputVector, long long nentries){
			double local_std=0;
			double diff=0;
			double mean=0;
			for (int i=0; i<nentries; ++i){
			    mean+=inputVector[i];
			} 
			mean/=nentries;
			for (int ientry=0; ientry<nentries; ++ientry){
                                std::cout << "mean: " << mean << std::endl;
				diff = (inputVector[ientry]-mean);
                                std::cout << "diff: " << diff << std::endl;
				local_std += diff*diff;
			}
			local_std /= nentries;
                        std::cout << "STD: " << local_std << std::endl;
			return sqrt(local_std);
		}
		
		void stdevStandardization(std::vector<double> &inputVector, long long nentries){
			double std = calcStd(inputVector, nentries);
			for (int ientry=0; ientry < nentries; ++ientry){
				inputVector[ientry] = inputVector[ientry]/std; 
			} 
                        std::cout << "Finished Stdev Standardization" << endl;
		}
};


// our distance calculation between two phase points
double calc_distance( int dim, double* phaseSpace_1, double* phaseSpace_2, bool verbose_outputDistCalc){
	double sum = 0;
	double diff=0;
        if(verbose_outputDistCalc){std::cout << "New event, new sum = " << sum << std::endl;}
	for (int i=0; i<dim; ++i){
		diff = phaseSpace_1[i]-phaseSpace_2[i];
		sum += diff*diff;
                if(verbose_outputDistCalc){
                    std::cout << "phasePoint1["<<i<<"]="<<phaseSpace_1[i]<<", phasePoint2["<<i<<"]="<<phaseSpace_2[i]<<" --- squared sum=" << diff*diff << "---- total so far="<<sum<<std::endl;
		}
	}
	return sum;
}

#endif
