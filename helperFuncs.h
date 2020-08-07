#ifndef GETSBWEIGHT_H
#define GETSBWEIGHT_H

#include <TMath.h>
#include <TPaveText.h>
#include <iostream>
#include <queue>
#include <chrono>

using namespace std;


///// **************************************************************
///// STEP0: DONT NEED TO MODIFY. RUN.PY WILL USE SED TO REPLACE THESE
///// **************************************************************
// !!!!!!
// NO SPACES BETWEEN THE = SIGNS. I USE SED TO REPLACE
// !!!!!!
bool verbose = true;
string rootFileLoc="/d/grid15/ln16/pi0eta/q-values/degALL_bcal_treeFlat_DSelector_UTweights.root";
string rootTreeName="degALL_bcal_tree_flat";
string fileTag="bcal";
string weightingScheme="as"; // "" or "as*bs"
string s_accWeight="AccWeight";
string s_sbWeight="weightSB";
string s_discrimVar="Meta";
string s_sideBandVar="Mpi0";
string s_uniquenessTracking="weighted"; // "default" or "weighted". Any neighbor is possible in the weighted setting

///// **************************************************************
///// STEP1: DEFINE SIGNAL/BKG DISTRIBUTIONS
///// **************************************************************

// We need to define a maximum value for the bernstein polynomial bkg. The domain is [0,1] so we have to scale the discriminating variable to be between 0,1. 
// Since we know the maximum value for Meta is < 1 then everything is fine. For larger masses we can divide by a larger value or just dividing by Ebeam_max = 12 is fine. 
// Without much calculation we know GlueX cannot produce anything with mass > 12 GeV. If you use other variables like then choose the appropriate maximum value
#define MAXVALUE 1

// -----------------------------
// Code for a single gaussian
double signal(double *x, double *par){
	return par[0]*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]));// + par[3]*exp(-0.5*((x[0]-par[4])/par[5])*((x[0]-par[4])/par[5]));

}
// Bernstein polynomial of degree 1
double background(double *x, double *par){
	return par[0]*x[0]/MAXVALUE+par[1]*(1-x[0]/MAXVALUE);
}
string namePar[numDOFbkg+numDOFsig] = {"bern01","bern11","amp","mass","sigma"};//,"ampRatio","sigmaRatio";
double fitFunc(double *x, double *par){
	return background(x,par)+signal(x,&par[numDOFbkg]);
}
int meanIndex=3; // used in getInitParams 
int widthIndex=4; // used in getInitParams
// -----------------------------
static const int numDOFsig = 3; // degrees of freedom for the signal distribution
static const int numDOFbkg = 2; //
// we need to tell "main" which variables need scaling when using the full reference distribution to the k nearest neighbors
std::vector<int> sigVarsNeedScaling = {2}; // these parameter indices are relative to fitFunc
std::vector<int> bkgVarsNeedScaling = {0,1};


// Code for a breit Wigner
//double signalBW(double* x, double* par) {
//    double arg1 = 14.0/22.0; // 2 over pi
//    double arg2 = par[2]*par[2]*par[1]*par[1]; //Gamma=par[1]  M=par[2]
//    double arg3 = ((x[0]*x[0]) - (par[1]*par[1]))*((x[0]*x[0]) - (par[1]*par[1]));
//    double arg4 = x[0]*x[0]*x[0]*x[0]*((par[2]*par[2])/(par[1]*par[1]));
//    
//    return par[0]*arg1*arg2/(arg3 + arg4);
//}

// ******* THIS FUNCTION ALREADY ACCOUNTS FOR THE USE OF AMP/WIDTH RATIOS
// DOuble gaus
//int numDOFsig = 5;
//double signal(double *x, double *par){
//	return par[0]/par[2]/TMath::Sqrt( 2*TMath::Pi() )*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]))
//	     + par[3]*par[0]/(par[4]*par[2])/TMath::Sqrt( 2*TMath::Pi() )*exp(-0.5*((x[0]-par[1])/(par[4]*par[2]))*((x[0]-par[1])/(par[4]*par[2])));
//
//}

struct parameterLimits{
    // this vector will contain the initialization parameters (variable initPar in main.h) from which we can use to set par limits if we want
    // These parameters would be a scaled version of the parameters taken from getInitParams.C, scaled down to kDim
    std::vector<double> initPars; 
    int numDOF = numDOFsig+numDOFbkg;
    double eventRatioSigToBkg;
    std::vector<double> lowerParLimits;
    std::vector<double> upperParLimits;

    void setupParLimits(){
        // We allow each parameter to vary based on the magnitude of the parameter times some scaleFactor.
        // Hopefully this is generic enough for different analyses but you can always set them directly here

        double scaleFactor = 3; // allow the parameters some flexibility
        double scaleSig = (1+1/eventRatioSigToBkg); // scale factor for the signal distribution to contain 100% of events
        double scaleBkg = (1+eventRatioSigToBkg); // scale factor for the bkg distribution to contain 100% of events
        // Since we do 3 different initializations, with 100%bkg, 50/50, 100% sig the parameters for the bkg all go to zero in the 100% signal case. 
        double abs_par = abs(initPars[0]*scaleBkg);
        lowerParLimits.push_back(0);
        upperParLimits.push_back(scaleBkg*initPars[0]+scaleFactor*abs_par);

        abs_par = abs(scaleBkg*initPars[1]);
        lowerParLimits.push_back(0);
        upperParLimits.push_back(scaleBkg*initPars[1]+scaleFactor*abs_par);

        abs_par = abs(scaleSig*initPars[2]);
        lowerParLimits.push_back(0);
        upperParLimits.push_back(scaleSig*initPars[2]+scaleFactor*abs_par);

        // since the mean and the width of gaussian is always + we can just multiply by a percentage
        lowerParLimits.push_back(initPars[3]*0.9);
        upperParLimits.push_back(initPars[3]*1.1);
        lowerParLimits.push_back(initPars[4]*0.70);
        upperParLimits.push_back(initPars[4]*1.30);
    }
    void printParLimits(){
        for (size_t i=0; i<numDOF; ++i){
            cout << "Par" << i << " has lower,upper limits: " << lowerParLimits[i] << ", " << upperParLimits[i] << endl;
        }
    }
};

std::vector<double> _SET_ParLimPercent = {300, 300, 300, 10, 75}; // These will be the percent deviations allowed for the parameters relative to initialization


///// **************************************************************
///// CLOSE STEP 1
///// **************************************************************


void getSBWeight(double discrimVar, double *sbWeight, std::string weightingScheme){
    // Updates the weight variable given a value for the discriminating variable
    // Determined by whether it is in the sideband region or the signal region.
    if (weightingScheme=="as*bs"){
        double sbRL = 0.09; // Right sideband left line
        double sbRR = 0.105; // Right sideband right line
        double sigL = 0.12;
        double sigR = 0.15;
        double sbLL = 0.165;
        double sbLR = 0.18;

        if ( discrimVar > sbRL && discrimVar < sbRR ) { *sbWeight = -1; } 
        else if ( discrimVar > sbLL && discrimVar < sbLR ) { *sbWeight = -1; } 
        else if ( discrimVar > sigL && discrimVar < sigR ) { *sbWeight = 1; } 
        else { *sbWeight = 0; }
    }
    else{
        *sbWeight=1;
    }
}






//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// DO NOT HAVE TO MODIFY ANYTHING BELOW THIS BLOCK
// DO NOT HAVE TO MODIFY ANYTHING BELOW THIS BLOCK

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// draws a textbox to output the fitted parameter values per Q-factor histogram
void drawText(Double_t *par, int dof, std::string tag, double qSigValue, double qBkgValue, double qTotValue){
    //TLatex parText;
    //parText.SetTextAlign(12);
    //parText.SetTextSize(0.01);
    //for (int iPar=0; iPar<dof; ++iPar){
    //    parText.DrawLatex(0.4,20,("par"+std::to_string(iPar)+":"+std::to_string(par[iPar])).c_str());
    //}
    //
    TPaveText *pt = new TPaveText(0.5,0.7,0.7,0.9,"NDC");
    for (int iPar=0; iPar<dof; ++iPar){
        pt->AddText((tag+std::to_string(iPar)+":"+std::to_string(par[iPar])).c_str());
    }
    pt->AddText(("qSigVal: "+std::to_string(qSigValue)).c_str());
    pt->AddText(("qBkgVal: "+std::to_string(qBkgValue)).c_str());
    pt->AddText(("qFitVal: "+std::to_string(qTotValue)).c_str());
    //pt->Paint("NDC");
    pt->Draw();
}

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


// Class to parse the string of phase space variables to consider
class parseVarString{
    private:
        int nVars=0;
        std::string delimiter=";";
        std::string token;
        size_t pos=0;
    public:
        std::vector<std::string> varStringSet;
        std::string varString;
        parseVarString (std::string inputString){
           varString =  inputString;
        }
        void parseString(){
            while ((pos = varString.find(delimiter)) != std::string::npos) {
                token = varString.substr(0, pos);
                varStringSet.push_back(token);
                //std::cout << token << std::endl;
                varString.erase(0, pos + delimiter.length());
                ++nVars;
            }
            varStringSet.push_back(varString);
            //std::cout << varString << std::endl;
        }
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











//void helperFuncs(){
//}



#endif
