#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <fstream>
#include <queue>
#include <ctime>
#include <chrono>
#include <TCanvas.h>
#include <TRandom.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>


//#include "Math/MinimizerOptions.h"

double peakWidth_eta[2] = {0.544,0.022};
double peakWidth_pi0[2] = {0.134,0.0069};

const int dim=7;
bool verbose2=false;
bool verbose_outputDistCalc=false;
TRandom rgen;

using namespace std;


double calc_distance( double phaseSpace_1[dim], double phaseSpace_2[dim] ){
	double sum = 0;
	double diff=0;
        if(verbose_outputDistCalc){cout << "New event, new sum = " << sum << endl;}
	for (int i=0; i<dim; ++i){
		diff = phaseSpace_1[i]-phaseSpace_2[i];
		sum += diff*diff;
                if(verbose_outputDistCalc){cout << "phasePoint1["<<i<<"]="<<phaseSpace_1[i]<<", phasePoint2["<<i<<"]="<<phaseSpace_2[i]<<" --- squared sum=" << diff*diff << "---- total so far="<<sum<<endl;}
	}
	return sum;
}


int numDOFbkg = 1;
Double_t background(Double_t *x, Double_t *par){
	return par[0];
	//return par[0]+par[1]*x[0];
	//return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
	//return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0];
	//return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0];
}

int numDOFsig = 3;
Double_t signal(Double_t *x, Double_t *par){
	return par[0]/par[2]/TMath::Sqrt( 2*TMath::Pi() )*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]));// + par[3]*exp(-0.5*((x[0]-par[1])/par[4])*((x[0]-par[1])/par[4]));
	//return par[0]*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]));// + par[3]*exp(-0.5*((x[0]-par[4])/par[5])*((x[0]-par[4])/par[5]));
	//return (x[0]-par[0])*(x[0]-par[0]);

}

Double_t fitFunc(Double_t *x, Double_t *par){
	return background(x,par)+signal(x,&par[numDOFbkg]);
}

double gausAmplitude( double xmin, double binsize, int binsToStep, int kDim, Double_t *par){
    // So this goal of this function is to calculate the Amplitude we need to scale a unit gaussian by to get the entries = k in kDim. 
    // In a perfect fit with some peak and width the gausssian overlays the bins heights perfectly. If we can evaluate the gaussian at various peaks
    // we can get the counts. Do a discrete sum over all the bin heights to get the count. 
    double discreteSum=0;
    double x = xmin;
    for (int iBin=0; iBin<binsToStep; ++iBin){
        discreteSum+=signal(&x, par);
        x+=binsize;
    }
    return (double)kDim/discreteSum;
}

class standardizeArray{
     // ERROR 1: not sure why this code did not work. It would sometimes give a max element that was not correct
     //double max_inputVector = *std::max(inputVector, inputVector+nentries);
     //double min_inputVector = *std::min(inputVector, inputVector+nentries);
     
     // ERROR 2: I tried to create a function that does this max,min calculation and returns an array. This was abit more difficult and quite buggy. Functions can return array points
     // but also would require static types or something like that to make it go outside the functions scope. This lead to wierd effects like calling the function twice might
     // pass the max, min from one function call to another. Safest is to just do it here. 
    private:
         long long _nentries;
         double _max_inputVector = DBL_MIN;
         double _min_inputVector = DBL_MAX;
         std::vector<double> _array;

     public:
        double mean; // initialization always calculates mean so safe to access. unlike std which has a useless value until stdevStandardization is run 
        standardizeArray (std::vector<double> inputVector, long long nentries){
            _nentries=nentries;
            _array.reserve(nentries);
            mean=0;
            for (int i=0; i<_nentries; ++i){
                _array.push_back(inputVector[i]);
                mean+=inputVector[i];
            } 
            mean/=_nentries;
            //cout << "nentries: " << _nentries << endl;
            //cout << "Mean: " << mean << endl;
        }
        
        void rangeStandardization(){
            for (int ientry=0; ientry<_nentries; ++ientry){
                if (_array[ientry] > _max_inputVector){
                    _max_inputVector = _array[ientry];
                }
                if (_array[ientry] < _min_inputVector){
                    _min_inputVector = _array[ientry];
                }
            }
            //cout << "Max,min: " << _max_inputVector << "," << _min_inputVector << endl;
	    for (int ientry=0; ientry<_nentries; ++ientry){
	    	_array[ientry] = (_array[ientry]-_min_inputVector)/(_max_inputVector-_min_inputVector);
	    }
	    //cout << "Max,min "+name+": " << max_inputVector << "," << min_inputVector << endl;
	    //cout << "	Finished standardizing "+name << endl;
        }

        double calcStd(){
            double local_std=0;
            double diff=0;
            for (int ientry=0; ientry<_nentries; ++ientry){
                //cout << "mean: " << mean << endl;
                diff = (_array[ientry]-mean);
                //cout << "diff: " << diff << endl;
                local_std += diff*diff;
            }
            local_std /= _nentries;
            //cout << "STD: " << local_std << endl;
            return sqrt(local_std);
        }

        void stdevStandardization(){
            double std = calcStd();
            for (int ientry=0; ientry < _nentries; ++ientry){
               _array[ientry] = _array[ientry]/std; 
            } 
        }

        std::vector<double> getVector(){
            return _array;
        }

        void printVector(){
            for (int ientry=0; ientry<_nentries; ++ientry){
                cout << _array[ientry] << ", ";
            }
            cout << endl;
        }

};


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
// distSort_kNN will setup our priority_que that keeps a maximum of kDim elements simply by popping and pushing data. 
class compareDist
{
    public:
        bool operator() (const pair<double,int>& p1, const pair<double,int>& p2)
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
            cout << "Priority Queue is set up to sort the {distance,index} pairs keeping " << kDim << " neighbors" << endl;
            _kDim = kDim;
        }

        std::priority_queue <pair<double,int>, vector<pair<double,int>>, compareDist > kNN;
        
        // depending on the size and the distanceis we will push, or pop+push
        void insertPair ( pair<double,int> newPair ){
            if ( kNN.size() >= _kDim ){
                // > for min heap and maybe < for max heap 
                if (kNN.top().first > newPair.first) {
                    //cout << "\nPOPPING OUT " << kNN.top().first << endl;
                    kNN.pop();
                    kNN.push( newPair );
                }
            } 
            else {      
                //cout << "\nPUSHING " << newPair.first << endl;
                kNN.push( newPair );
            }
        }

    private:
        UInt_t _kDim;
        pair<double,int> _pair;
};



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
            //cout << "mean: " << _mean << endl;
            for (int i=0; i<_timesCalled; ++i){
                _diff = (inputVector[i]-_mean);
                _sum += _diff*_diff;
                //cout << "sum: " << _sum << endl;
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

























#endif
