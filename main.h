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

const int dim=3;
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

int numDOFsig = 6;
Double_t signal(Double_t *x, Double_t *par){
	return par[0]*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2])) + par[3]*exp(-0.5*((x[0]-par[4])/par[5])*((x[0]-par[4])/par[5]));
	//return (x[0]-par[0])*(x[0]-par[0]);

}

Double_t fitFunc(Double_t *x, Double_t *par){
	return background(x,par)+signal(x,&par[numDOFbkg]);
}

void standardizeArray(double inputVector[], int nentries, string name){
        // ERROR 1: not sure why this code did not work. It would sometimes give a max element that was not correct
	//double max_inputVector = *std::max(inputVector, inputVector+nentries);
	//double min_inputVector = *std::min(inputVector, inputVector+nentries);
        
        // ERROR 2: I tried to create a function that does this max,min calculation and returns an array. This was abit more difficult and quite buggy. Functions can return array points
        // but also would require static types or something like that to make it go outside the functions scope. This lead to wierd effects like calling the function twice might
        // pass the max, min from one function call to another. Safest is to just do it here. 

        // This should be fool proof
        double max_inputVector = DBL_MIN;
        double min_inputVector = DBL_MAX;
        for (int ientry=0; ientry<nentries; ++ientry){
            if (inputVector[ientry] > max_inputVector){
                max_inputVector = inputVector[ientry];
            }
            else if (inputVector[ientry] < min_inputVector){
                min_inputVector = inputVector[ientry];
            }
        }
        if(verbose_outputDistCalc){
	    cout << "Max,min "+name+": " << max_inputVector << "," << min_inputVector << endl;
	    cout << "	Finished standardizing "+name << endl;
        }

	for (int ientry=0; ientry<nentries; ++ientry){
		inputVector[ientry] = (inputVector[ientry]-min_inputVector)/(max_inputVector-min_inputVector);
	}
}


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



//double elapsedTime(high_resolution_clock start) {
//	return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() 
//}



























#endif
