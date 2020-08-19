#ifndef GETINITPARAM_H
#define GETINITPARAM_H

#include "helperFuncs.h"

// Variables to load the data from your input tree
TTree *dataTree;
double var1;
double var2;
double Mpi0eta;
double AccWeight;
double sbWeight;


// Uniqueness tracking, the default way or using weighted by number of combos in an event 
bool isUniqueEtaB;
bool isUniquePi0B;
bool isUniquePi0EtaB;
double utWeight;

// Define your histograms
TH1F* massHistEta;
TH1F* massHistPi0;
TH1F* massHistPi0Eta;
TH2F* massHistPi0VsEta;
double weight; // weight to fill the histograms

// Variables related to fitting the discriminating variable
TF1* fit;
TF1* bkgFit;
TF1* sigFit;
TF2* fit2D;
TF2* bkgFit2D;
TF2* sigFit2D;
TF1* proj_var1;
TF1* proj_var2; 
double par[numDOFbkg+numDOFsig];
double par2D[numDOFbkg2+numDOFsig2];
double par2D_proj1[numDOFbkg2+numDOFsig2+3];
double par2D_proj2[numDOFbkg2+numDOFsig2+3];
std::vector<double> fitRange1;
std::vector<double> fitRange2;
		
// Define things to draw
TCanvas *allCanvases;


#endif
