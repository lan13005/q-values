#ifndef GETINITPARAM_H
#define GETINITPARAM_H

#include "helperFuncs.h"

// Variables to load the data from your input tree
TTree *dataTree;
double Mpi0eta;
double AccWeight;

// Uniqueness tracking, the default way or using weighted by number of combos in an event 
bool isUniqueEtaB;
bool isUniquePi0B;
bool isUniquePi0EtaB;
double utWeight;
double sbWeight;

// Define your histograms
TH1F* massHistEta;
TH1F* massHistPi0;
TH1F* massHistPi0Eta;
TH2F* massHistPi0VsEta;
double weight; // weight to fill the histograms

// Define things to draw
TCanvas *allCanvases;

// -----------------------------
// Code for a single gaussian
static const int numDOFsig = 3; // degrees of freedom for the signal distribution
static const int numDOFbkg = 2; //
string namePar[numDOFbkg+numDOFsig] = {"bern01","bern11","amp","mass","sigma"};//,"ampRatio","sigmaRatio";
int meanIndices[1]={3}; // used in getInitParams 
int widthIndices[1]={4}; // used in getInitParams
double signal(double *x, double *par){
	return par[0]*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]));// + par[3]*exp(-0.5*((x[0]-par[4])/par[5])*((x[0]-par[4])/par[5]));

}
// Bernstein polynomial of degree 1
double background(double *x, double *par){
	return par[0]*x[0]/MAXVALUE+par[1]*(1-x[0]/MAXVALUE);
}
double fitFunc(double *x, double *par){
	return background(x,par)+signal(x,&par[numDOFbkg]);
}
// -----------------------------
// we need to tell "main" which variables need scaling when using the full reference distribution to the k nearest neighbors
std::vector<int> sigVarsNeedScaling = {2}; // these parameter indices are relative to fitFunc
std::vector<int> bkgVarsNeedScaling = {0,1};
////////////////////////////////////////////
// -----------------------------
// Code for a 2d gaussian
static const int numDOFsig2=5;
static const int numDOFbkg2=4;
string namePar2[numDOFbkg2+numDOFsig2] = {"bernx01","bernx11","berny01","berny11","amp","massx","sigmax","massy","sigmay"};//,"ampRatio","sigmaRatio";
int meanIndices2[2]={5,7}; // used in getInitParams 
int widthIndices2[2]={6,8}; // used in getInitParams
Double_t signal2(Double_t *x, Double_t *par) {
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	Double_t r2 = Double_t((x[1]-par[3])/par[4]);
	return par[0]/2/TMath::Pi()/par[2]/par[4]*TMath::Exp(-0.5*(r1*r1+r2*r2));
}

Double_t background2(Double_t *x, Double_t *par) {
	return par[0]*x[0]/MAXVALUE+par[1]*(1-x[0]/MAXVALUE)+par[2]*x[1]/MAXVALUE+par[3]*(1-x[1]/MAXVALUE);
}
Double_t fitFunc2(Double_t *x, Double_t *par) {
	return background2(x,par)+signal2(x,&par[numDOFbkg2]);
}
Double_t fitFunc2_projectVar2(Double_t *x, Double_t *par){
        // integrates along var1 so it is a function of var2
        // 0-3=bernstein coeffs, 4=2D gauss amp, {5,6} = {var1 mean,std}, {7,8} = {var2 mean,std}
        // WE ALSO HAVE ONLY 1 VARIABLE SO THERE IS ONLY X[0] AND NO X[1]. WHEN PROJECTING 2D DISTRIBUTION TO 1D WE WILL ONLY HAVE X[0] FOR BOTH X,Y
	Double_t r2 = Double_t((x[0]-par[7])/par[8]);
        Double_t binWidth=binWidthPi0;
        Double_t min=binRangePi0[1];
        Double_t max=binRangePi0[2];
        Double_t integratedBernstein = 0.5*(par[0]-par[1])*(max*max-min*min)+(par[1]+par[3]+(par[2]-par[3])*x[0])*(max-min);
        Double_t gaus = par[4]/sqrt(2*TMath::Pi())/par[8]*TMath::Exp(-0.5*(r2*r2));
	return (integratedBernstein+gaus)/binWidth; // need to scale by binWidth of the integrated dimension
} 
Double_t fitFunc2_projectVar1(Double_t *x, Double_t *par){
        // integrates along var2 so it is a function of var1
        // 0-3=bernstein coeffs, 4=2D gauss amp, {5,6} = {var1 mean,std}, {7,8} = {var2 mean,std}
	Double_t r1 = Double_t((x[0]-par[5])/par[6]);
        Double_t binWidth=binWidthEta;
        Double_t min=binRangeEta[1];
        Double_t max=binRangeEta[2];
        Double_t integratedBernstein = 0.5*(par[2]-par[3])*(max*max-min*min)+(par[3]+par[1]+(par[0]-par[1])*x[0])*(max-min);
        Double_t gaus = par[4]/sqrt(2*TMath::Pi())/par[6]*TMath::Exp(-0.5*(r1*r1));
	return (integratedBernstein+gaus)/binWidth;
} 
//Double_t background2_projectVar1(Double_t *x, Double_t *par){
//        Double_t min=binRangeEta[1];
//        Double_t max=binRangeEta[2];
//        Double_t binWidth=binWidthEta;
//        Double_t integratedBernstein = 0.5*(par[2]-par[3])*(max*max-min*min)+(par[3]+par[1]+(par[0]-par[1])*x[0])*(max-min);
//        return integratedBernstein/binWidth;
//}
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

#endif
