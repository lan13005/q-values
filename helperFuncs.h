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
//bool verbose=true;
string rootFileLoc="/d/grid13/ln16/q-values-2/allMC_tree_ext.root";
string rootTreeName="degALL_acc_mEllipse_tree_flat";
string fileTag="all";
string runTag="";
string weightingScheme="as";
string s_accWeight="AccWeight";
string s_sbWeight="weightBS";
string s_discrimVar="Mpi0;Meta";
string s_utBranch="";
string s_mcprocessBranch="mcprocess";

///// **************************************************************
///// STEP0.5: DEFINE HISTOGRAM BIN PARAMETERS SO THERE IS CONSISTENCY
///// **************************************************************
std::vector<double> binRangeEta={200,0.34,0.8};
std::vector<double> binRangePi0={200,0.075,0.21};
//std::vector<double> binRangeEta={200,0.25,0.85};
//std::vector<double> binRangePi0={200,0.005,0.25};
double binWidthEta=(binRangeEta[2]-binRangeEta[1])/binRangeEta[0];
double binWidthPi0=(binRangePi0[2]-binRangePi0[1])/binRangePi0[0];

std::vector<double> fitRangeEta2={0.36,0.75};//{0.4,0.7};//{0.4,0.7};
std::vector<double> fitRangePi02={0.085,0.185};//{0.09,0.18};//{0.09,0.18};
//std::vector<double> fitRangeEta2={0.25,0.85};//{0.4,0.7};
//std::vector<double> fitRangePi02={0.005,0.25};//{0.09,0.18};

// used in getInitParams
std::vector<double> fitRange1=fitRangePi02;//{0.121,0.15};
std::vector<double> fitRange2=fitRangeEta2;//{0.516,0.58};


///// **************************************************************
///// STEP1: DEFINE SIGNAL/BKG DISTRIBUTIONS
///// **************************************************************

// We need to define a maximum value for the bernstein polynomial bkg. The domain is [0,1] so we have to scale the discriminating variable to be between 0,1. 
// Since we know the maximum value for Meta is < 1 then everything is fine. For larger masses we can divide by a larger value or just dividing by Ebeam_max = 12 is fine. 
// Without much calculation we know GlueX cannot produce anything with mass > 12 GeV. If you use other variables like then choose the appropriate maximum value
#define MAXVALUE 1

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
        parseVarString(){}; 
        void updateString(std::string inputString){
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









#endif
