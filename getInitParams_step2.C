#include "helperFuncs.h"

void getInitParams_step2(){
    //double binRange[3]={50,0.35,0.8};
    double binRange[3]={5000,0.35,0.8};
    double fitRange[3]={0.4,0.7};
    double binRange2[3]={50,0.05,0.25};
    double fitRange2[2]={0.1,0.17};
    int dof=numDOFbkg+numDOFsig;
    double binsize = (binRange[2]-binRange[1])/binRange[0];
    double binsize2 = (binRange2[2]-binRange2[1])/binRange2[0];
    //getParamInit(200,0, dof, true, fitRange[0], fitRange[1] ,binsize);
    //getParamInit(100,100, dof, true, fitRange[0], fitRange[1] ,binsize);
    //getParamInit(0,200, dof, true, fitRange[0], fitRange[1] ,binsize);
    //getParamInit(200,0, dof, false, fitRange2[0], fitRange2[1] ,binsize2);
    //getParamInit(100,100, dof, false, fitRange2[0], fitRange2[1] ,binsize2);
    //getParamInit(0,200, dof, false, fitRange2[0], fitRange2[1] ,binsize2);
    getParamInit paramInit = getParamInit ( 100, 6, 0.40, 0.7, binsize, 0.1, 0.17, binsize2);
    paramInit.loadData();
    std::vector<double> par0eta = {paramInit.getEta_par0(), paramInit.getEta_par0()/2, 0};
    std::vector<double> par1eta = {paramInit.getEta_par1(), paramInit.getEta_par1()/2, 0};
    std::vector<double> par2eta = {0, paramInit.getEta_par2()/2, paramInit.getEta_par2()};
    std::vector<double> par0pi0 = {paramInit.getPi0_par0(), paramInit.getPi0_par0()/2, 0};
    std::vector<double> par1pi0 = {paramInit.getPi0_par1(), paramInit.getPi0_par1()/2, 0};
    std::vector<double> par2pi0 = {0, paramInit.getPi0_par2()/2, paramInit.getPi0_par2()};
    std::vector<double> peakWidth_eta = {paramInit.getEta_peak(), paramInit.getEta_width()};
    std::vector<double> peakWidth_pi0 = {paramInit.getPi0_peak(), paramInit.getPi0_width()};

    std::vector< std::vector<double>  > allPars;
    std::vector <string> allParNames;
    allParNames.push_back("par0eta");
    allParNames.push_back("par1eta");
    allParNames.push_back("par2eta");
    allParNames.push_back("par0pi0");
    allParNames.push_back("par1pi0");
    allParNames.push_back("par2pi0");

    allPars.push_back(par0eta);
    allPars.push_back(par1eta);
    allPars.push_back(par2eta);
    allPars.push_back(par0pi0);
    allPars.push_back(par1pi0);
    allPars.push_back(par2pi0);

    cout << "\n\npeak width eta: " << peakWidth_eta[0] << " " << peakWidth_eta[1] << endl;
    cout << "peak width pi0: " << peakWidth_pi0[0] << " " << peakWidth_pi0[1] << endl;
    cout << "Max Bkg, 50/50 BkgSig, Max Sig" << endl;
    for(int iPar=0; iPar<allPars.size(); ++iPar){
        cout << allParNames[iPar] << ": " ;
        for(int iVal=0; iVal<allPars[iPar].size(); ++iVal){
            cout << allPars[iPar][iVal] << " ";
        }
        cout << endl;
    }


    //getParamInit paramInit_pi0 = getParamInit ( countSig,countBkg, dof, false, fitRange2[0], fitRange2[1] ,binsize2);
    //paramInit_eta

}
