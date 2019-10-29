#ifndef HELPER_H
#define HELPER_H

#include "main.h"

class getParamInit{
    private:
        std::string line;
        TF1* sigFit[2];
        TF1* bkgFit[2];
        std::vector<Double_t> _perFitParVals[2];
        std::vector<Double_t> _fullFitParVals[2];
        double _maxCount;
        int _dof;

    public:
        double scaleBkg[2];
        double scaleSig[2];
        double integralSig[2];
        double integralBkg[2];
        bool verbose=false;
        double fitMin[2];
        double fitMax[2];
        double binsize[2];
        std::vector<string> parNames[2];

        getParamInit ( int maxCount, int dof, double fitMin_eta, double fitMax_eta, double binsize_eta, double fitMin_pi0, double fitMax_pi0, double binsize_pi0){
            fitMin[0] = fitMin_eta;
            fitMin[1] = fitMin_pi0;
            fitMax[0] = fitMax_eta;
            fitMax[1] = fitMax_pi0;
            binsize[0] = binsize_eta;
            binsize[1] = binsize_pi0;
            _maxCount=maxCount;
            _dof=dof;

        }
        
        void loadData(){
            ifstream inFile[2];
            inFile[0].open("fitResults/etaFitNoAccSub.txt");
            inFile[1].open("fitResults/pi0FitNoAccSub.txt");
            
            for(int iFile=0; iFile<2; ++iFile){
                sigFit[iFile] = new TF1(("sigFit"+to_string(iFile)).c_str(),signalBW,fitMin[iFile],fitMax[iFile],numDOFsig);
                bkgFit[iFile] = new TF1(("bkgFit"+to_string(iFile)).c_str(),background,fitMin[iFile],fitMax[iFile],numDOFbkg);

                _fullFitParVals[iFile].reserve(_dof+1);
                _perFitParVals[iFile].reserve(_dof+1);
                parNames[iFile].reserve(_dof+1);
                while(std::getline(inFile[iFile], line)){
                    std::stringstream  lineStream(line);
                    string parName;
                    Double_t parVal;
                    lineStream >> parName;
                    lineStream >> parVal;
                    _fullFitParVals[iFile].push_back(parVal);    
                    parNames[iFile].push_back(parName);
                    cout << parName << " " << parVal << endl;
                }
                sigFit[iFile]->SetParameters(_fullFitParVals[iFile][3], _fullFitParVals[iFile][4],_fullFitParVals[iFile][5]);
                bkgFit[iFile]->SetParameters(_fullFitParVals[iFile][1], _fullFitParVals[iFile][2]);
                integralSig[iFile] = sigFit[iFile]->Integral(fitMin[iFile],fitMax[iFile]);
                integralBkg[iFile] = bkgFit[iFile]->Integral(fitMin[iFile],fitMax[iFile]);
                scaleBkg[iFile] = _maxCount/(integralBkg[iFile]/binsize[iFile]); 
                scaleSig[iFile] = _maxCount/(integralSig[iFile]/binsize[iFile]); 
		cout << "scaleBKG,SIG: " << scaleBkg[iFile] << ", " << scaleSig[iFile] << endl;
		cout << "_maxCount: " << _maxCount << endl;
		cout << "binsize[iFile]: " << binsize[iFile] << endl;
		cout << "integralBkg[iFile], integralSig[iFile]: " << integralBkg[iFile] << ", " << integralSig[iFile] << endl;
                if(verbose){cout << "Counts under curves for sig and bkg: " <<  integralSig[iFile]/binsize[iFile] << "," << integralBkg[iFile]/binsize[iFile] << endl;}
                _perFitParVals[iFile].push_back(_fullFitParVals[iFile][0]);
                _perFitParVals[iFile].push_back(_fullFitParVals[iFile][1]*scaleBkg[iFile]);
                _perFitParVals[iFile].push_back(_fullFitParVals[iFile][2]*scaleBkg[iFile]);
                _perFitParVals[iFile].push_back(_fullFitParVals[iFile][3]*scaleSig[iFile]);
                _perFitParVals[iFile].push_back(_fullFitParVals[iFile][4]);
                _perFitParVals[iFile].push_back(_fullFitParVals[iFile][5]);

                for (int iPar=0; iPar<numDOFbkg+numDOFsig+1; ++iPar){
                    if(verbose){cout << parNames[iFile][iPar] << " " << _fullFitParVals[iFile][iPar] <<" " << _perFitParVals[iFile][iPar] << endl;}
                }
                if(verbose){cout << "============================" << endl;}
            }
        }
        double getEta_par0(){ return _perFitParVals[0][1]; } 
        double getEta_par1(){ return _perFitParVals[0][2]; } 
        double getEta_par2(){ return _perFitParVals[0][3]; } 
        double getPi0_par0(){ return _perFitParVals[1][1]; } 
        double getPi0_par1(){ return _perFitParVals[1][2]; } 
        double getPi0_par2(){ return _perFitParVals[1][3]; } 
        double getEta_peak(){ return _perFitParVals[0][4]; }
        double getPi0_peak(){ return _perFitParVals[1][4]; }
        double getEta_width(){ return _perFitParVals[0][5]; }
        double getPi0_width(){ return _perFitParVals[1][5]; }
};
#endif
