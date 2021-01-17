#include "/d/grid13/ln16/q-values-2/helperFuncs.h"
#include "/d/grid13/ln16/q-values-2/auxilliary/drawPlots/drawPlots.C"

using namespace RooFit;
void checkGausFeedThrough(double initSigFrac){
    initSigFrac /= 10; // bash doesnt accept non integer increments in loops so we multiplied by 10

    RooRandom::randomGenerator()->SetSeed(1992);
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    // ---------------------------
    // RooFit Variables
    // ---------------------------
    // Some vars to initialize but not declare yet
    double sigPdfVal;
    double bkgPdfVal;
    double totPdfVal;
    double qvalue;
    
    // Fit parameters 
    double fittedMassX = 0.135059;
    double fittedMassY = 0.549354;
    double fittedSigmaX = 0.00612767;
    double fittedSigmaY = 0.0166494;
    double fittedBernA = 0.5;
    double fittedBernB = 0.5;
    double fittedBernC = 0.5;
    double fittedBernD = 0.5;
    double fittedBernE = 0.5;
    
    // For loading the data
    int iProcess=0;
    string sThread = to_string(iProcess);
    RooRealVar roo_Meta(("roo_Meta"+sThread).c_str(),"Mass GeV",fitRangeEta2[0],fitRangeEta2[1]);
    roo_Meta.setRange(("roo_fitRangeMeta"+sThread).c_str(),fitRangeEta2[0], fitRangeEta2[1]);
    RooRealVar roo_Mpi0(("roo_Mpi0"+sThread).c_str(),"Mass GeV",fitRangePi02[0],fitRangePi02[1]);
    roo_Mpi0.setRange(("roo_fitRangeMpi0"+sThread).c_str(),fitRangePi02[0], fitRangePi02[1]);
    RooRealVar roo_Weight(("roo_Weight"+sThread).c_str(), "Weight", -10, 10); // Weights can take a wide range
    roo_Meta.setBins(100);
    roo_Mpi0.setBins(100);
    
    RooRealVar peak_pi0(("peak_pi0"+sThread).c_str(),"peak_pi0",fittedMassX);//,fittedMassX,fittedMassX);
    RooRealVar width_pi0(("width_pi0"+sThread).c_str(),"width_pi0",fittedSigmaX,fittedSigmaX*0.5,fittedSigmaX*1.5);
    RooRealVar peak_eta(("peak_eta"+sThread).c_str(),"peak_eta",fittedMassY);//,fittedMassY,fittedMassY);
    RooRealVar width_eta(("width_eta"+sThread).c_str(),"width_eta",fittedSigmaY,fittedSigmaY*0.5,fittedSigmaY*1.5);
    
    RooRealVar bern_parA(("bern_parA"+sThread).c_str(),"bern_parA",fittedBernA,0,1);
    RooRealVar bern_parB(("bern_parB"+sThread).c_str(),"bern_parB",fittedBernB,0,1);
    RooRealVar bern_parC(("bern_parC"+sThread).c_str(),"bern_parC",fittedBernC,0,1);
    RooRealVar bern_parD(("bern_parD"+sThread).c_str(),"bern_parD",fittedBernD,0,1);
    //RooRealVar bern_parE(("bern_parE"+sThread).c_str(),"bern_parE",fittedBernE,0,1);///,fittedBernE,0,1);
    RooGenericPdf rooBkgX(("rooBkgX"+sThread).c_str(), "rooBkgX", ("bern_parA"+sThread+"*roo_Mpi0"+sThread+"+bern_parB"+sThread+"*(1-roo_Mpi0"+sThread+")").c_str(),RooArgSet(bern_parA,bern_parB,roo_Mpi0));
    //RooGenericPdf rooBkgY(("rooBkgY"+sThread).c_str(), "rooBkgY", ("bern_parC"+sThread+"*(1-roo_Meta"+sThread+")**2+bern_parD"+sThread+"*2*roo_Meta"+sThread+"*(1-roo_Meta"+    sThread+")+bern_parE"+sThread+"*(roo_Meta"+sThread+")**2").c_str(),RooArgSet(bern_parC,bern_parD,bern_parE,roo_Meta));
    RooGenericPdf rooBkgY(("rooBkgY"+sThread).c_str(), "rooBkgY", ("bern_parC"+sThread+"*roo_Meta"+sThread+"+bern_parD"+sThread+"*(1-roo_Meta"+sThread+")").c_str(),RooArgSet(bern_parC,bern_parD,roo_Meta));
    RooGaussian rooGausPi0_bkg(("rooGausPi0_bkg"+sThread).c_str(), "rooGausPi0_bkg", roo_Mpi0, peak_pi0, width_pi0);
    RooRealVar bkgPeakFrac(("bkgPeakFrac"+sThread).c_str(),"bkgPeakFrac",1,0,1);
    RooAddPdf rooBkgXplusPi0Peak(("rooBkgXplusPi0Peak"+sThread).c_str(), "rooBkgXplusPi0Peak", RooArgList(rooGausPi0_bkg,rooBkgX),RooArgSet(bkgPeakFrac));
    RooProdPdf rooBkg(("rooBkg"+sThread).c_str(),"rooBkg",RooArgList(rooBkgXplusPi0Peak,rooBkgY));

    RooGaussian rooGausPi0(("rooGausPi0"+sThread).c_str(), "rooGausPi0", roo_Mpi0, peak_pi0, width_pi0);
    RooGaussian rooGausEta(("rooGausEta"+sThread).c_str(), "rooGausEta", roo_Meta, peak_eta, width_eta);
    RooProdPdf rooGaus2D(("rooGaus2D"+sThread).c_str(), "rooGaus2D", RooArgSet(rooGausPi0,rooGausEta));

    RooRealVar nsig(("nsig"+sThread).c_str(),"nsig",0,1000);
    RooRealVar nbkg(("nbkg"+sThread).c_str(),"nbkg",0,1000);
    RooAddPdf rooSigPlusBkg(("rooSumPdf"+sThread).c_str(), "rooSumPdf", RooArgList(rooGaus2D,rooBkg),RooArgSet(nsig,nbkg));

    TCanvas* allCanvases = new TCanvas(("anyHists"+to_string(iProcess)).c_str(),"",1440,900);
    RooFitResult* roo_result;
    RooDataSet* rooData;


    int kDim;
    double sigFrac;
    double inputTotalSig;
    double outputTotalSig;
    double fitSigFrac=1;
    double initbkgPeakFrac=0;
    int totalEvnts=500;
    cout << "** totalEvnts kDim initSigFrac inputTotalSig outputTotalSig" << endl;
    std::vector<double> kDims={100,500,1000};
    for (auto kDim: kDims){
        inputTotalSig=0;
        outputTotalSig=0;
        for (int iter=0; iter<totalEvnts; ++iter){
            inputTotalSig += initSigFrac;
            allCanvases->Divide(3,2);
            peak_pi0.setVal(fittedMassX);
            peak_eta.setVal(fittedMassY);
            width_pi0.setVal(fittedSigmaX);
            width_eta.setVal(fittedSigmaY);
            bern_parA.setVal(fittedBernA);
            bern_parB.setVal(fittedBernB);
            bern_parC.setVal(fittedBernC);
            bern_parD.setVal(fittedBernD);
            bkgPeakFrac.setVal(initbkgPeakFrac);
            nsig.setVal(initSigFrac*kDim);
            nbkg.setVal((1-initSigFrac)*kDim);
            rooData=rooSigPlusBkg.generate(RooArgSet(roo_Mpi0,roo_Meta),kDim); 
            ////////////////////////////////////////////////////////////////////////////////////////////
            //  BEGIN 2D GAUS FIT
            ////////////////////////////////////////////////////////////////////////////////////////////
            nsig.setVal(fitSigFrac*kDim);
            nbkg.setVal((1-fitSigFrac)*kDim);
            roo_result = rooSigPlusBkg.fitTo(*rooData, Save(), Minos(kTRUE), Extended(kTRUE), RooFit::SumW2Error(true), BatchMode(kTRUE), PrintLevel(-1), PrintEvalErrors(-1));
            // setting parameters for bkg/sig and extracting q-value
            sigFrac = nsig.getVal()/(nsig.getVal()+nbkg.getVal());
            outputTotalSig += sigFrac;
            //cout << "sigfrac, nsig, nkg, ntot: " << sigFrac << ", "<< nsig.getVal() << ", " << nbkg.getVal() << ", " <<nsig.getVal()+nbkg.getVal() << endl;

            //-------------------
            // Doesnt do anything
            roo_Mpi0.setVal(0.5);
            roo_Meta.setVal(0.5);
            sigPdfVal = sigFrac*rooGaus2D.getVal(RooArgSet(roo_Mpi0,roo_Meta));
            bkgPdfVal = (1-sigFrac)*rooBkg.getVal(RooArgSet(roo_Mpi0,roo_Meta));
            totPdfVal = rooSigPlusBkg.getVal(RooArgSet(roo_Mpi0,roo_Meta));
            qvalue = sigPdfVal/(sigPdfVal+bkgPdfVal);
            //-------------------
                                        
            //drawPlots(&roo_Mpi0, &roo_Meta, 0.5, 0.5, kDim, &rooSigPlusBkg, &rooBkg, &rooGaus2D, rooData, &nsig, &nbkg, allCanvases, sThread);
            //allCanvases->SaveAs(("iter"+to_string(iter)+".png").c_str());
            //allCanvases->Clear();
        }
        cout << "** " << totalEvnts << " " << kDim << " " << initSigFrac << " " << inputTotalSig << " " << outputTotalSig << endl;
        
    }
}
