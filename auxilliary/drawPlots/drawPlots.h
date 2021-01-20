#ifndef DRAWPLOTS_H
#define DRAWPLOTS_H

#include "RooAbsRealLValue.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TLine.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"

void drawPlots(
        RooAbsRealLValue* x, RooAbsRealLValue* y, double valX, double valY, int nentries, 
        RooAbsPdf* model, RooAbsPdf* bkg, RooAbsPdf* sig, 
        RooDataSet* data, RooAbsRealLValue* nsig, RooAbsRealLValue* nbkg,
        TCanvas* c, TH1* model_hist, TH1* model_sig, TH1* model_bkg
        );

#endif
