#ifndef MAKEDIAGNOSTICHISTS_H
#define MAKEDIAGNOSTICHISTS_H

#include <ctime>
#include <math.h> 
#include "helperFuncs.h"

using namespace std;

// This will stack our histograms and make them pretty. 
void makeStackedHist(TH1F* tot, TH1F* sig, TH1F* bkg, TH1F* sig_sb, TH1F* bkg_sb, string name,string baseDir){
    	TCanvas *allCanvases = new TCanvas(name.c_str(),"",1440,900);
        allCanvases->Divide(2,2);
	TLegend* leg1 = new TLegend(0.8,0.8,1.0,1.0);
	TLegend* leg1_sb = new TLegend(0.8,0.8,1.0,1.0);
	TLegend* leg1_overlay = new TLegend(0.8,0.8,1.0,1.0);
	THStack* stackedHists = new THStack("stackedHists","");
	THStack* stackedHists_sb = new THStack("stackedHists_sb","");
	THStack* stackedHists_overlay = new THStack("stackedHists_overlay","");

	bkg->SetFillColorAlpha(kViolet-5,0.3);
	bkg->SetLineColorAlpha(kViolet-5,0);
	sig->SetLineColorAlpha(kRed+1,1);
	sig->SetLineWidth(2);

	bkg_sb->SetFillColorAlpha(kViolet-5,0.3);
	bkg_sb->SetLineColorAlpha(kViolet-5,0);
	sig_sb->SetLineColorAlpha(kBlue+1,1);
	sig_sb->SetLineWidth(2);

	leg1->AddEntry(bkg,"Q_Bkg","f");
	leg1->AddEntry(sig,"Q_Sig","l");
	leg1->AddEntry(tot,"Tot","l");
	leg1_sb->AddEntry(bkg_sb,"SB_Bkg","f");
	leg1_sb->AddEntry(sig_sb,"SB_Sig","l");
	leg1_sb->AddEntry(tot,"Tot","l");
	leg1_overlay->AddEntry(sig,"Q_Sig","l");
	leg1_overlay->AddEntry(sig_sb,"SB_Sig","l");

        allCanvases->cd(1);
	stackedHists->Add(tot,"HIST");
	stackedHists->Add(bkg,"HIST");
	stackedHists->Add(sig,"HIST");
	stackedHists->Draw("nostack");
	leg1->Draw();
	stackedHists->GetXaxis()->SetTitle(tot->GetXaxis()->GetTitle());
	stackedHists->GetYaxis()->SetTitle(tot->GetYaxis()->GetTitle());
	//string totTitle = tot->GetTitle();
	//stackedHists->SetTitle((totTitle+" overlay total and bkg").c_str());

        allCanvases->cd(2);
	stackedHists_sb->Add(tot,"HIST");
	stackedHists_sb->Add(bkg_sb,"HIST");
	stackedHists_sb->Add(sig_sb,"HIST");
	stackedHists_sb->Draw("nostack");
	leg1_sb->Draw();
	stackedHists_sb->GetXaxis()->SetTitle(tot->GetXaxis()->GetTitle());
	stackedHists_sb->GetYaxis()->SetTitle(tot->GetYaxis()->GetTitle());
	//string totTitle = tot->GetTitle();
	//stackedHists_sb->SetTitle((totTitle+" overlay total and bkg").c_str());

        allCanvases->cd(3);
        stackedHists_overlay->Add(sig,"HIST");
        stackedHists_overlay->Add(sig_sb,"HIST");
	stackedHists_overlay->Draw("nostack");
	stackedHists_overlay->GetXaxis()->SetTitle(tot->GetXaxis()->GetTitle());
	stackedHists_overlay->GetYaxis()->SetTitle(tot->GetYaxis()->GetTitle());
        leg1_overlay->Draw();
	//string totTitle = tot->GetTitle();
	//stackedHists_sb->SetTitle((totTitle+" ").c_str());

	allCanvases->SaveAs((baseDir+"/"+name+"_totBkg.png").c_str());
}


#endif


