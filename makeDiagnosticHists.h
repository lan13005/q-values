#ifndef MAKEDIAGNOSTICHISTS_H
#define MAKEDIAGNOSTICHISTS_H
using namespace std;

// This will stack our histograms and make them pretty. 
void makeStackedHist(TH1F* tot, TH1F* sig, TH1F* bkg, TH1F* sig_sb, TH1F* bkg_sb, string name,string baseDir){
    	TCanvas *allCanvases = new TCanvas(name.c_str(),"",1440,900);
	TLegend* leg1 = new TLegend(0.8,0.8,1.0,1.0);
	THStack* stackedHists = new THStack("stackedHists","");
	bkg->SetFillColorAlpha(kMagenta,0.3);
	bkg->SetLineColorAlpha(kMagenta,0);
	sig->SetLineColorAlpha(kRed+1,1);
	sig->SetLineWidth(2);

	bkg_sb->SetFillColorAlpha(kGreen,0.3);
	bkg_sb->SetLineColorAlpha(kGreen,0);
	sig_sb->SetLineColorAlpha(kCyan+1,1);
	sig_sb->SetLineWidth(2);

	leg1->AddEntry(bkg,"Q_Bkg","f");
	leg1->AddEntry(sig,"Q_Sig","l");
	leg1->AddEntry(bkg_sb,"SB_Bkg","f");
	leg1->AddEntry(sig_sb,"SB_Sig","l");
	leg1->AddEntry(tot,"Tot","l");

	stackedHists->Add(tot,"HIST");
	stackedHists->Add(bkg,"HIST");
	stackedHists->Add(sig,"HIST");
	stackedHists->Add(bkg_sb,"HIST");
	stackedHists->Add(sig_sb,"HIST");
	stackedHists->Draw("nostack");
	leg1->Draw();
	stackedHists->GetXaxis()->SetTitle(tot->GetXaxis()->GetTitle());
	stackedHists->GetYaxis()->SetTitle(tot->GetYaxis()->GetTitle());
	string totTitle = tot->GetTitle();
	stackedHists->SetTitle((totTitle+" overlay total and bkg").c_str());
	allCanvases->SaveAs((baseDir+"/"+name+"_totBkg.png").c_str());
}


#endif


