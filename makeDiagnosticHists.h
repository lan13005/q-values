#ifndef MAKEDIAGNOSTICHISTS_H
#define MAKEDIAGNOSTICHISTS_H
using namespace std;

// This will stack our histograms and make them pretty. 
void makeStackedHist(TH1F* tot, TH1F* sig, TH1F* bkg, string name,string baseDir){
    	TCanvas *allCanvases = new TCanvas(name.c_str(),"",1440,900);
	TLegend* leg1 = new TLegend(0.8,0.8,1.0,1.0);
	THStack* stackedHists = new THStack("stackedHists","");
	bkg->SetFillColorAlpha(kMagenta,0.5);
	bkg->SetLineColorAlpha(kMagenta,0);
	sig->SetLineColorAlpha(kGreen+2,1);
	sig->SetLineWidth(2);
	leg1->AddEntry(bkg,"Bkg","f");
	leg1->AddEntry(sig,"Sig","l");
	leg1->AddEntry(tot,"Tot","l");
	stackedHists->Add(tot,"HIST");
	stackedHists->Add(bkg,"HIST");
	stackedHists->Add(sig,"HIST");
	stackedHists->Draw("nostack");
	leg1->Draw();
	stackedHists->GetXaxis()->SetTitle(tot->GetXaxis()->GetTitle());
	stackedHists->GetYaxis()->SetTitle(tot->GetYaxis()->GetTitle());
	string totTitle = tot->GetTitle();
	stackedHists->SetTitle((totTitle+" overlay total and bkg").c_str());
	allCanvases->SaveAs((baseDir+"/"+name+"_totBkg.png").c_str());
}


#endif


