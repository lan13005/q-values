#ifndef MAKEDIAGNOSTICHISTS_H
#define MAKEDIAGNOSTICHISTS_H
using namespace std;

void makeStackedHist(TH1F* tot, TH1F* sig, TH1F* bkg, string name,string baseDir){
    	TCanvas *allCanvases = new TCanvas(name.c_str(),"",1440,900);
	stackedHists = new THStack("stackedHists","");
	bkg->SetFillColorAlpha(kMagenta,0.5);
	bkg->SetLineColorAlpha(kMagenta,0);
	stackedHists->Add(tot,"HIST");
	stackedHists->Add(bkg,"HIST");
	stackedHists->Draw("nostack");
	stackedHists->GetXaxis()->SetTitle(tot->GetXaxis()->GetTitle());
	stackedHists->GetYaxis()->SetTitle(tot->GetYaxis()->GetTitle());
	stackedHists->SetTitle((tot->GetTitle()+" overlay total and bkg");
	allCanvases->SaveAs((baseDir+"/"+name+"_totBkg.png").c_str());
	
	allCanvases->Clear();
	sig->SetLineColorAlpha(kGreen+2,1);
	sig->SetLineWidth(2);
	sig->Draw("HIST");
	sig->SetAxisRange( 1.05*(sig->GetMinimum()) , 1.05*(tot->GetMaximum()), "Y" ); 
	sig->SetTitle((sig->GetTitle()+" signal");
	allCanvases->SaveAs((baseDir+"/"+name+"_sig.png").c_str());

	//allCanvases->Clear();
	//sig->Draw("HIST");
	//allCanvases->SaveAs(("diagnosticPlots/"+name+"_sig.png").c_str());
}


#endif


