void convertROOTtoPNG(){
	const char* inDir = "histograms/bcal";
	
	char* dir = gSystem->ExpandPathName(inDir);
	void* dirp = gSystem->OpenDirectory(dir);
	
	const char* entry;
	const char* filename[100];
	Int_t n = 0;
	TString str;
	
	while((entry = (char*)gSystem->GetDirEntry(dirp))) {
		str = entry;
		if(str.EndsWith(".root"))
		 	filename[n++] = gSystem->ConcatFileName(dir, entry);
	}
	
	TFile* histFile;
	TH1F* hist;
	for (Int_t i = 0; i < n; i++){
	  	Printf("Opening file -> %s", filename[i]);
		histFile = TFile::Open(filename[i]);
		TIter keyList(histFile->GetListOfKeys());
		TKey *key;

   		while ((key = (TKey*)keyList())) {
   		   	TClass *cl = gROOT->GetClass(key->GetClassName());
   		   	if (cl->InheritsFrom("TCanvas")){
   		   		TCanvas* h = (TCanvas*)key->ReadObj();
				h->Draw();	
				string fileString = (string)filename[i];
				h->SaveAs((fileString.substr(0,fileString.size()-4)+"png").c_str());
			}
		}
		
		histFile->Close();	
	}
}
