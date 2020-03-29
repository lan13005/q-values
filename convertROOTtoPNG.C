void subdetector_convertROOTtoPNG(string inDir){
	cout << "Looking in directory: " << inDir << endl;
	char* dir = gSystem->ExpandPathName(inDir.c_str());
	void* dirp = gSystem->OpenDirectory(dir);
	
	const char* entry;
	static const int maxFiles = 500;
	const char* filename[maxFiles];
	Int_t n = 0;
	TString str;
	
	while((entry = (char*)gSystem->GetDirEntry(dirp))) {
		str = entry;
		if(str.EndsWith(".root")){
		 	filename[n++] = gSystem->ConcatFileName(dir, entry);
			cout << "Adding file: " << str << endl;
		}
		if (n>maxFiles) {break;}
	}
	cout << "------------------\nThere are " << n << " files added" << endl; 
	
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
   		   		TCanvas *h = (TCanvas*)key->ReadObj();
				Printf("Drawing canvas: %s", h->GetName());
				string fileString = (string)filename[i];
				h->Draw();
				h->SaveAs((fileString.substr(0,fileString.size()-4)+"png").c_str());
			}
		}
		histFile->Close();	
	}

	cout << "\n\n" << endl;
	cout << "rm " << inDir << "/*.root" << endl;
	gSystem->Exec(("rm "+inDir+"/*.root").c_str());
}

void convertROOTtoPNG(){
	subdetector_convertROOTtoPNG("histograms/fcal");
	subdetector_convertROOTtoPNG("histograms/bcal");
	subdetector_convertROOTtoPNG("histograms/split");
}
