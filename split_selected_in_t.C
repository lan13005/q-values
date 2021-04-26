#define split_selected_in_t_cxx
#include "split_selected_in_t.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// README!!!!!! IMPORTANT
// There are 3 connected files that we must use. split_selected_in_t.C/h and call_split_selected
// We use root -l call_split_selected to run the whole program. The file to be split is inside the
// header file and the values to split the  under is under the C file. There is some diagnostics
// plots that will be output so that we can check the t and mass distributions to make sure they look
// fine. 
// IMPORTANT TO NOTE: One thing to note is that event level values (i suspect) could be filled multiple times since this is a 
// tree, which is tened in combos. Thus the distributions might look weird but atleast we get an idea of
// what it should look like. Another things to note is that the  trees from the DSelector have allGeneralCuts applied
// to it whereas their comparisons like the pi0Mass and etaMass would not have the elliptical cut appleid. 
// So this program basically creates some files, trees, branches and reads the values in from the fchain loop
// and plots some diagnostic plots. The output of tihs program would be the split  that is ready to be
// input into amptools 



void split_selected_in_t::Loop()
{
//   In a ROOT session, you can do:
//      root> .L split_selected_in_t.C
//      root> split_selected_in_t t
//      root> t.GetEntry(12); // Fill t  members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TH1F *res_mass = new TH1F("h_resonance_mass", "resonance mass (GeV)", 175,0,3.5);
   TH2F *cosThetaVsResMass = new TH2F(",cosThetaVsResMass","CosTheta Vs Resonance Mass;resonance mass (GeV);cosTheta_Eta", 175,0,3.5,100,-1.00,1.00 );
   TH1F *psuedoScalar1_mass = new TH1F("h_psuedoScalar1_mass", "resonance mass (GeV)", 200,0.05,0.25);
   TH1F *psuedoScalar2_mass = new TH1F("h_psuedoScalar2_mass", "resonance mass (GeV)", 200,0.35,0.75);
   TH1F *tp_tot = new TH1F("h_tp_tot", "mandelstam tp (GeV)", 100,0,10);
   TH1F *tp_reg1 = new TH1F("h_tp_reg1", "mandelstam tp (GeV)", 100,0,10);
   TH1F *tp_reg2 = new TH1F("h_tp_reg2", "mandelstam tp (GeV)", 100,0,10);
   TH1F *tp_reg3 = new TH1F("h_tp_reg3", "mandelstam tp (GeV)", 100,0,10);
   TH1F *res_mass_t[3];
   TH2F *cosThetaVsResMass_t[3];
   TH1F *psuedoScalar1_mass_t[3];
   TH1F *psuedoScalar2_mass_t[3];
   string tag;
   for (int i=0; i<3; ++i){
      if (i==0) { tag = "tLT1"; }
      if (i==1) { tag = "tLT06"; }
      if (i==2) { tag = "tGT05LT1"; }
      res_mass_t[i] = new TH1F(("h_resonance_mass"+tag).c_str(), "resonance mass (GeV)", 175,0,3.5);
      cosThetaVsResMass_t[i]  = new TH2F(("cosThetaVsResMass"+tag).c_str(),"CosTheta Vs Resonance Mass;resonance mass (GeV);cosTheta_Eta", 175,0,3.5,100,-1.00,1.00 );
      psuedoScalar1_mass_t[i] = new TH1F(("h_psuedoScalar1_mass"+tag).c_str(), "resonance mass (GeV)", 200,0.05,0.25);
      psuedoScalar2_mass_t[i] = new TH1F(("h_psuedoScalar2_mass"+tag).c_str(), "resonance mass (GeV)", 200,0.35,0.75);
   }



   string fileNameTag = "deg000_data";
   outFile = new TFile((fileNameTag+"_pi0eta_2018_1_amptools.root").c_str(), "RECREATE");
   cout << "Making output tree file that is ready for amptools: " << (fileNameTag+"_pi0eta_2018_1_amptools.root").c_str() << endl;
   m_OutTree = new TTree("tree", "kin2");

   // Not so much diagnostic tree but rather allows me to grab the  for mandelstam_tp
   m_diagnostic = new TFile("holder.root", "RECREATE");
   diagnostic = new TTree("diagnostic","kin2");

   static size_t locNumFinalStateParticles = 3;
   n_LT1=0;
   n_LT06=0;
   n_GT05LT1=0;

// ************************** HERE WE DEFINE THE BRANCHES OF ALL THREE OUTPUT FILES ***************************
// ************************************************************************************************************
// ************************************************************************************************************
   m_OutTree->Branch("E_Beam", new float, "E_Beam/F");
   m_OutTree->Branch("Px_Beam", new float, "Px_Beam/F");
   m_OutTree->Branch("Py_Beam", new float, "Py_Beam/F");
   m_OutTree->Branch("Pz_Beam", new float, "Pz_Beam/F");
   m_OutTree->Branch("Target_Mass", new float, "Target_Mass/F");
   m_OutTree->Branch("NumFinalState", new int, "NumFinalState/I");
   m_OutTree->Branch("PID_FinalState", new int[locNumFinalStateParticles], "PID_FinalState[NumFinalState]/I");
   m_OutTree->Branch("E_FinalState", new float[locNumFinalStateParticles], "E_FinalState[NumFinalState]/F");
   m_OutTree->Branch("Px_FinalState", new float[locNumFinalStateParticles], "Px_FinalState[NumFinalState]/F");
   m_OutTree->Branch("Py_FinalState", new float[locNumFinalStateParticles], "Py_FinalState[NumFinalState]/F");
   m_OutTree->Branch("Pz_FinalState", new float[locNumFinalStateParticles], "Pz_FinalState[NumFinalState]/F");
   /// ****************** MAKE SURE THEY TYPES MATCH!!!!!!! CHECK THE TYPE BY LOOKING AT THE FLAT TREE USING PRINT()
   m_OutTree->Branch("Weight", new Float_t, "Weight/F");

   /// ****************** MAKE SURE THEY TYPES MATCH!!!!!!! CHECK THE TYPE BY LOOKING AT THE FLAT TREE USING PRINT()
   diagnostic->Branch("mandelstam_tp", new Double_t, "Mandelstam_tp/D");
// ************************************* SETTING BRANCH ADDRESSES *********************************************
// ************************************************************************************************************
// ************************************************************************************************************
   m_OutTree->SetBranchAddress("NumFinalState", &m_nPart);
   m_OutTree->SetBranchAddress("Target_Mass", &m_TargetMass);
   m_OutTree->SetBranchAddress("PID_FinalState", m_PID);
   m_OutTree->SetBranchAddress("E_FinalState", m_e);
   m_OutTree->SetBranchAddress("Px_FinalState", m_px);
   m_OutTree->SetBranchAddress("Py_FinalState", m_py);
   m_OutTree->SetBranchAddress("Pz_FinalState", m_pz);
   m_OutTree->SetBranchAddress("E_Beam", &m_eBeam);
   m_OutTree->SetBranchAddress("Px_Beam", &m_pxBeam);
   m_OutTree->SetBranchAddress("Py_Beam", &m_pyBeam);
   m_OutTree->SetBranchAddress("Pz_Beam", &m_pzBeam);
   m_OutTree->SetBranchAddress("Weight", &float_AccWeight);

   diagnostic->SetBranchAddress("mandelstam_tp", &mandelstam_tp);
// **************************************** FILLING THE VALUES! ***********************************************
// ************************************************************************************************************
// ************************************************************************************************************

// Define some constants across all events
   m_nPart = 3;
   m_TargetMass = 1*0.931494;          // Pb mass in GeV.
   m_PID[0] = 14; m_PID[1] = 7; m_PID[2] = 17;


   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) {
	cout << "ientry<0 ... breaking!" << endl;
	break;
      }
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      TLorentzVector resonance = *g1_p4_kin + *g2_p4_kin + *g3_p4_kin + *g4_p4_kin;
      TLorentzVector psuedoScalar_1 = *g1_p4_kin + *g2_p4_kin;
      TLorentzVector psuedoScalar_2 = *g3_p4_kin + *g4_p4_kin;
      // we will fill the proton, pi0, eta in this order since we already defined it that way in the PIDs section above. 
      m_e[0] = p_p4_kin->E();
      m_px[0] = p_p4_kin->Px();
      m_py[0] = p_p4_kin->Py();
      m_pz[0] = p_p4_kin->Pz();
      m_e[1] = psuedoScalar_1.E();
      m_px[1] = psuedoScalar_1.Px(); 
      m_py[1] = psuedoScalar_1.Py(); 
      m_pz[1] = psuedoScalar_1.Pz(); 
      m_e[2] =psuedoScalar_2.E(); 
      m_px[2]= psuedoScalar_2.Px(); 
      m_py[2]= psuedoScalar_2.Py(); 
      m_pz[2]= psuedoScalar_2.Pz(); 
      m_eBeam = beam_p4_kin->E();
      m_pxBeam = beam_p4_kin->Px();
      m_pyBeam = beam_p4_kin->Py();
      m_pzBeam = beam_p4_kin->Pz();
      m_tp = mandelstam_tp;

      Double_t resonanceMass = resonance.M();
      Double_t psuedoScalarMass_1 = psuedoScalar_1.M();
      Double_t psuedoScalarMass_2 = psuedoScalar_2.M();

      float_AccWeight = (Float_t)AccWeight;

      res_mass->Fill(resonanceMass,float_AccWeight);
      psuedoScalar1_mass->Fill(psuedoScalarMass_1,float_AccWeight);
      psuedoScalar2_mass->Fill(psuedoScalarMass_2,float_AccWeight);

	TLorentzVector dTargetP4 = TLorentzVector(TVector3(), m_TargetMass);
	TLorentzVector cm_vec = *beam_p4_kin+dTargetP4;
	TLorentzVector mixingPi0Eta_cm = psuedoScalar_1+psuedoScalar_2;
	TLorentzVector eta_cm = psuedoScalar_2;
        TLorentzVector beam_cm = *beam_p4_kin;
	mixingPi0Eta_cm.Boost(-cm_vec.BoostVector());
	eta_cm.Boost(-cm_vec.BoostVector());
	beam_cm.Boost(-cm_vec.BoostVector());
        TLorentzVector beam_res = beam_cm;
        TLorentzVector eta_res = eta_cm;
        eta_res.Boost(-mixingPi0Eta_cm.BoostVector());
        beam_res.Boost(-mixingPi0Eta_cm.BoostVector());
	TVector3 eta_res_unit = eta_res.Vect().Unit();	
        TVector3 z = beam_res.Vect().Unit();
        TVector3 y = mixingPi0Eta_cm.Vect().Cross(beam_cm.Vect()).Unit();
        TVector3 x = y.Cross(z).Unit();
	TVector3 angles_eta;
	angles_eta.SetXYZ ( eta_res_unit.Dot(x), eta_res_unit.Dot(y), eta_res_unit.Dot(z) );
	Double_t cosTheta_eta_GJ = angles_eta.CosTheta();
	cosThetaVsResMass->Fill(resonanceMass, cosTheta_eta_GJ, float_AccWeight);

      //if(float_AccWeight!=1){ cout << "accdouble_weightASBS != 1";} 
      tp_tot->Fill(m_tp,float_AccWeight);
      // Here we check what region of tp are we considering so we can od the amp analysis correctly split up
      if ((m_tp>0.5) && (m_tp<1)){
	n_GT05LT1++;
	tp_reg3->Fill(m_tp,float_AccWeight);
        res_mass_t[2]->Fill(resonanceMass, float_AccWeight);
        cosThetaVsResMass_t[2]->Fill(resonanceMass, cosTheta_eta_GJ, float_AccWeight);
        psuedoScalar1_mass_t[2]->Fill(psuedoScalarMass_1, float_AccWeight);
        psuedoScalar2_mass_t[2]->Fill(psuedoScalarMass_2, float_AccWeight);
      }
      if (m_tp<1){
	n_LT1++;
	tp_reg1->Fill(m_tp,float_AccWeight);
        res_mass_t[0]->Fill(resonanceMass, float_AccWeight);
        cosThetaVsResMass_t[0]->Fill(resonanceMass, cosTheta_eta_GJ, float_AccWeight);
        psuedoScalar1_mass_t[0]->Fill(psuedoScalarMass_1, float_AccWeight);
        psuedoScalar2_mass_t[0]->Fill(psuedoScalarMass_2, float_AccWeight);
      }
      if (m_tp<0.6){
	n_LT06++;
	tp_reg2->Fill(m_tp,float_AccWeight);
        res_mass_t[1]->Fill(resonanceMass, float_AccWeight);
        cosThetaVsResMass_t[1]->Fill(resonanceMass, cosTheta_eta_GJ, float_AccWeight);
        psuedoScalar1_mass_t[1]->Fill(psuedoScalarMass_1, float_AccWeight);
        psuedoScalar2_mass_t[1]->Fill(psuedoScalarMass_2, float_AccWeight);
      }
      m_OutTree->Fill();
   }
   // write out tree
   
   cout << "Number of events in tLT06: " << n_LT06 << endl;
   cout << "Number of events in tLT1: " << n_LT1 << endl;
   cout << "Number of events in tGT05LT1: " << n_GT05LT1 << endl;

   outFile->Write("tree");

   TFile *histFile = new TFile("AmpToolsDiagnosticHist.root", "RECREATE");
   cosThetaVsResMass->Write();
   res_mass->Write();
   psuedoScalar1_mass->Write();
   psuedoScalar2_mass->Write();
   for (int i=0; i<3; ++i){
      res_mass_t[i]->Write();
      cosThetaVsResMass_t[i]->Write(); 
      psuedoScalar1_mass_t[i]->Write(); 
      psuedoScalar2_mass_t[i]->Write();
   }
   tp_tot->Write();
   tp_reg1->Write();
   tp_reg2->Write();
   tp_reg3->Write();
   histFile->Close();
}
