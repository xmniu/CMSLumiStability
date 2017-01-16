//================================================================================================
//
// Select probes for muon efficiencies with Tag&Probe method
//
//  * outputs ROOT file with a TTree of probes
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include "TLorentzVector.h"               // 4-vector class

// structure for output ntuple
#include "../Utils/EffData.hh" 
#include "../Utils/DatasetBounds.hh"      // runnum, lumisec bounds for dividing datasets

#endif

//=== MAIN MACRO ================================================================================================= 

void selectProbesMuEff_Data_25invpb(const TString infilename,           // input ntuple
                       const TString outputDir, 	   // output directory
		       const Int_t   effType, 	           // type of efficiency to compute
		       const Bool_t  doGenMatch = kFALSE,  // match to generator leptons
		       const Bool_t  doWeighted = kFALSE,   // store events with weights
		       const Int_t ibound=-1,
		       const Int_t isoType=1
) {
  gBenchmark->Start("selectProbesMuEff");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
 
  const Double_t TAG_PT_CUT = 25;
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  enum { eHLTEff, eL1Eff, eSelEff, eTrkEff, eStaEff, eStaEff_iso, ePOGIDEff, ePOGIsoEff, eSelIsoTrkEff, eHLTSelStaEff, eHLTSelStaEff_iso };
  if(effType > eHLTSelStaEff_iso) {
    cout << "Invalid effType option! Exiting..." << endl;
    return;
  }
  
  enum { eMuMu2HLT=1, eMuMu1HLT1L1, eMuMu1HLT, eMuMuNoSel, eMuSta, eMuTrk };  // event category enum
  
  Double_t nProbes = 0;
  
  //
  // Set up output ntuple
  //
  gSystem->mkdir(outputDir,kTRUE);
  TFile *outFile = new TFile(outputDir+TString("/probes.root"),"RECREATE"); 
  TTree *outTree = new TTree("Events","Events");
  //EffData data;
  //outTree->Branch("Events",&data.mass,"mass/F:pt:eta:phi:weight:q/I:npv/i:npu:pass:runNum:lumiSec:evtNum");
  Float_t mass, pt, eta, phi;
  Double_t weight;
  Int_t q;
  UInt_t npv, npu, passes, runNum, lumiSec, evtNum;
  outTree->Branch("mass",   &mass,   "mass/F");
  outTree->Branch("pt",     &pt,     "pt/F");
  outTree->Branch("eta",    &eta,    "eta/F");
  outTree->Branch("phi",    &phi,    "phi/F");
  outTree->Branch("weight", &weight, "weight/D");
  outTree->Branch("q",      &q,      "q/I");
  outTree->Branch("npv",    &npv,    "npv/i");
  outTree->Branch("npu",    &npu,    "npu/i");
  outTree->Branch("pass",   &passes, "pass/i");
  outTree->Branch("runNum", &runNum, "runNum/i");
  outTree->Branch("lumiSec",&lumiSec,"lumiSec/i");
  outTree->Branch("evtNum", &evtNum, "evtNum/i");

  //
  // Declare output ntuple variables
  //
  //UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  //UInt_t  npv, npu;
  Float_t scale1fb;
  Float_t genWeight, PUWeight;
  Float_t met, metPhi, sumEt, u1, u2;
  Int_t   q1, q2;
  TLorentzVector *dilep=0, *lep1=0, *lep2=0;
  TLorentzVector *sta1=0, *sta2=0;
  Float_t trkIso1, trkIso2; 
  Float_t emIso1, emIso2;
  Float_t hadIso1, hadIso2;
  Float_t pfChIso1, pfChIso2;
  Float_t pfGamIso1, pfGamIso2;
  Float_t pfNeuIso1, pfNeuIso2;
  Float_t pfCombIso1, pfCombIso2;
  Float_t d01, dz1, d02, dz2;
  Float_t muNchi21,  muNchi22;
  UInt_t nPixHits1, nTkLayers1, nPixHits2, nTkLayers2;
  UInt_t nValidHits1, nMatch1, nValidHits2, nMatch2;
  UInt_t typeBits1, typeBits2;

  // Read input file and get the TTrees
  cout << "Processing " << infilename << "..." << endl;
  TFile *infile = new TFile(infilename);	 assert(infile);
  TTree *intree = (TTree*)infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("runNum",     &runNum);      // event run number
  intree->SetBranchAddress("lumiSec",    &lumiSec);     // event lumi section
  intree->SetBranchAddress("evtNum",     &evtNum);      // event number
  intree->SetBranchAddress("matchGen",   &matchGen);    // event has both leptons matched to MC Z->ll
  intree->SetBranchAddress("category",   &category);    // dilepton category
  intree->SetBranchAddress("npv",        &npv);	        // number of primary vertices
  intree->SetBranchAddress("npu",        &npu);	        // number of in-time PU events (MC)
  intree->SetBranchAddress("genWeight",  &genWeight);
  intree->SetBranchAddress("PUWeight",   &PUWeight);
  intree->SetBranchAddress("scale1fb",   &scale1fb);    // event weight per 1/fb (MC)
  intree->SetBranchAddress("met",        &met);	        // MET
  intree->SetBranchAddress("metPhi",     &metPhi);      // phi(MET)
  intree->SetBranchAddress("sumEt",      &sumEt);       // Sum ET
  intree->SetBranchAddress("u1",         &u1);	        // parallel component of recoil
  intree->SetBranchAddress("u2",         &u2);	        // perpendicular component of recoil
  intree->SetBranchAddress("q1",         &q1);	        // charge of tag lepton
  intree->SetBranchAddress("q2",         &q2);	        // charge of probe lepton
  intree->SetBranchAddress("dilep",      &dilep);       // dilepton 4-vector
  intree->SetBranchAddress("lep1",       &lep1);        // tag lepton 4-vector
  intree->SetBranchAddress("lep2",       &lep2);        // probe lepton 4-vector
  intree->SetBranchAddress("sta1",       &sta1);        // tag STA muon 4-vector
  intree->SetBranchAddress("sta2",       &sta2);        // probe STA muon 4-vector 
  intree->SetBranchAddress("trkIso1",    &trkIso1);       // track isolation of tag lepton
  intree->SetBranchAddress("trkIso2",    &trkIso2);       // track isolation of probe lepton
  intree->SetBranchAddress("emIso1",     &emIso1);        // ECAL isolation of tag lepton
  intree->SetBranchAddress("emIso2",     &emIso2);        // ECAL isolation of probe lepton
  intree->SetBranchAddress("hadIso1",    &hadIso1);       // HCAL isolation of tag lepton
  intree->SetBranchAddress("hadIso2",    &hadIso2);       // HCAL isolation of probe lepton
  intree->SetBranchAddress("pfChIso1",   &pfChIso1);      // PF charged hadron isolation of tag lepton
  intree->SetBranchAddress("pfChIso2",   &pfChIso2);      // PF charged hadron isolation of probe lepton
  intree->SetBranchAddress("pfGamIso1",  &pfGamIso1);     // PF photon isolation of tag lepton
  intree->SetBranchAddress("pfGamIso2",  &pfGamIso2);     // PF photon isolation of probe lepton
  intree->SetBranchAddress("pfNeuIso1",  &pfNeuIso1);     // PF neutral hadron isolation of tag lepton
  intree->SetBranchAddress("pfNeuIso2",  &pfNeuIso2);     // PF neutral hadron isolation of probe lepton
  intree->SetBranchAddress("pfCombIso1", &pfCombIso1);  // PF combined isolation of tag lepton
  intree->SetBranchAddress("pfCombIso2", &pfCombIso2);  // PF combined isolation of probe lepton    
  intree->SetBranchAddress("d01",        &d01); 	// transverse impact parameter of tag lepton
  intree->SetBranchAddress("d02",	 &d02); 	// transverse impact parameter of probe lepton      
  intree->SetBranchAddress("dz1",	 &dz1); 	// longitudinal impact parameter of tag lepton
  intree->SetBranchAddress("dz2",	 &dz2); 	// longitudinal impact parameter of probe lepton    
  intree->SetBranchAddress("muNchi21",	 &muNchi21);	// muon fit normalized chi^2 of tag lepton
  intree->SetBranchAddress("muNchi22",	 &muNchi22);	// muon fit normalized chi^2 of probe lepton	 
  intree->SetBranchAddress("nPixHits1",  &nPixHits1);   // number of pixel hits of tag muon
  intree->SetBranchAddress("nPixHits2",  &nPixHits2);   // number of pixel hits of probe muon
  intree->SetBranchAddress("nTkLayers1", &nTkLayers1);  // number of tracker layers of tag muon
  intree->SetBranchAddress("nTkLayers2", &nTkLayers2);  // number of tracker layers of probe muon
  intree->SetBranchAddress("nMatch1",	 &nMatch1);	// number of matched segments of tag muon
  intree->SetBranchAddress("nMatch2",	 &nMatch2);	// number of matched segments of probe muon    
  intree->SetBranchAddress("nValidHits1",&nValidHits1); // number of valid muon hits of tag muon
  intree->SetBranchAddress("nValidHits2",&nValidHits2); // number of valid muon hits of probe muon
  intree->SetBranchAddress("typeBits1",  &typeBits1);   // muon type of tag muon
  intree->SetBranchAddress("typeBits2",  &typeBits2);   // muon type of probe muon

  //
  // loop over events
  //
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    Bool_t doSelection = 0;

        if(ibound != NB-2){
                if(run_bound[ibound] != run_bound[ibound+1]){
                        if(runNum == run_bound[ibound] && lumiSec >= ls_bound[ibound]) doSelection = 1;
                        if(runNum > run_bound[ibound] && runNum < run_bound[ibound+1]) doSelection = 1;
                        if(runNum == run_bound[ibound+1] && lumiSec < ls_bound[ibound+1]) doSelection = 1;

                }else if(run_bound[ibound] == run_bound[ibound+1]){
                        if(runNum == run_bound[ibound] && lumiSec >= ls_bound[ibound] && lumiSec < ls_bound[ibound+1]) doSelection = 1;

                }else{
                        cout<<"ERROR! A"<<endl;
                        break;
                }
        }else if(ibound == NB-2){
                if(runNum == run_bound[ibound] && lumiSec >= ls_bound[ibound] && lumiSec <= ls_bound[ibound+1]) doSelection = 1;
        }else{
                cout<<"ERROR! B"<<endl;
                break;
        }

        if(doSelection != 1) continue;


    if(lep1->Pt() < TAG_PT_CUT) continue;
    
    // check GEN match if necessary
    if(doGenMatch && !matchGen) continue;
    
    Bool_t  pass=kFALSE;
    Float_t m=0;
    
    if(effType==eHLTEff) {
      //
      // probe = muon passing selection
      // pass  = matched to HLT
      // * MuMu2HLT event means a passing probe, MuMu1HLT(1L1) event means a failing probe
      //   all other categories do not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kFALSE; }
      else if(category==eMuMu1HLT)    { pass=kFALSE; }
      else if(category==eMuMuNoSel)   { continue; }
      else if(category==eMuSta)       { continue; }
      else                            { continue; }
      
      m = dilep->M();

    } else if(effType==eL1Eff) {
      //
      // probe = muon passing selection
      // pass  = matched to HLT
      // * MuMu2HLT, MuMu1HLT1L1 event means a passing probe, MuMu1HLT event means a failing probe
      //   all other categories do not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kFALSE; }
      else if(category==eMuMuNoSel)   { continue; }
      else if(category==eMuSta)       { continue; }
      else                            { continue; }
      
      m = dilep->M();
    
    } else if(effType==eSelEff) {
      //
      // probe = GLB muon
      // pass  = passing selection
      // * MuMu2HLT, MuMu1HLT(1L1) event means a passing probe, MuMuNoSel event means a failing probe,
      //   all other categories do not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kTRUE; }
      else if(category==eMuMuNoSel)   { pass=kFALSE; }
      else if(category==eMuSta)       { continue; }
      else                            { continue; }
      
      m = dilep->M();
    
    } else if(effType==eTrkEff) {
      //
      // probe = STA muon
      // pass  = is also a GLB muon
      // * MuMu2HLT, MuMu1HLT(1L1), MuMuNoSel event means a passing probe, MuSta event means a failing probe, 
      //   MuTrk event does not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kTRUE; }
      else if(category==eMuMuNoSel)   { pass=kTRUE; }
      else if(category==eMuSta)       { pass=kFALSE; }
      else                            { continue; }
      
      // compute mass using probe STA muon pT
      TLorentzVector tp = *lep1 + *sta2;
      m = tp.M();    
    
    } else if(effType==eStaEff) {
      //
      // probe = tracker track
      // pass  = is also a GLB muon
      // * MuMu2HLT, MuMu1HLT(1L1), MuMuNoSel event means a passing probe, MuTrk event means a failing probe, 
      //   MuSta event does not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kTRUE; }
      else if(category==eMuMuNoSel)   { pass=kTRUE; }
      else if(category==eMuSta)       { continue; }
      else                            { pass=kFALSE; }
      
      m = dilep->M();
    
    } else if(effType==eStaEff_iso) {
      //
      // probe = isolated tracker track
      // pass  = is also a GLB muon
      // * MuMu2HLT, MuMu1HLT(1L1), isolated MuMuNoSel event means a passing probe, isolated MuTrk event means a failing probe, 
      //   MuSta, non-isolated MuMuNoSel, and non-isolated MuTrk events do not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kTRUE; }
      else if(category==eMuMuNoSel)   { pass=kTRUE; }
      else if(category==eMuSta)       { continue; }
      else                            { if(pfCombIso2>0.12*(lep2->Pt())) continue; else pass=kFALSE; }
      
      m = dilep->M();
    } else if(effType==eSelIsoTrkEff) {

      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kTRUE; }
      else if(category==eMuMuNoSel)   { pass=kFALSE; }
      else if(category==eMuSta)       { pass=kFALSE; }
      else                            { continue; }
/*
      if     (category==eMuMu2HLT)    {
        if(isoType==0) pass=kTRUE;
        else if(isoType==1) pass=(pfCombIso2<0.15*lep2->Pt());
        else if(isoType==2) pass=(pfCombIso2<0.25*lep2->Pt());
        else if(isoType==3) pass=(trkIso2<0.05*lep2->Pt());
        else if(isoType==4) pass=(trkIso2<0.10*lep2->Pt());
      } 
      else if(category==eMuMu1HLT1L1) {
        if(isoType==0) pass=kTRUE;
        else if(isoType==1) pass=(pfCombIso2<0.15*lep2->Pt());
        else if(isoType==2) pass=(pfCombIso2<0.25*lep2->Pt());
        else if(isoType==3) pass=(trkIso2<0.05*lep2->Pt());
        else if(isoType==4) pass=(trkIso2<0.10*lep2->Pt());
     } 
      else if(category==eMuMu1HLT)    {
	if(isoType==0) pass=kTRUE;
	else if(isoType==1) pass=(pfCombIso2<0.15*lep2->Pt());
	else if(isoType==2) pass=(pfCombIso2<0.25*lep2->Pt());
	else if(isoType==3) pass=(trkIso2<0.05*lep2->Pt()); 
	else if(isoType==4) pass=(trkIso2<0.10*lep2->Pt());
      }
      else if(category==eMuMuNoSel)   { pass=kFALSE; }
      else if(category==eMuSta)       { pass=kFALSE; }
      else                            { continue; }
*/
      m = dilep->M();

    }
 else if(effType==eHLTSelStaEff_iso) {
      //
      // probe = isolated tracker track
      // pass  = is also a GLB muon
      // * MuMu2HLT, MuMu1HLT(1L1), isolated MuMuNoSel event means a passing probe, isolated MuTrk event means a failing probe, 
      //   MuSta, non-isolated MuMuNoSel, and non-isolated MuTrk events do not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kFALSE; }
      else if(category==eMuMu1HLT)    { pass=kFALSE; }
      else if(category==eMuMuNoSel)   { if(pfCombIso2>0.12*(lep2->Pt())) continue; else pass=kFALSE; }
      else if(category==eMuSta)       { continue; }
      else                            { if(pfCombIso2>0.12*(lep2->Pt())) continue; else pass=kFALSE; }
      
      m = dilep->M();    
    } else if(effType==eHLTSelStaEff) {
      //
      // probe = tracker track
      // pass = is also a GLB muon
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kFALSE; }
      else if(category==eMuMu1HLT)    { pass=kFALSE; }
      else if(category==eMuMuNoSel)   { pass=kTRUE; }
      else if(category==eMuSta)       { continue; }
      else                            { pass=kFALSE; }
      
      m = dilep->M();
    } else if(effType==ePOGIDEff) {
      //
      // probe = tracker track
      // pass  = passes "tight" muon ID
      // * 
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kTRUE; }
      else if(category==eMuMuNoSel)   { 
        
	pass=kTRUE; 
        if(nTkLayers2  < 6)   pass=kFALSE;
        if(nPixHits2   < 1)   pass=kFALSE;
        if(fabs(d02)   > 0.2) pass=kFALSE;
        if(fabs(dz2)   > 0.5) pass=kFALSE;
        if(muNchi22    > 10)  pass=kFALSE;
        if(nMatch2     < 2)   pass=kFALSE;
        if(nValidHits2 < 1)   pass=kFALSE;
        if(!(typeBits2 & 1))  pass=kFALSE;
        if(!(typeBits2 & 8))  pass=kFALSE;
      }
      else if(category==eMuSta)     { continue; }
      else                          { pass=kFALSE; }
      
      m = dilep->M();    
    
    } else if(effType==ePOGIsoEff) {
      //
      // probe = "tight" muon
      // pass  = passes isolation
      // * 
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kTRUE; }
      else if(category==eMuMuNoSel)   {
        
	if(nTkLayers2  < 6)   continue;
        if(nPixHits2   < 1)   continue;
        if(fabs(d02)   > 0.2) continue;
        if(fabs(dz2)   > 0.5) continue;
        if(muNchi22    > 10)  continue;
        if(nMatch2     < 2)   continue;
        if(nValidHits2 < 1)   continue;
        if(!(typeBits2 & 1))  continue;
        if(!(typeBits2 & 8))  continue;

	pass = (pfCombIso2 < 0.12*(lep2->Pt()));      
      }
      else if(category==eMuSta)     { continue; }
      else                          { continue; }
      
      m = dilep->M();    
    }
    
    nProbes += doWeighted ? genWeight*PUWeight/std::abs(genWeight) : 1;

    // Fill tree
    mass = m;
    pt	 = (effType==eTrkEff) ? sta2->Pt()  : lep2->Pt();
    eta	 = (effType==eTrkEff) ? sta2->Eta() : lep2->Eta();
    phi	 = (effType==eTrkEff) ? sta2->Phi() : lep2->Phi();
    weight  = doWeighted ? genWeight*PUWeight/std::abs(genWeight) : 1;
    q	    = q2;
    npv	    = npv;
    npu	    = npu;
    passes  = (pass) ? 1 : 0;
    runNum  = runNum;
    lumiSec = lumiSec;
    evtNum  = evtNum;
    outTree->Fill();
    
    if(category==eMuMu2HLT) {
      if(lep2->Pt() < TAG_PT_CUT) continue;

      nProbes += doWeighted ? genWeight*PUWeight/std::abs(genWeight) : 1;

      mass	   = m;
      pt	   = (effType==eTrkEff) ? sta1->Pt()  : lep1->Pt();
      eta	   = (effType==eTrkEff) ? sta1->Eta() : lep1->Eta();
      phi	   = (effType==eTrkEff) ? sta1->Phi() : lep1->Phi();
      weight  = doWeighted ? genWeight*PUWeight/std::abs(genWeight) : 1;
      q	   = q1;
      npv	   = npv;
      npu	   = npu;
      passes	   = 1;
      runNum  = runNum;
      lumiSec = lumiSec;
      evtNum  = evtNum;
      outTree->Fill();
    } 
  }  
  delete infile;
  infile=0, intree=0;	   


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;

  cout << " Number of probes selected: " << nProbes << endl;
  
  outFile->Write();
  outFile->Close();
  delete outFile;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectProbesMuEff"); 
}
