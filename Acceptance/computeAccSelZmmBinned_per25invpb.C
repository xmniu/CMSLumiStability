//================================================================================================
//
// Compute Z->mumu acceptance at full selection level
//
//  * outputs results summary text file
//
//  [!!!] propagation of efficiency scale factor uncertainties no yet implemented
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TH1D.h>                   // histogram class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "TLorentzVector.h"         // 4-vector class

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions

// helper class to handle efficiency tables
#include "../Utils/CEffUser1D.hh"
#include "../Utils/CEffUser2D.hh"
#endif

//=== MAIN MACRO ================================================================================================= 

void computeAccSelZmmBinned_per25invpb(const TString conf,      // input file
			    const TString purwDir,
			    const TString dteffDir,
                            const TString mceffDir,
                            const TString outputDir,  // output directory
			    const Int_t   doPU,
			    const Int_t   triType=0,
			    const Int_t   isoType=1,
			    const Int_t   ibound=-1
) {
  gBenchmark->Start("computeAccSelZmmBinned");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW   = 66;
  const Double_t MASS_HIGH  = 116;
  const Double_t PT_CUT     = 27;
  const Double_t ETA_CUT    = 2.4;
  const Double_t ETA_BARREL = 0.9;
  const Double_t ETA_ENDCAP = 2.4;
  const Double_t MUON_MASS  = 0.105658369;

  const Int_t BOSON_ID  = 23;
  const Int_t LEPTON_ID = 13;
  
  // efficiency files
  const TString dataHLTEffName_pos = dteffDir+"MuHLTEff/2MG/eff.root";
  const TString dataHLTEffName_neg = dteffDir+"MuHLTEff/2MG/eff.root";
  const TString zmmHLTEffName_pos  = mceffDir+"MuHLTEff/2CT/eff.root";
  const TString zmmHLTEffName_neg  = mceffDir+"MuHLTEff/2CT/eff.root";

  const TString dataSelEffName_pos = dteffDir+"MuSITEff/2MG/eff.root";
  const TString dataSelEffName_neg = dteffDir+"MuSITEff/2MG/eff.root";
  const TString zmmSelEffName_pos  = mceffDir+"MuSITEff/2CT/eff.root";
  const TString zmmSelEffName_neg  = mceffDir+"MuSITEff/2CT/eff.root";

  const TString dataTrkEffName_pos = dteffDir+"MuStaEff/2MG/eff.root";
  const TString dataTrkEffName_neg = dteffDir+"MuStaEff/2MG/eff.root";
  const TString zmmTrkEffName_pos  = mceffDir+"MuStaEff/2CT/eff.root";
  const TString zmmTrkEffName_neg  = mceffDir+"MuStaEff/2CT/eff.root";

  const TString dataStaEffName_pos = dteffDir+"MuStaEff/2MG/eff.root";
  const TString dataStaEffName_neg = dteffDir+"MuStaEff/2MG/eff.root";
  const TString zmmStaEffName_pos  = mceffDir+"MuStaEff/2CT/eff.root";
  const TString zmmStaEffName_neg  = mceffDir+"MuStaEff/2CT/eff.root";

  // load pileup reweighting file
  char rw_80X_file[200];
  sprintf(rw_80X_file, "%s/ReweightHistogram/pileup_rw_80X_%d.root", purwDir.Data(), ibound);
  TFile *f_rw_80X = TFile::Open(rw_80X_file, "read");
  TH1D *h_rw = (TH1D*) f_rw_80X->Get("h_rw");

 
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString> fnamev;  // file name per input file
  vector<TString> labelv;  // TLegend label per input file
  vector<Int_t>   colorv;  // plot color per input file
  vector<Int_t>   linev;   // plot line style per input file

  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    
    string fname;
    Int_t color, linesty;
    stringstream ss(line);
    ss >> fname >> color >> linesty;
    string label = line.substr(line.find('@')+1);
    fnamev.push_back(fname);
    labelv.push_back(label);
    colorv.push_back(color);
    linev.push_back(linesty);
  }
  ifs.close();

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);  

  TH2D *h=0;

  //
  // HLT efficiency
  //
  cout << "Loading trigger efficiencies..." << endl;
  
  TFile *dataHLTEffFile_pos = new TFile(dataHLTEffName_pos);
  CEffUser2D dataHLTEff_pos;
  dataHLTEff_pos.loadEff((TH2D*)dataHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataHLTEffFile_neg = new TFile(dataHLTEffName_neg);
  CEffUser2D dataHLTEff_neg;
  dataHLTEff_neg.loadEff((TH2D*)dataHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrhEtaPt"));
    
  TFile *zmmHLTEffFile_pos = new TFile(zmmHLTEffName_pos);
  CEffUser2D zmmHLTEff_pos;
  zmmHLTEff_pos.loadEff((TH2D*)zmmHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmHLTEffFile_neg = new TFile(zmmHLTEffName_neg);
  CEffUser2D zmmHLTEff_neg;
  zmmHLTEff_neg.loadEff((TH2D*)zmmHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmHLTEffFile_neg->Get("hErrhEtaPt"));

  h =(TH2D*)dataHLTEffFile_pos->Get("hEffEtaPt");
  TH2D *hHLTErr_pos = new TH2D("hHLTErr_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hHLTErr_neg = new TH2D("hHLTErr_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hHLTErrB_pos = new TH2D("hHLTErrB_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hHLTErrB_neg = new TH2D("hHLTErrB_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hHLTErrE_pos = new TH2D("hHLTErrE_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hHLTErrE_neg = new TH2D("hHLTErrE_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  
  //
  // Selection efficiency
  //
  cout << "Loading selection efficiencies..." << endl;
  
  TFile *dataSelEffFile_pos = new TFile(dataSelEffName_pos);
  CEffUser2D dataSelEff_pos;
  dataSelEff_pos.loadEff((TH2D*)dataSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataSelEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataSelEffFile_neg = new TFile(dataSelEffName_neg);
  CEffUser2D dataSelEff_neg;
  dataSelEff_neg.loadEff((TH2D*)dataSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataSelEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmSelEffFile_pos = new TFile(zmmSelEffName_pos);
  CEffUser2D zmmSelEff_pos;
  zmmSelEff_pos.loadEff((TH2D*)zmmSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmSelEffFile_pos->Get("hErrhEtaPt"));

  TFile *zmmSelEffFile_neg = new TFile(zmmSelEffName_neg);
  CEffUser2D zmmSelEff_neg;
  zmmSelEff_neg.loadEff((TH2D*)zmmSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmSelEffFile_neg->Get("hErrhEtaPt"));

  h =(TH2D*)dataSelEffFile_pos->Get("hEffEtaPt");
  TH2D *hSelErr_pos = new TH2D("hSelErr_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hSelErr_neg = new TH2D("hSelErr_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hSelErrB_pos = new TH2D("hSelErrB_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hSelErrB_neg = new TH2D("hSelErrB_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hSelErrE_pos = new TH2D("hSelErrE_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hSelErrE_neg = new TH2D("hSelErrE_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());

  //
  // Standalone efficiency
  //
  cout << "Loading standalone efficiencies..." << endl;
  
  TFile *dataStaEffFile_pos = new TFile(dataStaEffName_pos);
  CEffUser2D dataStaEff_pos;
  dataStaEff_pos.loadEff((TH2D*)dataStaEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataStaEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataStaEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataStaEffFile_neg = new TFile(dataStaEffName_neg);
  CEffUser2D dataStaEff_neg;
  dataStaEff_neg.loadEff((TH2D*)dataStaEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataStaEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataStaEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmStaEffFile_pos = new TFile(zmmStaEffName_pos);
  CEffUser2D zmmStaEff_pos;
  zmmStaEff_pos.loadEff((TH2D*)zmmStaEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmStaEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmStaEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmStaEffFile_neg = new TFile(zmmStaEffName_neg);
  CEffUser2D zmmStaEff_neg;
  zmmStaEff_neg.loadEff((TH2D*)zmmStaEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmStaEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmStaEffFile_neg->Get("hErrhEtaPt"));

  h =(TH2D*)dataStaEffFile_pos->Get("hEffEtaPt");
  TH2D *hStaErr_pos = new TH2D("hStaErr_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hStaErr_neg = new TH2D("hStaErr_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hStaErrB_pos = new TH2D("hStaErrB_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hStaErrB_neg = new TH2D("hStaErrB_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hStaErrE_pos = new TH2D("hStaErrE_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hStaErrE_neg = new TH2D("hStaErrE_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());

  //
  // Tracker efficiency
  //
  cout << "Loading track efficiencies..." << endl;
  
  TFile *dataTrkEffFile_pos = new TFile(dataTrkEffName_pos);
  CEffUser2D dataTrkEff_pos;
  dataTrkEff_pos.loadEff((TH2D*)dataTrkEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataTrkEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataTrkEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataTrkEffFile_neg = new TFile(dataTrkEffName_neg);
  CEffUser2D dataTrkEff_neg;
  dataTrkEff_neg.loadEff((TH2D*)dataTrkEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataTrkEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataTrkEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmTrkEffFile_pos = new TFile(zmmTrkEffName_pos);
  CEffUser2D zmmTrkEff_pos;
  zmmTrkEff_pos.loadEff((TH2D*)zmmTrkEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmTrkEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmTrkEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmTrkEffFile_neg = new TFile(zmmTrkEffName_neg);
  CEffUser2D zmmTrkEff_neg;
  zmmTrkEff_neg.loadEff((TH2D*)zmmTrkEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmTrkEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmTrkEffFile_neg->Get("hErrhEtaPt"));

  h =(TH2D*)dataTrkEffFile_pos->Get("hEffEtaPt");
  TH2D *hTrkErr_pos = new TH2D("hTrkErr_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hTrkErr_neg = new TH2D("hTrkErr_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());

  // Data structures to store info from TTrees
  baconhep::TEventInfo *info   = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr     = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *muonArr        = new TClonesArray("baconhep::TMuon");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;
   
  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv;
  vector<Double_t> nSelCorrv, nSelCorrVarv;
  vector<Double_t> accv, accCorrv;
  vector<Double_t> accErrv, accErrCorrv;

  vector<Double_t> nSelBBv;
  vector<Double_t> nSelBBCorrv, nSelBBCorrVarv;
  vector<Double_t> accBBv, accBBCorrv;
  vector<Double_t> accErrBBv, accErrBBCorrv;

  vector<Double_t> nSelEEv;
  vector<Double_t> nSelEECorrv, nSelEECorrVarv;
  vector<Double_t> accEEv, accEECorrv;
  vector<Double_t> accErrEEv, accErrEECorrv;

  const baconhep::TTrigger triggerMenu("/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_8_0_14/src/BaconAna/DataFormats/data/HLT_50nsGRun_80X");
  
  //
  // loop through files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << " ..." << endl;
    infile = TFile::Open(fnamev[ifile]); 
    assert(infile);
  
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
    eventTree->SetBranchAddress("Info",             &info); TBranch *infoBr = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("GenEvtInfo",        &gen); TBranch *genBr  = eventTree->GetBranch("GenEvtInfo");
    eventTree->SetBranchAddress("GenParticle",&genPartArr); TBranch* genPartBr = eventTree->GetBranch("GenParticle");
    eventTree->SetBranchAddress("Muon",          &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");   
    eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");

    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelCorrv.push_back(0);
    nSelCorrVarv.push_back(0);

    nSelBBv.push_back(0);
    nSelBBCorrv.push_back(0);
    nSelBBCorrVarv.push_back(0);

    nSelEEv.push_back(0);
    nSelEECorrv.push_back(0);
    nSelEECorrVarv.push_back(0);



    //
    // loop over events
    //      
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      genBr->GetEntry(ientry);
      genPartArr->Clear(); genPartBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);

      Int_t glepq1=-99;
      Int_t glepq2=-99;

      if (fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
      TLorentzVector *vec=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
      toolbox::fillGen(genPartArr, BOSON_ID, vec, lep1, lep2,&glepq1,&glepq2,1);
      if(vec->M()<MASS_LOW || vec->M()>MASS_HIGH) continue;
      delete vec; delete lep1; delete lep2;

      vertexArr->Clear();
      vertexBr->GetEntry(ientry);
      double npv  = vertexArr->GetEntries();
      Double_t weight=gen->weight;
      if(doPU>0) weight*=h_rw->GetBinContent(h_rw->FindBin(info->nPUmean));

      nEvtsv[ifile]+=weight;
      
      // trigger requirement               
      if (!isMuonTrigger(triggerMenu, info->triggerBits, triType)) continue;

      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
    
      muonArr->Clear();
      muonBr->GetEntry(ientry);

      for(Int_t i1=0; i1<muonArr->GetEntriesFast(); i1++) {
  	const baconhep::TMuon *mu1 = (baconhep::TMuon*)((*muonArr)[i1]);

        if(mu1->pt	  < PT_CUT)  continue;  // lepton pT cut
        if(fabs(mu1->eta) > ETA_CUT) continue;  // lepton |eta| cut
        if(!passMuonID(mu1, 0, isoType))	     continue;  // lepton selection
	
	TLorentzVector vMu1(0,0,0,0);
	vMu1.SetPtEtaPhiM(mu1->pt, mu1->eta, mu1->phi, MUON_MASS);

        for(Int_t i2=i1+1; i2<muonArr->GetEntriesFast(); i2++) {
          const baconhep::TMuon *mu2 = (baconhep::TMuon*)((*muonArr)[i2]);
        
          if(mu1->q == mu2->q)	       continue;  // opposite charge requirement
          if(mu2->pt        < PT_CUT)  continue;  // lepton pT cut
          if(fabs(mu2->eta) > ETA_CUT) continue;  // lepton |eta| cut
	  if(!passMuonID(mu2, 0, isoType))	       continue;  // lepton selection

          TLorentzVector vMu2(0,0,0,0);
	  vMu2.SetPtEtaPhiM(mu2->pt, mu2->eta, mu2->phi, MUON_MASS);  

          // trigger match
	  if(!isMuonTriggerObj(triggerMenu, mu1->hltMatchBits, kFALSE, triType) && !isMuonTriggerObj(triggerMenu, mu2->hltMatchBits, kFALSE, triType)) continue;
	  
	  // mass window
          TLorentzVector vDilep = vMu1 + vMu2;
          if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue;
          
          /******** We have a Z candidate! HURRAY! ********/

	  Bool_t isBB = (fabs(mu1->eta)<ETA_BARREL && fabs(mu2->eta)<ETA_BARREL) ? kTRUE : kFALSE;
	  Bool_t isEE = (fabs(mu1->eta)>ETA_BARREL && fabs(mu2->eta)>ETA_BARREL && fabs(mu1->eta)<ETA_ENDCAP && fabs(mu2->eta)<ETA_ENDCAP) ? kTRUE : kFALSE;

          Double_t effdata, effmc;
	  Double_t corr=1;
	  
	  effdata=1; effmc=1;    
          if(mu1->q>0) { 
            effdata *= (1.-dataHLTEff_pos.getEff(fabs(mu1->eta), mu1->pt)); 
            effmc   *= (1.-zmmHLTEff_pos.getEff(fabs(mu1->eta), mu1->pt)); 
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff(fabs(mu1->eta), mu1->pt)); 
            effmc   *= (1.-zmmHLTEff_neg.getEff(fabs(mu1->eta), mu1->pt)); 
          }
          if(mu2->q>0) {
            effdata *= (1.-dataHLTEff_pos.getEff(fabs(mu2->eta), mu2->pt)); 
            effmc   *= (1.-zmmHLTEff_pos.getEff(fabs(mu2->eta), mu2->pt));
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff(fabs(mu2->eta), mu2->pt)); 
            effmc   *= (1.-zmmHLTEff_neg.getEff(fabs(mu2->eta), mu2->pt));
          }
          effdata = 1.-effdata;
          effmc   = 1.-effmc;
          corr *= effdata/effmc;
    
          effdata=1; effmc=1;
          if(mu1->q>0) { 
            effdata *= dataSelEff_pos.getEff(fabs(mu1->eta), mu1->pt); 
            effmc   *= zmmSelEff_pos.getEff(fabs(mu1->eta), mu1->pt); 
          } else {
            effdata *= dataSelEff_neg.getEff(fabs(mu1->eta), mu1->pt); 
            effmc   *= zmmSelEff_neg.getEff(fabs(mu1->eta), mu1->pt); 
          }
          if(mu2->q>0) {
            effdata *= dataSelEff_pos.getEff(fabs(mu2->eta), mu2->pt); 
            effmc   *= zmmSelEff_pos.getEff(fabs(mu2->eta), mu2->pt);
          } else {
            effdata *= dataSelEff_neg.getEff(fabs(mu2->eta), mu2->pt); 
            effmc   *= zmmSelEff_neg.getEff(fabs(mu2->eta), mu2->pt);
          }
          corr *= effdata/effmc;
    
          effdata=1; effmc=1;
          if(mu1->q>0) { 
            effdata *= dataStaEff_pos.getEff(fabs(mu1->eta), mu1->pt); 
            effmc   *= zmmStaEff_pos.getEff(fabs(mu1->eta), mu1->pt); 
          } else {
            effdata *= dataStaEff_neg.getEff(fabs(mu1->eta), mu1->pt); 
            effmc   *= zmmStaEff_neg.getEff(fabs(mu1->eta), mu1->pt); 
          }
          if(mu2->q>0) {
            effdata *= dataStaEff_pos.getEff(fabs(mu2->eta), mu2->pt); 
            effmc   *= zmmStaEff_pos.getEff(fabs(mu2->eta), mu2->pt);
          } else {
            effdata *= dataStaEff_neg.getEff(fabs(mu2->eta), mu2->pt); 
            effmc   *= zmmStaEff_neg.getEff(fabs(mu2->eta), mu2->pt);
          }
          corr *= effdata/effmc;
    
          effdata=1; effmc=1;
          if(mu1->q>0) { 
            effdata *= dataTrkEff_pos.getEff(fabs(mu1->eta), mu1->pt); 
            effmc   *= zmmTrkEff_pos.getEff(fabs(mu1->eta), mu1->pt); 
          } else {
            effdata *= dataTrkEff_neg.getEff(fabs(mu1->eta), mu1->pt); 
            effmc   *= zmmTrkEff_neg.getEff(fabs(mu1->eta), mu1->pt); 
          }
          if(mu2->q>0) {
            effdata *= dataTrkEff_pos.getEff(fabs(mu2->eta), mu2->pt); 
            effmc   *= zmmTrkEff_pos.getEff(fabs(mu2->eta), mu2->pt);
          } else {
            effdata *= dataTrkEff_neg.getEff(fabs(mu2->eta), mu2->pt); 
            effmc   *= zmmTrkEff_neg.getEff(fabs(mu2->eta), mu2->pt);
          }
          //corr *= effdata/effmc;
	  
	  // scale factor uncertainties                                                                                                                                         
	  // TRACKER
          if(mu1->q>0) {
            Double_t effdata = dataTrkEff_pos.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errdata = TMath::Max(dataTrkEff_pos.getErrLow(fabs(mu1->eta), mu1->pt), dataTrkEff_pos.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t effmc   = zmmTrkEff_pos.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errmc   = TMath::Max(zmmTrkEff_pos.getErrLow(fabs(mu1->eta), mu1->pt), zmmTrkEff_pos.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t errTrk = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hTrkErr_pos->Fill(fabs(mu1->eta), mu1->pt, errTrk);
          } else {
            Double_t effdata = dataTrkEff_neg.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errdata = TMath::Max(dataTrkEff_neg.getErrLow(fabs(mu1->eta), mu1->pt), dataTrkEff_neg.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t effmc   = zmmTrkEff_neg.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errmc   = TMath::Max(zmmTrkEff_neg.getErrLow(fabs(mu1->eta), mu1->pt), zmmTrkEff_neg.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t errTrk = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hTrkErr_neg->Fill(fabs(mu1->eta), mu1->pt, errTrk);
          }

          if(mu2->q>0) {
            Double_t effdata = dataTrkEff_pos.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errdata = TMath::Max(dataTrkEff_pos.getErrLow(fabs(mu2->eta), mu2->pt), dataTrkEff_pos.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t effmc   = zmmTrkEff_pos.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errmc   = TMath::Max(zmmTrkEff_pos.getErrLow(fabs(mu2->eta), mu2->pt), zmmTrkEff_pos.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t errTrk = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    /*   if(fabs(mu2->eta)>1.2 && fabs(mu2->eta)<2.1) 
	      {
		errTrk=0.0013;
		}*/
            hTrkErr_pos->Fill(fabs(mu2->eta), mu2->pt, errTrk);
          } else {
            Double_t effdata = dataTrkEff_neg.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errdata = TMath::Max(dataTrkEff_neg.getErrLow(fabs(mu2->eta), mu2->pt), dataTrkEff_neg.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t effmc   = zmmTrkEff_neg.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errmc   = TMath::Max(zmmTrkEff_neg.getErrLow(fabs(mu2->eta), mu2->pt), zmmTrkEff_neg.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t errTrk = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    /*  if(fabs(mu2->eta)>1.2 && fabs(mu2->eta)<2.1) 
	      {
		errTrk=0.0013;
	      }
	      hTrkErr_neg->Fill(fabs(mu2->eta), mu2->pt, errTrk);*/
          }
	  // STANDALONE
          if(mu1->q>0) {
            Double_t effdata = dataStaEff_pos.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errdata = TMath::Max(dataStaEff_pos.getErrLow(fabs(mu1->eta), mu1->pt), dataStaEff_pos.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t effmc   = zmmStaEff_pos.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errmc   = TMath::Max(zmmStaEff_pos.getErrLow(fabs(mu1->eta), mu1->pt), zmmStaEff_pos.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hStaErr_pos->Fill(fabs(mu1->eta), mu1->pt, errSta);
            if(isBB) hStaErrB_pos->Fill(fabs(mu1->eta), mu1->pt, errSta);
            else if(isEE) hStaErrE_pos->Fill(fabs(mu1->eta), mu1->pt, errSta);
          } else {
            Double_t effdata = dataStaEff_neg.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errdata = TMath::Max(dataStaEff_neg.getErrLow(fabs(mu1->eta), mu1->pt), dataStaEff_neg.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t effmc   = zmmStaEff_neg.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errmc   = TMath::Max(zmmStaEff_neg.getErrLow(fabs(mu1->eta), mu1->pt), zmmStaEff_neg.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hStaErr_neg->Fill(fabs(mu1->eta), mu1->pt, errSta);
            if(isBB) hStaErrB_neg->Fill(fabs(mu1->eta), mu1->pt, errSta);
            else if(isEE) hStaErrE_neg->Fill(fabs(mu1->eta), mu1->pt, errSta);
          }

          if(mu2->q>0) {
            Double_t effdata = dataStaEff_pos.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errdata = TMath::Max(dataStaEff_pos.getErrLow(fabs(mu2->eta), mu2->pt), dataStaEff_pos.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t effmc   = zmmStaEff_pos.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errmc   = TMath::Max(zmmStaEff_pos.getErrLow(fabs(mu2->eta), mu2->pt), zmmStaEff_pos.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t errSta = ((effdata/effmc))*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hStaErr_pos->Fill(fabs(mu2->eta), mu2->pt, errSta);
            if(isBB) hStaErrB_pos->Fill(fabs(mu2->eta), mu2->pt, errSta);
            else if(isEE) hStaErrE_pos->Fill(fabs(mu2->eta), mu2->pt, errSta);
          } else {
            Double_t effdata = dataStaEff_neg.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errdata = TMath::Max(dataStaEff_neg.getErrLow(fabs(mu2->eta), mu2->pt), dataStaEff_neg.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t effmc   = zmmStaEff_neg.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errmc   = TMath::Max(zmmStaEff_neg.getErrLow(fabs(mu2->eta), mu2->pt), zmmStaEff_neg.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hStaErr_neg->Fill(fabs(mu2->eta), mu2->pt, errSta);
            if(isBB) hStaErrB_neg->Fill(fabs(mu2->eta), mu2->pt, errSta);
            else if(isEE) hStaErrE_neg->Fill(fabs(mu2->eta), mu2->pt, errSta);
	  }

	  // SELECTION
          if(mu1->q>0) {
            Double_t effdata = dataSelEff_pos.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errdata = TMath::Max(dataSelEff_pos.getErrLow(fabs(mu1->eta), mu1->pt), dataSelEff_pos.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t effmc   = zmmSelEff_pos.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errmc   = TMath::Max(zmmSelEff_pos.getErrLow(fabs(mu1->eta), mu1->pt), zmmSelEff_pos.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hSelErr_pos->Fill(fabs(mu1->eta), mu1->pt, errSel);
            if(isBB) hSelErrB_pos->Fill(fabs(mu1->eta), mu1->pt, errSel);
            else if(isEE) hSelErrE_pos->Fill(fabs(mu1->eta), mu1->pt, errSel);
          } else {
            Double_t effdata = dataSelEff_neg.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errdata = TMath::Max(dataSelEff_neg.getErrLow(fabs(mu1->eta), mu1->pt), dataSelEff_neg.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t effmc   = zmmSelEff_neg.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errmc   = TMath::Max(zmmSelEff_neg.getErrLow(fabs(mu1->eta), mu1->pt), zmmSelEff_neg.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hSelErr_neg->Fill(fabs(mu1->eta), mu1->pt, errSel);
            if(isBB) hSelErrB_neg->Fill(fabs(mu1->eta), mu1->pt, errSel);
            else if(isEE) hSelErrE_neg->Fill(fabs(mu1->eta), mu1->pt, errSel);
          }

          if(mu2->q>0) {
            Double_t effdata = dataSelEff_pos.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errdata = TMath::Max(dataSelEff_pos.getErrLow(fabs(mu2->eta), mu2->pt), dataSelEff_pos.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t effmc   = zmmSelEff_pos.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errmc   = TMath::Max(zmmSelEff_pos.getErrLow(fabs(mu2->eta), mu2->pt), zmmSelEff_pos.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hSelErr_pos->Fill(fabs(mu2->eta), mu2->pt, errSel);
            if(isBB) hSelErrB_pos->Fill(fabs(mu2->eta), mu2->pt, errSel);
            else if(isEE) hSelErrE_pos->Fill(fabs(mu2->eta), mu2->pt, errSel);
          } else {
            Double_t effdata = dataSelEff_neg.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errdata = TMath::Max(dataSelEff_neg.getErrLow(fabs(mu2->eta), mu2->pt), dataSelEff_neg.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t effmc   = zmmSelEff_neg.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errmc   = TMath::Max(zmmSelEff_neg.getErrLow(fabs(mu2->eta), mu2->pt), zmmSelEff_neg.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hSelErr_neg->Fill(fabs(mu2->eta), mu2->pt, errSel);
            if(isBB) hSelErrB_neg->Fill(fabs(mu2->eta), mu2->pt, errSel);
            else if(isEE) hSelErrE_neg->Fill(fabs(mu2->eta), mu2->pt, errSel);
	  }

	  //HLT
          if(mu1->q>0) {
            Double_t effdata = dataHLTEff_pos.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errdata = TMath::Max(dataHLTEff_pos.getErrLow(fabs(mu1->eta), mu1->pt), dataHLTEff_pos.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t effmc   = zmmHLTEff_pos.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errmc   = TMath::Max(zmmHLTEff_pos.getErrLow(fabs(mu1->eta), mu1->pt), zmmHLTEff_pos.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hHLTErr_pos->Fill(fabs(mu1->eta), mu1->pt, errHLT);
            if(isBB) hHLTErrB_pos->Fill(fabs(mu1->eta), mu1->pt, errHLT);
            else if(isEE) hHLTErrE_pos->Fill(fabs(mu1->eta), mu1->pt, errHLT);
          } else {
            Double_t effdata = dataHLTEff_neg.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errdata = TMath::Max(dataHLTEff_neg.getErrLow(fabs(mu1->eta), mu1->pt), dataHLTEff_neg.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t effmc   = zmmHLTEff_neg.getEff(fabs(mu1->eta), mu1->pt);
            Double_t errmc   = TMath::Max(zmmHLTEff_neg.getErrLow(fabs(mu1->eta), mu1->pt), zmmHLTEff_neg.getErrHigh(fabs(mu1->eta), mu1->pt));
            Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hHLTErr_neg->Fill(fabs(mu1->eta), mu1->pt, errHLT);
            if(isBB) hHLTErrB_neg->Fill(fabs(mu1->eta), mu1->pt, errHLT);
            else if(isEE) hHLTErrE_neg->Fill(fabs(mu1->eta), mu1->pt, errHLT);
          }

          if(mu2->q>0) {
            Double_t effdata = dataHLTEff_pos.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errdata = TMath::Max(dataHLTEff_pos.getErrLow(fabs(mu2->eta), mu2->pt), dataHLTEff_pos.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t effmc   = zmmHLTEff_pos.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errmc   = TMath::Max(zmmHLTEff_pos.getErrLow(fabs(mu2->eta), mu2->pt), zmmHLTEff_pos.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hHLTErr_pos->Fill(fabs(mu2->eta), mu2->pt, errHLT);
            if(isBB) hHLTErrB_pos->Fill(fabs(mu2->eta), mu2->pt, errHLT);
            else if(isEE) hHLTErrE_pos->Fill(fabs(mu2->eta), mu2->pt, errHLT);
          } else {
            Double_t effdata = dataHLTEff_neg.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errdata = TMath::Max(dataHLTEff_neg.getErrLow(fabs(mu2->eta), mu2->pt), dataHLTEff_neg.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t effmc   = zmmHLTEff_neg.getEff(fabs(mu2->eta), mu2->pt);
            Double_t errmc   = TMath::Max(zmmHLTEff_neg.getErrLow(fabs(mu2->eta), mu2->pt), zmmHLTEff_neg.getErrHigh(fabs(mu2->eta), mu2->pt));
            Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hHLTErr_neg->Fill(fabs(mu2->eta), mu2->pt, errHLT);
            if(isBB) hHLTErrB_neg->Fill(fabs(mu2->eta), mu2->pt, errHLT);
            else if(isEE) hHLTErrE_neg->Fill(fabs(mu2->eta), mu2->pt, errHLT);
          }
	  
	  nSelv[ifile]    +=weight;
	  nSelCorrv[ifile]+=weight*corr;
	  nSelCorrVarv[ifile]+=weight*weight*corr*corr;

	  if(isBB){
            nSelBBv[ifile]+=weight;
            nSelBBCorrv[ifile]+=weight*corr;
            nSelBBCorrVarv[ifile]+=weight*weight*corr*corr;
	  }else if(isEE){
            nSelEEv[ifile]+=weight;
            nSelEECorrv[ifile]+=weight*corr;
            nSelEECorrVarv[ifile]+=weight*weight*corr*corr;
	  }
        }
      }      
    }

    Double_t var=0, varBB=0, varEE=0;
    for(Int_t iy=0; iy<=hHLTErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hHLTErr_pos->GetNbinsX(); ix++) {
        Double_t err;
	err=hHLTErr_pos->GetBinContent(ix,iy); var+=err*err;
	err=hHLTErrB_pos->GetBinContent(ix,iy); varBB+=err*err;
        err=hHLTErrE_pos->GetBinContent(ix,iy); varEE+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hHLTErr_neg->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hHLTErr_neg->GetNbinsX(); ix++) {
        Double_t err;
        err=hHLTErr_neg->GetBinContent(ix,iy); var+=err*err;
        err=hHLTErrB_neg->GetBinContent(ix,iy); varBB+=err*err;
        err=hHLTErrE_neg->GetBinContent(ix,iy); varEE+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hSelErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hSelErr_pos->GetNbinsX(); ix++) {
        Double_t err;
        err=hSelErr_pos->GetBinContent(ix,iy); var+=err*err;
        err=hSelErrB_pos->GetBinContent(ix,iy); varBB+=err*err;
        err=hSelErrE_pos->GetBinContent(ix,iy); varEE+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hSelErr_neg->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hSelErr_neg->GetNbinsX(); ix++) {
        Double_t err;
        err=hSelErr_neg->GetBinContent(ix,iy); var+=err*err;
        err=hSelErrB_neg->GetBinContent(ix,iy); varBB+=err*err;
        err=hSelErrE_neg->GetBinContent(ix,iy); varEE+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hTrkErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hTrkErr_pos->GetNbinsX(); ix++) {
        Double_t err=hTrkErr_pos->GetBinContent(ix,iy);
        //var+=err*err;
	var+=0.0;
      }
    }
    for(Int_t iy=0; iy<=hTrkErr_neg->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hTrkErr_neg->GetNbinsX(); ix++) {
        Double_t err=hTrkErr_neg->GetBinContent(ix,iy);
	var+=0.0;
        //var+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hStaErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hStaErr_pos->GetNbinsX(); ix++) {
        Double_t err;
        err=hStaErr_pos->GetBinContent(ix,iy); var+=err*err;
        err=hStaErrB_pos->GetBinContent(ix,iy); varBB+=err*err;
        err=hStaErrE_pos->GetBinContent(ix,iy); varEE+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hStaErr_neg->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hStaErr_neg->GetNbinsX(); ix++) {
	Double_t err;
        err=hStaErr_neg->GetBinContent(ix,iy); var+=err*err;
        err=hStaErrB_neg->GetBinContent(ix,iy); varBB+=err*err;
        err=hStaErrE_neg->GetBinContent(ix,iy); varEE+=err*err;
      }
    }

    nSelCorrVarv[ifile]+=var;
    nSelBBCorrVarv[ifile]+=varBB;
    nSelEECorrVarv[ifile]+=varEE;


    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);     accErrv.push_back(accv[ifile]*sqrt((1.+accv[ifile])/nEvtsv[ifile]));
    accCorrv.push_back(nSelCorrv[ifile]/nEvtsv[ifile]); accErrCorrv.push_back(accCorrv[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrv[ifile]*nSelCorrv[ifile]) + 1./nEvtsv[ifile]));

    accBBv.push_back(nSelBBv[ifile]/nEvtsv[ifile]); accErrBBv.push_back(sqrt(accBBv[ifile]*(1.+accBBv[ifile])/nEvtsv[ifile]));
    accBBCorrv.push_back(nSelBBCorrv[ifile]/nEvtsv[ifile]); accErrBBCorrv.push_back(accBBCorrv[ifile]*sqrt(nSelBBCorrVarv[ifile]/nSelBBCorrv[ifile]/nSelBBCorrv[ifile] + 1./nEvtsv[ifile]));

    accEEv.push_back(nSelEEv[ifile]/nEvtsv[ifile]); accErrEEv.push_back(sqrt(accEEv[ifile]*(1.+accEEv[ifile])/nEvtsv[ifile]));
    accEECorrv.push_back(nSelEECorrv[ifile]/nEvtsv[ifile]); accErrEECorrv.push_back(accEECorrv[ifile]*sqrt(nSelEECorrVarv[ifile]/nSelEECorrv[ifile]/nSelEECorrv[ifile] + 1./nEvtsv[ifile]));
    
    delete infile;
    infile=0, eventTree=0;  
  }  
  delete info;
  delete gen;
  delete muonArr;

    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================    
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " Z -> mu mu" << endl;
  cout << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "   ================================================" << endl;
    cout << "    Label: " << labelv[ifile] << endl;
    cout << "     File: " << fnamev[ifile] << endl;
    cout << endl;
    cout << "    *** Acceptance ***" << endl;
    cout << "       AllNominal: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    cout << "   AllSFCorrected: " << accCorrv[ifile] << " +/- " << accErrCorrv[ifile] << endl;
    cout << "    BarrelNominal: " << setw(12) << nSelBBv[ifile]   << " / " << nEvtsv[ifile] << " = " << accBBv[ifile]   << " +/- " << accErrBBv[ifile] << endl;
    cout << "BarrelSFCorrected: " << accBBCorrv[ifile] << " +/- " << accErrBBCorrv[ifile] << endl;
    cout << "    EndcapNominal: " << setw(12) << nSelEEv[ifile]   << " / " << nEvtsv[ifile] << " = " << accEEv[ifile]   << " +/- " << accErrEEv[ifile] << endl;
    cout << "EndcapSFCorrected: " << accEECorrv[ifile] << " +/- " << accErrEECorrv[ifile] << endl;
    cout << endl;
  }
  
  char txtfname[300];
  sprintf(txtfname,"%s/binned.txt",outputDir.Data());
  ofstream txtfile;
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << " Z -> mu mu" << endl;
  txtfile << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  txtfile << "  pT > " << PT_CUT << endl;
  txtfile << "  |eta| < " << ETA_CUT << endl;
  txtfile << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile << "   ================================================" << endl;
    txtfile << "    Label: " << labelv[ifile] << endl;
    txtfile << "     File: " << fnamev[ifile] << endl;
    txtfile << endl;
    txtfile << "    *** Acceptance ***" << endl;
    txtfile << "       AllNominal: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    txtfile << "   AllSFCorrected: " << accCorrv[ifile] << " +/- " << accErrCorrv[ifile] << endl;
    txtfile << "    BarrelNominal: " << setw(12) << nSelBBv[ifile]   << " / " << nEvtsv[ifile] << " = " << accBBv[ifile]   << " +/- " << accErrBBv[ifile] << endl;
    txtfile << "BarrelSFCorrected: " << accBBCorrv[ifile] << " +/- " << accErrBBCorrv[ifile] << endl;
    txtfile << "    EndcapNominal: " << setw(12) << nSelEEv[ifile]   << " / " << nEvtsv[ifile] << " = " << accEEv[ifile]   << " +/- " << accErrEEv[ifile] << endl;
    txtfile << "EndcapSFCorrected: " << accEECorrv[ifile] << " +/- " << accErrEECorrv[ifile] << endl;
    txtfile << endl;
  }
  txtfile.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccSelZmmBinned"); 
}
