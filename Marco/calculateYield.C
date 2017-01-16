#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include "TLorentzVector.h"         // 4-vector class
#include "TH1D.h"
#include "TRandom.h"
#include "../Utils/DatasetBounds.hh"      // runnum, lumisec bounds for dividing datasets


#endif

void calculateYield( const TString inputFile,
		     const TString outputDir,
                     const TString outputName,
		     const Int_t ibound=-1)
{

const Double_t MASS_LOW   = 66;
const Double_t MASS_HIGH  = 116;

//const Double_t PT_CUT     = 30;
//const Double_t ETA_CUT    = 2.5;
//const Double_t ETA_BARREL = 1.4442;
//const Double_t ETA_ENDCAP = 1.566;

const Double_t PT_CUT     = 27;
const Double_t ETA_CUT    = 2.4;
const Double_t ETA_BARREL = 0.9;
const Double_t ETA_ENDCAP = 0.9;

UInt_t runNum, lumiSec, evtNum;
UInt_t  category;
Float_t pfCombIso1, pfCombIso2;
TLorentzVector *dilep=0, *lep1=0, *lep2=0;

TFile *infile = new TFile(inputFile);         assert(infile);
TTree *intree = (TTree*)infile->Get("Events"); assert(intree);

intree->SetBranchAddress("runNum",     &runNum);      // event run number
intree->SetBranchAddress("lumiSec",    &lumiSec);     // event lumi section
intree->SetBranchAddress("evtNum",     &evtNum);      // event number
intree->SetBranchAddress("category",   &category);    // dilepton category
intree->SetBranchAddress("dilep",      &dilep);       // dilepton 4-vector
intree->SetBranchAddress("lep1",       &lep1);        // tag lepton 4-vector
intree->SetBranchAddress("lep2",       &lep2);        // probe lepton 4-vector
intree->SetBranchAddress("pfCombIso1", &pfCombIso1);  // PF combined isolation of tag lepton
intree->SetBranchAddress("pfCombIso2", &pfCombIso2);  // PF combined isolation of probe lepton   

Double_t Yield = 0;
Double_t YieldBB = 0;
Double_t YieldEE = 0;
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

	if(!(category==1 || category==2 || category==3)) continue;
//	if(pfCombIso1 > 0.15 * lep1->Pt()) cout<<"Not expected..."<<endl;
//	if(pfCombIso1 > 0.15 * lep1->Pt() || pfCombIso2 > 0.15 * lep2->Pt()) continue;
	if(!(dilep->M()<MASS_HIGH && dilep->M()>MASS_LOW)) continue;
	if(!(lep1->Pt()>PT_CUT && lep2->Pt()>PT_CUT && fabs(lep1->Eta())<ETA_CUT && fabs(lep2->Eta())<ETA_CUT)) continue;
	Yield += 1.;
	if(fabs(lep1->Eta())<ETA_BARREL && fabs(lep2->Eta())<ETA_BARREL) YieldBB += 1.;
	else if(fabs(lep1->Eta())>ETA_ENDCAP && fabs(lep2->Eta())>ETA_ENDCAP) YieldEE += 1.;
}	
        char outputfile[100];
        sprintf(outputfile,"%s/%s.txt",outputDir.Data(), outputName.Data());
        ofstream yieldfile;
        yieldfile.open(outputfile);
        yieldfile<<Yield<<endl;
        yieldfile.close();

        char outputBBfile[100];
        sprintf(outputBBfile,"%s/BB%s.txt",outputDir.Data(), outputName.Data());
        ofstream yieldBBfile;
        yieldBBfile.open(outputBBfile);
        yieldBBfile<<YieldBB<<endl;
        yieldBBfile.close();

        char outputEEfile[100];
        sprintf(outputEEfile,"%s/EE%s.txt",outputDir.Data(), outputName.Data());
        ofstream yieldEEfile;
        yieldEEfile.open(outputEEfile);
        yieldEEfile<<YieldEE<<endl;
        yieldEEfile.close();
}
