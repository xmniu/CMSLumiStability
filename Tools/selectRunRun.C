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
#include "../Utils/DatasetBounds.hh"// runnum, lumisec bounds for dividing datasets
//#include "DatasetBounds.hh"
#endif

void selectRunRun( const TString inputFile,
		   const TString outputFile,
		   const Int_t ibound=-1)
{

TFile *outFile = new TFile(outputFile,"RECREATE");
TTree *outTree = new TTree("Events","Events");
UInt_t runNum, lumiSec, evtNum;
outTree->Branch("runNum", &runNum, "runNum/i");
outTree->Branch("lumiSec",&lumiSec,"lumiSec/i");
outTree->Branch("evtNum", &evtNum, "evtNum/i");

TFile *infile = new TFile(inputFile);         assert(infile);
TTree *intree = (TTree*)infile->Get("Events"); assert(intree);

intree->SetBranchAddress("runNum",     &runNum);      // event run number
intree->SetBranchAddress("lumiSec",    &lumiSec);     // event lumi section
intree->SetBranchAddress("evtNum",     &evtNum);      // event number

for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
	intree->GetEntry(ientry);
	Bool_t doSelection = 0;
	//int ibound = -1;

	//for(int i = 0; i < 55; i++){
	//ibound = ibound_array[i];
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
	//}

        if(doSelection != 1) continue;

	runNum  = runNum;
	lumiSec = lumiSec;
	evtNum  = evtNum;

	outTree->Fill();
}	

outFile->Write();
outFile->Close();
}
