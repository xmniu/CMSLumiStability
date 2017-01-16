#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include <TLegend.h>
#include "TCanvas.h"
#include "TStyle.h"
#include "TFractionFitter.h"
#include <string>
#include <vector>
#include <math.h>
#include <TPaveStats.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

void make_f_rw_80X(
		   const TString outputDir,
		   const Int_t ibound
){

	TFile *f_mc_pu = TFile::Open("pileup_dis_80X_MC.root", "read");
	TH1D  *h_mc_pu = (TH1D*) f_mc_pu->Get("pileup");
	h_mc_pu->Scale(1./h_mc_pu->Integral());

	char data_pu_file[100];
	sprintf(data_pu_file,"%s/PileupHistogram/pileup_%d.root", outputDir.Data(), ibound);
	TFile *f_data_pu = TFile::Open(data_pu_file, "read");
	TH1D  *h_data_pu = (TH1D*) f_data_pu->Get("pileup");
	h_data_pu->Scale(1./h_data_pu->Integral());


	TH1D *h_rw = new TH1D("h_rw","",60,0,60);

	for(int ibin = 0; ibin < 60; ibin++){
		if(h_mc_pu->GetBinContent(ibin+1) == 0.) h_rw->SetBinContent(ibin+1, 0.);
		else h_rw->SetBinContent(ibin+1, h_data_pu->GetBinContent(h_data_pu->FindBin(ibin+0.5))/h_mc_pu->GetBinContent(h_mc_pu->FindBin(ibin+0.5)));
	}

        char pu_rw_file[100];
        sprintf(pu_rw_file,"%s/ReweightHistogram/pileup_rw_80X_%d.root", outputDir.Data(), ibound);
	TFile *f = new TFile(pu_rw_file, "RECREATE");
	h_rw->Write();
	f->Write();
	f->Close();
}
