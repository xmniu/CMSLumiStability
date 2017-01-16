#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TString.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#endif

// Main macro function
//--------------------------------------------------------------------------------------------------
void check_f_rw_80X(
		    const TString outputDir,
		    const Int_t ibound
){

        TFile *f_mc_pu = TFile::Open("pileup_dis_80X_MC.root", "read");
        TH1D  *h_mc_pu = (TH1D*) f_mc_pu->Get("pileup");
        h_mc_pu->Scale(1./h_mc_pu->Integral());
	h_mc_pu->SetStats(0);

        char data_pu_file[100];
        sprintf(data_pu_file,"%s/PileupHistogram/pileup_%d.root", outputDir.Data(), ibound);
        TFile *f_data_pu = TFile::Open(data_pu_file, "read");
        TH1D  *h_data_pu = (TH1D*) f_data_pu->Get("pileup");
        h_data_pu->Scale(1./h_data_pu->Integral());
	float ymax = (h_data_pu->GetMaximum() > h_mc_pu->GetMaximum()) ? h_data_pu->GetMaximum() : h_mc_pu->GetMaximum(); 

        char pu_rw_file[100];
        sprintf(pu_rw_file,"%s/ReweightHistogram/pileup_rw_80X_%d.root", outputDir.Data(), ibound);
	TFile *f_pu_rw = TFile::Open(pu_rw_file, "read");
	TH1D  *h_pu_rw = (TH1D*) f_pu_rw->Get("h_rw");
	h_pu_rw->SetStats(0);

	TH1D *h_mc_pu_aft = new TH1D("h_mc_pu_aft","",60,0,60);
	for(int ibin = 0; ibin < 60; ibin++){
		h_mc_pu_aft->SetBinContent(ibin, h_mc_pu->GetBinContent(ibin)*h_pu_rw->GetBinContent(ibin));
	}

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  h_data_pu->GetXaxis()->SetTitle("nPUmean");
  h_data_pu->GetXaxis()->SetTitleSize(0.05);
  h_data_pu->GetYaxis()->SetRangeUser(0.,ymax+0.01);
  h_data_pu->SetMarkerColor(kBlack);
  h_data_pu->SetMarkerStyle(20);
//  h_data_pu->SetMarkerSize();
  h_data_pu->Draw("P");
  h_mc_pu->SetLineColor(kBlue);
  h_mc_pu->SetLineWidth(2);
  h_mc_pu->Draw("same");
  h_mc_pu_aft->SetLineColor(kRed);
  h_mc_pu_aft->SetLineWidth(2);
  h_mc_pu_aft->Draw("same");

  TLegend *l = new TLegend(0.6,0.5,0.9,0.7);
  l->AddEntry(h_data_pu,"sub-dataset of 2016","p");
  l->AddEntry(h_mc_pu,"MC before reweight","l");
  l->AddEntry(h_mc_pu_aft,"MC after reweight","l");
  l->Draw();

  char check_plot[100];
  sprintf(check_plot,"%s/CheckReweight/rw_check_%d.png", outputDir.Data(), ibound);
  c1->SaveAs(check_plot);
}
