#include <iostream>
#include <TLegend.h>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include "TH1D.h"
#include "TH2D.h"
#include <THStack.h>
#include "TProfile.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFractionFitter.h"
#include <string>
#include <vector>
#include <math.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPave.h>
#include <TPaveStats.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TString.h>
#include "TGraphErrors.h"
#include "PileupLumiArray.h"

using namespace std;

void DrawFiducialZRate_Time(
                        const TString  lepType,  //"Mu" or "Ele"
                        const TString  dirLabel, //for example, "Fill5198"
			const Int_t    lepRegion,//0 for all, 1 for central, 2 for forward, 3 for 3 curves
			const Int_t    varType,  //0 for Z rates, 1 for Z counts, 2 for fiducial Z xsec, 3 for total Z xsec
                        const Double_t lumLabel, //lumonosity number on plot
                        const Double_t ymin=0.,  //y axis range for plot
                        const Double_t ymax=0.   //y axis range for plot
){
        TCanvas *c = new TCanvas("c","", 800, 600);
	TH1D *h = new TH1D("h","", ND, Time_bin_array);

        Float_t y[ND] = {};
        Float_t ey[ND] = {};

	float bkgpercent = (lepType=="Mu"? 0.0105: 0.0064);

	float genacc = 0.;
	if(lepRegion == 0) genacc = 0.350312;
	if(lepRegion == 1) genacc = 0.0864849;

	for(int ibin = 0; ibin < ND; ibin++){
		if(lumi[ibin] == 0.) continue;

		char yieldfile[200];
	        char accfile[200];
	        char errfile[200];
	        sprintf(yieldfile,"Yield%s%s/%syield_%d.txt", lepType.Data(), dirLabel.Data(), ((lepRegion == 0) ? "" : "BB"), ibin+D0);
	        sprintf(accfile,"Efficiency%s%s/%sacc_%d.txt",lepType.Data(), dirLabel.Data(), ((lepRegion == 0) ? "" : "BB"), ibin+D0);
	        sprintf(errfile,"Efficiency%s%s/%serr_%d.txt",lepType.Data(), dirLabel.Data(), ((lepRegion == 0) ? "" : "BB"), ibin+D0);
	        std::ifstream ydfile(yieldfile);
	        std::ifstream acfile(accfile);
		if(!acfile.good()) continue;
	        std::ifstream erfile(errfile);
	        float yield, acc, err = 0;
	        while(ydfile >> yield){}
	        while(acfile >> acc){}
	        while(erfile >> err){}

	        if(lepType=="Mu"){
			if(varType == 0) y[ibin] = yield * (1-bkgpercent)/(acc*Time_Corr[ibin]/genacc); //Z rates
			if(varType == 1) y[ibin] = yield * (1-bkgpercent)/(acc/genacc);		        //Z counts 
			if(varType == 2) y[ibin] = yield * (1-bkgpercent)/(acc*lumi[ibin]/genacc);      //Z fiducial xsecs
			if(varType == 3) y[ibin] = yield * (1-bkgpercent)/(acc*lumi[ibin]);		//Z total xsecs
		}

	        ey[ibin] = y[ibin] * sqrt((sqrt(yield)/yield)*(sqrt(yield)/yield) + err*err);// + 0.013*0.013 + 0.014*0.014 +0.027*0.027);
                h->SetBinContent(ibin+1, y[ibin]);
                h->SetBinError(ibin+1, ey[ibin]);

//		cout<<Fill_bound[ibin]<<", "<<Time_bound[ibin]<<", "<<Time_bound[ibin+1]<<", "<<y[ibin]<<", "<<lumi[ibin]/Time_Corr[ibin]<<", "<<lumi[ibin]*Time_Real[ibin]/Time_Corr[ibin]<<", "<<yield * (1-bkgpercent)*Time_Real[ibin]/(Time_Corr[ibin]*acc/genacc)<<endl;
		cout<<Fill_bound[ibin]<<", "<<Time_bound[ibin]<<", "<<Time_bound[ibin+1]<<", "<<y[ibin]<<", "<<lumi[ibin]/Time_Corr[ibin]<<", "<<lumi[ibin]<<", "<<yield * (1-bkgpercent)/(acc/genacc)<<endl;
	}		

        char lumitext[100];
        sprintf(lumitext,"%.1f fb^{-1}  at  #sqrt{s} = 13 TeV",lumLabel);
        TPaveText *lumitb = new TPaveText(0.52,0.92,0.84,0.99,"NDC");
        lumitb->SetBorderSize(0);
        lumitb->SetFillColor(0);
        lumitb->AddText(lumitext);

	TPaveText *cmstb = new TPaveText(0.15,0.83,0.44,0.89,"NDC");
	cmstb->SetBorderSize(0);
	cmstb->SetFillColor(0);
	cmstb->AddText("CMS Preliminary");

	h->SetTitle("");
	h->SetStats(0);
	//for(int ibin = 0; ibin < ND; ibin++) h->GetXaxis()->SetBinLabel(ibin+1, Time_bound[ibin]);
	h->GetXaxis()->SetLabelSize(0.02);
        h->GetXaxis()->SetTitle("real time [min]");
	if(varType == 0) {h->GetYaxis()->SetTitle("Z rates [s-1]");			h->GetYaxis()->SetRangeUser(0,10);}
	if(varType == 1) {h->GetYaxis()->SetTitle("Z counts");				h->GetYaxis()->SetRangeUser(5500,7500);}
	if(varType == 2) {h->GetYaxis()->SetTitle("Z fiducial cross section [pb]");	h->GetYaxis()->SetRangeUser(550,750);}
	if(varType == 3) {h->GetYaxis()->SetTitle("Z total cross section [pb]");	h->GetYaxis()->SetRangeUser(1600,2100);}
        h->GetXaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetTitleSize(0.05);
        h->SetMarkerStyle(20);
        h->SetLineColor(1);
        h->SetLineWidth(2);
        h->Draw();

        lumitb->Draw("same");
        cmstb->Draw("same");

	char plotname[200];
	if(varType == 0) sprintf(plotname,"%s_%s_%s.pdf", dirLabel.Data(), "ZRates", ((lepRegion == 0) ? "All" : "Central"));
	if(varType == 1) sprintf(plotname,"%s_%s_%s.pdf", dirLabel.Data(), "ZCounts", ((lepRegion == 0) ? "All" : "Central")); 
	if(varType == 2) sprintf(plotname,"%s_%s_%s.pdf", dirLabel.Data(), "ZFiducialXSec", ((lepRegion == 0) ? "All" : "Central"));
	if(varType == 3) sprintf(plotname,"%s_%s_%s.pdf", dirLabel.Data(), "ZTotalXSec", ((lepRegion == 0) ? "All" : "Central"));
	
	c->SaveAs(plotname);

	delete h;
	delete c;
}
