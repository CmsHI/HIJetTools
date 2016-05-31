#include "TCanvas.h"
#include "TFile.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TString.h"
#include "TList.h"

#include "plotUtils.C"

const Double_t eps = 0.000001;

void plotRhoEtaSlices(TString str = "RandomConesPF.root", TString tag = "HYDJETMB") {

  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);


  TFile *f = new TFile(str.Data());

  TH3F *fh3RhoCentEtaBin   = dynamic_cast<TH3F*>(f->Get("fh3RhoCentEtaBin"));
  TH3F *fh3RhoMCentEtaBin  = dynamic_cast<TH3F*>(f->Get("fh3RhoMCentEtaBin"));
  
  // const int netabins = 7;
  // double etaMin[netabins] = {-5.,-3.,-2.1,-1.3,1.3,2.1,3.};
  // double etaMax[netabins] = {-3.,-2.1,-1.3,1.3,2.1,3.,5.};
  const int netabins = 9;
  double etaMin[netabins] = {-5.,-3.,-2.,-1.5,-1.,1.,1.5,2.1,3.};
  double etaMax[netabins] = {-3.,-2.,-1.5,-1.,1.,1.5,2.1,3.,5.};

  TH2D *h2RhoCent[netabins];
  TH2D *h2RhoMCent[netabins];
  TH2D *h2RhoJetCent[netabins];
  for(int i = 1; i<=netabins; ++i) {
    fh3RhoCentEtaBin->GetZaxis()->SetRange(i,i);
    h2RhoCent[i-1] = dynamic_cast<TH2D*>(fh3RhoCentEtaBin->Project3D("yx"));
    h2RhoCent[i-1]->SetName(Form("hRhoCentEta%d",i));
    h2RhoCent[i-1]->SetTitle(Form("hRhoCentEta%d",i));

    fh3RhoMCentEtaBin->GetZaxis()->SetRange(i,i);
    h2RhoMCent[i-1] = dynamic_cast<TH2D*>(fh3RhoMCentEtaBin->Project3D("yx"));
    h2RhoMCent[i-1]->SetName(Form("hRhoMCentEta%d",i));
    h2RhoMCent[i-1]->SetTitle(Form("hRhoMCentEta%d",i));
  }

  TCanvas *c1 = new TCanvas("c1","c1: rho",900,840);
  c1->Divide(3,3);
  for(int i = 1; i<=netabins; ++i) {
    c1->cd(i);
    TH1F *fr1 = DrawFrame(0.,100.,0.,500.,"centrality (%)","#rho GeV/c");
    gPad->SetLogz();
    h2RhoCent[i-1]->Draw("colz same");

    DrawLatex(0.5,0.77,Form("%.1f<#eta<%.1f",etaMin[i-1],etaMax[i-1]));
  }

  c1->SaveAs("RhoCentEtaBins.png");
  c1->SaveAs("RhoCentEtaBins.pdf");
  
  TCanvas *c2 = new TCanvas("c2","c2: rho_m",900,840);
  c2->Divide(3,3);
  for(int i = 1; i<=netabins; ++i) {
    c2->cd(i);
    TH1F *fr2 = DrawFrame(0.,100.,0.,2.,"centrality (%)","#rho_{m} GeV/c");
    gPad->SetLogz();
    h2RhoMCent[i-1]->Draw("colz same");

    DrawLatex(0.5,0.77,Form("%.1f<#eta<%.1f",etaMin[i-1],etaMax[i-1]));
  }

  c2->SaveAs("RhoMCentEtaBins.png");
  c2->SaveAs("RhoMCentEtaBins.pdf");
  
}
