#include "plotUtils.C"
#include "fitUtils.C"

void plotRandomCones(TString str = "output.root", TString tag = "", double etamin = 1.5, double etamax = 1.5, int run = -1) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TString tagSave = tag;
  tagSave.ReplaceAll(" ","");
  
  TFile *f = new TFile(str.Data());

  TH3D *h3CentPtEta = dynamic_cast<TH3D*>(f->Get("h2CentPtRCEta"));
  TH3D *h3CentPtEtaVS = dynamic_cast<TH3D*>(f->Get("h2CentPtRCEtaVS"));

  TH3D *h3EHFPtEta = dynamic_cast<TH3D*>(f->Get("h2EHFPtRCEta"));
  TH3D *h3EHFPtEtaVS = dynamic_cast<TH3D*>(f->Get("h2EHFPtRCEtaVS"));
  
  int minbinEta, maxbinEta;
  minbinEta = h3CentPtEta->GetZaxis()->FindBin(etamin+0.000001);
  maxbinEta = h3CentPtEta->GetZaxis()->FindBin(etamax-0.000001);
  h3CentPtEta->GetZaxis()->SetRange(minbinEta,maxbinEta);
  TH2D *h2CentPt = dynamic_cast<TH2D*>(h3CentPtEta->Project3D("yx"));
  minbinEta = h3CentPtEtaVS->GetZaxis()->FindBin(etamin+0.000001);
  maxbinEta = h3CentPtEtaVS->GetZaxis()->FindBin(etamax-0.000001);
  h3CentPtEtaVS->GetZaxis()->SetRange(minbinEta,maxbinEta);
  TH2D *h2CentPtVS = dynamic_cast<TH2D*>(h3CentPtEtaVS->Project3D("yx"));

  minbinEta = h3EHFPtEta->GetZaxis()->FindBin(etamin+0.000001);
  maxbinEta = h3EHFPtEta->GetZaxis()->FindBin(etamax-0.000001);
  h3EHFPtEta->GetZaxis()->SetRange(minbinEta,maxbinEta);
  TH2D *h2EHFPt = dynamic_cast<TH2D*>(h3EHFPtEta->Project3D("yx"));
  minbinEta = h3EHFPtEtaVS->GetZaxis()->FindBin(etamin+0.000001);
  maxbinEta = h3EHFPtEtaVS->GetZaxis()->FindBin(etamax-0.000001);
  h3EHFPtEtaVS->GetZaxis()->SetRange(minbinEta,maxbinEta);
  TH2D *h2EHFPtVS = dynamic_cast<TH2D*>(h3EHFPtEtaVS->Project3D("yx"));

  TCanvas *c1 = new TCanvas("c1","c1",500,400);
  gPad->SetLogz();
  TH1F *f1 = DrawFrame(0.,100.,-50.,100.,"centrality (%)","p_{T}^{RC}",0.1);
  h2CentPt->Draw("colz same");

  TCanvas *c11 = new TCanvas("c11","c11",500,400);
  gPad->SetLogz();
  TH1F *f11 = DrawFrame(0.,6000.,-50.,200.,"sum E_{HF}","p_{T}^{RC}",0.1);
  h2EHFPt->Rebin2D(5,5);
  h2EHFPt->Draw("colz same");

  TCanvas *c2 = new TCanvas("c2","c2",500,400);
  gPad->SetLogz();
  TH1F *f2 = DrawFrame(0.,100.,-80.,100.,"centrality (%)","p_{T}^{RC}",0.1);
  h2CentPtVS->Draw("colz same");

  DrawLatex(0.25,0.87,tag.Data());
  if(run>0)   DrawLatex(0.25,0.81,Form("run: %d",run));

  TCanvas *c22 = new TCanvas("c22","c22",500,400);
  gPad->SetLogz();
  TH1F *f22 = DrawFrame(0.,6000.,-80.,100.,"sum E_{HF}","p_{T}^{RC}",0.1);
  h2EHFPtVS->Rebin2D(5,5);
  h2EHFPtVS->Draw("colz same");
  
  c1->SaveAs(Form("RCPtVsCentNoSub%s.png",tagSave.Data()));
  c2->SaveAs(Form("RCPtVsCentSub%s.png",tagSave.Data()));
  c1->SaveAs(Form("RCPtVsCentNoSub%s.pdf",tagSave.Data()));
  c2->SaveAs(Form("RCPtVsCentSub%s.pdf",tagSave.Data()));
  c1->SaveAs(Form("RCPtVsCentNoSub%s.eps",tagSave.Data()));
  c2->SaveAs(Form("RCPtVsCentSub%s.eps",tagSave.Data()));

  c11->SaveAs(Form("RCPtVsEHFNoSub%s.png",tagSave.Data()));
  c22->SaveAs(Form("RCPtVsEHFSub%s.png",tagSave.Data()));
  c11->SaveAs(Form("RCPtVsEHFNoSub%s.pdf",tagSave.Data()));
  c22->SaveAs(Form("RCPtVsEHFSub%s.pdf",tagSave.Data()));
  c11->SaveAs(Form("RCPtVsEHFNoSub%s.eps",tagSave.Data()));
  c22->SaveAs(Form("RCPtVsEHFSub%s.eps",tagSave.Data()));

  const int nCent = 5;
  double centMin[nCent] = {0.,10.,30.,50.,80.};
  double centMax[nCent] = {10.,30.,50.,80.,100.};

  TH1D *h1RC[nCent];
  int binMin, binMax;
  for(int i = 0; i<nCent; ++i) {
    binMin = h2CentPtVS->GetXaxis()->FindBin(centMin[i]+0.000001);
    binMax = h2CentPtVS->GetXaxis()->FindBin(centMax[i]-0.000001);
    h1RC[i] = dynamic_cast<TH1D*>(h2CentPtVS->ProjectionY(Form("h1RC_%d",i),binMin,binMax));
  }

  TH1D *h1RCNoSub[nCent];
   for(int i = 0; i<nCent; ++i) {
    binMin = h2CentPt->GetXaxis()->FindBin(centMin[i]+0.000001);
    binMax = h2CentPt->GetXaxis()->FindBin(centMax[i]-0.000001);
    h1RCNoSub[i] = dynamic_cast<TH1D*>(h2CentPt->ProjectionY(Form("h1RCNoSub_%d",i),binMin,binMax));
  }

  TF1 *fitLHS[nCent];
  double meanLHS[nCent] = {0.};
  double meanLHSerr[nCent] = {0.};
  
  TCanvas *c3 = new TCanvas("c3","c3",500,400);
  TH1F *fr3 = DrawFrame(-60.,150.,3e-5,2.,"p_{T}^{RC} - #rho A (GeV/c)","probability/bin");
  gPad->SetLogy();
  TLegend *leg3 = CreateLegend(0.62,0.9,0.43,0.93,"",0.045);
  for(int i = 0; i<nCent; ++i) {
    if(i<2) h1RC[i]->Rebin(5);
    h1RC[i]->SetLineColor(GetColor(i));
    h1RC[i]->SetMarkerColor(GetColor(i));
    h1RC[i]->SetMarkerStyle(GetMarker(i));
    h1RC[i]->Scale(1./h1RC[i]->Integral());
    h1RC[i]->Draw("same");

    //if(i==0) {
    //fitLHS[i] = fitLHSGaus(h1RC[i]);
    TF1 *f1tmp = fitLHSGaus(h1RC[i],2.,0.5);
    fitLHS[i] = dynamic_cast<TF1*>(f1tmp->Clone(Form("fitLHS%d",i)));
      fitLHS[i]->SetName(Form("fitLHS%d",i));
      fitLHS[i]->SetTitle(Form("fitLHS%d",i));
      fitLHS[i]->SetLineColor(GetColor(i));
      if(i==0) fitLHS[i]->Draw("same");
      
      meanLHS[i] = fitLHS[i]->GetParameter(1);
      meanLHSerr[i] = fitLHS[i]->GetParError(1);
      //}
    TString strleg = Form("%.0f-%0.f",centMin[i],centMax[i]);
    strleg+="%";
    leg3->AddEntry(h1RC[i],Form("#splitline{%s}{#mu^{LHS}=%.1f #sigma=%.1f}",strleg.Data(),meanLHS[i],h1RC[i]->GetRMS()),"pl");
  }
  leg3->Draw();

  for(int i = 0; i<nCent; ++i)
    fitLHS[i]->Draw("same");

  
  DrawLatex(0.25,0.88,tag.Data(),0.05);
  if(run>0)   DrawLatex(0.25,0.88,Form("run: %d",run),0.05);

  c3->SaveAs(Form("RCPt1DCentSub%s.png",tagSave.Data()));
  c3->SaveAs(Form("RCPt1DCentSub%s.pdf",tagSave.Data()));
  c3->SaveAs(Form("RCPt1DCentSub%s.eps",tagSave.Data()));

  // return;

  //no sub
  TF1 *fitLHSNoSub[nCent];
  double meanLHSNoSub[nCent] = {0.};
  
  TCanvas *c4 = new TCanvas("c4","c4",500,400);
  TH1F *fr4 = DrawFrame(0.,300.,3e-5,2.,"p_{T}^{RC}","probability/bin");
  gPad->SetLogy();
  TLegend *leg4 = CreateLegend(0.62,0.9,0.4,0.92,"",0.05);
  for(int i = 0; i<nCent; ++i) {
    h1RCNoSub[i]->Rebin(5);
    h1RCNoSub[i]->SetLineColor(GetColor(i));
    h1RCNoSub[i]->SetMarkerColor(GetColor(i));
    h1RCNoSub[i]->SetMarkerStyle(GetMarker(i));
    h1RCNoSub[i]->Scale(1./h1RCNoSub[i]->Integral());
    h1RCNoSub[i]->Draw("same");

    //if(i==0) {
      fitLHSNoSub[i] = fitLHSGaus(h1RCNoSub[i]);
      fitLHSNoSub[i]->SetName(Form("fitLHSNoSub%d",i));
      fitLHSNoSub[i]->SetTitle(Form("fitLHSNoSub%d",i));
      fitLHSNoSub[i]->SetLineColor(GetColor(i));
      if(i==0) fitLHSNoSub[i]->Draw("same");
      meanLHSNoSub[i] = fitLHSNoSub[i]->GetParameter(1);
      //}
    TString strleg = Form("%.0f-%0.f",centMin[i],centMax[i]);
    strleg+="%";
    leg4->AddEntry(h1RCNoSub[i],Form("#splitline{%s}{#mu^{LHS}=%.1f #sigma=%.1f}",strleg.Data(),meanLHSNoSub[i],h1RCNoSub[i]->GetRMS()),"pl");
  }
  leg4->Draw();

  DrawLatex(0.25,0.82,tag.Data(),0.05);
  if(run>0)   DrawLatex(0.25,0.88,Form("run: %d",run),0.05);

  c4->SaveAs(Form("RCPt1DCentNoSub%s.png",tagSave.Data()));
  c4->SaveAs(Form("RCPt1DCentNoSub%s.pdf",tagSave.Data()));
  c4->SaveAs(Form("RCPt1DCentNoSub%s.eps",tagSave.Data()));

  for(int i = 0; i<nCent; ++i) {
    Printf("%.0f-%.0f : mu=%.1f +/- %.1f",centMin[i],centMax[i],meanLHS[i],meanLHSerr[i]);
  }
  
}
