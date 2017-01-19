
//Plotting macro to simply sum all the stuff in the deriveDijetResponse macro
//K. Jung, Raghav Kunnawalkam Elayavalli, December 2015


#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TRandom.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"


using namespace std;

void fillRelHist (TH1D* inHist, TH1F* outHist, int j){
	
	if(inHist->GetEntries()>0) outHist->SetBinContent(j+1,(1+inHist->GetMean())/(1-inHist->GetMean()));
	if(inHist->GetEntries()>0) outHist->SetBinError(j+1,inHist->GetRMS()*(1./TMath::Sqrt(inHist->GetEntries()>0)));
	else outHist->SetBinContent(j+1,0);
}

void fillMPFHist (TH1D* inHist, TH1F* outHist, int j){
	
	if(inHist->GetEntries()>0) outHist->SetBinContent(j+1,inHist->GetMean());
	if(inHist->GetEntries()>0) outHist->SetBinError(j+1,inHist->GetRMS()*(1./TMath::Sqrt(inHist->GetEntries()>0)));
	else outHist->SetBinContent(j+1,0);
}

void fillRelRatioHist (TH1D* inHist, TH1D* inHistMC, TH1F* outHist, int j){
	
	float relratio = 0;
	float aerr = 0;
	float amcerr = 0;
	float amean = 0;
	float amcmean = 0;
	
	float data_fit_mean = 0;
	float data_fit_err = 0;
	float mc_fit_mean = 0;
	float mc_fit_err = 0;
	inHist->Fit("gaus", "Q L M", "");
	inHistMC->Fit("gaus", "Q L M", "");
	
	data_fit_mean = inHist->GetFunction("gaus")->GetParameter(1);
	data_fit_err = inHist->GetFunction("gaus")->GetParError(1);
	mc_fit_mean = inHistMC->GetFunction("gaus")->GetParameter(1);
	mc_fit_err = inHistMC->GetFunction("gaus")->GetParError(1);
/*
	if(inHist->GetEntries()>0 && inHistMC->GetEntries()>0 ){
	relratio =  (1+inHistMC->GetMean())*(1-inHist->GetMean())/((1-inHistMC->GetMean())*(1+inHist->GetMean()));
	aerr = inHist->GetRMS()*(1./TMath::Sqrt(inHist->Integral()));
	amcerr = inHistMC->GetRMS()*(1./TMath::Sqrt(inHistMC->Integral()));
	amean = (1+inHist->GetMean())/(1-inHist->GetMean());
	amcmean = (1+inHistMC->GetMean())/(1-inHistMC->GetMean());
	outHist->SetBinContent(j+1,relratio);
	outHist->SetBinError(j+1,relratio*sqrt(pow(aerr/amean,2.)+pow(amcerr/amcmean,2.)));
*/
	if(inHist->GetEntries()>0 && inHistMC->GetEntries()>0 ){
	relratio =  (1+mc_fit_mean)*(1-data_fit_mean)/((1-mc_fit_mean)*(1+data_fit_mean));
	//aerr = data_fit_err*(1./TMath::Sqrt(inHist->Integral()));
	//amcerr = mc_fit_err*(1./TMath::Sqrt(inHistMC->Integral()));
	aerr = data_fit_err;
	amcerr = mc_fit_err;
	
	amean = (1+data_fit_mean)/(1-data_fit_mean);
	amcmean = (1+mc_fit_mean)/(1-mc_fit_mean);
	outHist->SetBinContent(j+1,relratio);
	outHist->SetBinError(j+1,relratio*sqrt(pow(aerr/amean,2.)+pow(amcerr/amcmean,2.)));
	
	
	}
	else outHist->SetBinContent(j+1,0);
}

void fillMPFRatioHist (TH1D* inHist, TH1D* inHistMC, TH1F* outHist, int j){
	
	float mpfratio = 0;
	float berr = 0;
	float bmcerr = 0;
	float bmean = 0;
	float bmcmean = 0;
	
	float data_fit_mean = 0;
	float data_fit_err = 0;
	float mc_fit_mean = 0;
	float mc_fit_err = 0;
	inHist->Fit("gaus", "Q L M", "");
	inHistMC->Fit("gaus", "Q L M", "");
	
	data_fit_mean = inHist->GetFunction("gaus")->GetParameter(1);
	data_fit_err = inHist->GetFunction("gaus")->GetParError(1);
	mc_fit_mean = inHistMC->GetFunction("gaus")->GetParameter(1);
	mc_fit_err = inHistMC->GetFunction("gaus")->GetParError(1);
	
/*
	if(inHist->GetEntries()>0 && inHistMC->GetEntries()>0 ){
	mpfratio =  (inHistMC->GetMean())/((inHist->GetMean()));
	berr = inHist->GetRMS()*(1./TMath::Sqrt(inHist->Integral()));
	bmcerr = inHistMC->GetRMS()*(1./TMath::Sqrt(inHistMC->Integral()));
	bmean = inHist->GetMean();
	bmcmean = inHistMC->GetMean();
	outHist->SetBinContent(j+1,mpfratio);
	outHist->SetBinError(j+1,mpfratio*sqrt(pow(berr/bmean,2.)+pow(bmcerr/bmcmean,2.)));
*/
	if(inHist->GetEntries()>0 && inHistMC->GetEntries()>0 ){
	mpfratio =  (mc_fit_mean)/(data_fit_mean);
	//berr = data_fit_err*(1./TMath::Sqrt(inHist->Integral()));
	//bmcerr = mc_fit_err*(1./TMath::Sqrt(inHistMC->Integral()));
	berr = data_fit_err;
	bmcerr = mc_fit_err;

	bmean = data_fit_mean;
	bmcmean = mc_fit_mean;
	outHist->SetBinContent(j+1,mpfratio);
	outHist->SetBinError(j+1,mpfratio*sqrt(pow(berr/bmean,2.)+pow(bmcerr/bmcmean,2.)));

	}
	else outHist->SetBinContent(j+1,0);
}

void setHistParameters (TH1F* inHist, std::string name, std::string type, float min, float max, double lower_bin, double higher_bin, std::string bin_type){
	if (type == "MPF") inHist->SetMarkerColor(8);
	if (type == "MPF1") inHist->SetMarkerColor(1);
	if (type == "MPF2") inHist->SetMarkerColor(2);
	if (type == "MPF3") inHist->SetMarkerColor(3);
	if (type == "Pt")  inHist->SetMarkerColor(9);
	if (type == "MPF_data")  inHist->SetMarkerColor(1);
	if (type == "MPF_mc")  inHist->SetMarkerColor(2);
	std::string bin_exp = "";
	std::string bin_units = "";
	if (bin_type == "pt"){
	bin_exp = "p^{avg}_{T}";
	bin_units = "GeV";
	}
	if (bin_type == "eta"){
	bin_exp = "|#eta|";
	bin_units = "";
	}		
    	inHist->SetMarkerStyle(21);
    	inHist->SetLineStyle(2);
    	inHist->SetLineColor(1);
    	inHist->SetTitle(Form("%s %s, %g<%s<%g %s", name.c_str(), type.c_str(), lower_bin, bin_exp.c_str(), higher_bin, bin_units.c_str()));
    	inHist->SetMaximum(max);
    	inHist->SetMinimum(min);
    	/*
    	if (lower_bin == 2.853){
    	inHist->SetMaximum(0.8);
    	inHist->SetMinimum(0.6);
    					}
     	if (lower_bin == 3.139){
    	inHist->SetMaximum(0.7);
    	inHist->SetMinimum(0.5);
    					}
     	if (lower_bin == 2.650){
    	inHist->SetMaximum(1.2);
    	inHist->SetMinimum(1.0);
    					}
    	*/
    	inHist->Write();
	}

void setHistParameters_fit (TH1F* inHist, TF1* fit, std::string name, std::string type, float min, float max, double lower_bin, double higher_bin, std::string bin_type){
	std::string bin_exp = "";
	std::string bin_units = "";	
	if (type == "MPF") {
	inHist->SetMarkerColor(8);
	fit->SetLineColor(8);
						}
	if (type == "Pt")  {
	inHist->SetMarkerColor(9);
	fit->SetLineColor(9);	
						}
	if (bin_type == "pt"){
	bin_exp = "p^{avg}_{T}";
	bin_units = "GeV";	
	}
	if (bin_type == "eta"){
	bin_exp = "|#eta|";
	bin_units = "";	
	}		
    	inHist->SetMarkerStyle(21);
    	inHist->SetLineStyle(2);
    	inHist->SetLineColor(1);
    	inHist->SetTitle(Form("%s %s, %g<%s<%g %s", name.c_str(), type.c_str(), lower_bin, bin_exp.c_str(), higher_bin, bin_units.c_str()));
    	inHist->SetMaximum(max);
    	inHist->SetMinimum(min);
    	/*
    	if (lower_bin == 2.853){
    	inHist->SetMaximum(0.8);
    	inHist->SetMinimum(0.6);
    					}
     	if (lower_bin == 3.139){
    	inHist->SetMaximum(0.7);
    	inHist->SetMinimum(0.5);
    					}
     	if (lower_bin == 2.650){
    	inHist->SetMaximum(1.2);
    	inHist->SetMinimum(1.0);
    					}
    	*/
    	inHist->Write();
	    fit->Write();
	}	
	
void drawMultiCanvas(std::string figure, TH1F** inHist1, TH1F** inHist2, std::string name1, std::string name2, int div1, int div2, int nbins, std::string xaxis, std::string yaxis, int xlog, std::string prefix){


TCanvas *canv = new TCanvas("canvas","",1700,700);
TLegend *l = new TLegend(0.4,0.83,0.7,0.93);
l->AddEntry(inHist1[1],name1.c_str());
l->AddEntry(inHist2[1],name2.c_str());
canv->Divide(div1,div2,0,0);
for(int i=0; i<nbins; i++){
	canv->cd(i+1);
	if (xlog == 1) gPad->SetLogx();
	inHist1[i]->GetYaxis()->SetTitle(yaxis.c_str());
    inHist1[i]->GetXaxis()->SetTitle(xaxis.c_str());
    inHist1[i]->Draw();
    inHist2[i]->Draw("same");
    l->Draw("");
  }
canv->SaveAs(Form("plots/%s_%s.gif",prefix.c_str(), figure.c_str()));
delete canv;
}

void sumResponse(std::string filename1="data_1512.root", std::string filename2="pythia_1512.root", std::string prefix = "1512_pythia"){
  double xbins_pt[] = {55.0,75.0,95.0,120.0,150.0,200.0,300.0}; 
  //double xbins_pt[] = {25.0,35.0,55.0,75.0,95.0,120.0,150.0,300.0}; 
  const int nbins_pt = sizeof(xbins_pt)/sizeof(double)-1;

  //double xbins_eta[] = {0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 3.139, 5.191};
  //double xbins_eta[] = {0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 4.013, 5.191};
  double xbins_eta[] = {0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964};
  
  const int nbins_eta = sizeof(xbins_eta)/sizeof(double)-1;

  double xbins_alpha[] = {0.1,0.15,0.2,0.25,0.3,0.35,0.4};
  const int nbins_alpha = sizeof(xbins_alpha)/sizeof(double)-1;

  TH1F *hRelResponse[nbins_pt];
  TH1F *hRelResponseMC[nbins_pt];
  TH1F *hMPFResponse[nbins_pt];
  TH1F *hMPFResponseMC[nbins_pt];

  TH1F *hMPFResponsePt[nbins_eta];
  TH1F *hMPFResponseMCPt[nbins_eta];
  TH1F *hRelResponsePt[nbins_eta];
  TH1F *hRelResponseMCPt[nbins_eta];

  TH1F *hMPFRatio[nbins_eta];
  TH1F *hRelRatio[nbins_eta];
  TH1F *hMPFRatioPt[nbins_pt];
  TH1F *hMPFRatioPt1[nbins_pt];
  TH1F *hMPFCRatioPt[nbins_pt];
  TH1F *hMPFCRatioPt2[nbins_pt];
  TH1F *hMPFCRatioPt3[nbins_pt];  
  TH1F *hRelRatioPt[nbins_pt];
  
  TH1F *hMPFRatioPt_cfit;
  TH1F *hRelRatioPt_cfit;
  TH1F *hMPFRatioPt_lfit;
  TH1F *hRelRatioPt_lfit;

  TH1D *h_avgA[nbins_pt][nbins_eta];
  TH1D *h_avgB[nbins_pt][nbins_eta];
  TH1D *h_avgBMC[nbins_pt][nbins_eta];
  TH1D *h_avgAMC[nbins_pt][nbins_eta];

  TH1D *h_avgA_alpha[nbins_pt][nbins_eta][nbins_alpha];
  TH1D *h_avgB_alpha[nbins_pt][nbins_eta][nbins_alpha];

  TH1D *h_avgAMC_alpha[nbins_pt][nbins_eta][nbins_alpha];
  TH1D *h_avgBMC_alpha[nbins_pt][nbins_eta][nbins_alpha];

  TH1F *hMPFRatio_alpha[nbins_pt][nbins_eta];
  TH1F *hRelRatio_alpha[nbins_pt][nbins_eta];
  
  TH1F *kFSR_Rel[nbins_pt];
  TH1F *kFSR_MPF[nbins_pt];
  
    
  TH1F *kFSR_Rel_Avg;
  TH1F *kFSR_MPF_Avg;

  TFile *fin = new TFile(filename1.c_str());
  TFile *finmc = new TFile(filename2.c_str());
// Get Histograms and create new ones
    kFSR_Rel_Avg = new TH1F("kFSR_Rel_Avg","",nbins_eta,xbins_eta);
	kFSR_MPF_Avg = new TH1F("kFSR_MPF_Avg","",nbins_eta,xbins_eta);
	
	hMPFRatioPt_cfit = new TH1F("hMPFRatioPt_cfit","",nbins_eta,xbins_eta);
    hRelRatioPt_cfit = new TH1F("hRelRatioPt_cfit","",nbins_eta,xbins_eta);
    hMPFRatioPt_lfit = new TH1F("hMPFRatioPt_lfit","",nbins_eta,xbins_eta);
    hRelRatioPt_lfit = new TH1F("hRelRatioPt_lfit","",nbins_eta,xbins_eta);

  for(int i=0; i<nbins_pt; i++){
    hRelResponse[i] = new TH1F(Form("hRelResponse_pt%d",i),"",nbins_eta,xbins_eta);
    hRelResponseMC[i] = new TH1F(Form("hRelResponseMC_pt%d",i),"",nbins_eta,xbins_eta);
    hMPFResponse[i] = new TH1F(Form("hMPFResponse_pt%d",i),"",nbins_eta,xbins_eta);
    hMPFResponseMC[i] = new TH1F(Form("hMPFResponseMC_pt%d",i),"",nbins_eta,xbins_eta);
    hMPFRatioPt[i] = new TH1F(Form("hMPFRatio_pt%d",i),"",nbins_eta,xbins_eta);
    hMPFRatioPt1[i] = new TH1F(Form("hMPFRatio_pt1%d",i),"",nbins_eta,xbins_eta);
    hMPFCRatioPt[i] = new TH1F(Form("hMPFCRatio_pt%d",i),"",nbins_eta,xbins_eta);
    hMPFCRatioPt2[i] = new TH1F(Form("hMPFCRatio_pt2%d",i),"",nbins_eta,xbins_eta);
    hMPFCRatioPt3[i] = new TH1F(Form("hMPFCRatio_pt3%d",i),"",nbins_eta,xbins_eta);    
    hRelRatioPt[i] = new TH1F(Form("hRelRatio_pt%d",i),"",nbins_eta,xbins_eta);
        
    kFSR_Rel[i] = new TH1F(Form("kFSR_Rel_pt%d",i),"",nbins_eta,xbins_eta);
	kFSR_MPF[i] = new TH1F(Form("kFSR_MPF_pt%d",i),"",nbins_eta,xbins_eta);

    for(int j=0; j<nbins_eta; j++){

	h_avgA[i][j]=(TH1D*)fin->Get(Form("h_avgA_%d_%d",i,j));
	h_avgB[i][j]=(TH1D*)fin->Get(Form("h_avgB_%d_%d",i,j));

	h_avgBMC[i][j]=(TH1D*)finmc->Get(Form("h_avgB_%d_%d",i,j));
	h_avgAMC[i][j]=(TH1D*)finmc->Get(Form("h_avgA_%d_%d",i,j));

	hRelRatio_alpha[i][j] = new TH1F(Form("hRelRatio_alpha%d_%d",i,j),"",nbins_alpha,xbins_alpha);
	hMPFRatio_alpha[i][j] = new TH1F(Form("hMPFRatio_alpha%d_%d",i,j),"",nbins_alpha,xbins_alpha);
      
        for(int k=0; k<nbins_alpha; k++){
      h_avgA_alpha[i][j][k]=(TH1D*)fin->Get(Form("h_avgA_alpha_%d_%d_%d",i,j,k));
      h_avgB_alpha[i][j][k]=(TH1D*)fin->Get(Form("h_avgB_alpha_%d_%d_%d",i,j,k));

      h_avgAMC_alpha[i][j][k]=(TH1D*)finmc->Get(Form("h_avgA_alpha_%d_%d_%d",i,j,k));
      h_avgBMC_alpha[i][j][k]=(TH1D*)finmc->Get(Form("h_avgB_alpha_%d_%d_%d",i,j,k));
     }
    }
  }
  
  for(int j=0; j<nbins_eta; j++){

    hMPFResponsePt[j] = new TH1F(Form("hMPFResponse_eta%d",j),"",nbins_pt,xbins_pt);
    hMPFResponseMCPt[j] = new TH1F(Form("hMPFResponseMC_eta%d",j),"",nbins_pt,xbins_pt);

    hRelResponsePt[j] = new TH1F(Form("hRelResponse_eta%d",j),"",nbins_pt,xbins_pt);
    hRelResponseMCPt[j] = new TH1F(Form("hRelResponseMC_eta%d",j),"",nbins_pt,xbins_pt);

    hMPFRatio[j] = new TH1F(Form("hMPFRatio_eta%d",j),"",nbins_pt,xbins_pt);
    hRelRatio[j] = new TH1F(Form("hRelRatio_eta%d",j),"",nbins_pt,xbins_pt);
				}


//Filling histograms
for(int i=0; i<nbins_pt; i++){
    for(int j=0; j<nbins_eta; j++){
	fillRelHist(h_avgA[i][j], hRelResponse[i], j);
	fillRelHist(h_avgAMC[i][j], hRelResponseMC[i], j);
	fillMPFHist(h_avgB[i][j], hMPFResponse[i], j);
	fillMPFHist(h_avgBMC[i][j],  hMPFResponseMC[i], j);
	fillRelRatioHist(h_avgA[i][j], h_avgAMC[i][j], hRelRatioPt[i], j);
	fillMPFRatioHist(h_avgB[i][j], h_avgBMC[i][j], hMPFRatioPt[i], j);
    }
  }
  

	  

	  for(int i=0; i<nbins_pt; i++){
	  	  for(int j=0; j<nbins_eta; j++){
	fillRelHist(h_avgA[i][j], hRelResponsePt[j], i);
	fillRelHist(h_avgAMC[i][j], hRelResponseMCPt[j], i);
	fillMPFHist(h_avgB[i][j], hMPFResponsePt[j], i);
	fillMPFHist(h_avgBMC[i][j], hMPFResponseMCPt[j], i);
	fillRelRatioHist(h_avgA[i][j], h_avgAMC[i][j], hRelRatio[j], i);
	fillMPFRatioHist(h_avgB[i][j], h_avgBMC[i][j], hMPFRatio[j], i);
		for(int k=0; k<nbins_alpha; k++){
		fillRelRatioHist(h_avgA_alpha[i][j][k], h_avgAMC_alpha[i][j][k], hRelRatio_alpha[i][j], k);
		fillMPFRatioHist(h_avgB_alpha[i][j][k], h_avgBMC_alpha[i][j][k], hMPFRatio_alpha[i][j], k);
		}
	}
}

//Doing the fits
TF1* mpffit[nbins_eta];
TF1* relfit[nbins_eta];

TF1* mpffit_log[nbins_eta];
TF1* relfit_log[nbins_eta];

for(int j=0; j<nbins_eta; j++){
   mpffit[j] = new TF1(Form("mpffit%d",j),"[0]",25,300);
   relfit[j] = new TF1(Form("relfit%d",j),"[0]",25,300);
   
   mpffit_log[j] = new TF1(Form("mpffit_log%d",j),"[0]+[1]*TMath::Log(x)",25,300);
   relfit_log[j] = new TF1(Form("relfit_log%d",j),"[0]+[1]*TMath::Log(x)",25,300);
   mpffit_log[j]->SetParameter(0,1.0);
   mpffit_log[j]->SetParameter(1,0.1);
   relfit_log[j]->SetParameter(0,1.0);
   relfit_log[j]->SetParameter(1,0.1);
   	//	}
   
   relfit[j]->SetLineColor(9);
   mpffit[j]->SetLineColor(8);
   hMPFRatio[j]->Fit(mpffit[j], "R");
   hRelRatio[j]->Fit(relfit[j], "R");
   relfit_log[j]->SetLineColor(9);
   mpffit_log[j]->SetLineColor(8);
   relfit_log[j]->SetLineStyle(4);
   mpffit_log[j]->SetLineStyle(4);
   hMPFRatio[j]->Fit(mpffit_log[j], "R+");
   hRelRatio[j]->Fit(relfit_log[j], "R+");   
				}

TF1* mpffit_alpha[nbins_pt][nbins_eta];
TF1* relfit_alpha[nbins_pt][nbins_eta];
for(int i=0; i<nbins_pt; i++){
	for(int j=0; j<nbins_eta; j++){

	mpffit_alpha[i][j] = new TF1(Form("mpffit_alpha%d_%d",i,j),"[0]+[1]*x",0,0.4);
	mpffit_alpha[i][j]->SetLineColor(8);
	mpffit_alpha[i][j]->SetParameter(0,1.0);
	mpffit_alpha[i][j]->SetParameter(1,0.01);
	relfit_alpha[i][j] = new TF1(Form("relfit_alpha%d_%d",i,j),"[0]+[1]*x",0,0.4);
	relfit_alpha[i][j]->SetParameter(0,1.0);
	relfit_alpha[i][j]->SetParameter(1,0.01);
	relfit_alpha[i][j]->SetLineColor(9);

   hMPFRatio_alpha[i][j]->Fit(mpffit_alpha[i][j], "R");
   hRelRatio_alpha[i][j]->Fit(relfit_alpha[i][j], "R");
				}
			}

for(int i=0; i<nbins_pt; i++){
	for(int j=0; j<nbins_eta; j++){
	kFSR_Rel[i]->SetBinContent(j+1,relfit_alpha[i][j]->Eval(0)/relfit_alpha[i][j]->Eval(0.2));
	kFSR_Rel[i]->SetBinError(j+1,relfit_alpha[i][j]->GetParError(0));
	kFSR_MPF[i]->SetBinContent(j+1,mpffit_alpha[i][j]->Eval(0)/mpffit_alpha[i][j]->Eval(0.2));
	kFSR_MPF[i]->SetBinError(j+1,mpffit_alpha[i][j]->GetParError(0));
	}
}

double relfit_avg[nbins_eta];
double mpffit_avg[nbins_eta];
	for(int j=0; j<nbins_eta; j++){
		for(int i=0; i<nbins_pt; i++){
	relfit_avg[j]+=relfit_alpha[i][j]->Eval(0)/relfit_alpha[i][j]->Eval(0.2);
	mpffit_avg[j]+=mpffit_alpha[i][j]->Eval(0)/mpffit_alpha[i][j]->Eval(0.2);
	}
}
	for(int i=0; i<nbins_eta; i++){
		kFSR_Rel_Avg->SetBinContent(i+1,relfit_avg[i]/nbins_pt);
		kFSR_Rel_Avg->SetBinError(i+1,0.01);		
		kFSR_MPF_Avg->SetBinContent(i+1,mpffit_avg[i]/nbins_pt);
		kFSR_MPF_Avg->SetBinError(i+1,0.01);	
									}

TF1* mpffit_alpha_avg;
TF1* relfit_alpha_avg;
	mpffit_alpha_avg = new TF1("mpffit_alpha_avg","[0]",xbins_eta[0],xbins_eta[nbins_eta]);
	mpffit_alpha_avg->SetLineColor(8);
	relfit_alpha_avg = new TF1("relfit_alpha_avg","[0]+[1]*TMath::CosH(x)/(1+[2]*TMath::CosH(x))",xbins_eta[0],xbins_eta[nbins_eta-1]);
	relfit_alpha_avg->SetParameter(0,0.99);
	relfit_alpha_avg->SetParameter(1,2.2e-03);
	relfit_alpha_avg->SetParameter(2,1.6e-02);
	relfit_alpha_avg->SetLineColor(9);
	kFSR_MPF_Avg->Fit(mpffit_alpha_avg,"R");
	kFSR_Rel_Avg->Fit(relfit_alpha_avg,"R");

 for(int i=0; i<nbins_eta; i++){
hMPFRatioPt_cfit->SetBinContent(i+1,mpffit[i]->GetParameter(0));
hMPFRatioPt_cfit->SetBinError(i+1,mpffit[i]->GetParError(0));
hRelRatioPt_cfit->SetBinContent(i+1,relfit[i]->GetParameter(0));
hRelRatioPt_cfit->SetBinError(i+1,relfit[i]->GetParError(0));

//hMPFRatioPt_lfit->SetBinContent(i+1,mpffit_log[i]->Eval(177));
hMPFRatioPt_lfit->SetBinContent(i+1,mpffit_log[i]->Eval(122));
hMPFRatioPt_lfit->SetBinError(i+1,mpffit_log[i]->GetParError(1));
hRelRatioPt_lfit->SetBinContent(i+1,relfit_log[i]->Eval(122));
//hRelRatioPt_lfit->SetBinContent(i+1,relfit_log[i]->Eval(177));
hRelRatioPt_lfit->SetBinError(i+1,relfit_log[i]->GetParError(1));
								}

TFile *fout = new TFile(Form("%s_response_fits.root",prefix.c_str()),"recreate");
fout->cd();

setHistParameters_fit(kFSR_Rel_Avg,relfit_alpha_avg,"kFSR Avg","Pt",0.9,1.1,xbins_pt[0],xbins_pt[1],"pt avg");
setHistParameters_fit(kFSR_MPF_Avg,mpffit_alpha_avg,"kFSR Avg","MPF",0.9,1.1,xbins_pt[0],xbins_pt[1],"pt avg");

//setHistParameters(kFSR_Rel_Avg,"kFSR Avg","Pt",0.9,1.1,xbins_pt[0],xbins_pt[1],"pt avg");
//setHistParameters(kFSR_MPF_Avg,"kFSR Avg","MPF",0.9,1.1,xbins_pt[0],xbins_pt[1],"pt avg");

setHistParameters(hRelRatioPt_cfit,"Response Ratio constant fit","Pt",0.8,1.2,xbins_pt[0],xbins_pt[1],"pt avg");
setHistParameters(hRelRatioPt_lfit,"Response Ratio loglin fit","Pt",0.8,1.2,xbins_pt[0],xbins_pt[1],"pt avg");  
setHistParameters(hMPFRatioPt_cfit,"Response Ratio constant fit","MPF",0.8,1.2,xbins_pt[0],xbins_pt[1],"pt avg");
setHistParameters(hMPFRatioPt_lfit,"Response Ratio loglin fit","MPF",0.8,1.2,xbins_pt[0],xbins_pt[1],"pt avg");  


for(int i=0; i<nbins_pt; i++){
	for(int j=0; j<nbins_eta; j++){
		if (j > 11){
//setHistParameters_fit(hRelRatio_alpha[i][j],relfit_alpha[i][j],"Response Ratio loglin fit","Pt",0.5,0.9,xbins_eta[j],xbins_eta[j+1],"eta");
//setHistParameters_fit(hMPFRatio_alpha[i][j],mpffit_alpha[i][j],"Response Ratio loglin fit","MPF",0.5,0.9,xbins_eta[j],xbins_eta[j+1],"eta");     	    
    	    }
		else{
setHistParameters_fit(hRelRatio_alpha[i][j],relfit_alpha[i][j],"Response Ratio loglin fit","Pt",0.8,1.2,xbins_eta[j],xbins_eta[j+1],"eta");
setHistParameters_fit(hMPFRatio_alpha[i][j],mpffit_alpha[i][j],"Response Ratio loglin fit","MPF",0.8,1.2,xbins_eta[j],xbins_eta[j+1],"eta");   
			}  	
    }
}
    
  for(int i=0; i<nbins_pt; i++){
setHistParameters(hRelResponse[i],"Response Data","Pt",0.9,1.1,xbins_pt[i],xbins_pt[i+1],"pt");
setHistParameters(hRelResponseMC[i],"Response MC","Pt",0.9,1.1,xbins_pt[i],xbins_pt[i+1],"pt");
setHistParameters(hMPFResponse[i],"Response Data","MPF_data",0.9,1.1,xbins_pt[i],xbins_pt[i+1],"pt");
setHistParameters(hMPFResponseMC[i],"Response MC","MPF_mc",0.9,1.1,xbins_pt[i],xbins_pt[i+1],"pt");
setHistParameters(hRelRatioPt[i],"Response Ratio","Pt",0.6,1.2,xbins_pt[i],xbins_pt[i+1],"pt");
setHistParameters(hMPFRatioPt[i],"Response Ratio","MPF",0.6,1.2,xbins_pt[i],xbins_pt[i+1],"pt");
setHistParameters(kFSR_Rel[i],"kFSR","Pt",0.8,1.2,xbins_pt[i],xbins_pt[i+1],"pt");
setHistParameters(kFSR_MPF[i],"kFSR","MPF",0.8,1.2,xbins_pt[i],xbins_pt[i+1],"pt");
  								}

  for(int i=0; i<nbins_eta; i++){
setHistParameters(hRelResponsePt[i],"Response Data","Pt",0.9,1.1,xbins_eta[i],xbins_eta[i+1],"eta");  
setHistParameters(hMPFResponsePt[i],"Response Data","MPF",0.9,1.1,xbins_eta[i],xbins_eta[i+1],"eta"); 
setHistParameters(hRelResponseMCPt[i],"Response MC","Pt",0.9,1.1,xbins_eta[i],xbins_eta[i+1],"eta");  
setHistParameters(hMPFResponseMCPt[i],"Response MC","MPF",0.9,1.1,xbins_eta[i],xbins_eta[i+1],"eta"); 
setHistParameters_fit(hRelRatio[i],relfit[i],"Response Ratio","Pt",0.9,1.1,xbins_eta[i],xbins_eta[i+1],"eta");  
setHistParameters_fit(hRelRatio[i],relfit_log[i],"Response Ratio","Pt",0.9,1.1,xbins_eta[i],xbins_eta[i+1],"eta"); 
setHistParameters_fit(hMPFRatio[i],mpffit[i],"Response Ratio","MPF",0.9,1.1,xbins_eta[i],xbins_eta[i+1],"eta"); 
setHistParameters_fit(hMPFRatio[i],mpffit_log[i],"Response Ratio","MPF",0.9,1.1,xbins_eta[i],xbins_eta[i+1],"eta"); 
								}


cout<<"going to draw histograms"<<endl;
gStyle->SetOptStat(0);

drawMultiCanvas("response_ratio_vs_pt", hRelRatio, hMPFRatio, "Pt", "MPF", 6, 2, nbins_eta, "P_{T} avg, GeV", "Rmc/Rdata",1,prefix);
drawMultiCanvas("response_ratio_vs_eta", hRelRatioPt, hMPFRatioPt, "Pt", "MPF", 4, 1, nbins_pt, "|#eta|", "Rmc/Rdata",0,prefix);
drawMultiCanvas("kfsr_vs_eta", kFSR_Rel, kFSR_MPF, "Pt", "MPF", 4, 1, nbins_pt, "|#eta|", "K_fsr",0,prefix);

TCanvas *c7 = new TCanvas("c7","",1700,700);
TLegend *l7 = new TLegend(0.4,0.8,0.7,0.9);
    l7->AddEntry(hMPFRatio_alpha[1][0],"MPF");
    l7->AddEntry(hRelRatio_alpha[1][0],"Pt");
    l7->SetTextSize(0.045); 
c7->Divide(6,2,0,0);
  for(int j=0; j<nbins_eta; j++){
c7->cd(j+1);
    mpffit_alpha[0][j]->SetLineColor(8);
    relfit_alpha[0][j]->SetLineColor(9);
    hMPFRatio_alpha[0][j]->SetMarkerColor(8);
    hMPFRatio_alpha[0][j]->SetMarkerStyle(21);
    hMPFRatio_alpha[0][j]->SetTitle(Form("%g < |#eta| < %g ",xbins_eta[j],xbins_eta[j+1]));
    hMPFRatio_alpha[0][j]->GetYaxis()->SetTitle("Rmc/Rdata");
    hMPFRatio_alpha[0][j]->GetXaxis()->SetTitle("cut on #alpha");
    hRelRatio_alpha[0][j]->SetMarkerColor(9);
    hRelRatio_alpha[0][j]->SetMarkerStyle(21);
	hMPFRatio_alpha[0][j]->Draw();
	hRelRatio_alpha[0][j]->Draw("same");
l7->Draw("");
				}
c7->SaveAs(Form("plots/%s_rel_mpf_alpha_1.gif",prefix.c_str()));
  
  TCanvas *c8 = new TCanvas("c8","",1700,700);
TLegend *l8 = new TLegend(0.4,0.8,0.7,0.9);
    l8->AddEntry(hMPFRatio_alpha[2][0],"MPF");
    l8->AddEntry(hRelRatio_alpha[2][0],"Pt");
    l8->SetTextSize(0.045); 
c8->Divide(6,2,0,0);
  for(int j=0; j<nbins_eta; j++){
c8->cd(j+1);
    mpffit_alpha[1][j]->SetLineColor(8);
    relfit_alpha[1][j]->SetLineColor(9);
    hMPFRatio_alpha[1][j]->SetMarkerColor(8);
    hMPFRatio_alpha[1][j]->SetMarkerStyle(21);
    hMPFRatio_alpha[1][j]->SetTitle(Form("%g < |#eta| < %g ",xbins_eta[j],xbins_eta[j+1]));
    hMPFRatio_alpha[1][j]->GetYaxis()->SetTitle("Rmc/Rdata");
    hMPFRatio_alpha[1][j]->GetXaxis()->SetTitle("cut on #alpha");
    hRelRatio_alpha[1][j]->SetMarkerColor(9);
    hRelRatio_alpha[1][j]->SetMarkerStyle(21);
	hMPFRatio_alpha[1][j]->Draw();
	hRelRatio_alpha[1][j]->Draw("same");
l8->Draw("");
				}
c8->SaveAs(Form("plots/%s_rel_mpf_alpha_2.gif",prefix.c_str()));

  TCanvas *c9 = new TCanvas("c9","",1700,700);
TLegend *l9 = new TLegend(0.4,0.8,0.7,0.9);
    l9->AddEntry(hMPFRatio_alpha[3][0],"MPF");
    l9->AddEntry(hRelRatio_alpha[3][0],"Pt");
    l9->SetTextSize(0.045); 
c9->Divide(6,2,0,0);
  for(int j=0; j<nbins_eta; j++){
c9->cd(j+1);
    mpffit_alpha[2][j]->SetLineColor(8);
    relfit_alpha[2][j]->SetLineColor(9);
    hMPFRatio_alpha[2][j]->SetMarkerColor(8);
    hMPFRatio_alpha[2][j]->SetMarkerStyle(21);
    hMPFRatio_alpha[2][j]->SetTitle(Form("%g < |#eta| < %g ",xbins_eta[j],xbins_eta[j+1]));
    hMPFRatio_alpha[2][j]->GetYaxis()->SetTitle("Rmc/Rdata");
    hMPFRatio_alpha[2][j]->GetXaxis()->SetTitle("cut on #alpha");
    hRelRatio_alpha[2][j]->SetMarkerColor(9);
    hRelRatio_alpha[2][j]->SetMarkerStyle(21);
	hMPFRatio_alpha[2][j]->Draw();
	hRelRatio_alpha[2][j]->Draw("same");
l9->Draw("");
				}
c9->SaveAs(Form("plots/%s_rel_mpf_alpha_3.gif",prefix.c_str()));

  TCanvas *c12 = new TCanvas("c12","",1700,700);
TLegend *l12 = new TLegend(0.4,0.8,0.7,0.9);
    l12->AddEntry(hMPFRatio_alpha[4][0],"MPF");
    l12->AddEntry(hRelRatio_alpha[4][0],"Pt");
    l12->SetTextSize(0.045); 
c12->Divide(6,2,0,0);
  for(int j=0; j<nbins_eta; j++){
c12->cd(j+1);
    mpffit_alpha[3][j]->SetLineColor(8);
    relfit_alpha[3][j]->SetLineColor(9);
    hMPFRatio_alpha[3][j]->SetMarkerColor(8);
    hMPFRatio_alpha[3][j]->SetMarkerStyle(21);
    hMPFRatio_alpha[3][j]->SetTitle(Form("%g < |#eta| < %g ",xbins_eta[j],xbins_eta[j+1]));
    hMPFRatio_alpha[3][j]->GetYaxis()->SetTitle("Rmc/Rdata");
    hMPFRatio_alpha[3][j]->GetXaxis()->SetTitle("cut on #alpha");
    hRelRatio_alpha[3][j]->SetMarkerColor(9);
    hRelRatio_alpha[3][j]->SetMarkerStyle(21);
	hMPFRatio_alpha[3][j]->Draw();
	hRelRatio_alpha[3][j]->Draw("same");
l12->Draw("");
				}
c12->SaveAs(Form("plots/%s_rel_mpf_alpha_4.gif",prefix.c_str()));

  TCanvas *c13 = new TCanvas("c13","",1700,700);
TLegend *l13 = new TLegend(0.4,0.8,0.7,0.9);
    l13->AddEntry(hMPFRatio_alpha[5][0],"MPF");
    l13->AddEntry(hRelRatio_alpha[5][0],"Pt");
    l13->SetTextSize(0.045); 
c13->Divide(6,2,0,0);
  for(int j=0; j<nbins_eta; j++){
c13->cd(j+1);
    mpffit_alpha[4][j]->SetLineColor(8);
    relfit_alpha[4][j]->SetLineColor(9);
    hMPFRatio_alpha[4][j]->SetMarkerColor(8);
    hMPFRatio_alpha[4][j]->SetMarkerStyle(21);
    hMPFRatio_alpha[4][j]->SetTitle(Form("%g < |#eta| < %g ",xbins_eta[j],xbins_eta[j+1]));
    hMPFRatio_alpha[4][j]->GetYaxis()->SetTitle("Rmc/Rdata");
    hMPFRatio_alpha[4][j]->GetXaxis()->SetTitle("cut on #alpha");
    hRelRatio_alpha[4][j]->SetMarkerColor(9);
    hRelRatio_alpha[4][j]->SetMarkerStyle(21);
	hMPFRatio_alpha[4][j]->Draw();
	hRelRatio_alpha[4][j]->Draw("same");
l13->Draw("");
				}
c13->SaveAs(Form("plots/%s_rel_mpf_alpha_5.gif",prefix.c_str()));

TCanvas *c14 = new TCanvas("c14","",700,700);
TLegend *l14 = new TLegend(0.4,0.8,0.7,0.9);
    l14->AddEntry(kFSR_MPF_Avg,"MPF");
    l14->AddEntry(kFSR_Rel_Avg,"Pt");
    l14->SetTextSize(0.045);
    mpffit_alpha_avg->SetLineColor(8);
    relfit_alpha_avg->SetLineColor(9);
	kFSR_Rel_Avg->SetMarkerColor(9);
    kFSR_Rel_Avg->SetMarkerStyle(21);
    kFSR_Rel_Avg->SetTitle("Kfsr Avg");
    kFSR_Rel_Avg->GetYaxis()->SetTitle("Kfsr");
    kFSR_Rel_Avg->GetXaxis()->SetTitle("|#eta|");

    
    kFSR_MPF_Avg->SetMarkerColor(8);
    kFSR_MPF_Avg->SetMarkerStyle(21);
    kFSR_MPF_Avg->SetTitle("Kfsr MPF Avg");

   	kFSR_Rel_Avg->Draw();
	kFSR_MPF_Avg->Draw("same");	
	l14->Draw(""); 
c14->SaveAs(Form("plots/%s_kfsr_avg.gif",prefix.c_str())); 
 
TCanvas *c6 = new TCanvas("c6","",700,700);
TLegend *l6 = new TLegend(0.15,0.75,0.45,0.88);

	hRelRatioPt_cfit->SetMarkerColor(9);
    hRelRatioPt_cfit->SetMarkerStyle(21);
    hRelRatioPt_cfit->SetTitle("Response ratio");
    hRelRatioPt_cfit->GetYaxis()->SetTitle("R_{mc}/R_{data}");
    hRelRatioPt_cfit->GetXaxis()->SetTitle("|#eta|");

    
    hMPFRatioPt_cfit->SetMarkerColor(8);
    hMPFRatioPt_cfit->SetMarkerStyle(21);
    hMPFRatioPt_cfit->SetTitle("Response MPF Avg");

    hMPFRatioPt_lfit->SetMarkerColor(8);
    hMPFRatioPt_lfit->SetMarkerStyle(25);
    hMPFRatioPt_lfit->SetTitle("Response MPF Avg");
    
    hRelRatioPt_lfit->SetMarkerColor(9);
    hRelRatioPt_lfit->SetMarkerStyle(25);
    hRelRatioPt_lfit->SetTitle("Response Rel Avg"); 
       
    l6->AddEntry(hMPFRatioPt_cfit,"MPF Flat");
    l6->AddEntry(hRelRatioPt_cfit,"Pt Flat");
    l6->AddEntry(hMPFRatioPt_lfit,"MPF Loglin");
    l6->AddEntry(hRelRatioPt_lfit,"Pt Loglin");    
    l6->SetTextSize(0.03);
   	hRelRatioPt_cfit->Draw();
	hMPFRatioPt_cfit->Draw("same");	
   	hRelRatioPt_lfit->Draw("same");
	hMPFRatioPt_lfit->Draw("same");		
	l6->Draw("");  
c6->SaveAs(Form("plots/%s_avg_response.gif",prefix.c_str()));

  fout->Close();


}
