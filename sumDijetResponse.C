
//Plotting macro to simply sum all the stuff in the deriveDijetResponse macro
//K. Jung, December 2015

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"

using namespace std;

void sumDijetResponse(std::string filename){

	bool doDraw = true;

	const int nbins_pt = 5;
	double xbins_pt[nbins_pt+1] = {40,60,80,110,200,1000};
	const int nbins_eta = 20;
	double xbins_eta[nbins_eta+1];
	for(int i=0; i<=nbins_eta; i++){
		xbins_eta[i] = -5+(10./(double)nbins_eta)*i;
	}

	TH1F *hRelResponse[nbins_pt];
	TH1F *hMPFResponse[nbins_pt];

	TH1F *hAbsPhoResponse[nbins_pt];
	TH1F *hMPFAbsPhoResponse[nbins_pt];

	TH1F *avgAHisto[nbins_pt][nbins_eta];
	TH1F *avgBHisto[nbins_pt][nbins_eta];
	TH1F *hAvgAbsPhoResponse[nbins_pt][nbins_eta];

	TFile *fin = new TFile(filename.c_str());

	for(int i=0; i<nbins_pt; i++){
		hRelResponse[i] = new TH1F(Form("hRelResponse_pt%d",i),"",nbins_eta,xbins_eta);
		hMPFResponse[i] = new TH1F(Form("hMPFResponse_pt%d",i),"",nbins_eta,xbins_eta);
		hAbsPhoResponse[i] = new TH1F(Form("hAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);
		hMPFAbsPhoResponse[i] = new TH1F(Form("hMPFAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);
		for(int j=0; j<nbins_eta; j++){
			avgAHisto[i][j] = (TH1F*)fin->Get(Form("avgAHisto_pt%d_eta%d",i,j))->Clone(Form("avgAHisto_pt%d_eta%d",i,j));
			avgBHisto[i][j] = (TH1F*)fin->Get(Form("avgBHisto_pt%d_eta%d",i,j))->Clone(Form("avgBHisto_pt%d_eta%d",i,j));
			//hAvgAbsPhoResponse[i][j] = (TH1F*)fin->Get(Form("hAvgAbsPhoResponse_pt%d_eta%d",i,j))->Clone(Form("hAvgAbsPhoResponse_pt%d_eta%d",i,j));
		}
	}

	for(int i=0; i<nbins_pt; i++){
		for(int j=0; j<nbins_eta; j++){

		    hRelResponse[i]->SetBinContent(j+1,(1+avgAHisto[i][j]->GetMean())/(1-avgAHisto[i][j]->GetMean()));
			hRelResponse[i]->SetBinError(j+1,hRelResponse[i]->GetBinContent(j+1)*(1./sqrt(avgAHisto[i][j]->GetEntries())));

			hMPFResponse[i]->SetBinContent(j+1,(1+avgBHisto[i][j]->GetMean())/(1-avgBHisto[i][j]->GetMean()));
			hMPFResponse[i]->SetBinError(j+1,hMPFResponse[i]->GetBinContent(j+1)*(1./sqrt(avgBHisto[i][j]->GetEntries())));

			//hMPFAbsPhoResponse[i]->SetBinContent(j+1,hAvgAbsPhoResponse[i][j]->GetMean());
			//hMPFAbsPhoResponse[i]->SetBinError(j+1,hMPFAbsPhoResponse[i]->GetBinContent(j+1)*(1./sqrt(hAvgAbsPhoResponse[i][j]->GetEntries())));

			int totEntriesA=0, totEntriesB=0, totEntriesAbs=0;
			for(int j=0; j<nbins_eta; j++){ 
				totEntriesA+=avgAHisto[i][j]->GetEntries(); 
				totEntriesB+=avgBHisto[i][j]->GetEntries();
				//totEntriesAbs+=hAvgAbsPhoResponse[i][j]->GetEntries();
			}
			hRelResponse[i]->SetEntries(totEntriesA);
			hMPFResponse[i]->SetEntries(totEntriesB);
			//hMPFAbsPhoResponse[i]->SetEntries(totEntriesAbs);
		}
	}
	
	int color[5] = {1,2,4,8,20};
	TFile *fout = new TFile("simplePlots.root","recreate");
	fout->cd();
	for(int i=0; i<nbins_pt; i++){
		hRelResponse[i]->SetMarkerColor(color[i]);
		hRelResponse[i]->SetLineColor(color[i]);
		hRelResponse[i]->SetTitle(Form("Rel, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
		hRelResponse[i]->Write();

		hMPFResponse[i]->SetMarkerColor(color[i]);
		hMPFResponse[i]->SetLineColor(color[i]);
		hMPFResponse[i]->SetMarkerStyle(21);
		hMPFResponse[i]->SetLineStyle(2);
		hRelResponse[i]->SetTitle(Form("MPF, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
		hMPFResponse[i]->Write();

		hMPFAbsPhoResponse[i]->SetMarkerColor(color[i]);
		hMPFAbsPhoResponse[i]->SetLineColor(color[i]);
		hMPFAbsPhoResponse[i]->SetTitle(Form("MPF Abs, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
		hMPFAbsPhoResponse[i]->Write();
	}

	fout->Close();

	if(doDraw){
		TCanvas *c1 = new TCanvas("c1","",600,600);
		c1->cd();
		TLegend *l1 = new TLegend(0.4,0.8,0.7,0.9);
		hMPFResponse[0]->Draw();
		l1->AddEntry(hMPFResponse[0],Form("%g<p_{T}<%g",xbins_pt[0],xbins_pt[1]));
		for(int i=1; i<nbins_pt; i++){
			hMPFResponse[i]->Draw("same");
			l1->AddEntry(hMPFResponse[i],Form("%g<p_{T}<%g",xbins_pt[i],xbins_pt[i+1]));
		}
		l1->Draw("Same");

		TCanvas *c2 = new TCanvas("c2","",600,600);
		c2->cd();
		hMPFAbsPhoResponse[0]->Draw();
		TLegend *l2 = new TLegend(0.4,0.8,0.7,0.9);
		l2->AddEntry(hMPFAbsPhoResponse[0],Form("%g<p_{T}<%g",xbins_pt[0],xbins_pt[1]));
		for(int i=1; i<nbins_pt; i++){
			hMPFAbsPhoResponse[i]->Draw("same");
			l2->AddEntry(hMPFAbsPhoResponse[i],Form("%g<p_{T}<%g",xbins_pt[i],xbins_pt[i+1]));
		}
		l2->Draw("Same");

	}

}