#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
#include "TLorentzVector.h"

using namespace std;

int findEtaBin(float eta, int nbins_eta, double* xbins_eta){

	for(int i=0; i<nbins_eta; i++){
		if(eta<xbins_eta[i+1]) return i;
	}
	return nbins_eta-1;
}

int findPtBin(float pt, int nbins_pt, double* xbins_pt){

	for(int i=0; i<nbins_pt; i++){
		if(pt<xbins_pt[i+1]) return i;
	}
	return nbins_pt-1;
}


void deriveDijetResponse(int startfile=0, int endfile=1){

	const int nbins_pt = 3;
	double xbins_pt[nbins_pt+1] = {60,110,200,1000};
	const int nbins_eta = 20;
	double xbins_eta[nbins_eta+1];
	for(int i=0; i<=nbins_eta; i++){
		xbins_eta[i] = -5+(10./(double)nbins_eta)*i;
	}
	
	TH1F *hRelResponse[nbins_pt];
	TH1F *hMPFResponse[nbins_pt];
	TH1F *avgAHisto[nbins_pt][nbins_eta];
	TH1F *avgBHisto[nbins_pt][nbins_eta];
	for(int i=0; i<nbins_pt; i++){
		hRelResponse[i] = new TH1F(Form("hRelResponse_pt%d",i),"",nbins_eta,xbins_eta);
		hMPFResponse[i] = new TH1F(Form("hMPFResponse%d",i),"",nbins_eta,xbins_eta);
		for(int j=0; j<nbins_eta; j++){
			avgAHisto[i][j] = new TH1F(Form("avgAHisto_pt%d_eta%d",i,j),"",30,-2,2);
			avgBHisto[i][j] = new TH1F(Form("avgBHisto_pt%d_eta%d",i,j),"",30,-2,2);
		}
	}

	double avgA[nbins_pt][nbins_eta];
	int nEntries[nbins_pt][nbins_eta];
	double avgB[nbins_pt][nbins_eta];
	int nEntriesB[nbins_pt][nbins_eta];
	for(int i=0; i<nbins_pt; i++){
		for(int j=0; j<nbins_eta; j++){
			avgA[i][j]=0;
			nEntries[i][j]=0;
			avgB[i][j]=0;
			nEntriesB[i][j]=0;
		}
	}

	std::string infile_Forest = "ppHighPtJet80_list.txt";
	std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
	std::string filename_Forest;
	TChain *jtTree = new TChain("ak4PFJetAnalyzer/t");
	TChain *pfTree = new TChain("pfcandAnalyzer/pfTree");

	for(int ifile = 0;ifile<startfile;ifile++){
		instr_Forest>>filename_Forest;
	}
	for(int ifile = startfile; ifile<endfile; ++ifile){
		
		instr_Forest>>filename_Forest;
		cout << "opening file " << filename_Forest << endl;
		TFile *fin = TFile::Open(filename_Forest.c_str());
		TTree *jtTree = (TTree*)fin->Get("ak4PFJetAnalyzer/t");
		TTree *pfTree = (TTree*)fin->Get("pfcandAnalyzer/pfTree");

		float pt_F[1000], eta_F[1000], phi_F[1000], rawpt_F[1000], m_F[1000];
		float pfPt[10000], pfEta[10000], pfPhi[10000];
		int nref, nPFpart, pfId[10000];

		jtTree->SetBranchAddress("jtpt",pt_F);
		jtTree->SetBranchAddress("jteta",eta_F);
		jtTree->SetBranchAddress("jtphi",phi_F);
		jtTree->SetBranchAddress("jtm",m_F);
		jtTree->SetBranchAddress("rawpt",rawpt_F);
		jtTree->SetBranchAddress("nref",&nref);
		pfTree->SetBranchAddress("nPFpart",&nPFpart);
		pfTree->SetBranchAddress("pfPt",pfPt);
		pfTree->SetBranchAddress("pfEta",pfEta);
		pfTree->SetBranchAddress("pfPhi",pfPhi);
		pfTree->SetBranchAddress("pfId",pfId);

		cout << "entries: "<< jtTree->GetEntries() << endl;
		int totalEntries = jtTree->GetEntries();

		for(int ientry=0; ientry<totalEntries; ientry++){
			jtTree->GetEntry(ientry);
			pfTree->GetEntry(ientry);

			if(ientry && ientry%1000000==0) cout << "entry: "<< ientry << endl;

		//Doing standard Rrel method
			if(nref>1){
				int rJet, pJet;
				if(abs(eta_F[0])<1.3 && abs(eta_F[1])<1.3){
					pJet = rand()%2;
					rJet = (pJet+1)%2;
				}
				else if(abs(eta_F[0])<1.3 && abs(eta_F[1])<5.){
					pJet = 1; rJet = 0;
				}
				else if(abs(eta_F[1])<1.3 && abs(eta_F[0])<5.){
					pJet = 0; rJet = 1;
				}
				else continue;
				double avgPt = 0.5*(pt_F[rJet]+pt_F[pJet]);
				if(nref>2){
					if(pt_F[2]/avgPt > 0.2){
						continue;
					}
				}
				int etaBin = findEtaBin(eta_F[pJet],nbins_eta,xbins_eta);
				int ptBin = findPtBin(pt_F[pJet],nbins_pt,xbins_pt);
	             //if(etaBin>xbins_eta[nbins_eta] || ptBin>xbins_pt[nbins_pt]) cout << "OH NO! Bin Mismatch!" << endl;
				avgA[ptBin][etaBin] += (pt_F[pJet]-pt_F[rJet])/avgPt;
				avgAHisto[ptBin][etaBin]->Fill((pt_F[pJet]-pt_F[rJet])/avgPt);
				nEntries[ptBin][etaBin]++;
			}

	    //Starting MPF Method
			if(nref>1){
				double avgPt = 0.5*(pt_F[0]+pt_F[1]);
				if(nref>2){
					if(pt_F[2]/avgPt > 0.2){
						continue;
					}
				}
	    	//start by summing the pf candidates
				TLorentzVector missEt(0,0,0,0);
				for(int i=0; i<nPFpart; i++){
					TLorentzVector pfp;
					double pfMass=0;
					if(abs(pfEta[i])<1.3){
						if(pfId[i]==1 || pfId[i]==4){ pfMass = 0.1395702; }
						pfp.SetPtEtaPhiE(pfPt[i],pfEta[i],pfPhi[i],pfMass);
						missEt += pfp;
					}
				}
	    	//include jet residuals
				for(int i=0; i<nref; i++){
					if(abs(eta_F[i])<1.3){
						TLorentzVector jt(pt_F[i]-rawpt_F[i],eta_F[i],phi_F[i],m_F[i]);
						missEt += jt;
					}
				}
	        //rotate to get *missing* ET from residual ET
				missEt*=-1;

				for(int i=0; i<nref; i++){
					int etaBin = findEtaBin(eta_F[i],nbins_eta,xbins_eta);
					int ptBin = findPtBin(pt_F[i],nbins_pt,xbins_pt);
					TLorentzVector jetVec(pt_F[i],eta_F[i],phi_F[i],m_F[i]);
					avgB[ptBin][etaBin] += jetVec.Dot(missEt) / ( 2*avgPt * jetVec.Mag());
					avgBHisto[ptBin][etaBin]->Fill(jetVec.Dot(missEt) / ( 2*avgPt * jetVec.Mag()));
					nEntriesB[ptBin][etaBin]++;
				}
			}
		} //close entries
	} //close file loop

	for(int i=0; i<nbins_pt; i++){
		for(int j=0; j<nbins_eta; j++){
			avgA[i][j] = avgA[i][j]/(double)nEntries[i][j];
			avgB[i][j] = avgB[i][j]/(double)nEntriesB[i][j];
			if(nEntries[i][j]) hRelResponse[i]->SetBinContent(j+1,(1+avgA[i][j])/(1-avgA[i][j]));
			if(nEntries[i][j]) hRelResponse[i]->SetBinError(j+1,hRelResponse[i]->GetBinContent(j+1)*(1./sqrt(nEntries[i][j])));
			else hRelResponse[i]->SetBinContent(j+1,0);

			if(nEntriesB[i][j]) hMPFResponse[i]->SetBinContent(j+1,(1+avgB[i][j])/(1-avgB[i][j]));
			if(nEntriesB[i][j]) hMPFResponse[i]->SetBinError(j+1,hMPFResponse[i]->GetBinContent(j+1)*(1./sqrt(nEntriesB[i][j])));
			else hMPFResponse[i]->SetBinContent(j+1,0);

			int totEntriesA=0, totEntriesB=0;
			for(int j=0; j<nbins_eta; j++){ 
				totEntriesA+=nEntries[i][j]; 
				totEntriesB+=nEntriesB[i][j]; 
			}
			hRelResponse[i]->SetEntries(totEntriesA);
			hMPFResponse[i]->SetEntries(totEntriesB);
		}
	}
	
	int color[4] = {1,2,4,8};
	TFile *fout = new TFile("relDijetResponse.root","recreate");
	fout->cd();
	for(int i=0; i<nbins_pt; i++){
		hRelResponse[i]->Write();
		hRelResponse[i]->SetMarkerColor(color[i]);
		hRelResponse[i]->SetLineColor(color[i]);
		hRelResponse[i]->Scale(1./hRelResponse[i]->Integral());

		hMPFResponse[i]->Write();
		hMPFResponse[i]->SetMarkerColor(color[i]);
		hMPFResponse[i]->SetLineColor(color[i]);
		hMPFResponse[i]->Scale(1./hMPFResponse[i]->Integral());
		for(int j=0; j<nbins_eta; j++){
			avgAHisto[i][j]->Write();
			avgBHisto[i][j]->Write();
		}
	}

	fout->Close();

}