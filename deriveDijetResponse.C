#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>

using namespace std;

int findEtaBin(float eta, int nbins_eta, double* xbins_eta){

	for(int i=0; i<nbins_eta; i++){
		if(abs(eta)<xbins_eta[i+1]) return i;
	}
	return nbins_eta-1;
}

int findPtBin(float pt, int nbins_pt, double* xbins_pt){

	for(int i=0; i<nbins_pt; i++){
		if(pt<xbins_pt[i+1]) return i;
	}
	return nbins_pt-1;
}


void deriveDijetResponse(int startfile=0, int endfile=5){

	const double eta_cut_min = -3.;
	const double eta_cut_max = 3.;

	const int nbins_pt = 4;
	double xbins_pt[nbins_pt+1] = {30,90,180,360,600};
	const int nbins_eta = 20;
	double xbins_eta[nbins_eta+1];
	for(int i=0; i<=nbins_eta; i++){
		xbins_eta[i] = (5./(double)nbins_eta)*i;
	}
	
	TH1F *hRelResponse[nbins_pt];
	TH1F *avgAHisto[nbins_pt][nbins_eta];
	for(int i=0; i<nbins_pt; i++){
		hRelResponse[i] = new TH1F(Form("hRelResponse_pt%d",i),"",nbins_eta,xbins_eta);
		for(int j=0; j<nbins_eta; j++){
			avgAHisto[i][j] = new TH1F(Form("avgAHisto_pt%d_eta%d",i,j),"",30,-2,2);
		}
	}

	double avgA[nbins_pt][nbins_eta];
	int nEntries[nbins_pt][nbins_eta];
	for(int i=0; i<nbins_pt; i++){
		for(int j=0; j<nbins_eta; j++){
			avgA[i][j]=0;
			nEntries[i][j]=0;
		}
	}

	std::string infile_Forest = "ppHighPtJet80_list.txt";
	std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
	std::string filename_Forest;
	TChain *jtTree = new TChain("ak4PFJetAnalyzer/t");

	for(int ifile = 0;ifile<startfile;ifile++){
		instr_Forest>>filename_Forest;
	}
	for(int ifile = startfile; ifile<endfile; ++ifile){
		instr_Forest>>filename_Forest;
		cout << "opening file " << filename_Forest << endl;
		jtTree->Add(filename_Forest.c_str());
	}

	float pt_F[1000], eta_F[1000], phi_F[1000], rawpt_F[1000];
	int nref;

	jtTree->SetBranchAddress("jtpt",pt_F);
	jtTree->SetBranchAddress("jteta",eta_F);
	jtTree->SetBranchAddress("jtphi",phi_F);
	jtTree->SetBranchAddress("rawpt",rawpt_F);
	jtTree->SetBranchAddress("nref",&nref);

	cout << "entries: "<< jtTree->GetEntries() << endl;
	int totalEntries = jtTree->GetEntries();
	//first need to calculate all the <A> values in pt,eta bins
	for(int ientry=0; ientry<totalEntries; ientry++){
		jtTree->GetEntry(ientry);

		if(nref>2){
			int rJet, pJet;
			if(abs(eta_F[0])<1.3 && abs(eta_F[1])<1.3){
				pJet = rand()%2;
				rJet = (pJet+1)%2;
			}
			else if(abs(eta_F[0])<1.3 && abs(eta_F[1])<3.){
				pJet = 1; rJet = 0;
			}
			else if(abs(eta_F[1])<1.3 && abs(eta_F[0])<3.){
				pJet = 0; rJet = 1;
			}
			else continue;
			double avgPt = 0.5*(pt_F[rJet]+pt_F[pJet]);
			if(pt_F[2]/avgPt<0.2){
	             int etaBin = findEtaBin(eta_F[pJet],nbins_eta,xbins_eta);
	             int ptBin = findPtBin(pt_F[pJet],nbins_pt,xbins_pt);
	             //if(etaBin>xbins_eta[nbins_eta] || ptBin>xbins_pt[nbins_pt]) cout << "OH NO! Bin Mismatch!" << endl;
	             avgA[ptBin][etaBin] += (pt_F[pJet]-pt_F[rJet])/avgPt;
	             avgAHisto[ptBin][etaBin]->Fill((pt_F[pJet]-pt_F[rJet])/avgPt);
	             nEntries[ptBin][etaBin]++;
	        }
	    }
	}
	cout << "filled <A>" << endl;

	for(int i=0; i<nbins_pt; i++){
		for(int j=0; j<nbins_eta; j++){
			if(nEntries[i][j]){
				avgA[i][j] = avgA[i][j]/(double)nEntries[i][j];
				hRelResponse[i]->SetBinContent(j+1,(1+avgA[i][j])/(1-avgA[i][j]));
				hRelResponse[i]->SetBinError(j+1,hRelResponse[i]->GetBinContent(j+1)*(1./sqrt(nEntries[i][j])));
				int totEntries=0;
				for(int j=0; j<nbins_eta; j++){ totEntries+=nEntries[i][j]; }
					hRelResponse[i]->SetEntries(totEntries);
			}
			else{
				hRelResponse[i]->SetBinContent(j+1,0);
			}
		}
	}

	int color[4] = {1,2,4,8};
	TFile *fout = new TFile("relDijetResponse_5Files.root","recreate");
	fout->cd();
	for(int i=0; i<nbins_pt; i++){
		hRelResponse[i]->Write();
		hRelResponse[i]->SetMarkerColor(color[i]);
		hRelResponse[i]->SetLineColor(color[i]);
		hRelResponse[i]->Scale(1./hRelResponse[i]->Integral());
		for(int j=0; j<nbins_eta; j++){
			avgAHisto[i][j]->Write();
		}
	}

	fout->Close();

}
