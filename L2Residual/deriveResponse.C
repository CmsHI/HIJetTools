
/*****************
 * Original disciption
 * Macro for deriving the dijet relative and absolute response from data (or MC)
 * Uses JME-13-004 as a reference - "A" and "B" formulae are defined there
 * Kurt Jung, Raghav Kunnawalkam Elayavalli, December 2015
 
 
 * Functionality was extened in order to be able to run on all kind of datasets
 * Dropped L3 residual corrections
 * Different MPF method, alpha extrapolation
 * Stepan Obraztsov, December 2016
 *******************/

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include "TLorentzVector.h"
#include <TF1.h>

using namespace std;

int findAlphaBin(float alpha, int nbins_alpha, double* xbins_alpha){

  for(int i=0; i<nbins_alpha; i++){
    if(alpha<xbins_alpha[i+1]) return i;
  }
  return nbins_alpha-1;
}

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


TLorentzVector findMissEtJets(int nref, float* pt_F, float* rawpt_F, float* eta_F, float* phi_F, float* m_F){
double sum_ex = 0.0;
double sum_ey = 0.0;
  TLorentzVector missEtJets(0,0,0,0);

  for(int i=0; i<nref; i++){
  if (pt_F[i]<10) continue; 
  if (rawpt_F[i]<10 && rawpt_F[i]>5){
  sum_ex += rawpt_F[i]*cos(phi_F[i]);
  sum_ey += rawpt_F[i]*sin(phi_F[i]);  
  }
  //else {
  if (rawpt_F[i]>10){
	sum_ex += pt_F[i]*cos(phi_F[i]);
	sum_ey += pt_F[i]*sin(phi_F[i]);  
  }
  }
  double met = sqrt( sum_ex*sum_ex + sum_ey*sum_ey );
  missEtJets.SetPxPyPzE( -sum_ex, -sum_ey, 0.0, met);

  return missEtJets;
}

void deriveResponse(int startfile = 0,
			int endfile = 1,
			int mc = 0,
			int vec = 0,
			std::string infile_Forest = "pp_mc_forests.txt",
			std::string mode = "minbias",
			int radius = 3,
			std::string out = "test.root")
{
  TH1::SetDefaultSumw2(); 
  
  const double minJetPt = 10;

  //double xbins_pt[] = {55.0,75.0,95.0,120.0,150.0,300.0};
  double xbins_pt[] = {25.0,35.0,55.0,75.0,95.0,120.0,150.0,300.0}; 
  const int nbins_pt = sizeof(xbins_pt)/sizeof(double)-1;
  double xbins_eta[] = {0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 3.139, 5.191};
  const int nbins_eta = sizeof(xbins_eta)/sizeof(double)-1;
  double xbins_alpha[] = {0.1,0.15,0.2,0.25,0.3,0.35,0.4};
  const int nbins_alpha = sizeof(xbins_alpha)/sizeof(double)-1;


  TH1D *h_avgA[nbins_pt][nbins_eta];
  TH1D *h_avgB[nbins_pt][nbins_eta];

  TH1D *h_avgA_alpha[nbins_pt][nbins_eta][nbins_alpha];
  TH1D *h_avgB_alpha[nbins_pt][nbins_eta][nbins_alpha];

  
  for(int i=0; i<nbins_pt; i++){
     for(int j=0; j<nbins_eta; j++){
      h_avgA[i][j]=new TH1D(Form("h_avgA_%d_%d",i,j),"",200, -2, 2);;
      h_avgB[i][j]=new TH1D(Form("h_avgB_%d_%d",i,j),"",200, -1, 3);;
    }
  }
  
  
  for(int i=0; i<nbins_pt; i++){
     for(int j=0; j<nbins_eta; j++){
	for(int k=0; k<nbins_alpha; k++){
      h_avgA_alpha[i][j][k]=new TH1D(Form("h_avgA_alpha_%d_%d_%d",i,j,k),"",200, -2, 2);;
      h_avgB_alpha[i][j][k]=new TH1D(Form("h_avgB_alpha_%d_%d_%d",i,j,k),"",200, -1, 3);;
    }
  }
}
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;

  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  int totalEntries = 0;
  for(int ifile = startfile; ifile<endfile; ++ifile){
    	
    instr_Forest>>filename_Forest;
    cout << "opening file " << filename_Forest << endl;
    TFile *fin = TFile::Open(filename_Forest.c_str());
    TTree *jtTree = (TTree*)fin->Get(Form("ak%dPFJetAnalyzer/t",radius));
    TTree *pfTree = (TTree*)fin->Get("pfcandAnalyzer/pfTree");
    TTree *phoTree = (TTree*)fin->Get("ggHiNtuplizer/EventTree");
    TTree *hltTree = (TTree*)fin->Get("hltanalysis/HltTree");
    TTree *skimTree = (TTree*)fin->Get("skimanalysis/HltTree");

    float pt_F[1000], eta_F[1000], phi_F[1000], rawpt_F[1000], m_F[1000], eSum[1000], phoSum[1000];

    vector<float> *phoEt=0, *phoEta=0, *phoPhi=0; 
    int nref, nPFpart, nPho;
    
    int pfId[10000];
    float pfPt[10000], pfEta[10000], pfPhi[10000], pfEnergy[10000];
    vector<float> *pfPt_vec=0, *pfEta_vec=0, *pfPhi_vec=0, *pfEnergy_vec=0;
    vector<int> *pfId_vec=0;    
    if (vec == 0){

    pfTree->SetBranchAddress("pfPt",pfPt);
    pfTree->SetBranchAddress("pfEta",pfEta);
    pfTree->SetBranchAddress("pfPhi",pfPhi);
    pfTree->SetBranchAddress("pfId",pfId);
    pfTree->SetBranchAddress("pfEnergy",pfEnergy); 
    			 }
    else {
    pfTree->SetBranchAddress("pfPt",&pfPt_vec);
    pfTree->SetBranchAddress("pfEta",&pfEta_vec);
    pfTree->SetBranchAddress("pfPhi",&pfPhi_vec);
    pfTree->SetBranchAddress("pfId",&pfId_vec);
    pfTree->SetBranchAddress("pfEnergy",&pfEnergy_vec);
    	 }
    int HLT_AK4PFJet40_Eta5p1_v1, HLT_AK4PFJet60_Eta5p1_v1, HLT_AK4PFJet80_Eta5p1_v1;
    int HBHENoiseFilterResultRun2Loose, pPAprimaryVertexFilter, pBeamScrapingFilter;
	float pthat;
	if (mc == 1){
    jtTree->SetBranchAddress("pthat",&pthat);
    			}
    jtTree->SetBranchAddress("jtpt",pt_F);
    jtTree->SetBranchAddress("jteta",eta_F);
    jtTree->SetBranchAddress("jtphi",phi_F);
    jtTree->SetBranchAddress("jtm",m_F);
    jtTree->SetBranchAddress("rawpt",rawpt_F);
    jtTree->SetBranchAddress("nref",&nref);
    jtTree->SetBranchAddress("eSum",&eSum);
    jtTree->SetBranchAddress("photonSum",&phoSum);
    pfTree->SetBranchAddress("nPFpart",&nPFpart);
    phoTree->SetBranchAddress("nPho",&nPho);
    phoTree->SetBranchAddress("phoEt",&phoEt);
    phoTree->SetBranchAddress("phoEta",&phoEta);
    phoTree->SetBranchAddress("phoPhi",&phoPhi);
    if (mc == 0){
    hltTree->SetBranchAddress("HLT_AK4PFJet40_Eta5p1_v1", &HLT_AK4PFJet40_Eta5p1_v1);
    hltTree->SetBranchAddress("HLT_AK4PFJet60_Eta5p1_v1", &HLT_AK4PFJet60_Eta5p1_v1);
    hltTree->SetBranchAddress("HLT_AK4PFJet80_Eta5p1_v1", &HLT_AK4PFJet80_Eta5p1_v1);
				}
    skimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose);
    skimTree->SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter);
    skimTree->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter);

    cout << "file entries: "<< jtTree->GetEntries() << endl;

    //for(int ientry=0; ientry<jtTree->GetEntries(); ientry++){
    for(int ientry=0; ientry<100; ientry++){
      jtTree->GetEntry(ientry);
      pfTree->GetEntry(ientry);
      phoTree->GetEntry(ientry);
      hltTree->GetEntry(ientry);
      totalEntries++;
      if(nPFpart > 10000 || nref > 1000) cout << " warning! nPF: "<< nPFpart << " and njet: "<< nref << endl;
      if(totalEntries && totalEntries%100000==0) cout << "entry: "<< ientry << endl;

    if(nref>1){
	int rJet, pJet;
	if(abs(eta_F[0])>5.191 || abs(eta_F[1])>5.191) continue;
	if(abs(eta_F[0])<1.3 && abs(eta_F[1])<1.3){
	  pJet = rand()%2;
	  rJet = (pJet+1)%2;
	}
	else if(abs(eta_F[0])<1.3 && abs(eta_F[1])>=1.3){
	  pJet = 1; rJet = 0;
	}
	else if(abs(eta_F[1])<1.3 && abs(eta_F[0])>=1.3){
	  pJet = 0; rJet = 1;
	}
	else continue;
	double avgPt = 0.5*(pt_F[rJet]+pt_F[pJet]);
	

	if(pt_F[rJet]<minJetPt) continue;
	double dphi = abs(phi_F[rJet]-phi_F[pJet]);
	if(dphi>(2*3.14159)) dphi-=(2*3.14159);
	if (dphi>3.14159) dphi=2*3.14159-dphi;
	if(dphi < 2.7) continue;
	int etaBin = findEtaBin(eta_F[pJet],nbins_eta,xbins_eta);
	int ptBin = findPtBin(avgPt,nbins_pt,xbins_pt);
	if (mc == 0){
	if (mode=="jet80" && HLT_AK4PFJet80_Eta5p1_v1 == 0) continue;
	if (mode=="lower" && HLT_AK4PFJet60_Eta5p1_v1 == 0 && HLT_AK4PFJet40_Eta5p1_v1 == 0) continue;
	if (mode=="jet80" && HLT_AK4PFJet80_Eta5p1_v1 == 1 && avgPt<95) continue;
	if (mode=="lower" && HLT_AK4PFJet60_Eta5p1_v1 == 1 && (avgPt<75 || avgPt>=95)) continue;
	if (mode=="lower" && HLT_AK4PFJet40_Eta5p1_v1 == 1 && (avgPt<55 || avgPt>=75)) continue;
	if (mode=="minbias" && (avgPt<25 || avgPt>=55)) continue;
				}
	else {
    if (mode=="pthat30" && (avgPt>=55 || avgPt<25 )) continue;
    if (mode=="pthat50" && (avgPt>=75 || avgPt<55 )) continue;
    if (mode=="pthat80" && (avgPt>=120 || avgPt<75 )) continue;
    if (mode=="pthat120" && avgPt<120 ) continue;    
		  }
    double w;
    if (mode=="herwig100") w = 0.000000052;
    else if (mode=="herwig200") w = 0.000000002;
    else w = 1.0;
	//Starting MPF Method
	TLorentzVector missEtJets = findMissEtJets(nref, pt_F, rawpt_F, eta_F, phi_F, m_F);	
	double dphi2 = abs(missEtJets.Phi()-phi_F[rJet]);
	if(dphi2>(2*3.14159)) dphi2-=(2*3.14159);
	if (dphi2>3.14159) dphi2=2*3.14159-dphi2;
	if(nref>2){
		for (int alphaBin = 0; alphaBin<nbins_alpha; alphaBin++){
			if (pt_F[2]/avgPt<xbins_alpha[alphaBin+1]){
			
			h_avgA_alpha[ptBin][etaBin][alphaBin]->Fill(0.5*(pt_F[pJet]-pt_F[rJet])/avgPt,w);
			
			h_avgB_alpha[ptBin][etaBin][alphaBin]->Fill((1 + missEtJets.Et()*cos(dphi2)/pt_F[rJet]),w);
				}
			}
		}
	if(nref>2) {
	  if(pt_F[2]/avgPt > 0.2){
	    continue;
	 		 				 }
			   }
	h_avgA[ptBin][etaBin]->Fill(0.5*(pt_F[pJet]-pt_F[rJet])/avgPt,w);
	h_avgB[ptBin][etaBin]->Fill((1 + missEtJets.Et()*cos(dphi2)/pt_F[rJet]),w);
      }
    } //close entries
    fin->Close();
	}
  TFile *fout;
  fout = new TFile(Form("%s",out.c_str()),"recreate");
  fout->cd();
  
  for(int i=0; i<nbins_pt; i++){

    for(int j=0; j<nbins_eta; j++){
      h_avgA[i][j]->Write();
      h_avgB[i][j]->Write();
    }
  }
  for(int i=0; i<nbins_pt; i++){
     for(int j=0; j<nbins_eta; j++){
	for(int k=0; k<nbins_alpha; k++){
      h_avgA_alpha[i][j][k]->Write();
      h_avgB_alpha[i][j][k]->Write();
	}
	}
	}


  fout->Close();

}
int main(int argc, char *argv[])
{ 
  deriveResponse(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),argv[5],argv[6]);
  return 0;
}
