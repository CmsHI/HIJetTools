#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH3D.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TVector2.h"
#include "TRandom3.h"

#include <iostream>

#include "ForestPFsVector.h"
#include "ForestFJRho.h"

Double_t deltaR(const Double_t phi1, const Double_t phi2, const Double_t eta1, const Double_t eta2);

void runPFRandomConeRho(TString str = "root://eoscms//eos/cms/store/cmst3/group/hintt/mverweij/CS/MC/velicanu-Hydjet_Quenched_MinBias_5020GeV_758p2_RECODEBUG_v0_forest_csjet_v1/0.root", Int_t nRCPerEvent = 1, int maxEvents = -1) {

  double minEta = -1.5;
  double maxEta = 1.5;
  double minPhi = -TMath::Pi();
  double maxPhi = TMath::Pi();
  double radiusRC = 0.4;
  
  TChain *fChain = new TChain("hiEvtAnalyzer/HiTree");
  TChain *pfTree = new TChain("pfcandAnalyzer/pfTree");
  TChain *rhoTree = new TChain("hiFJRhoAnalyzer/t");
  if(!rhoTree) Printf("Couldn't find rhoTree");
  TChain *skimTree = new TChain("skimanalysis/HltTree");
  TChain *hltTree = new TChain("hltanalysis/HltTree");
  //  TFile *f = TFile::Open(str.Data());
  fChain->Add(str.Data());
  pfTree->Add(str.Data());
  rhoTree->Add(str.Data());
  skimTree->Add(str.Data());
  hltTree->Add(str.Data());
  fChain->AddFriend(pfTree);
  fChain->AddFriend(rhoTree);
  fChain->AddFriend(skimTree);
  fChain->AddFriend(hltTree);
  
  if(!fChain) {
    Printf("Couldn't find fChain. Aborting!");
    return;
  }

  Int_t MinBiasTriggerBit = 1;
  Int_t phfCoincFilter = 1;
  
  Int_t           hiBin;
  TBranch        *b_hiBin;   //!
  fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
//  fChain->SetBranchAddress("HLT_HIL1MinimumBiasHF1ANDExpress_v1",&MinBiasTriggerBit);i
  fChain->SetBranchAddress("HLT_HIL1MinimumBiasHF1AND_v1",&MinBiasTriggerBit);
  fChain->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
  
  ForestPFsVector                    fPFs;              //!PFs in tree
  if (fChain->GetBranch("nPFpart"))
    fChain->SetBranchAddress("nPFpart", &fPFs.nPFpart, &fPFs.b_nPFpart);
  if (fChain->GetBranch("pfId"))
    fChain->SetBranchAddress("pfId", &fPFs.pfId, &fPFs.b_pfId);
  if (fChain->GetBranch("pfPt"))
    fChain->SetBranchAddress("pfPt", &fPFs.pfPt, &fPFs.b_pfPt);
  if (fChain->GetBranch("pfEta"))
    fChain->SetBranchAddress("pfEta", &fPFs.pfEta, &fPFs.b_pfEta);
  if (fChain->GetBranch("pfPhi"))
    fChain->SetBranchAddress("pfPhi", &fPFs.pfPhi, &fPFs.b_pfPhi);

  ForestFJRho fFJRho; //!rho in tree
  fChain->SetBranchStatus("eta*",1);
  if (fChain->GetBranch("etaMin"))
    fChain->SetBranchAddress("etaMin", &fFJRho.etaMin, &fFJRho.b_etaMin);
  else
    Printf("etaMin not found in fChain");
  fChain->SetBranchAddress("etaMax", &fFJRho.etaMax, &fFJRho.b_etaMax);
  fChain->SetBranchStatus("rho*",1);
  fChain->SetBranchAddress("rho", &fFJRho.rho, &fFJRho.b_rho);
  fChain->SetBranchAddress("rhom", &fFJRho.rhom, &fFJRho.b_rhom);
  
  //Printf("nentries: %d",(Int_t)fChain->GetEntries());

  TList *fOutput =  new TList();
  TH1::SetDefaultSumw2();
  
  TH3D *h2CentPtRCEta = new TH3D("h2CentPtRCEta","h2CentPtRCEta;centrality;p_{T,RC};#eta",100,0,100,300,-100.,200.,60,-6,6);
  fOutput->Add(h2CentPtRCEta);
  TH3D *h2CentPtRCEtaVS = new TH3D("h2CentPtRCEtaVS","h2CentPtRCEtaVS;centrality;p_{T,RC}-#rho A;#eta",100,0,100,300,-100.,200.,60,-6,6);
  fOutput->Add(h2CentPtRCEtaVS);

  TH3D *h2MultPtRCEta = new TH3D("h2MultPtRCEta","h2MultPtRCEta;centrality;p_{T,RC};#eta",3000,0,6000,200,-10.,100.,60,-6,6);
  fOutput->Add(h2MultPtRCEta);
  TH3D *h2MultPtRCEtaVS = new TH3D("h2MultPtRCEtaVS","h2MultPtRCEtaVS;centrality;p_{T,RC}-#rho A;#eta",3000,0,6000,200,-100.,100.,60,-6,6);
  fOutput->Add(h2MultPtRCEtaVS);
  Printf("histos defined");

  TH2F *fh2RhoCent = new TH2F("fh2RhoCent","fh2RhoCent;centrality;#rho",100,0,100,500,0,500);
  fOutput->Add(fh2RhoCent);

  TH2F *fh2RhoMCent = new TH2F("fh2RhoMCent","fh2RhoMCent;centrality;#rho_{m}",100,0,100,500,0,5);
  fOutput->Add(fh2RhoMCent);

  Int_t startEntry = 0;
  Int_t lastEntry = fChain->GetEntries();//100;
  Printf("events in chain: %d",lastEntry);
  if(maxEvents>0 && maxEvents<lastEntry)
    lastEntry = maxEvents;
  Printf("lastEntry: %d",lastEntry);

  TRandom3 *rnd = new TRandom3();
 
  for (int j=startEntry; j<lastEntry; j++) {
    fChain->GetEntry(j);
    if(j%10000==0)
      std::cout << "entry: "<< j << std::endl;

    if(!MinBiasTriggerBit) continue;
    if(!phfCoincFilter) continue;
    
    double etaRC = rnd->Rndm() * (maxEta - minEta) + minEta;
    double phiRC = rnd->Rndm() * (maxPhi - minPhi) + minPhi;
    
    float rhoCur = -1.;
    float rhomCur = -1.;
    for(unsigned int ieta = 0; ieta<fFJRho.etaMin->size(); ++ieta) {
      if(etaRC>=fFJRho.etaMin->at(ieta) && etaRC<fFJRho.etaMax->at(ieta)) {
        rhoCur = fFJRho.rho->at(ieta);
        rhomCur = fFJRho.rhom->at(ieta);
      }
    }

    double ptRC = 0.;
    Int_t pfCount = 0;
    for(Int_t i = 0; i<fPFs.nPFpart; i++) {
      double pt = fPFs.pfPt->at(i);
      double phi = fPFs.pfPhi->at(i);
      double eta = fPFs.pfEta->at(i);

      double dr = deltaR(phi,phiRC,eta,etaRC);
      if(dr<radiusRC) {
        ptRC+=pt;
      }
      if(std::abs(eta)<2.) pfCount++;
    }

    double ptRCSub = ptRC - rhoCur*TMath::Pi()*radiusRC*radiusRC;
    
    Double_t cent = (Double_t)hiBin/2.;
    h2CentPtRCEta->Fill(cent,ptRC,etaRC);
    h2CentPtRCEtaVS->Fill(cent,ptRCSub,etaRC);

    h2MultPtRCEta->Fill(pfCount,ptRC,etaRC);
    h2MultPtRCEtaVS->Fill(pfCount,ptRCSub,etaRC);

    fh2RhoCent->Fill(cent,rhoCur);
    fh2RhoMCent->Fill(cent,rhomCur);
  }

  if(rnd) delete rnd;

  TFile *fout = new TFile("RandomConesPF.root","RECREATE");
  fOutput->Write();
  fout->Write();
  fout->Close();
}

Double_t deltaR(const Double_t phi1, const Double_t phi2, const Double_t eta1, const Double_t eta2) {
  //calculate distance
  Double_t dPhi = phi1 - phi2;
  Double_t dEta = eta1 - eta2;
  dPhi = TVector2::Phi_mpi_pi(dPhi);
  return TMath::Sqrt(dPhi * dPhi + dEta * dEta);
}
