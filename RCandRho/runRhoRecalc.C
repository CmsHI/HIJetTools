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

void runRhoRecalc(TString str = "root://eoscms//eos/cms/store/cmst3/group/hintt/mverweij/CS/MC/velicanu-Hydjet_Quenched_MinBias_5020GeV_758p2_RECODEBUG_v0_forest_csjet_v1/0.root", Int_t nRCPerEvent = 1, int maxEvents = -1, int useRhoCorr = 0) {

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
  Int_t CaloJet100TriggerBit = 1;
  Int_t HBHENoiseFilterResult = 1;
  Int_t pprimaryVertexFilter = 1;
  Int_t pcollisionEventSelection = 1;
 
  Int_t           hiBin;
  float           hiHF;
  TBranch        *b_hiBin;   //!
  TBranch        *b_hiHF;   //!
  fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
  fChain->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
//  fChain->SetBranchAddress("HLT_HIL1MinimumBiasHF1ANDExpress_v1",&MinBiasTriggerBit);i
  fChain->SetBranchAddress("HLT_HIL1MinimumBiasHF1AND_v1",&MinBiasTriggerBit);
  fChain->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1",&CaloJet100TriggerBit);
  fChain->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
  fChain->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);
  fChain->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
  fChain->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  
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
  if(useRhoCorr) {
  fChain->SetBranchAddress("rhoCorr", &fFJRho.rho, &fFJRho.b_rho);
  fChain->SetBranchAddress("rhomCorr", &fFJRho.rhom, &fFJRho.b_rhom);
  } else {
  fChain->SetBranchAddress("rho", &fFJRho.rho, &fFJRho.b_rho);
  fChain->SetBranchAddress("rhom", &fFJRho.rhom, &fFJRho.b_rhom);
  }
  fChain->SetBranchAddress("ptJets", &fFJRho.ptJets, &fFJRho.b_ptJets);
  fChain->SetBranchAddress("etaJets", &fFJRho.etaJets, &fFJRho.b_etaJets);
  fChain->SetBranchAddress("areaJets", &fFJRho.areaJets, &fFJRho.b_areaJets);
  //Printf("nentries: %d",(Int_t)fChain->GetEntries());

  Int_t fgkNCentBins = 100;
  Float_t kMinCent   = 0.;
  Float_t kMaxCent   = 100;
  Double_t *binsCent = new Double_t[fgkNCentBins+1];
  for(Int_t i=0; i<=fgkNCentBins; i++) binsCent[i]=(Double_t)kMinCent + (kMaxCent-kMinCent)/fgkNCentBins*(Double_t)i ;

  Int_t fgkNMultBins = 300;
  Float_t kMinMult   = 0.;
  Float_t kMaxMult   = 6000;
  Double_t *binsMult = new Double_t[fgkNMultBins+1];
  for(Int_t i=0; i<=fgkNMultBins; i++) binsMult[i]=(Double_t)kMinMult + (kMaxMult-kMinMult)/fgkNMultBins*(Double_t)i ;

  Int_t fgkNEHFBins = 300;
  Float_t kMinEHF   = 0.;
  Float_t kMaxEHF   = 6000;
  Double_t *binsEHF = new Double_t[fgkNEHFBins+1];
  for(Int_t i=0; i<=fgkNEHFBins; i++) binsEHF[i]=(Double_t)kMinEHF + (kMaxEHF-kMinEHF)/fgkNEHFBins*(Double_t)i ;

  Int_t fgkNRhoBins = 500;
  Float_t kMinRho   = 0.;
  Float_t kMaxRho   = 500;
  Double_t *binsRho = new Double_t[fgkNRhoBins+1];
  for(Int_t i=0; i<=fgkNRhoBins; i++) binsRho[i]=(Double_t)kMinRho + (kMaxRho-kMinRho)/fgkNRhoBins*(Double_t)i ;

  Int_t fgkNRhoMBins = 500;
  Float_t kMinRhoM   = 0.;
  Float_t kMaxRhoM   = 5;
  Double_t *binsRhoM = new Double_t[fgkNRhoMBins+1];
  for(Int_t i=0; i<=fgkNRhoMBins; i++) binsRhoM[i]=(Double_t)kMinRhoM + (kMaxRhoM-kMinRhoM)/fgkNRhoMBins*(Double_t)i ;

  //WARNING: eta edges hard-coded. If we change this in the rho producer will not correspond anymore to the eta map
  std::map<int,double> fMapEtaRanges;   //eta ranges
  fMapEtaRanges[1] = -5.;
  fMapEtaRanges[2] = -3.;
  fMapEtaRanges[3] = -2.;//-2.1;
  fMapEtaRanges[4] = -1.5;//-1.3;
  fMapEtaRanges[5] = -1.;// 1.3;
  fMapEtaRanges[6] =  1.;//2.1;
  fMapEtaRanges[7] =  1.5;
  fMapEtaRanges[8] =  2.;
  fMapEtaRanges[9] =  3.;
  fMapEtaRanges[10] =  5.;

  Int_t fgkNEtaBins = (Int_t)fMapEtaRanges.size()-1;
  Printf("fgkNEtaBins: %d",fgkNEtaBins);
  Double_t *binsEta = new Double_t[fgkNEtaBins+1];
  for(Int_t i=0; i<=fgkNEtaBins; i++)
    binsEta[i]=fMapEtaRanges.at(i+1);

  // return;
  
  TList *fOutput =  new TList();
  TH1::SetDefaultSumw2();

  TH1F *fh1EventSel = new TH1F("fh1EventSel","fh1EventSel;sel step",10,0.,10.);
  fOutput->Add(fh1EventSel);
  
  TH3F *fh3RhoCentEtaBin = new TH3F("fh3RhoCentEtaBin","fh3RhoCentEtaBin;centrality;#rho;#eta",fgkNCentBins,binsCent,fgkNRhoBins,binsRho,fgkNEtaBins,binsEta);
  fOutput->Add(fh3RhoCentEtaBin);

  TH3F *fh3RhoRecalcCentEtaBin = new TH3F("fh3RhoRecalcCentEtaBin","fh3RhoRecalcCentEtaBin;centrality;#rho;#eta",fgkNCentBins,binsCent,fgkNRhoBins,binsRho,fgkNEtaBins,binsEta);
  fOutput->Add(fh3RhoRecalcCentEtaBin);


  TH3F *fh3RhoMCentEtaBin = new TH3F("fh3RhoMCentEtaBin","fh3RhoMCentEtaBin;centrality;#rho;#eta",fgkNCentBins,binsCent,fgkNRhoMBins,binsRhoM,fgkNEtaBins,binsEta);
  fOutput->Add(fh3RhoMCentEtaBin);

  TH3F *fh3RhoEHFEtaBin = new TH3F("fh3RhoEHFEtaBin","fh3RhoEHFEtaBin;sum E_{HF};#rho;#eta",fgkNEHFBins,binsEHF,fgkNRhoBins,binsRho,fgkNEtaBins,binsEta);
  fOutput->Add(fh3RhoEHFEtaBin);

  TH3F *fh3RhoMEHFEtaBin = new TH3F("fh3RhoMEHFEtaBin","fh3RhoMEHFEtaBin;sum E_{HF};#rho;#eta",fgkNEHFBins,binsEHF,fgkNRhoMBins,binsRhoM,fgkNEtaBins,binsEta);
  fOutput->Add(fh3RhoMEHFEtaBin);

  delete [] binsCent;
  delete [] binsRho;
  delete [] binsRhoM;
  delete [] binsEta;

  std::vector<float> rhoRecalc;

  Int_t startEntry = 0;
  Int_t lastEntry = fChain->GetEntries();//100;
  Printf("events in chain: %d",lastEntry);
  if(maxEvents>0 && maxEvents<lastEntry)
    lastEntry = maxEvents;
  Printf("lastEntry: %d",lastEntry);

  for (int j=startEntry; j<lastEntry; j++) {
    rhoRecalc.clear();
    fChain->GetEntry(j);
    if(j%1000==0)
      std::cout << "entry: "<< j << std::endl;

    double selStep = 0.5;
    fh1EventSel->Fill(selStep);
    if(!MinBiasTriggerBit) continue;
    selStep+=1.;
    fh1EventSel->Fill(selStep);
    //if(!CaloJet100TriggerBit) continue;
    if(!phfCoincFilter) continue;
    selStep+=1.;
    fh1EventSel->Fill(selStep);
    if(!HBHENoiseFilterResult) continue;
    selStep+=1.;
    fh1EventSel->Fill(selStep);
    if(!pprimaryVertexFilter) continue;
    selStep+=1.;
    fh1EventSel->Fill(selStep);
    if(!pcollisionEventSelection) continue;
    selStep+=1.;
    fh1EventSel->Fill(selStep);
    
    Double_t cent = (Double_t)hiBin/2.;

    //store rho and rhom in eta slices used in derivation
    for(unsigned int ieta = 0; ieta<fFJRho.etaMin->size(); ++ieta) {
      fh3RhoCentEtaBin->Fill(cent,fFJRho.rho->at(ieta),fFJRho.etaMin->at(ieta) + 0.5*(fFJRho.etaMax->at(ieta)-fFJRho.etaMin->at(ieta)));
      fh3RhoMCentEtaBin->Fill(cent,fFJRho.rhom->at(ieta),fFJRho.etaMin->at(ieta) + 0.5*(fFJRho.etaMax->at(ieta)-fFJRho.etaMin->at(ieta)));

      fh3RhoEHFEtaBin->Fill(hiHF,fFJRho.rho->at(ieta),fFJRho.etaMin->at(ieta) + 0.5*(fFJRho.etaMax->at(ieta)-fFJRho.etaMin->at(ieta)));
      fh3RhoMEHFEtaBin->Fill(hiHF,fFJRho.rhom->at(ieta),fFJRho.etaMin->at(ieta) + 0.5*(fFJRho.etaMax->at(ieta)-fFJRho.etaMin->at(ieta)));
    }

    //Recalculate rho
    int nExcl_ = 2;
    double ptMinExcl_ = 20.;
    double etaMaxExcl_ = 2.;
    int nExcl2_ = 1;
    double ptMinExcl2_ = 20.;
    double etaMaxExcl2_ = 3.;
    
    static double rhoVec[999];
    static double rhomVec[999];
    static double etaVec[999];
    int nacc = 0;
    unsigned int njetsEx = 0;
    unsigned int njetsEx2 = 0;
    for(unsigned int ij = 0; ij<fFJRho.ptJets->size(); ++ij) {

      // //excluce leading kt clusters from rho calculation
      // if(njetsEx<nExcl_ && fabs(fFJRho.etaJets->at(ij))<etaMaxExcl_ && fFJRho.ptJets->at(ij)>ptMinExcl_) {
      //   njetsEx++;
      //   continue;
      // }
      // if(njetsEx2<nExcl2_ && fabs(fFJRho.etaJets->at(ij))<etaMaxExcl2_ && fabs(fFJRho.etaJets->at(ij))>etaMaxExcl_ && fFJRho.ptJets->at(ij)>ptMinExcl2_) {
      //   njetsEx2++;
      //   continue;
      // }

      if(fFJRho.areaJets->at(ij)>0.) {
        rhoVec[nacc] = fFJRho.ptJets->at(ij)/fFJRho.areaJets->at(ij);
        etaVec[nacc] = fFJRho.etaJets->at(ij);
        ++nacc;
      }
    }

    double radius = 0.2; //distance kt clusters needs to be from edge
    for(unsigned int ieta = 0; ieta<fFJRho.etaMin->size(); ++ieta) {
      static double rhoVecCur[999] = {0.};
      static double rhomVecCur[999]= {0.};

      double etaMin = fFJRho.etaMin->at(ieta)+radius;
      double etaMax = fFJRho.etaMax->at(ieta)-radius;
/*
      if(etaMax<1.) {
        etaMin = -1.3+radius;
        etaMax =  1.3-radius;
      }     
 */
      int    naccCur    = 0 ;
      double rhoCurSum  = 0.;
      double rhomCurSum = 0.;
      for(int i = 0; i<nacc; i++) {
        if(etaVec[i]>=etaMin && etaVec[i]<etaMax) {
          rhoVecCur[naccCur] = rhoVec[i];
          rhoCurSum += rhoVec[i];
          ++naccCur;
        }//eta selection
      }//accepted jet loop
      
      if(naccCur>0) {
        double rhoCur = TMath::Median(naccCur, rhoVecCur);
        rhoRecalc.push_back(rhoCur);
//        fh3RhoRecalcCentEtaBin->Fill(cent,rhoCur,etaMin + 0.5*(etaMax - etaMin));
        fh3RhoRecalcCentEtaBin->Fill(cent,rhoCur,fFJRho.etaMin->at(ieta) + 0.5*(fFJRho.etaMax->at(ieta) - fFJRho.etaMin->at(ieta)));
        //        Printf("ieta: %d eta:%f-%f rhoRecalc: %f  rhoOrig: %f",ieta,etaMin,etaMax,rhoCur,fFJRho.rho->at(ieta));
      }
    }//eta ranges

  }


  TFile *fout = new TFile("RhoRecalc.root","RECREATE");
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
