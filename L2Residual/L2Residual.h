#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TF1.h"

using namespace std;
class L2Residual
{
 private:
  int neta;
  double etacut;
  double lower_pt_cut;
  double higher_pt_cut;
  int radius;
  TFile *correction_file;
  TF1 *relfit[100];
  TF1 *mpffit[100];
  
  
  TF1 *relfit_log[100];
  TF1 *mpffit_log[100];
  
  TF1 *relfit_alpha_avg;
  TF1 *mpffit_alpha_avg;    
  //TF1 *relfit_alpha[100];
  //TF1 *mpffit_alpha[100];
  double bins_eta[100];
  int nbins_eta;

  public:
  
   void reset()
   { 
    for(int ieta=0;ieta<100;ieta++){
    relfit[ieta] = NULL;
    }
    correction_file = NULL;
   } 
  
  L2Residual(int radius=3)
  {
  reset();


  double xbins_eta[] = {0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 3.139, 5.191};

  nbins_eta = sizeof(xbins_eta)/sizeof(double)-1;

  for(int ieta = 0; ieta < nbins_eta; ieta++){
	 bins_eta[ieta] = xbins_eta[ieta];
	}

  lower_pt_cut = 30;
  higher_pt_cut = 300;
  etacut = 5.191;
   
  //correction_file = new TFile("5tev_alpha_fits_2505_5.root");
  correction_file = new TFile("5tev_alpha_fits_2207.root");  

  for(int i=0; i<nbins_eta; i++){	
	relfit[i] = (TF1*)correction_file->Get(Form("relfit%d",i));
	mpffit[i] = (TF1*)correction_file->Get(Form("mpffit%d",i));
	
	relfit_log[i] = (TF1*)correction_file->Get(Form("relfit_log%d",i));
	mpffit_log[i] = (TF1*)correction_file->Get(Form("mpffit_log%d",i));	
	
	//relfit_alpha[i] = (TF1*)correction_file->Get(Form("relfit_alpha%d",i));
	//mpffit_alpha[i] = (TF1*)correction_file->Get(Form("mpffit_alpha%d",i));
				} 
				
	relfit_alpha_avg = (TF1*)correction_file->Get("relfit_alpha_avg");				
	mpffit_alpha_avg = (TF1*)correction_file->Get("mpffit_alpha_avg");	
  }
  

  double relcor(double jetpt, double jeteta, int isConst, int isFSR)
  {
   double correction = 1;
   if( abs(jeteta)> ((double)etacut)) return jetpt*correction;
   int etaindex = 0;
   for(int i=0; i<nbins_eta; i++){
    if(abs(jeteta)<bins_eta[i]) continue;
    else etaindex = i;
				 				 }
  // if(jetpt < lower_pt_cut) return jetpt*relfit[etaindex]->Eval(lower_pt_cut);
  // if(jetpt > higher_pt_cut) return jetpt*relfit[etaindex]->Eval(higher_pt_cut);

   if(jetpt < lower_pt_cut) jetpt = lower_pt_cut;
   if(jetpt > higher_pt_cut) jetpt = higher_pt_cut;
   
  if (isConst == 1) correction *= relfit[etaindex]->Eval(jetpt);
  else correction *= relfit_log[etaindex]->Eval(jetpt);
  
  if (isFSR == 1) correction *= relfit_alpha_avg->Eval(bins_eta[etaindex]);
  
  return jetpt*correction;
  }


  double mpfcor(double jetpt, double jeteta, int isConst, int isFSR)
  {
   double correction = 1;
   if( abs(jeteta)> ((double)etacut)) return jetpt*correction;
   int etaindex = 0;
   for(int i=0; i<nbins_eta; i++){
    if(abs(jeteta)<bins_eta[i]) continue;
    else etaindex = i;
								 }
   //if(jetpt < lower_pt_cut) return jetpt*mpffit[etaindex]->Eval(lower_pt_cut);
   //if(jetpt > higher_pt_cut) return jetpt*mpffit[etaindex]->Eval(higher_pt_cut);


   //return jetpt*mpffit[etaindex]->Eval(jetpt);
   
   //if(jetpt < lower_pt_cut) jetpt = lower_pt_cut;
   //if(jetpt > higher_pt_cut) jetpt = higher_pt_cut;
   
  if (isConst == 1) correction *= mpffit[etaindex]->Eval(jetpt);
  else correction *= mpffit_log[etaindex]->Eval(jetpt);
  
  if (isFSR == 1) correction *= mpffit_alpha_avg->Eval(bins_eta[etaindex]);
  
  //return jetpt*correction;
  return jetpt*1.05;
  }

  double relfactor(double jetpt, double jeteta)
  {
   double correction = 1;
   if( abs(jeteta)> ((double)etacut)) return correction;
   int etaindex = 0;
   for(int i=0; i<nbins_eta; i++){
    if(abs(jeteta)<bins_eta[i]) continue;
    else etaindex = i;
				 }
   if(jetpt < lower_pt_cut) return relfit[etaindex]->Eval(lower_pt_cut);
   if(jetpt > higher_pt_cut) return relfit[etaindex]->Eval(higher_pt_cut);

   return relfit[etaindex]->Eval(jetpt);
  }

double mpffactor(double jetpt, double jeteta)
  {
   double correction = 1;
   if( abs(jeteta)> ((double)etacut)) return correction;
   int etaindex = 0;
   for(int i=0; i<nbins_eta; i++){
    if(abs(jeteta)<bins_eta[i]) continue;
    else etaindex = i;
				 }
   if(jetpt < lower_pt_cut) return mpffit[etaindex]->Eval(lower_pt_cut);
   if(jetpt > higher_pt_cut) return mpffit[etaindex]->Eval(higher_pt_cut);

   return mpffit[etaindex]->Eval(jetpt);
  }

};
