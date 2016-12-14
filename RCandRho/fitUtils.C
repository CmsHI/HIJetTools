#include "TF1.h"
#include "TH1.h"

TF1 *fg1 = NULL;
TF1 *fitTruncatedGaus(TH1 *h, Double_t n = 3) {

  if(fg1) delete fg1;
  fg1 = new TF1("g1","gaus",h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
  h->Fit(fg1,"R0","",h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
  Double_t sigma = fg1->GetParameter(2);
  Double_t mean  = fg1->GetParameter(1);
  Printf("--1 fit mean: %f sigma: %f",mean,sigma);
  h->Fit(fg1,"R0","",mean-n*sigma,mean+n*sigma);

  sigma = fg1->GetParameter(2);
  mean  = fg1->GetParameter(1);
  Printf("--2 fit mean: %f sigma: %f",mean,sigma);
  h->Fit(fg1,"R0","",mean-n*sigma,mean+n*sigma);

  sigma = fg1->GetParameter(2);
  mean  = fg1->GetParameter(1);
  Printf("--3 fit mean: %f sigma: %f",mean,sigma);

  return fg1;
}

TF1 *fitLHSGaus(TH1 *h, double nsig = 3.,double nsig2 = 0.5) {

  int maxIter = 10;

  if(fg1) delete fg1;

  double meanold = h->GetMean();
  double mean = h->GetMean();
  double sigma = h->GetRMS();
  fg1 = new TF1(h->GetName(),"gaus",-nsig*h->GetRMS(),nsig2*h->GetRMS());
  for(int i = 0; i<maxIter; ++i) {
    h->Fit(fg1,"R0","",mean-nsig*sigma,mean+nsig2*sigma);
    meanold = mean;
    sigma = fg1->GetParameter(2);
    mean  = fg1->GetParameter(1);
    if(fabs(mean-meanold)<0.001) return fg1;
  }
  return fg1;
}
