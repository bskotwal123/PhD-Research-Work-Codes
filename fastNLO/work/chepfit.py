#!/usr/bin/env python
#-----------------------------------------------------------------------------
# File:        chebfit.py
# Description: Fit Chebyshev polynomials to m(gg) signal
# Created:     Fri Jul 21, 2012
# Author:      Harrison B. Prosper
# $Revision:$
#-----------------------------------------------------------------------------
import os, sys, re
from ROOT import *
# ---------------------------------------------------------------------------
N = 5
class Ansatz(ROOT.Math.IGenFunction):
    def __init__(self):
        pass
    
  Ansatz() {}
  double DoEval(double x) const { return ansatz(x); }
  IBaseFunctionOneDim* Clone() const
  {
    return new Ansatz();
  }
};
Ansatz g;

class Ansatz2 : public ROOT::Math::IGenFunction
{
public:
  Ansatz2() {}
  virtual ~Ansatz2() {}
  double DoEval(double x) const 
  { 
    return scale*ansatz(x)*pow(x/100, N); 
  }
  IBaseFunctionOneDim* Clone() const
  {
    return new Ansatz2();
  }
};
Ansatz2 g2;

double xmin=10;
double xmax=1000;
int nord = 136;
int nord2= 8;

ROOT::Math::Chebyshev cheb(g, xmin, xmax, nord);

ROOT::Math::Chebyshev cheb2(g2, xmin, xmax, nord2);


double func1(double* x, double* p)
{
  return ansatz(x[0]);
}

double func2(double* x, double* p)
{
  return cheb(x[0]);
}

double fun1(double* x, double* p)
{
  return scale*ansatz(x[0])*pow(x[0]/100, N); 
}

double fun2(double* x, double* p)
{
  return cheb2(x[0]);
}
#-----------------------------------------------------------------------------
void chebfit()
{
  gROOT->ProcessLine(".L histutil.C+");
  setStyle();

  TFIle* hfile = new TFile("spectrum_chebyshev.root", "recreate");

  TCanvas* cspec = kit::canvas("spectrum_chebyshev");

  TF1 f1("spectrum1", func1, xmin, xmax);
  f1.SetMinimum(1.e-1);
  f1.GetHistogram()->GetXaxis()->SetTitle("p_{T} [GeV/c^{2}]");
  f1.GetHistogram()->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb/GeV/c^{2}]");

  TF1 f2("spectrum2", func2, xmin, xmax);
  f2.SetLineColor(kBlue);

  TH1F* h = (TH1F*)f2.GetHistogram()->Clone();
  h->Divide(f1.GetHistogram());
  h->SetMinimum(0);
  h->SetMaximum(2);
  h->GetXaxis()->SetTitle("p_{T} [GeV/c^{2}]");
  h->GetYaxis()->SetTitle("Chebyshev(p_{T})/f(p_{T})");

  cspec->cd();
  cspec->SetLogy();
  cspec->SetLogx();
  f1.Draw();
  f2.Draw("same");
  cspec->Update();

  string ans;
  cout << "enter 1>>> ";
  cin >> ans;

  TF1 ff1("spect1", fun1, xmin, xmax);
  ff1.GetHistogram()->GetXaxis()->SetTitle("p_{T} [GeV/c^{2}]");
  ff1.GetHistogram()->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb/GeV/c^{2}]");

  TF1 ff2("spect2", fun2, xmin, xmax);
  ff2.SetLineColor(kBlue);

  cspec->cd();
  cspec->SetLogy(0);
  cspec->SetLogx(0);
  ff1.Draw();
  ff2.Draw("same");
  cspec->Update();

  cout << "enter 2>>> ";
  cin >> ans;

  TCanvas* cratio = canvas("spectrum_chebyshev_ratio", 10);
  cratio->cd();
  h->Draw();
  cratio->Update();

  hfile->cd();
  cspec->Write();
  f1.Write();
  f2.Write();

  cratio->Write();
  h->Write();
  hfile->close();
}
