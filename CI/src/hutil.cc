// ---------------------------------------------------------------------------
// file: hutil.cc
// HBP 2012 - 2014
// Updated: 11-October-2014 HBP - add non-perturbative correction
// ---------------------------------------------------------------------------
#include <iostream>
#include "TFile.h"
#include "TClass.h"
#include "TList.h"
#include "TKey.h"
#include "TH1D.h"
#include "hutil.h"
// ---------------------------------------------------------------------------
using namespace std;
vector<double> hutil::contents(TH1* hist)
{
  vector<double> c(hist->GetNbinsX());
  for(int i=0; i < hist->GetNbinsX(); i++)
      c[i] = hist->GetBinContent(i+1);
  return c;
}

vector<double> hutil::binlowedges(TH1* hist)
{
  vector<double> c(hist->GetNbinsX());
  for(int i=0; i < hist->GetNbinsX(); i++)
    c[i] = hist->GetBinLowEdge(i+1);
    return c;
}

vector<double> hutil::binhighedges(TH1* hist)
{
  vector<double> c(hist->GetNbinsX());
  for(int i=0; i < hist->GetNbinsX(); i++)
    c[i] = hist->GetBinLowEdge(i+1) + hist->GetBinWidth(i+1);
  return c;
}

vector<double> hutil::bincenters(TH1* hist)
{
  vector<double> c(hist->GetNbinsX());
  for(int i=0; i < hist->GetNbinsX(); i++)
    c[i] = hist->GetBinLowEdge(i+1) + 0.5*hist->GetBinWidth(i+1);
  return c;
}

vector<double> hutil::errors(TH1* hist)
{
  vector<double> c(hist->GetNbinsX());
  for(int i=0; i < hist->GetNbinsX(); i++)
    c[i] = hist->GetBinError(i+1);
  return c;
}

vector<double> hutil::cdf(TH1* hist)
{
  vector<double> c(hist->GetNbinsX());
  c[0] = hist->GetBinContent(1);
  for(int i=1; i < hist->GetNbinsX(); i++)
    c[i] = c[i-1]+hist->GetBinContent(i+1);
  return c;
}

void hutil::setContents(TH1* hist, vector<double>& c)
{
  int n = hist->GetNbinsX();
  int nbin = n < (int)c.size() ? n : c.size();
  for(int i=0; i < nbin; i++) hist->SetBinContent(i+1, c[i]);
}

void hutil::setErrors(TH1* hist, vector<double>& err)
{
  int n = hist->GetNbinsX();
  int nbin = n < (int)err.size() ? n : err.size();
  for(int i=0; i < nbin; i++) hist->SetBinError(i+1, err[i]);
}

void hutil::divideByWidth(TH1* hist)
{
  for(int i=0; i < hist->GetNbinsX(); i++)
    {
      double w = hist->GetBinWidth(i+1);
      double c = hist->GetBinContent(i+1);
      double e = hist->GetBinError(i+1);
      hist->SetBinContent(i+1, c/w);
      hist->SetBinError(i+1, e/w);
    }
}

vector<TH1*> hutil::histograms(TFile* hfile, string)
{
  vector<TH1*> hist;
  TIter nextkey(hfile->GetListOfKeys());
  hfile->cd();
  while ( TKey* key = (TKey*)nextkey() )
    {
      TObject* o = key->ReadObj();
      if ( o->InheritsFrom("TH1") )
          hist.push_back((TH1*)o);
    }
  return hist;
}

vector<string> hutil::histogramNames(TFile* hfile)
{
  vector<string> hist;
  TIter nextkey(hfile->GetListOfKeys());
  hfile->cd();
  while ( TKey* key = (TKey*)nextkey() )
    {
      TObject* o = key->ReadObj();
      if ( o->InheritsFrom("TH1") )
	hist.push_back(((TH1*)o)->GetName());
    }
  return hist;
}

void hutil::error(string program, string message)
{
  cout << "** " << program << " ** " << message << endl;
  exit(0);
}
