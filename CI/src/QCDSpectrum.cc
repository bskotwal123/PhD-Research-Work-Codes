// ---------------------------------------------------------------------------
// Read a QCD histogram
// HBP 2014
// ---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "QCDSpectrum.h"
#include "TFile.h"
#include "TClass.h"
#include "TKey.h"
// ---------------------------------------------------------------------------
ClassImp(QCDSpectrum)
// ---------------------------------------------------------------------------
using namespace std;

namespace {
void error(string program, string message)
{
  cout << "** " << program << " ** " << message << endl;
  exit(0);
}
};

int QCDSpectrum::histid=0;
			 
QCDSpectrum::QCDSpectrum(const  char* _name, const char* _title,
			 string hdir, string hname, string fname,
			 string prefix)
  : RooAbsReal(_name, _title),
    histdir_(hdir),
    histname_(hname),
    filename_(fname),
    y(vector<double>()),
    pt(vector<double>()),
    xsec(0)
{
  get(hdir, hname, y);
  
  // create a histogram to receive spectrum
  string title(hdir+hname);
  if ( prefix != "" ) title = prefix + title;

  QCDSpectrum::histid++;
  char namen[80];
  sprintf(namen, "hQCD%4.4d", QCDSpectrum::histid);
  
  xsec = new TH1D(namen, "", pt.size() - 1, &pt[0]);
  xsec->SetTitle(title.c_str());
  xsec->GetXaxis()->SetTitle("Jet p_{T} (GeV)");
  for(int c=0; c < (int)pt.size()-1; c++) xsec->SetBinError(c+1, 0);
}

QCDSpectrum::QCDSpectrum(const QCDSpectrum& other, const  char* newname)
  : RooAbsReal(other, newname),
    histdir_(other.histdir_),
    histname_(other.histname_),
    filename_(other.filename_),
    y(other.y),
    pt(other.pt),
    xsec(other.xsec)
{
}

TH1D* QCDSpectrum::operator()()
{
  for(int c=0; c < (int)pt.size()-1; c++) xsec->SetBinContent(c+1, y[c]);
  return xsec;
}

double QCDSpectrum::operator()(int index) const
{
  if ( index < 0 ) return -1;
  if ( index > (int)pt.size()-1 ) return -2;
  return y[index];
}

double QCDSpectrum::evaluate() const { return 1; }

void QCDSpectrum::get(string hdir, 
		      string hname, 
		      vector<double>& h)
{
  // open root file
  char filename[80];
  sprintf(filename, "%s/%s", hdir.c_str(), filename_.c_str());
  TFile hfile(filename);
  if ( ! hfile.IsOpen() )
    error("QCDSpectrum", 
		 string("can't open ") + string(filename));

  // get histogram
  TH1D* hh = (TH1D*)hfile.Get(hname.c_str());
  if ( hh == 0 )
    error("QCDSpectrum", 
		 string("can't get histogram ") + hname);

  // fill internal buffer
  int nbins = hh->GetNbinsX();
  for(int ii=0; ii < nbins; ii++)
    {
      h.push_back(hh->GetBinContent(ii+1));
      pt.push_back(hh->GetBinLowEdge(ii+1));
    }
  pt.push_back(pt.back()+hh->GetBinWidth(nbins));
}

int QCDSpectrum::count(TFile* hfile, string prefix)
{
  hfile->cd();
  int n = 0;
  char name[80];
  sprintf(name, "%s%3.3d", prefix.c_str(), n);
  QCDSpectrum* o = (QCDSpectrum*)hfile->Get(name);
  n++;
  while ( o )
    {
      sprintf(name, "%s%3.3d", prefix.c_str(), n);
      o = (QCDSpectrum*)hfile->Get(name);
      if ( !o ) break;
      n++;
    }
  return n;
}
