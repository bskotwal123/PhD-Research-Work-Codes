// ---------------------------------------------------------------------------
// Read coefficients histograms c and compute CI cross
// section. HBP 2014
// ---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "CISpectrum.h"
#include "CIXsection.h"
#include "TFile.h"
#include "TClass.h"
#include "TKey.h"
// ---------------------------------------------------------------------------
ClassImp(CISpectrum)
// ---------------------------------------------------------------------------
using namespace std;

namespace {
void error(string program, string message)
{
  cout << "** " << program << " ** " << message << endl;
  exit(0);
}
  
  //#include "CIXsection.cc"
};

int CISpectrum::histid=0;

CISpectrum::CISpectrum(const char* _name, const char* _title,
		       string hdir, string hname, string fname, string prefix)
  : RooAbsReal(_name, _title),
    histdir_(hdir),
    histname_(hname),
    filename_(fname),

    bi(vector<vector<double> >()),
    aig(vector<vector<double> >()),
    ai(vector<vector<double> >()),

    bij(vector<vector<double> >()),
    aijg(vector<vector<double> >()),
    aij(vector<vector<double> >()),

    bi4(vector<vector<double> >()),
    ai4g(vector<vector<double> >()),
    ai4(vector<vector<double> >()),

    pt(vector<double>()),
    xsec(0),
    which(0)
{
  char name[2048];
  if ( fname == "" )
    {
      sprintf(name, "%s/ai1.root", hdir.c_str());
      init(name, hname);

      get(hdir, hname, "bi",  bi);
      get(hdir, hname, "aig", aig);
      get(hdir, hname, "ai",  ai);

      get(hdir, hname, "bij", bij);
      get(hdir, hname, "aijg",aijg);
      get(hdir, hname, "aij", aij);

      get(hdir, hname, "bi4", bi4);
      get(hdir, hname, "ai4g",ai4g);
      get(hdir, hname, "ai4", ai4);
    }
  else
    {
      sprintf(name, "%s/%s", hdir.c_str(), fname.c_str());
      init(string(name), hname);

      string number = hname.substr(1, hname.size()-1);
      sprintf(name, "b%s", number.c_str());
      get(hdir, string(name), fname, bi, -1.0/25);

      sprintf(name, "a%s", number.c_str());
      get(hdir, string(name), fname, bij, 1.0/625);
    }


  // create a histogram to receive spectrum
  string title(hdir+hname);
  if ( prefix != "" ) title = prefix + title;

  CISpectrum::histid++;
  char namen[80];
  sprintf(namen, "hCI%4.4d", CISpectrum::histid);
  
  xsec = new TH1D(namen, "", pt.size() - 1, &pt[0]);
  xsec->SetTitle(title.c_str());
  xsec->GetXaxis()->SetTitle("Jet p_{T} (GeV)");
  for(int c=0; c < (int)pt.size()-1; c++) xsec->SetBinError(c+1, 0);
}

CISpectrum::CISpectrum(const CISpectrum& other, const char* newname)
  : RooAbsReal(other, newname),
    histdir_(other.histdir_),
    histname_(other.histname_),
    filename_(other.filename_),

    bi(other.bi),
    aig(other.aig),
    ai(other.ai),

    bij(other.bij),
    aijg(other.aijg),
    aij(other.aij),

    bi4(other.bi4),
    ai4g(other.ai4g),
    ai4(other.ai4),

    pt(other.pt),
    xsec(other.xsec),
    which(other.which)
{
}

TH1D* CISpectrum::operator()(double lambda, vector<double>& kappa)
{
  for(int c=0; c < xsec->GetNbinsX(); c++)
    {
      vector<double> xsecs = CIXsection::evaluate(lambda, kappa,
						  bi[c],  aig[c],  ai[c], 
						  bij[c], aijg[c], aij[c],
						  bi4[c], ai4g[c], ai4[c]);
      if      ( which == 0 )
	xsec->SetBinContent(c+1, xsecs[0]+xsecs[1]);
      else if ( which == 1 )
	xsec->SetBinContent(c+1, xsecs[0]);
      else if ( which == 2 )
	xsec->SetBinContent(c+1, xsecs[1]);
      else if ( which == 3 )
	xsec->SetBinContent(c+1, xsecs[2]);
      else if ( which == 4 )
	xsec->SetBinContent(c+1, xsecs[3]);
      else if ( which == -1 )
	xsec->SetBinContent(c+1, xsecs[2]+xsecs[3]);
      else
	xsec->SetBinContent(c+1, xsecs[0]+xsecs[1]);
    }
  return xsec;
}

double CISpectrum::operator()(double lambda, vector<double>& kappa,
			      int c) const
{
  if ( c < 0 ) return -999;
  if ( c > xsec->GetNbinsX()-1 ) return -999;
  
  vector<double> xsecs = CIXsection::evaluate(lambda, kappa,
					      bi[c],  aig[c],  ai[c], 
					      bij[c], aijg[c], aij[c],
					      bi4[c], ai4g[c], ai4[c]);
  if ( which < 0 )
    return xsecs[2] + xsecs[3];
  else
    return xsecs[0] + xsecs[1];
}

double CISpectrum::evaluate() const { return 1; }

void CISpectrum::init(string rootfile, 
		      string hname)
{
  // open root file
  TFile hfile(rootfile.c_str());
  if ( ! hfile.IsOpen() )
    error("CISpectrum", string("can't open ") + rootfile);

  // get histogram
  TH1D* hh = (TH1D*)hfile.Get(hname.c_str());
  if ( hh == 0 )
    error("CISpectrum", 
		 string("can't get histogram ") + hname);
  int nbins = hh->GetNbinsX();

  for(int ii=0; ii < nbins; ii++)
    pt.push_back(hh->GetBinLowEdge(ii+1));
  pt.push_back(pt.back()+hh->GetBinWidth(nbins));

  // initialize internal buffers
  for(int ii=0; ii < nbins; ii++)
    {
      bi.push_back(vector<double>(6, 0));
      aig.push_back(vector<double>(6, 0));
      ai.push_back(vector<double>(6, 0));

      bij.push_back(vector<double>(9, 0));
      aijg.push_back(vector<double>(9, 0));
      aij.push_back(vector<double>(9, 0));

      bi4.push_back(vector<double>(4, 0));
      ai4g.push_back(vector<double>(4, 0));
      ai4.push_back(vector<double>(4, 0));
    }
}

void CISpectrum::get(string hdir, 
		     string hname, 
		     string fname, 
		     vector<vector<double> >& h,
		     double scale)
{
  int nh = h[0].size();
  bool is7TeV = fabs(scale) > 0;

  // for 7 TeV data, we have only one set of coefficients
  if ( is7TeV ) nh = 1;

  for(int c=0; c < nh; c++)
    { 
      // open root file
      char filename[2048];
      if ( is7TeV )
	sprintf(filename, "%s/%s", hdir.c_str(), fname.c_str());
      else
	sprintf(filename, "%s/%s%d.root", hdir.c_str(), fname.c_str(), c);

      TFile hfile(filename);
      if ( ! hfile.IsOpen() )
	error("CISpectrum", 
		     string("can't open ") + string(filename));

      // get histogram
      TH1D* hh = (TH1D*)hfile.Get(hname.c_str());
      if ( hh == 0 )
	error("CISpectrum", 
		     string("can't get histogram ") + hname);
      int nbins = hh->GetNbinsX();

      // fill internal buffer
      if ( is7TeV )
	for(int ii=0; ii < nbins; ii++)
	  h[ii][c] = hh->GetBinContent(ii+1) * scale;
      else
	for(int ii=0; ii < nbins; ii++)    
	  h[ii][c] = hh->GetBinContent(ii+1);
    }
}


int CISpectrum::count(TFile* hfile, string prefix)
{
  hfile->cd();
  int n = 0;
  char name[80];
  sprintf(name, "%s%3.3d", prefix.c_str(), n);
  CISpectrum* o = (CISpectrum*)hfile->Get(name);
  n++;
  while ( o )
    {
      sprintf(name, "%s%3.3d", prefix.c_str(), n);
      o = (CISpectrum*)hfile->Get(name);
      if ( !o ) break;
      n++;
    }
  return n;
}
