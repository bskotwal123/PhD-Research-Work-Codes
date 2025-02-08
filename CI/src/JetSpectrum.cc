// ----------------------------------------------------------------------------
// file: JetSpectrum.cc
// Read a spectrum from a root file and created an interpolated function from
// it. 
// HBP 2012 - 2014
// Updated: 11-October-2014 HBP - add non-perturbative correction
// ----------------------------------------------------------------------------
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>

#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"

#include "Math/Interpolator.h"
#include "Math/WrappedFunction.h"
#include "Math/Integrator.h"

#include "JetSpectrum.h"
#include "hutil.h"
// ----------------------------------------------------------------------------
using namespace std;
using namespace ROOT::Math;

const double ABSTOL=1.e-9;
const double RELTOL=1.e-6;

// ----------------------------------------------------------------------------
// Inclusive jet spectrum in |y| < 0.5
// ----------------------------------------------------------------------------
JetSpectrum::JetSpectrum(string filename,
		   string histname,
		   bool positive_,
		   bool applyNPC_,
		   TH1F* hEWKC_)
  : ptlo(std::vector<double>()),
    pthi(std::vector<double>()),
    ptcn(std::vector<double>()),
    xsection(std::vector<double>()),
    _positive(positive_),
    applyNPC(applyNPC_),
    hEWKC(hEWKC_),
    nullHist(false)
{
  //cout << "==> reading: " << filename << endl;
  TFile hfile(filename.c_str());
  if ( ! hfile.IsOpen() )
    hutil::error("JetSpectrum", "can't open file " + filename);

  TH1D* h = (TH1D*)hfile.Get(histname.c_str());
  if ( h == 0 )
    hutil::error("JetSpectrum", "can't find histogram " + histname);
  
  ptlo = hutil::binlowedges(h);
  pthi = hutil::binhighedges(h);
  ptcn = hutil::bincenters(h);
  xsection = hutil::contents(h);

  // check for a null histogram
  nullHist = h->Integral() == 0;

  // apply non-perturbative correction
  if ( nullHist )
    interp = 0;
  else
    {
      if ( applyNPC )
	{
	  for(unsigned int c=0; c < ptcn.size(); c++)
	    xsection[c] = NPC( ptcn[c] ) * xsection[c];
	}
      if ( applyEWKC )
	{
	  for(unsigned int c=0; c < ptcn.size(); c++)
	    xsection[c] = EWKC( ptcn[c] ) * xsection[c];
	}
      if ( _positive )
	for(unsigned int c=0; c < ptcn.size(); c++)
	  xsection[c] = log(xsection[c]);
      interp = new Interpolator(ptcn, xsection, Interpolation::kLINEAR);
    }
}

double JetSpectrum::NPC(double pt)
{
  // from Sanmay Ganguly (Nov. 2014)
  double A =  1.003; // +/-  0.27   (old = 1.003)
  double B = 77.374; // +/-295.0    (old = 46.774)
  double n =  1.385; // +/-  8.5    (old = 1.253)
  return A + B/pow(pt, n);
}

double JetSpectrum::EWKC(double pt)
{
  if ( hEWKC )
    return hEWKC->Interpolate(pt);
  else
    return 1;
}


double JetSpectrum::operator()(double pt)
{
  if ( pt < ptcn.front() ) return 0;
  if ( pt > ptcn.back() )  return 0;  
  if ( nullHist ) return 0;
  if ( _positive )
    return exp(interp->Eval(pt));
  else
    return interp->Eval(pt);
}

double JetSpectrum::operator()(double pTlow, double pThigh)
{
  if ( pTlow < ptcn.front() ) return 0;
  if ( pThigh > ptcn.back() )  return 0;
  if ( nullHist ) return 0;
  ROOT::Math::WrappedMemFunction<JetSpectrum, double (JetSpectrum::*)(double)>
    fn(*this, &JetSpectrum::operator());
  ROOT::Math::Integrator ifn(fn);
  ifn.SetAbsTolerance(ABSTOL);
  ifn.SetRelTolerance(RELTOL);
  return ifn.Integral(pTlow, pThigh);
}
