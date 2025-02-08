// ---------------------------------------------------------------------------
// file: JetSpectrumSmeared.cc
// apply jet response function.
// HBP 2012 - 2014
// Updated: 11-October-2014 HBP - add non-perturbative correction
// ---------------------------------------------------------------------------
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

#include "JetSpectrumSmeared.h"
#include "JetSpectrum.h"
#include "JetCorrectionUncertainty.h"
// ---------------------------------------------------------------------------
using namespace std;
using namespace ROOT::Math;

const double ABSTOL=1.e-8;
const double RELTOL=1.e-5;
// ---------------------------------------------------------------------------
// Inclusive jet spectrum in |y| < 0.5
//----------------------------------------------------------------------------
JetSpectrumSmeared::JetSpectrumSmeared(JetSpectrum* spectrum_,
				       JetCorrectionUncertainty* JESunc_,
				       double JERunc_,
				       double x_, double y_,
				       double pTmin_, double pTmax_,
				       int npT_)
  : spectrum(spectrum_), 
    JESunc(JESunc_),
    JERunc(JERunc_),
    x(x_), y(y_),
    pTmin(pTmin_), pTmax(pTmax_), npT(npT_),
    pT(std::vector<double>(npT+1,0)),
    xsection(std::vector<double>(npT+1,0)),
    interp(0)
{
  if ( spectrum->null() )
    interp = 0;
  else
    {
      double pTstep = (pTmax-pTmin)/npT;
      for(int c=0; c <= npT; c++)
	{
	  pT[c] = pTmin + c*pTstep;
	  if ( JESunc )
	    xsection[c] = applySmearing_(pT[c]);
	  else
	    xsection[c] = (*spectrum)(pT[c]);
	  if ( spectrum->positive() )
	    xsection[c] = log(xsection[c]);
	  
	}
      interp = new Interpolator(pT, xsection, Interpolation::kLINEAR);
    }
}

JetSpectrumSmeared::JetSpectrumSmeared(const JetSpectrumSmeared& o)
  : spectrum(o.spectrum), 
    JESunc(o.JESunc),
    JERunc(o.JERunc),
    x(o.x),
    y(o.y),
    pTmin(o.pTmin),
    pTmax(o.pTmax),
    npT(o.npT),
    pT(o.pT),
    xsection(o.xsection),
    interp(o.interp)
{
}

double JetSpectrumSmeared::operator()(double pt)
{
  if ( pt < pT.front() ) return 0;
  if ( pt > pT.back() )  return 0;
  if ( spectrum->null() ) return 0;
  if ( spectrum->positive() )
    return exp(interp->Eval(pt));
  else
    return interp->Eval(pt);
}

double JetSpectrumSmeared::operator()(double pTlow, double pThigh)
{
  if ( pTlow  < pT.front() ) return 0;
  if ( pThigh > pT.back() )  return 0;
  if ( spectrum->null() ) return 0;
  ROOT::Math::WrappedMemFunction<JetSpectrumSmeared, 
				 double (JetSpectrumSmeared::*)(double)>
    fn(*this, &JetSpectrumSmeared::operator());
  ROOT::Math::Integrator ifn(fn);
  ifn.SetAbsTolerance(ABSTOL);
  ifn.SetRelTolerance(RELTOL);  
  return ifn.Integral(pTlow, pThigh);
}

// Convolution of response function with NLO spectrum;
double JetSpectrumSmeared::applySmearing_(double pTreco)
{
  pTreco_ = pTreco; // NB: cache reco-level pT
  if ( spectrum->null() ) return 0;

  double offset = 8 * sigmapT(pTreco);
  double ptmin = TMath::Max(20.0, pTreco - offset);
  double ptmax = TMath::Min(3000.0, pTreco + offset);

  ROOT::Math::WrappedMemFunction<JetSpectrumSmeared, 
				 double (JetSpectrumSmeared::*)(double)>
    fn(*this, &JetSpectrumSmeared::integrand_);
  ROOT::Math::Integrator ifn(fn);
  ifn.SetAbsTolerance(ABSTOL);
  ifn.SetRelTolerance(RELTOL);  
  return ifn.Integral(ptmin, ptmax);
}
 
double JetSpectrumSmeared::integrand_(double pT)
{
  return response(pTreco_, pT) * (*spectrum)(pT);
}

// Jet response function
double JetSpectrumSmeared::response(double pTreco, double pT)
{
  JESunc->setJetPt(pTreco);
  JESunc->setJetEta(0.0);
  double X = TMath::Max(1.e-3, 1.0 + x * JESunc->getUncertainty(false));
  double Y = TMath::Max(1.e-3, 1.0 + y * JERunc);
  return TMath::Gaus(pTreco/X, pT, Y*sigmapT(pT), kTRUE); 
}

double JetSpectrumSmeared::sigmapT(double pT)
{    
  // August 15, 2014, Sanmay Ganguly
  // https://indico.cern.ch/event/335191/contribution/1/material/slides/0.pdf
  double C = 0.029; //0.031;
  double S = 0.984; //0.949;
  double N = 5.7936;//6.130;
  double CdataMC = 1.052;// 1.12;
  return pT*CdataMC*sqrt(N*N/(pT*pT) + S*S/pT + C*C); 
}
