//---------------------------------------------------------------------------
// File: RooMultinomial.cc
// Description:
// Compute Multinominal pdf (N choose n) theta^n (1 - theta)^(N-n)
// Created: 11-Mar-2013 Supriya Jain and Harrison B. Prosper
//          CERN
// Updated: 19-Mar-2013 HBP - fix bug in constructor. See how RooArgLists
//                      are handled in RooMultiVarGaussian
//          20-Mar-2014 HBP - use RooArgSets instead of lists
//          10-Jun-2014 HBP - fix bug in evaluate (need to update _total
//                      because counts may have changed)
//$Revision: 1.8 $
//---------------------------------------------------------------------------
#include <iostream>
#include <limits>
#include "TMath.h"
#include "RooFit.h"
#include "RooFormulaVar.h"
#include "RooMultinomial.h"
#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"
//---------------------------------------------------------------------------
ClassImp(RooMultinomial)

namespace {
ROOT::Math::Random<ROOT::Math::GSLRngMT>* gslrandom = 
  new ROOT::Math::Random<ROOT::Math::GSLRngMT>();
};
//---------------------------------------------------------------------------
using namespace std;

// RooMultinomial::RooMultinomial(const char *name, const char *title,
//   			       RooArgList& counts,  RooArgList& thetas)
//   : RooAbsPdf(name, title), 
//     _counts("counts", "observables", this, kTRUE, kFALSE),
//     _thetas("thetas", "parameters", this, kTRUE, kFALSE),
//     _nbins(counts.getSize()),
//     _smallest(numeric_limits<double>::denorm_min())
// {  
//   _counts.add(counts);
//   _thetas.add(thetas);
//   _total = 0.0;
//   string formula(counts[0].GetName());
//   for(int ibin=1; ibin < _nbins; ++ibin)
//     _total = dynamic_cast<RooRealVar*>(_counts.at(ibin))->getVal();
// }

RooMultinomial::RooMultinomial(const char* name, const char* title,
  			       RooArgSet& counts,  RooArgSet& thetas)
  : RooAbsPdf(name, title), 
    _counts("counts", "observables", this, kTRUE, kFALSE),
    _thetas("thetas", "parameters", this, kTRUE, kFALSE),
    _nbins(counts.getSize()),
    _smallest(numeric_limits<double>::denorm_min())
{  
  _counts.add(counts);
  _thetas.add(thetas);
  _total = 0.0;
  for(int ibin=0; ibin < _nbins; ++ibin)
    _total = dynamic_cast<RooRealVar*>(_counts.at(ibin))->getVal();
}

//---------------------------------------------------------------------------
RooMultinomial::RooMultinomial(const RooMultinomial& other, 
			       const char* name) 
  : RooAbsPdf(other, name), 
    _counts("counts",  this, other._counts),
    _thetas("thetas",  this, other._thetas),
    _total(other._total),
    _nbins(other._nbins),
    _smallest(other._smallest)
{
}
//---------------------------------------------------------------------------
void RooMultinomial::setTotal(int total)
{
  _total = total;
}

double RooMultinomial::evaluate() const
{
  std::vector<double> n(_nbins);
  std::vector<double> p(_nbins);
  // be sure to update total count
  double sum = 0.0;
  double total = 0.0;
  for(int ibin=0; ibin < _nbins; ++ibin)
    {
      p[ibin] = ((RooAbsReal*)_thetas.at(ibin))->getVal();  
      sum += p[ibin];
      n[ibin] = ((RooAbsReal*)_counts.at(ibin))->getVal();
      total += n[ibin];
    }
 
  // normalize thetas
  for(int ibin=0; ibin < _nbins; ++ibin) p[ibin] /= sum;

  // compute multinomial

  long double u = TMath::LnGamma(total+1);
  if ( u != u ) return _smallest;

  long double y = 0.0;
  long double z = 0.0;
  for(unsigned int i=0; i < p.size(); i++)
    {
      long double lng = TMath::LnGamma(n[i]+1);
      if ( lng != lng ) return _smallest;
      y += lng;
      long double nlnp = n[i] * log(p[i]);
      if ( nlnp != nlnp ) return _smallest;
      z += nlnp;
    }
  long double value = z + u - y;
  if ( value != value) return _smallest;
  
  value = exp(value);

  if (value != value) return _smallest;
  return (double)value; 
}

//_____________________________________________________________________________
Int_t RooMultinomial::getGenerator(const RooArgSet& directVars, 
				   RooArgSet& generateVars, Bool_t /*staticInitOK*/) const
{
  // Advertise internal generator
  if (matchArgs(directVars,generateVars,_counts)) 
    return 1;
  else
    return 0;
}

//_____________________________________________________________________________
void RooMultinomial::generateEvent(Int_t code)
{
  // Implement internal generator using gslrandom->Multinomial
  if (code==1)
    {
      std::vector<double> n(_nbins);
      std::vector<double> p(_nbins);
      double sum = 0.0;
      for(int ibin=0; ibin < _nbins; ++ibin)
      	{
      	  p[ibin] = ((RooAbsReal*)_thetas.at(ibin))->getVal();
      	  sum += p[ibin];
      	} 

      for(int ibin=0; ibin < _nbins; ++ibin) p[ibin] /= sum;

      assert(gslrandom != 0);
      vector<unsigned int> N = gslrandom->Multinomial(_total, p);
      
      // set RooFit parameters
      for(int ibin=0; ibin < _nbins; ++ibin)
      	((RooRealVar*)_counts.at(ibin))->setVal(N[ibin]);
    }
  return;
}

//_____________________________________________________________________________
Int_t RooMultinomial::getAnalyticalIntegral(RooArgSet& allVars,
					    RooArgSet& analVars, const char* /*rangeName*/) const 
{
  if (matchArgs(allVars,analVars,_counts)) 
    return 1;
  else
    return 0;
}


//_____________________________________________________________________________
Double_t RooMultinomial::analyticalIntegral(Int_t code, 
					    const char* /*rangeName*/) const 
{
  if(code==1)
    {
      // already normalized!
      return 1;
    } 
  else
    return 0; 
}
