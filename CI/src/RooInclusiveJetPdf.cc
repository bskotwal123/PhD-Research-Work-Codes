//---------------------------------------------------------------------------
// File: RooInclusiveJetPdf.cc
// Description: compute inclusive jet pdf
// Created: 11-Nov-2014 Harrison B. Prosper
//---------------------------------------------------------------------------
#include <vector>
#include <iostream>
#include <limits>
#include "TMath.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooInclusiveJetPdf.h"
#include "TMath.h"
#include "TRandom3.h"
#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"
//---------------------------------------------------------------------------
ClassImp(RooInclusiveJetPdf)
//---------------------------------------------------------------------------
using namespace std;
namespace {
  ROOT::Math::Random<ROOT::Math::GSLRngMT>* gslrandom = 
    new ROOT::Math::Random<ROOT::Math::GSLRngMT>();
  TRandom3 rand3;
};

RooInclusiveJetPdf::RooInclusiveJetPdf(const char *name, const char *title,
				       RooArgSet&   _count,
				       RooAbsReal&  _lambda,
				       RooArgSet&   _kappa)
  : RooAbsPdf(name, title),
    count("count",   "N",       this, kTRUE, kFALSE), 
    lambda("lambda", "#lambda", this, _lambda),
    kappa("kappa",   "#kappa",  this, kTRUE, kFALSE),
    qcd(vector<QCDSpectrum>()),
    ci(vector<CISpectrum>()),
    smallest(numeric_limits<double>::denorm_min()),
    number(-1),
    useasimov(false),
    asimov(vector<double>(_count.getSize(), 0)),
    qcdxsect(vector<double>(_count.getSize(), 0)),
    xsection(vector<double>(_count.getSize(), 0)),
    firstbin(0),
    lastbin(_count.getSize()-1),
    useLO(false),
    maxsize(0),
    nbootstrap(1),
    useinterpolation(false),
    interp(0)
{  
  count.add(_count);
  kappa.add(_kappa);
}

RooInclusiveJetPdf::RooInclusiveJetPdf(const RooInclusiveJetPdf& other, 
				       const char* name) 
  : RooAbsPdf(other, name),
    count("count",   this, other.count), 
    lambda("lambda", this, other.lambda),
    kappa("kappa",   this, other.kappa),
    qcd(other.qcd),
    ci(other.ci),
    smallest(other.smallest),
    number(other.number),
    useasimov(other.useasimov),
    asimov(other.asimov),
    qcdxsect(other.qcdxsect),
    xsection(other.xsection),
    firstbin(other.firstbin),
    lastbin(other.lastbin),
    useLO(other.useLO),
    maxsize(other.maxsize),
    nbootstrap(other.nbootstrap),
    useinterpolation(other.useinterpolation),
    interp(other.interp)
{}

void RooInclusiveJetPdf::add(QCDSpectrum& _qcd, CISpectrum&  _ci)
{
  qcd.push_back(_qcd);
  ci.push_back(_ci);
  maxsize = qcd.size()-1;
}

QCDSpectrum* RooInclusiveJetPdf::QCD(int c)
{
  if ( c < 0 ) return 0;
  if ( c >= (int)qcd.size() ) return 0;
  return &(qcd[c]);
}

CISpectrum* RooInclusiveJetPdf::CI(int c)
{
  if ( c < 0 ) return 0;
  if ( c >= (int)ci.size() ) return 0;
  return &(ci[c]);
}

vector<double>& RooInclusiveJetPdf::crossSection(int c)
{
  for(size_t ii=0; ii < xsection.size(); ii++)  xsection[ii] = 0;
  if ( c < 0 ) return xsection;
  if ( c >= (int)ci.size() ) return xsection;

  // CI parameters
  double l = (double)lambda;
  std::vector<double> k(6);
  for(int ii=0; ii < 6; ii++)
    k[ii] = dynamic_cast<RooRealVar*>(&kappa[ii])->getVal();

  if ( useLO )
    ci[c].set(-1);
  else
    ci[c].set(0);
  for(size_t ii=0; ii < xsection.size(); ii++)
    {
      qcdxsect[ii] = qcd[c](ii);
      xsection[ii] = qcdxsect[ii] + ci[c](l, k, ii);
    }
  return xsection;
}

void RooInclusiveJetPdf::setNumber(int which)
{
  number = which;
  if ( number > (int)qcd.size()-1 ) number =-1;
}

void RooInclusiveJetPdf::setBinRange(int first, int last)
{
  if ( first < 0 ) first = 0;
  if ( first > count.getSize()-1 ) first = count.getSize()-1;
  if ( last  < 0 ) last  = count.getSize()-1;
  if ( last  > count.getSize()-1 ) last  = count.getSize()-1;
  if ( first > last ) first = last;
  firstbin = first;
  lastbin  = last;
  cout << "RooInclusiveJetPdf: bin range["
       << firstbin << "..." << lastbin << "]"
       << endl;
}

void RooInclusiveJetPdf::setAsimov(bool yes, double lumi, double l)
{
  useasimov = yes;
  if ( ! useasimov ) return;

  std::vector<double> k(6);
  if ( l > 0 )
    for(int c=0; c < 6; c++)
      k[c] = dynamic_cast<RooRealVar*>(&kappa[c])->getVal();
  
  // Note: the first element is presumed to be the
  // nominal cross section.
  double scale = 1000*lumi;
  for(size_t ii=0; ii < asimov.size(); ii++)
    {
      double xsec = qcd[0](ii);
      if ( l > 0 )
	{
	  xsec += ci[0](l, k, ii);
	}
      asimov[ii] = scale * xsec;
    }
}

void RooInclusiveJetPdf::setInterpolate(bool yes)
{
  useinterpolation = false;
  if ( ! yes ) return;
  
  int npts=50;
  vector<double> x(npts+1);
  vector<double> y(npts+1);
  double xmin = dynamic_cast<RooRealVar*>(&lambda)->getMin();
  double xmax = dynamic_cast<RooRealVar*>(&lambda)->getMax();
  double xstep= (xmax-xmin)/npts;
  for(int c=0; c <= npts; c++)
    {
      x[c] = xmin + c * xstep;
      dynamic_cast<RooRealVar*>(&lambda)->setVal(x[c]);
      y[c] = evaluate();
    }

  if ( interp ) delete interp;
  interp = new ROOT::Math::Interpolator(x, y);
  useinterpolation = true;
}

double RooInclusiveJetPdf::evaluate() const
{
  double l = (double)lambda;
  if ( useinterpolation ) return interp->Eval(l);
  
  // CI parameters
  std::vector<double> k(6);
  for(int c=0; c < 6; c++)
    k[c] = dynamic_cast<RooRealVar*>(&kappa[c])->getVal();
  
  // Counts
  long double y  = 0;
  
  // firstbin, lastbin determine the range of bins to use in
  // calculation of likelihood
  
  int nbins = lastbin-firstbin+1;
  vector<double> n(nbins, 0);
  vector<double> p(nbins, 0);
  
  if ( useasimov )
    {
      int jj=0;
      for(int ii=firstbin; ii <= lastbin; ii++)
	{
	  n[jj] = asimov[ii];
	  jj++;
	}
    }
  else
    {
      int jj = 0;
      for(int ii=firstbin; ii <= lastbin; ii++)
	{
	  n[jj] = dynamic_cast<RooRealVar*>(&count[ii])->getVal();
	  jj++;
	}
    }

  if ( number > -1 )
    {
      // if ( useLO )
      // 	ci[number].set(-1);
      // else
      // 	ci[number].set(0);
      
      int jj = 0;
      for(int ii=firstbin; ii <= lastbin; ii++)
	{
	  p[jj] = qcd[number](ii) + ci[number](l, k, ii);
	  jj++;
	}
      y = RooInclusiveJetPdf::multinomial(n, p);
    }
  else
    {
      // Note: the first element is presumed to be the
      // nominal cross section so we skip it.
      // nbootstrap is the number of bootstrap samples
      for(int B=0; B < nbootstrap; B++)
	{
	  for(int kk=1; kk < maxsize; kk++)
	    {
	      int c = kk;
	      if ( nbootstrap > 1 ) c = rand3.Integer(maxsize-2)+1;
	      
	      // if ( useLO )
	      // 	ci[c].set(-1);
	      // else
	      // 	ci[c].set(0);
      
	      int jj=0;
	      for(int ii=firstbin; ii <= lastbin; ii++)
		{
		  p[jj] = qcd[c](ii) + ci[c](l, k, ii);
		  jj++;
		}
	      y += RooInclusiveJetPdf::multinomial(n, p);
	    }
	}
      y /= (maxsize * nbootstrap);
    }
  if ( y != y || y <= 0 ) 
    return 0;
  else
    return (double)y;
}

long double RooInclusiveJetPdf::multinomial(vector<double>& n, 
					    vector<double>& p)
{
  double sum = 0;
  double total  = 0;
  for(size_t c=0; c < n.size(); c++)
    {
      sum += p[c];
      total += n[c];
    }
  long double zplus  = 0.0;
  long double zminus = 0.0;
  for(size_t c=0; c < p.size(); c++)
    {
      long double nlnp  = 0;
      if ( n[c] > 0 )
	{
	  long double x = p[c]/sum;
	  long double y = n[c]/total;
	  nlnp = n[c] * log(x/y);
	}
      if ( nlnp < 0 )
	zminus += -nlnp;
      else
	zplus  +=  nlnp;
    }
  long double z = zplus - zminus;
  long double q = exp(z);
  return q;
}


// long double RooInclusiveJetPdf::multinomial(vector<double>& n, 
// 				     vector<double>& p)
// {
//   double sum = 0;
//   double total = 0;
//   for(size_t c=0; c < n.size(); c++)
//     {
//       total += n[c];
//       sum   += p[c];
//     }
//   long double u = TMath::LnGamma(total+1);
//   if ( u != u ) return 0;

//   long double y = 0.0;
//   long double z = 0.0;
//   for(size_t i=0; i < p.size(); i++)
//     {
//       long double lng = TMath::LnGamma(n[i]+1);
//       if ( lng != lng ) return 0;
//       y += lng;

//       long double nlnp = n[i] * log(p[i]/sum);
//       if ( nlnp != nlnp ) return 0;
//       z += nlnp;
//     }
//   long double q = z + u - y;
//   return exp(q);
// }

int RooInclusiveJetPdf::getAnalyticalIntegral(RooArgSet& allVars,
					      RooArgSet& analVars, 
					      const char* /*rangeName*/)
  const 
{
  if (matchArgs(allVars, analVars, count)) 
    return 1;
  else
    return 0;
}

double RooInclusiveJetPdf::analyticalIntegral(int code, 
					      const char* /*rangeName*/)
  const 
{
  if(code==1)
    {
      // already normalized!
      return 1;
    } 
  else
    return 0; 
}
