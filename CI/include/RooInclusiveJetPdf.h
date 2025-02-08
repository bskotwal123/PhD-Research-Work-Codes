#ifndef ROOINCLUSIVEJETPDF_H
#define ROOINCLUSIVEJETPDF_H
//---------------------------------------------------------------------------
// Description: model pdf of inclusive jet spectrum averaged
//              over nuisance parameters
// Created: 11-Nov-2014 Harrison B. Prosper
//          09-Feb-2015 HBP add setNumber to allowing picking individual
//                      likelihoods from the ensemble of likelihoods.
//---------------------------------------------------------------------------
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "QCDSpectrum.h"
#include "CISpectrum.h"
#include "Math/Interpolator.h"
//---------------------------------------------------------------------------
class RooInclusiveJetPdf : public RooAbsPdf
{
public:
  RooInclusiveJetPdf() {}
  RooInclusiveJetPdf(const char *name, const char *title,
		     RooArgSet&   _count,
		     RooAbsReal&  _lambda,
		     RooArgSet&   _kappa);

  RooInclusiveJetPdf(const RooInclusiveJetPdf& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const 
  { return new RooInclusiveJetPdf(*this, newname); }
  
  void add(QCDSpectrum& _qcd, CISpectrum& _ci);

  virtual ~RooInclusiveJetPdf() {}

  int getAnalyticalIntegral(RooArgSet& allVars,
  			    RooArgSet& analVars,
  			    const char* /*rangeName*/) const;

  double analyticalIntegral(int code,
  			    const char* /*rangeName*/) const;

  QCDSpectrum* QCD(int c);

  CISpectrum*  CI(int c);

  void setBinRange(int first=-1, int last=-1);
  
  size_t size() { return maxsize; }
  void setSize(int n=-1)
  {
    int mxsize = qcd.size()-1;
    maxsize = n < 1 ? mxsize : n < mxsize ? n : mxsize;
  }
  void setBootstrapSize(int n=1) { nbootstrap = n; }
  void setNumber(int which=0);
  void setAsimov(bool yes=true, double lumi=19.71, double l=0);
  void setInterpolate(bool yes=true);
  void setLO(bool yes=true) { useLO = yes; }
  
  std::vector<double>& Asimov() { return asimov; }
  std::vector<double>& crossSection(int n);
  std::vector<double>& qcdXsection() { return qcdxsect; }
  
  static long double multinomial(std::vector<double>& n,
				 std::vector<double>& p);
		    
 protected:
  double evaluate() const;

  RooListProxy count;
  RooRealProxy lambda;
  RooListProxy kappa;

 private:

  std::vector<QCDSpectrum> qcd;
  std::vector<CISpectrum>  ci;
 
  double smallest;
  int number;
  bool useasimov;
  std::vector<double> asimov;
  std::vector<double> qcdxsect;
  std::vector<double> xsection;
  int firstbin;
  int lastbin;
  bool useLO;
  int maxsize;
  int nbootstrap;
  bool useinterpolation;
  mutable ROOT::Math::Interpolator* interp;
  
  ClassDef(RooInclusiveJetPdf,1)  
};
//---------------------------------------------------------------------------
#endif
