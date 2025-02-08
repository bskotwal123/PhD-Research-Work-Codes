//---------------------------------------------------------------------------
#ifndef RooMultinomial_h
#define RooMultinomial_h
//---------------------------------------------------------------------------
// Description:
// Compute Multinominal pdf
// Created: 11-Mar-2013 Supriya Jain and Harrison B. Prosper
//          CERN
//          10-Jun-2014 HBP - compute total count using a RooFormulaVar
//---------------------------------------------------------------------------
#include <vector>
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooFormulaVar.h"
//---------------------------------------------------------------------------
class RooMultinomial : public RooAbsPdf
{
public:
  RooMultinomial() {}
  /* RooMultinomial(const char *name, const char *title, */
  /* 		 RooArgList& counts, RooArgList& thetas); */

  RooMultinomial(const char* name, const char* title,
		 RooArgSet& counts, RooArgSet& thetas);

  RooMultinomial(const RooMultinomial& other,
		 const char* name = 0);

  virtual TObject* clone(const char* newname) const 
  { return new RooMultinomial(*this, newname); }
  
  virtual ~RooMultinomial() {}

  Int_t getAnalyticalIntegral(RooArgSet& allVars, 
			      RooArgSet& analVars, 
			      const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;
  
  Int_t getGenerator(const RooArgSet& directVars, 
		     RooArgSet &generateVars, 
		     Bool_t staticInitOK=kTRUE) const;
  void generateEvent(Int_t code);

  void setTotal(int total);

protected:

  double evaluate() const;

  // Observables
  RooListProxy  _counts;
  // Parameters
  RooListProxy  _thetas;

 private:
  int _total;
  int _nbins;
  double _smallest;
  
  ClassDef(RooMultinomial,1)  
};
//---------------------------------------------------------------------------
#endif
