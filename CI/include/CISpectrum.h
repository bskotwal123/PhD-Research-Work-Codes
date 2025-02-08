#ifndef CISPECTRUM_H
#define CISPECTRUM_H
// ---------------------------------------------------------------------------
#include <fstream>
#include <vector>
#include <string>
#include "RooAbsReal.h"
#include "TH1D.h"
#include "TFile.h"

class CISpectrum : public RooAbsReal
{
 public:
 CISpectrum() : RooAbsReal() {}
  CISpectrum(const char* namen, const char* title,
	     std::string _histdir, std::string _histname,
	     std::string _filename="",
	     std::string _prefix="");
  CISpectrum(const CISpectrum& other, const char* newname=0);
  ~CISpectrum() {}

  TObject* clone(const char* newname) const
  {
    return new CISpectrum(*this, newname);
  }
 
  TH1D*  operator()(double lambda, std::vector<double>& kappa);
  double operator()(double lambda,
		    std::vector<double>& kappa, int index) const;
  void set(int which_=0) { which = which_; }

  std::vector<double> get()
    {
      std::vector<double> xsec(4);
      xsec[0] = blo;
      xsec[1] = alo;
      xsec[2] = bnlo;
      xsec[3] = anlo;
      return xsec;
    }
  
  std::vector<double> pT() { return pt; }
  std::string dirname()    { return histdir_; }
  std::string histname()   { return histname_; }

  static int count(TFile* hfile, std::string prefix="CI");
  static int histid;

 protected:
  double evaluate() const;
    
 private:
  std::string histdir_;
  std::string histname_;
  std::string filename_;

  std::vector<std::vector<double> > bi;  
  std::vector<std::vector<double> > aig;  
  std::vector<std::vector<double> > ai; 
			
  std::vector<std::vector<double> > bij; 
  std::vector<std::vector<double> > aijg; 
  std::vector<std::vector<double> > aij; 
  
  std::vector<std::vector<double> > bi4; 
  std::vector<std::vector<double> > ai4g; 
  std::vector<std::vector<double> > ai4;

  std::vector<double> pt;
  TH1D* xsec;

  void init(std::string rootfile, 
	    std::string hname);

  void get(std::string histdir, 
	   std::string histname, 
	   std::string name, 
	   std::vector<std::vector<double> >& h,
	   double scale=0);

  int which;
  double blo, alo, bnlo, anlo;
  ClassDef(CISpectrum,1)
};


#endif
