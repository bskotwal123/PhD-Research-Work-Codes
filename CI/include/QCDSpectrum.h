#ifndef QCDSPECTRUM_H
#define QCDSPECTRUM_H
// ---------------------------------------------------------------------------
#include <fstream>
#include <vector>
#include <string>
#include "RooAbsReal.h"
#include "TH1D.h"
#include "TFile.h"

class QCDSpectrum : public RooAbsReal
{
 public:
 QCDSpectrum() : RooAbsReal() {}
  QCDSpectrum(const char* namen, const char* title,
	      std::string _histdir, std::string _histname, 
	      std::string _filename="qcd.root",
	      std::string _prefix="");
  QCDSpectrum(const QCDSpectrum& other, const char* newname=0);
  ~QCDSpectrum() {}

  TObject* clone(const char* newname) const
  {
    return new QCDSpectrum(*this, newname);
  }
  
  TH1D*  operator()();
  double operator()(int index) const;
  std::vector<double> pT() { return pt; }
  std::string dirname() { return histdir_; }
  std::string histname() { return histname_; }

  static int count(TFile* hfile, std::string prefix="QCD");
  static int histid;

 protected:
  double evaluate() const;
  
 private:
  std::string histdir_;
  std::string histname_;
  std::string filename_;

  std::vector<double> y;
  std::vector<double> pt;

  TH1D* xsec;

  void get(std::string histdir, 
	   std::string histname, 
	   std::vector<double>& h);

  ClassDef(QCDSpectrum,1)
};


#endif
