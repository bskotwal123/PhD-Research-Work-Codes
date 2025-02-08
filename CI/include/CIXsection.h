#ifndef CIXSECTION_H
#define CIXSECTION_H
// ---------------------------------------------------------------------------
#include <fstream>
#include <vector>
#include <string>
#include "TDirectory.h"
#include "TH1D.h"

class CIXsection
{
public:
  CIXsection() {}
  CIXsection(std::string filename);
  ~CIXsection() {}

  double mu0() { return _mu0; }

  double muf(int index) 
  { 
    return (index > -1) && (index < 9) ? scaf[index] : -1;
  }

  double mur(int index) 
  { 
    return (index > -1) && (index < 9) ? scar[index] : -1;
  }

  double LOxsec(int index) 
  { 
    return (index > -1) && (index < 9) ? loxsec[index] : -1;
  }

  double NLOxsec(int index) 
  { 
    return (index > -1) && (index < 9) ? nloxsec[index] : -1;
  }

  double bi (int index, int ii, int which=1);
  double ai (int index, int ii, int which=1);
  double bij(int index, int ii, int which=1);
  double aij(int index, int ii, int which=1);
  double bi4(int index, int ii, int which=1);
  double ai4(int index, int ii, int which=1);
  double aig (int index, int ii, int which=1);
  double aijg(int index, int ii, int which=1);
  double ai4g(int index, int ii, int which=1);

  // lambda is in TeV^-2 on input
  std::vector<double>&
    operator()(double lambda, std::vector<double>& kappa, int which=1);

  static void evaluate(TH1D* xsec,
		       double lambda, std::vector<double>& kappa,
		       std::vector<TH1D*>& bi,  
		       std::vector<TH1D*>& aig,  
		       std::vector<TH1D*>& ai, 
			
		       std::vector<TH1D*>& bij, 
		       std::vector<TH1D*>& aijg, 
		       std::vector<TH1D*>& aij, 

		       std::vector<TH1D*>& bi4, 
		       std::vector<TH1D*>& ai4g, 
		       std::vector<TH1D*>& ai4);

  static std::vector<double>
    evaluate(double lambda, std::vector<double>& kappa,
	     const std::vector<double>& bi,  
	     const std::vector<double>& aig,  
	     const std::vector<double>& ai, 
	     
	     const std::vector<double>& bij, 
	     const std::vector<double>& aijg, 
	     const std::vector<double>& aij, 
	     
	     const std::vector<double>& bi4, 
	     const std::vector<double>& ai4g, 
	     const std::vector<double>& ai4);

private:
  std::string _filename;
  std::vector<std::vector<double> >  loinf,  lobi1,  lobi2; 
  std::vector<std::vector<double> > nloinf, nlobi1, nlobi2;
  std::vector<double> loxsec, nloxsec;
  std::vector<double> scaf, scar, coe, null;

  double _mu0;
  double rpar, chid, chiu, mad, mau, norm1, norm2;
  int pdfmember, tcall, ncount; 
  std::string pdfname;

  void read(std::ifstream& fin, 
	    std::vector<double>& x, int n, 
	    std::string message="");
  inline
    double compute(std::vector<double>& inf,
		   std::vector<double>& bi1,
		   std::vector<double>& bi2,
		   std::vector<double>& la);
};

#endif
