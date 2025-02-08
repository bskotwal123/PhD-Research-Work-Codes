// ---------------------------------------------------------------------------
// Read coefficients computed using ciconv (from CIJET) and compute CI cross
// section. HBP 2014
// ---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "CIXsection.h"
#include "TFile.h"
#include "hutil.h"
// ---------------------------------------------------------------------------
using namespace std;

CIXsection::CIXsection(std::string filename)
  : _filename(filename),
    loinf (std::vector<std::vector<double> >(9, std::vector<double>(12))), 
    lobi1 (std::vector<std::vector<double> >(9, std::vector<double>(18))), 
    lobi2 (std::vector<std::vector<double> >(9, std::vector<double>(8))),
    
    nloinf(std::vector<std::vector<double> >(9, std::vector<double>(12))), 
    nlobi1(std::vector<std::vector<double> >(9, std::vector<double>(18))), 
    nlobi2(std::vector<std::vector<double> >(9, std::vector<double>(8))),

    loxsec(std::vector<double>(9,0)),
    nloxsec(std::vector<double>(9,0)),

    scaf(std::vector<double>(9,0)),
    scar(std::vector<double>(9,0)),
    coe(std::vector<double>(10,0)),

    null(std::vector<double>())
{
  ifstream fin(_filename.c_str());
  if ( !fin.good() )
    {
      std::cout << "" << std::endl;
      exit(0);
    }

  string line, datum;
  // check file size
  tcall = 0;
  ncount = 0;
  {
    getline(fin, line);
    getline(fin, line);
    stringstream sin(line);
    sin >> pdfname >> pdfmember;
    ncount += 2;
    while (getline(fin, line)) ncount++;
    fin.close();
    tcall = ncount/68;
    if (  (ncount % 68 != 0) || (tcall == 0) )
      {
	cout << "wrong file format!" << endl;
	exit(0);
      }
  }

  // reopen and read coefficients
  fin.open(_filename.c_str());
  getline(fin, line);
  getline(fin, line);
  {
    stringstream sin(line);
    sin >> pdfname >> pdfmember;
  }
  getline(fin, line);
  getline(fin, line);
  {
    stringstream sin(line);
    sin >> chid >> chiu >> mad >> mau >> _mu0;
  }

  // Read coefficients
  getline(fin, line);
  for (int i = 0; i < 9; i++)
    {
      getline(fin, line);
      {
	stringstream sin(line);
	sin >> scaf[i] >> scar[i];   // muf, mur
      }

      // LO
      read(fin, loinf[i], 12);
      read(fin, lobi1[i], 18);
      read(fin, lobi2[i],  8);

      // NLO
      read(fin, nloinf[i], 12);
      read(fin, nlobi1[i], 18);
      read(fin, nlobi2[i],  8);
    }
}

double CIXsection::bi(int index, int k, int which)
{
  if (index < 0 ) return 0;
  if (index > 8 ) return 0;
  if (k < 0) return 0;
  if (k > 5) return 0;
  vector<double>& in = which == 0 ? loinf[index] : nloinf[index];
  return in[2*k];
}

double CIXsection::ai(int index, int k, int which)
{
  if (index < 0 ) return 0;
  if (index > 8 ) return 0;
  if (k < 0) return 0;
  if (k > 5) return 0;
  vector<double>& in = which == 0 ? loinf[index] : nloinf[index];
  return in[2*k+1];
}

double CIXsection::aig(int index, int k, int which)
{
  if (index < 0 ) return 0;
  if (index > 8 ) return 0;
  if (k < 0) return 0;
  if (k > 5) return 0;
  vector<double>& in = which == 0 ? loinf[index] : nloinf[index];
  return -in[2*k+1]*log(_mu0);
}


double CIXsection::bij(int index, int k, int which)
{
  if (index < 0 ) return 0;
  if (index > 8 ) return 0;
  if (k < 0) return 0;
  if (k > 8) return 0;
  vector<double>& in = which == 0 ? lobi1[index] : nlobi1[index];
  return in[2*k];
}

double CIXsection::aij(int index, int k, int which)
{
  if (index < 0 ) return 0;
  if (index > 8 ) return 0;
  if (k < 0) return 0;
  if (k > 8) return 0;
  vector<double>& in = which == 0 ? lobi1[index] : nlobi1[index];
  return in[2*k+1];
}

double CIXsection::aijg(int index, int k, int which)
{
  if (index < 0 ) return 0;
  if (index > 8 ) return 0;
  if (k < 0) return 0;
  if (k > 8) return 0;
  vector<double>& in = which == 0 ? lobi1[index] : nlobi1[index];
  return -in[2*k+1]*log(_mu0);
}

double CIXsection::bi4(int index, int k, int which)
{
  if (index < 0 ) return 0;
  if (index > 8 ) return 0;
  if (k < 0) return 0;
  if (k > 3) return 0;
  vector<double>& in = which == 0 ? lobi2[index] : nlobi2[index];
  return in[2*k];
}

double CIXsection::ai4(int index, int k, int which)
{
  if (index < 0 ) return 0;
  if (index > 8 ) return 0;
  if (k < 0) return 0;
  if (k > 3) return 0;
  vector<double>& in = which == 0 ? lobi2[index] : nlobi2[index];
  return in[2*k+1];
}

double CIXsection::ai4g(int index, int k, int which)
{
  if (index < 0 ) return 0;
  if (index > 8 ) return 0;
  if (k < 0) return 0;
  if (k > 3) return 0;
  vector<double>& in = which == 0 ? lobi2[index] : nlobi2[index];
  return -in[2*k+1]*log(_mu0);
}


// lambda is in TeV^-2 on input
vector<double>&
CIXsection::operator()(double lambda, std::vector<double>& kappa, int which)
{
  // double Lambda = 1000.0/sqrt(lambda + 1.e-200);
  // double ll = Lambda/5000;
  // norm = ll*ll;
  // norm2= norm*norm;
  // rpar = log(Lambda/_mu0);

  norm1 = 25*lambda;
  norm2 = norm1*norm1;
  rpar  = log(1.e3/_mu0) - 0.5 * log(lambda+1.e-300);  

  if ( which == 0 )
    {
      for (int i = 0; i < 9; i++)
	loxsec[i] = compute(loinf[i],  lobi1[i],  lobi2[i], kappa);
      return loxsec;
    }
  else
    {
      for (int i = 0; i < 9; i++)
	nloxsec[i] = compute(nloinf[i], nlobi1[i], nlobi2[i], kappa);
      return nloxsec;
    }
}

void 
CIXsection::read(ifstream& fin, vector<double>& x, int n, string message)
{
  string line;
  getline(fin, line);
  if ( message != "" )
    {
      cout << line << endl;
      cout << message << endl;
    }
  stringstream sin(line);
  for(int jj=0; jj < n; jj++) 
    {
      string datum;
      sin >> datum;
      int j = datum.find("D");
      if ( j > -1 ) datum.replace(j,1,"e");
      stringstream in(datum);
      in >> x[jj];
    }
}

inline
double 
CIXsection::compute(vector<double>& inf,
		    vector<double>& bi1,
		    vector<double>& bi2,
		    vector<double>& la)
{
  double xsec = 0.0;
  for(int k=0; k < 6; k++)
    xsec += la[k]*(inf[2*k] + inf[2*k+1]*rpar) * norm1;

  for(int k=0; k < 10; k++) coe[k]=0;

  for(int k=0; k < 3; k++)
    {
      coe[  3*k] = la[2*k]  * la[2*k];
      coe[1+3*k] = la[2*k+1]* la[2*k+1];
      coe[2+3*k] = la[2*k]  * la[2*k+1];
    }

  for (int k = 0; k < 9; k++)
    xsec += coe[k]*(bi1[2*k] + bi1[2*k+1]*rpar) * norm2;

  for(int k=0; k < 10; k++) coe[k]=0;

  coe[0] = la[0]*la[3];
  coe[1] = la[1]*la[3];
  coe[2] = la[4]*la[3];
  coe[3] = la[5]*la[3];

  for (int k = 0; k < 4; k++)
    xsec += coe[k]*(bi2[2*k] + bi2[2*k+1]*rpar) * norm2;
  return xsec;
}


void CIXsection::evaluate(TH1D* xsec,
			  double lambda, std::vector<double>& la,
			  vector<TH1D*>& bi,
			  vector<TH1D*>& aig,
			  vector<TH1D*>& ai,
			      
			  vector<TH1D*>& bij,
			  vector<TH1D*>& aijg,
			  vector<TH1D*>& aij,
 
			  vector<TH1D*>& bi4, 
			  vector<TH1D*>& ai4g,
			  vector<TH1D*>& ai4)
{
  // double Lambda = 1000.0/sqrt(lambda+1.e-100);
  // double f  = log(Lambda);
  // double ll = Lambda/5000;
  // double norm  = 1.0/(ll*ll);
  // double norm2 = norm*norm;

  double sqrtlambda = sqrt(lambda);
  double norm1 = 25*lambda;
  double norm2 = norm1*norm1;
  double x0 = log(1000.0);
  double f1 = log(1.0/pow(sqrtlambda+1.e-300, norm1));
  double f2 = log(1.0/pow(sqrtlambda+1.e-300, norm2));
  
  xsec->Reset();
  for(int k=0; k < 6; k++)
    *xsec = *xsec + la[k]*(*bi[k] + *aig[k] + *ai[k]*x0)*norm1
      + la[k]*(*ai[k]*f1);

  vector<double> coe(10,0);
  for(int k=0; k < 3; k++)
    {
      coe[  3*k] = la[2*k]  * la[2*k];
      coe[1+3*k] = la[2*k+1]* la[2*k+1];
      coe[2+3*k] = la[2*k]  * la[2*k+1];
    }

  for (int k = 0; k < 9; k++)
    *xsec = *xsec + coe[k]*(*bij[k] + *aijg[k] + *aij[k]*x0)*norm2
      + coe[k]*(*aij[k]*f2);

  for(int k=0; k < 10; k++) coe[k]=0;

  coe[0] = la[0]*la[3];
  coe[1] = la[1]*la[3];
  coe[2] = la[4]*la[3];
  coe[3] = la[5]*la[3];

  for (int k = 0; k < 4; k++)
    *xsec = *xsec + coe[k]*(*bi4[k] + *ai4g[k] + *ai4[k]*x0)*norm2
      + coe[k]*(*ai4[k]*f2);
}

vector<double> CIXsection::evaluate(double lambda, std::vector<double>& kappa,
				    const vector<double>& bi,
				    const vector<double>& aig,
				    const vector<double>& ai,
				    
				    const vector<double>& bij,
				    const vector<double>& aijg,
				    const vector<double>& aij,
 
				    const vector<double>& bi4, 
				    const vector<double>& ai4g,
				    const vector<double>& ai4)
{
  // double Lambda = 1000.0/sqrt(lambda+1.e-100);
  // double f  = log(Lambda);
  // double ll = Lambda/5000;
  // double norm  = 1.0/(ll*ll);
  // double norm2 = norm*norm;

  double sqrtlambda = sqrt(lambda);
  double norm1 = 25*lambda;
  double norm2 = norm1*norm1;
  double x0 = log(1000.0);
  // these terms should go to zero as lambda => 0 and, therefore,
  // norm1 and norm2 => 0
  double f1 = log(1.0/pow(sqrtlambda+1.e-300, norm1));
  double f2 = log(1.0/pow(sqrtlambda+1.e-300, norm2));
  
  double blo  = 0;
  double bnlo = 0;
  for(int k=0; k < 6; k++)
    {
      double lo = kappa[k]*bi[k]*norm1;
      blo  += lo;
      bnlo += lo + kappa[k]*(aig[k] + ai[k]*x0)*norm1 + kappa[k]*ai[k]*f1;
    }
  
  vector<double> coe(10,0);
  for(int k=0; k < 3; k++)
    {
      coe[  3*k] = kappa[2*k]  * kappa[2*k];
      coe[1+3*k] = kappa[2*k+1]* kappa[2*k+1];
      coe[2+3*k] = kappa[2*k]  * kappa[2*k+1];
    }
  
  double alo  = 0;
  double anlo = 0;  
  for (int k = 0; k < 9; k++)
    {
      double lo = coe[k]*bij[k]*norm2;
      alo  += lo;
      anlo += lo + coe[k]*(aijg[k] + aij[k]*x0)*norm2 + coe[k]*aij[k]*f2;
    }
  for(int k=0; k < 10; k++) coe[k]=0;
  coe[0] = kappa[0]*kappa[3];
  coe[1] = kappa[1]*kappa[3];
  coe[2] = kappa[4]*kappa[3];
  coe[3] = kappa[5]*kappa[3];

  for (int k = 0; k < 4; k++)
    {
      double lo = coe[k]*bi4[k]*norm2;
      alo  += lo;
      anlo += lo + coe[k]*(ai4g[k] + ai4[k]*x0)*norm2 + coe[k]*ai4[k]*f2;
    }
  vector<double> xsecs(4);
  xsecs[0] = bnlo;
  xsecs[1] = anlo;
  xsecs[2] = blo;
  xsecs[3] = alo;
  return xsecs;
}

