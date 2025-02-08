#ifndef JETSPECTRUM_H
#define JETSPECTRUM _H
// ----------------------------------------------------------------------------
// file: JetSpectrum.h
// HBP 2012 - 2014
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
#include "TH1F.h"
#include "TH1D.h"
#include "Math/Interpolator.h"
// --------------------------------------------------------------------------
// Inclusive jet spectrum in |y| < 0.5
// --------------------------------------------------------------------------
class JetSpectrum
{
 public:
  JetSpectrum() {}

  JetSpectrum(std::string filename,
	      std::string histname,
	      bool positive_=true,
	      bool applyNPC_=true,
	      TH1F* hEWKC_=0);

  ~JetSpectrum() {}
  
  double operator()(double pT);
  double operator()(double pTlow, double pThigh);
  double NPC(double pT);
  double EWKC(double pT);
  bool   positive() {return _positive;}
  std::vector<double> ptMin() {return ptlo;}
  std::vector<double> ptMax() {return pthi;}
  std::vector<double> ptCenter() {return ptcn;}
  std::vector<double> crossSection() {return xsection;}
  
  bool null() {return nullHist;}

private:
  std::vector<double> ptlo;
  std::vector<double> pthi;
  std::vector<double> ptcn;
  std::vector<double> xsection;
  ROOT::Math::Interpolator* interp;
  bool _positive;
  bool applyNPC;
  bool applyEWKC;
  bool nullHist;
  TH1F*  hEWKC;
};

#endif
