#ifndef JETSPECTRUMSMEARED_H
#define JETSPECTRUMSMEARED_H
// ---------------------------------------------------------------------------
// file: JetSpectrumSmeared.h
// HBP 2012 - 2014
// ---------------------------------------------------------------------------
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include "TMath.h"
#include "TH1D.h"
#include "Math/Interpolator.h"
// --------------------=-----------------------------------------------------
// Inclusive jet spectrum in |y| < 0.5
// --------------------------------------------------------------------------
class JetSpectrum;
class JetCorrectionUncertainty;

class JetSpectrumSmeared
{
 public:
  JetSpectrumSmeared() {}

  JetSpectrumSmeared(JetSpectrum* spectrum_,
		     JetCorrectionUncertainty* JESunc_=0,
		     double JERunc_=0.1,
		     double x_=0, double y_=0,
		     double pTmin_=500.0,
		     double pTmax_=2500.0,
		     int npT_=20);

  JetSpectrumSmeared(const JetSpectrumSmeared&);

  ~JetSpectrumSmeared() {}
  
  double operator()(double pT);
  double operator()(double pTlow, double pThigh);
  double response(double pTreco, double pT);
  double sigmapT(double pT);
  void setxy(double x_, double y_) { x = x_; y = y_; }

private:
  JetSpectrum* spectrum;
  JetCorrectionUncertainty* JESunc;
  double JERunc;
  double x, y;
  double pTmin, pTmax;
  int npT;
  std::vector<double> pT;
  std::vector<double> xsection;
  ROOT::Math::Interpolator* interp;

  double pTreco_;
  double applySmearing_(double pT);
  double integrand_(double pT);
};

#endif
