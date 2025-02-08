#ifndef HUTIL_H
#define HUTIL _H
// ----------------------------------------------------------------------------
// file: hutil.h
// HBP 2014
// ----------------------------------------------------------------------------
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "QCDSpectrum.h"
#include "CISpectrum.h"
// ----------------------------------------------------------------------------
struct hutil
{
  static std::vector<double> contents(TH1* hist);
  static std::vector<double> binlowedges(TH1* hist);
  static std::vector<double> binhighedges(TH1* hist);
  static std::vector<double> bincenters(TH1* hist);
  static std::vector<double> errors(TH1* hist);
  static std::vector<double> cdf(TH1* hist);
  static void setContents(TH1* hist, std::vector<double>& c);
  static void setErrors(TH1* hist, std::vector<double>& err);
  static void divideByWidth(TH1* hist);
  static void error(std::string program, std::string message);
  static std::vector<TH1*> histograms(TFile* hfile, std::string hdir="");
  static std::vector<std::string> histogramNames(TFile* hfile);
};

#endif
