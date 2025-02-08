#!/usr/bin/env python
#-----------------------------------------------------------------------------
# loadSpectra.py
# read spectra and cache them in a simple server
# created 12-Oct-2014 Harrison B. Prosper - a rewrite
#-----------------------------------------------------------------------------
import os, sys, re, optparse, histutil
from math import *
from string import *
from array import array
from time import sleep
from histutil import *
from ROOT import *
#-----------------------------------------------------------------------------
LHAPATH = os.environ['LHAPDF_DATA_PATH']
LAMBDA  = [10.0, 15.0, 20.0]
COLOR   = [kRed, kOrange+1, kGreen+1]
#-----------------------------------------------------------------------------
def main():
    print "\n\t\t\t=== createWorkspace.py ==="
    gSystem.Load("libCI.so")
    setStyle()
    
    QCDdir = '../fastNLO/CT10/099'
    CIdir  = replace(QCDdir, 'fastNLO', 'fastCI')
    print QCDdir
    print CIdir
    
    # -------------------------------------                
    # get list of histograms
    # -------------------------------------            
    rootfile = '%s/qcd.root' % QCDdir
    if not os.path.exists(rootfile):
        hutil.error("checkSpectra.py",
                    "can't find rootfile %s" % rootfile)
    hfile = TFile(rootfile)
    if not hfile.IsOpen():
        hutil.error("checkSpectra.py",
                    "can't open rootfile %s" % rootfile)
    histnames = hutil.histogramNames(hfile)
    hfile.Close()
        
    QCD = QCDSpectrum(QCDdir, histnames[0]) 
    CI  = CISpectrum(CIdir, histnames[0])

    # --------------------------------------------------------
    filename = 'fig_check_spectra'
    xmin = 500
    xmax = 2500
    ymin = -0.00004
    ymax =  0.00002
    xx = array('d'); xx.append(xmin); xx.append(xmax)
    yy = array('d'); yy.append(0); yy.append(0)
    horiz = TGraph(2, xx, yy)
    horiz.SetLineWidth(2)
    horiz.SetLineColor(kBlue)
    horiz.GetHistogram().SetAxisRange(xmin, xmax, "X")
    horiz.GetHistogram().SetAxisRange(ymin, ymax, "Y")
    horiz.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    horiz.GetYaxis().SetTitle('(QCD+CI)/QCD')
    
    cratio = TCanvas(filename,
                     filename,
                     10, 10, 500, 500)
    cratio.cd()
    gPad.SetLogy(kFALSE)
    horiz.Draw('ac')
   
    # Set up model
    kappa = vector('double')(6, 0)
    kappa[0] = 1
    hist = []
    option = 'c'
    for ii, L in enumerate(LAMBDA):
        l = 1.0/L**2
        color = COLOR[ii]
        
        hci = CI(l, kappa); #hutil.divideByWidth(hci)  # convert to a density
        hci.SetAxisRange(xmin, xmax, "X")
        hci.SetLineColor(color)
        hci.SetLineWidth(2)         
        #h   = QCD();        hutil.divideByWidth(h)    # convert to a density
        #h.Add(hci)
        #h.Divide(h)
        hist.append(hci.Clone())
        #h.SetLineColor(color)
        #h.SetLineWidth(2)
        cratio.cd()
        hist[-1].Draw('c same')
        option = 'c same'
        cratio.Update()
    cratio.Update()
    gApplication.Run()
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
