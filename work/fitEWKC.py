#!/usr/bin/env python
import os, sys, re
from string import *
from histutil import *
from array import array
from time import sleep
from ROOT import *
#-------------------------------------------------------------------------------
INPFILENAME = '../corrections/ElectroWeakCorrection.root'

def testf(x):
    return x**2 + x**3 + x**4
#-----------------------------------------------------------------------------
def main():
    setStyle()
    gSystem.Load('libCI.so')
    gROOT.ProcessLine('.L Function.cc+')
    
    
    hfile = TFile(INPFILENAME)
    if not hfile.IsOpen():
        print "** can't open file %s" % INPFILENAME
        sys.exit(0)

    h = hfile.Get("IncJet_Ewk_Y0")
    print h.Interpolate(50)
    print h.Interpolate(3000)
    print h.GetBinLowEdge(1)
    print h.GetBinLowEdge(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())
    
    
    f = HistFunction(h)

    N = 5
    xmin  =   50.0
    xmax  = 2500.0
    fcheb = ROOT.Math.ChebyshevApprox(f, xmin, xmax, N)
        
    #hutil.divideByWidth(h)
    
    h.GetYaxis().SetTitle('dC_{EWK}/dp_{T}  (1/GeV)')
    h.GetYaxis().SetTitleOffset(1.6)
    h.SetNdivisions(505, "Y")
                
    h.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    h.GetXaxis().SetTitleOffset(1.2)
    h.SetNdivisions(505, "X")
        
    c = TCanvas('figs/ElectroWeakCorrection', "correction", 10, 10, 500, 500)
    #gPad.SetLogy()
    #gPad.SetLogx()
    c.cd()
    h.Draw('cmore ')
    fcheb.Draw('same')
    scribe = addTitle('CMS Preliminary  #surds=8TeV CI Search L=19.7/fb', 0.035)
    c.Update()
    c.SaveAs('.pdf')
    sleep(10)
    
try:    
    main()
except KeyboardInterrupt:
    print "ciao!"
    
