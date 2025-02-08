#!/usr/bin/env python
import os, sys, re
from string import *
from histutil import *
from array import array
from time import sleep
from ROOT import *
#-------------------------------------------------------------------------------
INPFILENAME = '../corrections/ElectroWeakCorrection.root'

def funNP(x, p):
    # from Sanmay Ganguly (Nov. 2014)
    A =  1.003
    B = 77.374
    n =  1.385
    pt = x[0]
    return A + B/pt**n    
#-----------------------------------------------------------------------------
def main():

    gSystem.Load('libCI.so')

    A =  1.003
    B = 77.374
    n =  1.385    
    fNP = TF1('NPcor', "%f + %f/pow(x[0], %f)" % (A, B, n), 500, 2500)
    fNP.SetLineColor(kRed)
    fNP.SetLineWidth(2)

    setStyle()        
    hfile = TFile(INPFILENAME)
    if not hfile.IsOpen():
        print "** can't open file %s" % INPFILENAME
        sys.exit(0)

    htmp1 = mkhist1('htmp1', 'Jet p_{T} (GeV)',
                    'Correction Factor', 50, 500, 2500)
    htmp1.SetMinimum(0.9)
    htmp1.SetMaximum(1.2)
    htmp1.SetNdivisions(505, "Y")
    
    h = hfile.Get("IncJet_Ewk_Y0")


    h.GetYaxis().SetTitle('Corrections')
    h.GetYaxis().SetTitleOffset(1.6)
    h.SetNdivisions(505, "Y")
                
    h.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    h.GetXaxis().SetTitleOffset(1.2)
    h.SetNdivisions(505, "X")
    h.SetMinimum(0)
    h.SetMaximum(1.5)
    h.SetLineColor(kBlue)
    h.SetLineWidth(2)
    h.SetAxisRange(500, 2500)

    legend = mklegend(0.24, 0.70, 0.36, 0.18)
    legend.AddEntry(h, 'EWK correction', 'l')    
    legend.AddEntry(fNP, 'NP correction', 'l')

    
    c1 = TCanvas('figs/Corrections', "Corrections", 10, 10, 500, 500)
    c1.cd()
    htmp1.Draw()
    fNP.Draw('c same')
    h.Draw('c same')
    legend.Draw()
    scribe = addTitle('CMS Preliminary  #surds=8TeV CI Search', 0.035)
    c1.Update()
    c1.SaveAs('.pdf')

    sleep(10)
#-----------------------------------------------------------------------------
try:    
    main()
except KeyboardInterrupt:
    print "ciao!"
    
