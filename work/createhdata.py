#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Make a histogram of the jet count per bin (19.71/fb) and simultaneously
# create a table data.tex containing the observed counts
# Fall 2014 HBP
#-------------------------------------------------------------------------------
import os, sys, re
from string import *
from histutil import *
from array import array
from time import sleep
from ROOT import *
#-------------------------------------------------------------------------------
INPFILENAME = '../data/JetYield_Winter14_V8_JEC.txt'
OUTFILENAME = '../data/data_8TeV_L19.71_V8.root'
def main():
    inp = open(INPFILENAME)
    table = '''
\\begin{table}[htp]
  \\caption{Jet yield for each $p_\\textrm{T}$ bin and $|y| < 0.5$.}
  \\label{tab:yield}
  \\medskip
  \\centering
  \\begin{tabular}{|r|rr|r|}
  \\hline
  bin\t& $p_\\textrm{T,min}$\t& $p_\\textrm{T,max}$\t& jet yield \\\\ \\hline
 '''
    pt = array('d')
    y  = []
    ptmin = 0.0
    ptmax = 0.0
    ii = 0
    while 1:
        t = split(inp.readline())
        if len(t) == 0: continue
        if t[0] != 'Bin': continue
        ptmin = atof(t[7])
        if ptmin < 507: continue

        ptmax = atof(t[9][:-1])
        count = atof(t[12])
        pt.append(ptmin)
        y.append(count)
        
        ii += 1
        record = '%4d\t&%8.0f\t&%8.0f\t&%10.0f\t\\\\ \\hline' % \
          (ii, ptmin, ptmax, count)
        print record
        table += '%s\n' % record
        if ptmax > 2400: break
    pt.append(ptmax)
    
    table += '\\end{tabular}\n\\end{table}\n'
    open('data.tex','w').write(table)

    print

    # ---------------------------------------------------------------
    # Make and plot histogram
    # ---------------------------------------------------------------
    setStyle() # set standard graphics style (see python/histutil.py)
    hfile = TFile(OUTFILENAME, 'recreate')

    hdata = TH1D('hdata', '', len(y), pt)
    hdata.GetYaxis().SetTitle('count / bin')
    hdata.GetYaxis().SetTitleOffset(1.6)
    hdata.SetNdivisions(505, "Y")
                
    hdata.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    hdata.GetXaxis().SetTitleOffset(1.2)
    hdata.SetNdivisions(505, "X")
    
    for ii in xrange(len(y)):
        hdata.SetBinContent(ii+1, y[ii])
        hdata.SetBinError(ii+1, sqrt(y[ii]))
        
    cdata = TCanvas('figs/fig_data_8TeV_L19710_V8',
                    "observed", 10, 10, 500, 500)
    cdata.cd()
    gPad.SetLogy()
    hdata.Draw('ep')
    # add a header (see python/histutil.py)
    scribe = addTitle('CMS Preliminary  #surds=8TeV CI Search L=19.71/fb',
                      0.035)
    cdata.Update()
    cdata.SaveAs('.pdf')
    sleep(5)
    hfile.Write()
    hfile.Close()
#------------------------------------------------------------------------------
try:    
    main()
except KeyboardInterrupt:
    print "ciao!"
    
