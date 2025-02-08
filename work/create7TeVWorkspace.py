#!/usr/bin/env python
#-----------------------------------------------------------------------------
# create7TeVWorkspace.py
# read spectra and create a workspace (7 TeV data)
# created 12-Oct-2014 Harrison B. Prosper - a rewrite
#-----------------------------------------------------------------------------
import os, sys, re, optparse, histutil
from math import *
from string import *
from glob import glob
from array import array
from time import sleep
from random import shuffle
from ROOT import *
#-----------------------------------------------------------------------------
LHAPATH = os.environ['LHAPDF_DATA_PATH']
MAXHIST = 500
#-----------------------------------------------------------------------------
def check(o, message):
    if o == None:
        print "** %s **" % message
        sys.exit(0)
#-----------------------------------------------------------------------------
def main():
    print "\n\t\t\t=== create7TeVWorkspace.py ==="
    gSystem.Load("libCI.so")

    PDFset   = 'CTEQ6.6'
    filename = '%s_JESJERPDF_workspace.root' % PDFset
    dirname  = '../data'
    fname    = 'CI7TeVhistograms.root'
    rootfile = '%s/%s' % (dirname, fname)
    
    # -------------------------------------                
    # get list of histograms
    # -------------------------------------            
    if not os.path.exists(rootfile):
        hutil.error("create7TeVWorkspace.py",
                    "can't find rootfile %s" % rootfile)
    hfile = TFile(rootfile)
    if not hfile.IsOpen():
        hutil.error("create7TeVWorkspace.py",
                    "can't open rootfile %s" % rootfile)
    histnames = hutil.histogramNames(hfile)
    hfile.Close()

    # -------------------------------------                
    # make a random collection of histograms
    # -------------------------------------                
    histcollection = []
    for histname in histnames:
            if histname[0] == 'c': 
                histcollection.append(histname)
    
    print "="*80
    print " number of spectra: %d" % len(histcollection)
    print " first spectrum:      ", histcollection[0]
    print " last  spectrum:      ", histcollection[-1]
    print "="*80
    
    # --------------------------------------
    # create RooFit workspace
    # --------------------------------------
    print "\t==> create workspace..."
    ws = RooWorkspace("CI")
    # Suppress info messages
    RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)    

    # -------------------------------------
    # create model parameters
    # Note: l = 1/Lambda^2
    # -------------------------------------    
    ws.factory('lambda[1.e-20, 1.e-20, 0.015]')
    l = ws.var("lambda")
    l.setBins(50)    
    # note upper case "S"; yes, annoying!
    l.SetTitle('#lambda = 1/#Lambda^{2} (TeV^{-2})')

    # create kapppa parameters
    record = []
    for ii in xrange(6):
        name = 'kappa%d' % ii
        ws.factory('%s[0,-1,1]' % name)
        record.append(name)
    ws.var('kappa0').setVal(-1) # default: LL model
    record = joinfields(record, ',')
    ws.defineSet('kappaset', record)
        
    # -------------------------------------            
    # create data parameters and
    # import data histogram into workspace
    # -------------------------------------
    hdfile = TFile(rootfile)
    hdata  = hdfile.Get('data').Clone('hdata')
    check(hdata, 'data not found')
    hdata.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    hdata.GetYaxis().SetTitle('count / bin')
    getattr(ws, 'import')(hdata)
    
    nbins = hdata.GetNbinsX()
    record = []
    for ii in xrange(nbins):
        name = 'N%2.2d' % ii
        c = hdata.GetBinContent(ii+1)
        j = int(c + 6*sqrt(c))/10000
        cmax = (j+1)*10000
        variable = '%s[%10.0f,0,%10.0f]' % (name, c, cmax)
        ws.factory(variable)
        record.append (name)
        print "%5.0f\t%s" % (hdata.GetBinLowEdge(ii+1), variable)
    record = joinfields(record, ',')
    ws.defineSet('Nset', record)    
    total = hdata.Integral()
    print "\t==> observed count: %d" % int(total)

    
    # make a RooDataSet from the counts
    data = RooDataSet('data',
                      'counts per bin',
                      ws.set('Nset'))
    data.add(ws.set('Nset'))
    getattr(ws, 'import')(data)

    # -------------------------------------            
    # create model
    # -------------------------------------
    model = RooQCDCIPdf('model', 'p(D|#lambda, #kappa)',
                        ws.set('Nset'),
                        ws.var('lambda'),
                        ws.set('kappaset'))

    print "\t==> saving spectra..."
    hfile = TFile(filename, "recreate")
    qcdspectrum = []
    cispectrum  = []
    xsec = 0.0
    for index, histname in enumerate(histcollection):
        indx   = '%3.3d' % index
        prefix = 'hQCD%s' % indx
        qcdspectrum.append( QCDSpectrum(dirname, histname, fname, prefix) )
        hfile.WriteObject(qcdspectrum[-1], "QCD%s" % indx)

        prefix = 'hCI%s' % indx
        cispectrum.append( CISpectrum(dirname, histname, fname, prefix) )
        hfile.WriteObject(cispectrum[-1], "CI%s" % indx)

        model.add(qcdspectrum[-1], cispectrum[-1])
        xsec += qcdspectrum[-1]().Integral()
        
        if index % 100 == 0: print index
    xsec /= len(histcollection)
    print "\t==> total xsection: %10.1f pb" % xsec
    scale = total / xsec
    print "\t==> data /<theory>: %10.1f/pb" % scale
    getattr(ws, 'import')(model)

    # -------------------------------------                
    print "="*80    
    #ws.Print()
    print "\t==> writing workspace to file: %s" % filename
    hfile.Close()
    ws.writeToFile(filename, kFALSE)
    print "\tdone!"    
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
