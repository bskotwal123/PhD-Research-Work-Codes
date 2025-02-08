#!/usr/bin/env python
#-----------------------------------------------------------------
# File:        createWorkspace.py
# Description: create 8 TeV CI workspace
# Created:     02-Mar-2014 HBP
# Modified:    20-May-2014 Greg Myers
#              22-May-2014 HBP add 1.e-100 to logs to avoid a NAN
#                          when l = 0
#              09-Jun-2014 Greg - fixed order to histograms in main
#              18-Jun-2014 Greg - added various debugging methods
#              11-Nov-2014 HBP - rework
#-----------------------------------------------------------------
import os,sys,re,math
from array import array
from time import sleep
from histutil import *
from string import *
from ROOT import *
#-----------------------------------------------------------------
# Get histogram contents
def getContents(h):
    xbins = h.GetNbinsX()
    c = []
    for ii in xrange(xbins):
        c.append(h.GetBinContent(ii+1))
    return c
#-----------------------------------------------------------------
# decode name[range] into (name, [range])
def decodeParameter(param):
    from string import split
    t = split(param, '[')
    if len(t) == 1:
        return (param, None)
    else:
        return (t[0], '[%s' % t[1])
#-----------------------------------------------------------------
# get list of variables from an expression
expformat = re.compile('[1-9][.]?e[-+]?[0-9]+')
words = re.compile('[a-zA-Z][a-zA-Z0-9_]*[\[\]\.0-9,]*')
def nameList(record):
    from string import joinfields
    record = expformat.sub("&", record)
    # hide min and max first
    record = replace(record, 'min', '!')
    record = replace(record, 'max', '@')
    record = replace(record, 'log', '#')
    record = replace(record, 'exp', '$')
    record = replace(record, 'sqrt','%')
    vlist = words.findall(record)
    vl = []
    for x in vlist:
        if not (x in vl): vl.append(x)
    return vl
#-----------------------------------------------------------------
# strip away ranges
striprange = re.compile('\[.+?\]')
def stripRange(expression):
    return striprange.sub("", expression)
#-----------------------------------------------------------------
def setValues(ws, name, h):
    nbins = h.GetNbinsX()
    for ii in xrange(nbins):
        cmd = '%s_%2.2d' % (name, ii)
        if ws.var(cmd) == None:
            print "** can't find: %s" % cmd
            sys.exit(0)
        x = h.GetBinContent(ii+1)
        ws.var(cmd).setVal(x)
#-----------------------------------------------------------------
def main():
    filenames = sys.argv[1:]
    if len(filenames) == 0:
        print '''
    Usage:
        ./createWorkspace.py <spectrum-filenames>
        '''
        sys.exit(0)
        
    lumi0  = 19.34 # 1/pb
    energy = '8TeV'
        
    # --------------------------------------
    # load various codes needed for the
    # calculatioms
    # --------------------------------------
    gSystem.Load('libCI.so')

    # set up some standard graphics style
    setStyle()
    
    # --------------------------------------
    # create a RooFit workspace for
    # the model to be constructed.
    # --------------------------------------
    ws = RooWorkspace("CI")
    
    # Suppress info messages
    RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)
    
    print "="*80
    print "create %s workspace" % energy
    print "integrated luminosity: %6.2f" % lumi0
    print
    
    # -------------------------------------
    # create model parameters of interest
    # (l = 1/Lambda^2)
    # -------------------------------------    
    ws.factory('l[0, 0, 0.015]')
    l = ws.var("l")
    l.setBins(100)    
    # note upper case "S"; yes, annoying!
    l.SetTitle('#lambda = 1/#Lambda^{2} (TeV^{-2})')

    # create kapppa parameters
    for ii in xrange(6):
        ws.factory('kappa%d[0,-1,1]' % ii)

    # -------------------------------------            
    # create data parameters
    # -------------------------------------
    hdfile = TFile('../data/data_8TeV_L19.34.root')
    hdata  = hdfile.Get('hdata')
    hdata.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    hdata.GetYaxis().SetTitle('count / bin')
    nbins = hdata.GetNbinsX()
    record = []
    for ii in xrange(nbins):
        name = 'N%2.2d' % ii
        c = hdata.GetBinContent(ii+1)
        j = int(c + 6*sqrt(c))/1000
        mxc = (j+1)*1000
        ws.factory('%s[%f,0,%f]' % (name, c, mxc))
        record.append (name)
    record = joinfields(record, ',')
    ws.defineSet('setN', record)    
    # --------------------------------------
    # make a RooDataSet from the counts
    # --------------------------------------
    data = RooDataSet('data',
                      'counts per bin',
                      ws.set('setN'))
    data.add(ws.set('setN'))
    getattr(ws, 'import')(data)

    # -------------------------------------            
    # load spectra
    # -------------------------------------
    print "\nloading spectra..."
    hfile = []
    for filename in filenames:
        hfile.append( TFile(filename) )
        nspectra = min(QCDSpectrum.count(hfile[-1]),
                       CISpectrum.count(hfile[-1]))
        for ii in xrange(nspectra):
            jj = ii + 1
            qcd = hfile[-1].Get("QCD%3.3d" % jj)
            print qcd.GetName()
            if ii > 3: sys.exit(0)
            continue
            
            qcd = QCDSpectrum.get(hfile[-1], "QCD%3.3d" % jj)
            if qcd == None: break

            ci = CISpectrum.get(hfile[-1], "CI%3.3d" % jj)
            if ci == None: break
            

            getattr(ws, 'import')(ci)
                                        
            if jj % 100 == 0:
                print "%4d %s %s" % (jj, qcd.dirname(), qcd.histname())
                print "%5s %s %s" % ('', ci.dirname(), ci.histname())
    ws.Print()
    sys.exit(0)


    # --------------------------------------
    # create multinomial model
    # --------------------------------------
    counts = RooArgList(ws.set('setN')); counts.setName('counts')
    sigmas = RooArgList(ws.set('setsigma')); sigmas.setName('sigmas')
    pdf = RooMultinomial("pdf",
                         "p(N_{0},..,N_{19}|#theta_{0},..,#theta_{19})",
                         counts,
                         sigmas)
 

    filename = 'CI%sworkspace.root' % prefix

    print '\n\n=== writing workspace to %s ===' % filename
    ws.writeToFile(filename)
    print 'done!'
#----------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print
    print "ciao!"
