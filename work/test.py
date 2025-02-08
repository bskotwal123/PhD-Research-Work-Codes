#!/usr/bin/env python
#-----------------------------------------------------------------
# File:        test.py
#-----------------------------------------------------------------
import os,sys,re,math
from array import array
from time import sleep
from histutil import *
from string import *
from ROOT import *
#-----------------------------------------------------------------
def check(o, message):
    if o == None:
        print "** %s **" % message
        sys.exit(0)
#-----------------------------------------------------------------        
def main():
    gSystem.Load('libCI.so')
    setStyle()

    print "="*80
    print "\t==> create workspace..."
    fname = 'CTEQ6.6_JESJERPDF_workspace.root'
    hfile = TFile(fname)
    if not hfile.IsOpen():
        print "can't open %s" % fname
        sys.exit()
        
    ws = hfile.Get("CI")
    # Suppress info messages
    RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)    
    model = ws.pdf('model')
    
    rname = '../data/CI7TeVhistograms.root'
    rfile = TFile(rname)
    if not rfile.IsOpen():
        print "can't open %s" % rname
        sys.exit()

    count = model.size()
    Lambda = 10.0
    l  = 1.0/Lambda**2
    ll = l*l
    kappa = vector('double')(6, 0)
    kappa[0] = -1

    count = 10
    for ii in xrange(count):
        qcd  = model.QCD(ii)
        ci   = model.CI(ii)
        name = qcd.histname()
        numb = name[1:]
        hqcd = rfile.Get(name)
        check(hqcd, "can't find %s" % name)
        
        h = rfile.Get('b%s' % numb).Clone()
        check(h, "can't find b%s" % numb)
        
        ha= rfile.Get('a%s' % numb).Clone()
        check(ha, "can't find a%s" % numb)

        h.Scale(l)
        ha.Scale(ll)
        h.Add(ha)

        print "%10s\t%10.3f %10.3f" % (name,
                                       h.Integral(),
                                       ci(l, kappa).Integral())

#----------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print
    print "ciao!"
