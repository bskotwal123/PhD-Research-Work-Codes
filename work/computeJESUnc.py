#!/usr/bin/env python
#------------------------------------------------------------------------------
# test of JES uncertainty code
#------------------------------------------------------------------------------
import os, sys, re
from ROOT import *
#------------------------------------------------------------------------------
def main():
    gSystem.Load("../lib/libCI.so")
    filename= '../corrections/Winter14_V5_DATA_UncertaintySources_AK7PF.txt'
    jetcor = JetCorrectorParameters(filename, "TotalNoTime")
    jetunc = JetCorrectionUncertainty(jetcor)
    
    ptmin = 500.0
    ptmax =2500.0
    nstep = 25
    ptstep= (ptmax-ptmin)/nstep
    eta   = 0.0
    print "%10s\t%12s" % ("pT", "d(pT)/pT(%)")
    for ii in xrange(nstep+1):
        pt = ptmin + ii*ptstep
        
        jetunc.setJetPt(pt)
        jetunc.setJetEta(eta)
        unc = jetunc.getUncertainty(False)

        print "%10.1f\t%12.2f" % (pt, unc*100)
#------------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "ciao!"
