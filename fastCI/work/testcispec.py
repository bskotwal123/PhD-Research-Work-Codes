#!/usr/bin/env python
#-------------------------------------------------------------------------------
import os, sys, re
from string import *
from glob import glob
from ROOT import *
#-------------------------------------------------------------------------------
HISTNAMES = '''
_0.500_0.500
_0.500_1.000
_1.000_0.500
_1.000_1.000
_1.000_2.000
_2.000_1.000
_2.000_2.000
'''
HISTNAMES = split(strip(HISTNAMES))
#-------------------------------------------------------------------------------
def main():
    gSystem.Load('libCI.so')
    
    cinlo= []
    print '\t...loading NLO CI histograms...'
    for member in xrange(1, 2):
        histdir = '../CT10/%3.3d' % member
        swatch = TStopwatch()
        for postfix in HISTNAMES:
            hname = 'nlo%s' % postfix
            cinlo.append( CISpectrum(histdir, hname) )
        print "%s\ttime/constructor: %10.3fs" % \
            (histdir, swatch.RealTime()/(2*len(HISTNAMES)))

    CI = cinlo[-1]

    Lambda = 10.0
    l = 1.0/Lambda**2
    kappa = vector('double')(6, 0)
    kappa[0] =-1
    
    h  = CI(l, kappa)
#------------------------------------------------------------------------------ 
main()

