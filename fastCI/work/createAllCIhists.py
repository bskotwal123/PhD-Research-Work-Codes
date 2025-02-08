#!/usr/bin/env python
#-----------------------------------------------------------------------------
# create inclusive jet CI histograms
# created 13-Oct-2014 Harrison B. Prosper
#-----------------------------------------------------------------------------
import os, sys, re
from math import *
from histutil import *
from string import *
from glob import glob
from array import array
from time import sleep
from ROOT import *
#-----------------------------------------------------------------------------
def main():
    print "\n\t<=== makeAllCIhists.py ===>"
    for PDFset in ['MSTW', 'NNPDF']:
        cmd = '\n./mkCIhists.py %s' % PDFset
        print cmd
        os.system(cmd)                   
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
