#!/usr/bin/env python
#-----------------------------------------------------------------------------
# Run fast NLO
# created 09-Oct-2012 Harrison B. Prosper
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
TABLE7TeV = 'fnl2332y0.tab'
TABLE8TeV = 'fnl3332y0.tab'
LHAPATH   = os.environ['LHAPDF_DATA_PATH']

getxsect = re.compile('The scale factor chosen here is:[.\s]+(?=\==)',
              re.M)

def main():
    print "\n\t<=== runfastNLO.py ===>"
    setStyle()

    argv = sys.argv[1:]
    argc = len(argv)
    if argc < 1:
        print '''
    ./runfastnlo.py PDFset [PDFindex=all]
        '''
        sys.exit(0)

    table  = TABLE8TeV
    PDFset = argv[0]
    PDFsetfilename = '%s/%s' % (LHAPATH, PDFset)

    if not os.path.exists(PDFsetfilename):
        print "** can't find %s" % PDFset
        sys.exit(0)

    if   PDFset[:2] == 'CT':
        pdfsetdir = '../CT10'
    elif PDFset[:2] == 'MS':
        pdfsetdir = '../MSTW'
    elif PDFset[:2] == 'NN':
        pdfsetdir = '../NNPDF'
    else:
        print "wrong PDFset %s" % PDFset
        sys.exit(0)

    if argc > 1:
        PDFindexMin = atoi(argv[1])
        PDFindexMax = PDFindexMin
    else:
        PDFindexMin =   0
        PDFindexMax = 500

    numberOfScales = 7        # Number of renormalization/factorization scales
    evolutionCode  = 'LHAPDF'
    
    print
    print "\t==> PDFset:            %s" % PDFset
    print "\t==> PDFset index(min): %d" % PDFindexMin
    print "\t==> PDFset index(max): %d" % PDFindexMax

    for index in xrange(PDFindexMin, PDFindexMax+1):
        dirname = '%s/%3.3d' % (pdfsetdir, index)
        cmd = 'mkdir -p %s' % dirname
        os.system(cmd)
        
        cmd = 'fnlo-tk-cppread %s %s %d %s %d' % (table,
                                                  PDFset,
                                                  numberOfScales,
                                                  evolutionCode,
                                                  index)
        print "\n%s" % cmd
        records = os.popen(cmd).readlines()

        ii = 0
        while ii < len(records):
            record = strip(records[ii])
            if record != "My Cross Sections":
                ii += 1
                continue
            
            ii += 1
            record = records[ii]
            mur, muf = split(record)[-2:]; mur = mur[:-1]
            outfilename = '%s/qcd_%s_%s.txt' % (dirname, mur, muf)
            print "\t%s" % outfilename
            out = open(outfilename, 'w')
            rec = '%7s %7s %20s %20s' % ('pT-low', 'pT-high',
                                         'LO-xsection', 'NLO-xsection')
            out.write('%s\n' % rec)
            ii += 3
            for jj in xrange(48):
                ii += 1
                record = records[ii]
                t = split(record)[-6:-1]
                t = (t[0],t[1],t[-2],t[-1])
                rec = '%7s %7s %20s %20s' % t
                out.write('%s\n' % rec)
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
