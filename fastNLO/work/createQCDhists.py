#!/usr/bin/env python
#-----------------------------------------------------------------------------
# create inclusive jet QCD histograms
# created 09-Oct-2014 Harrison B. Prosper
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
LHAPATH = os.environ['LHAPDF_DATA_PATH']
YMIN    = 1.0e-10
YMAX    = 1.0e+8
#-----------------------------------------------------------------------------
def getBins():    
    dirname = '../CT10/001'
    txtfilenames = glob('%s/*.txt' % dirname)
    txtfilenames.sort()
    filename = txtfilenames[0]
    records = map(lambda x: map(atof, x),
                  map(split, open(filename).readlines()[1:]))
    pt = array('d')
    for ptlow, pthigh, lo, nlo in records: pt.append(ptlow)
    ptlow, pthigh, lo, nlo = records[-1]
    pt.append(pthigh)
    print "\t==> pT_min: %5.1f\tpT_max: %6.1f" % (pt[0], pt[-1])
    return pt

def makeHist(txtfilename, pt, which=1):
    records = map(lambda x: map(atof, x),
                  map(split, open(txtfilename).readlines()[1:]))    
    name = split(txtfilename, 'qcd')[-1]
    if which == 0:
        hname = "lo%s" % replace(name, '.txt', '')
    else:
        hname = "nlo%s" % replace(name, '.txt', '')
    
    h = TH1D(hname, '', len(pt)-1, pt)
    h.GetXaxis().SetTitle('p_{T} (GeV)')
    h.SetNdivisions(504, "Y")
    h.GetYaxis().SetTitle('d^{2}#sigma /dp_{T}|y| (pb/GeV)')
    h.SetMinimum(YMIN)
    h.SetMaximum(YMAX)
    h.SetLineWidth(1)
    
    for ii, (ptlo, pthi, lo, nlo) in enumerate(records):
        if which == 0:
            h.SetBinContent(ii+1, lo)
        else:
            h.SetBinContent(ii+1, nlo)
        h.SetBinError(ii+1, 0)
    return h
#-----------------------------------------------------------------------------
def main():
    print "\n\t<=== createQCDhists.py ===>"
    setStyle()

    argv = sys.argv[1:]
    argc = len(argv)
    if argc < 1:
        print '''
    ./createQCDhists.py PDFset [PDFindex=all]
        '''
        sys.exit(0)

    PDFset = argv[0]
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

    print
    print "\t==> PDFset:            %s" % PDFset
    print "\t==> PDFset index(min): %d" % PDFindexMin
    print "\t==> PDFset index(max): %d" % PDFindexMax

    pt = getBins()
    hfile = []
    for index in xrange(PDFindexMin, PDFindexMax+1):
        dirname = '%s/%3.3d' % (pdfsetdir, index)
        hfilename = '%s/qcd.root' % dirname
        print "\t%s" % hfilename
        hfile.append( TFile(hfilename, 'recreate') )
        txtfilenames = glob('%s/*.txt' % dirname)
        txtfilenames.sort()
        hnlo = []
        hlo  = []
        for txtfilename in txtfilenames:
            hfile[-1].cd()
            hnlo.append( makeHist(txtfilename, pt, 1) )
        for txtfilename in txtfilenames:
            hfile[-1].cd()
            hlo.append( makeHist(txtfilename, pt, 0) )
        hfile[-1].Write()
        hfile[-1].Close()
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
