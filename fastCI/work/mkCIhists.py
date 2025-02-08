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
LHAPATH = os.environ['LHAPDF_DATA_PATH']
YMIN    = 1.0e-10
YMAX    = 1.0e+8
getbin  = re.compile('(?<=ci).+(?=.txt)')
gSystem.Load("libCI.so")

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
MURMUF = [0, 3, 1, 4, 7, 5, 8]
#-----------------------------------------------------------------------------
def getBins():    
    dirname = '../CT10/000'
    txtfilenames = glob('%s/*.txt' % dirname)
    txtfilenames.sort()

    pt = array('d')
    for filename in txtfilenames:
        pt.append( atof(getbin.findall(filename)[0]) )
    pt.append(3500)
    return pt

def makeHist(hname, name, pt, data):
    h = TH1D(hname, '', len(pt)-1, pt)
    h.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    h.SetNdivisions(504, "Y")
    h.GetYaxis().SetTitle('d^{2}%s/dp_{T}dy (pb/GeV)' % name)
    #h.SetMinimum(YMIN)
    #h.SetMaximum(YMAX)
    h.SetLineWidth(1)
    
    for ii, d in enumerate(data):
        width = pt[ii+1]-pt[ii]
        d = d / (width*0.5)
        h.SetBinContent(ii+1, d)
        h.SetBinError(ii+1, 0)
    return h
#-----------------------------------------------------------------------------
def makeHistograms(dirname, name, number, pt, xsection):
    y = [0]*len(xsection)
    for ii in xrange(number):
        namen = '%s%d' % (name, ii)
        hfilename = '%s/%s.root' % (dirname, namen)
        #print "\t%s" % hfilename
        hfile = TFile(hfilename, 'recreate')
        hist = []
        for which, prefix in [(1, 'nlo'),
                              (0, 'lo')]: # NLO, LO
            for jj, jjj in enumerate(MURMUF):
                hname = '%s%s' % (prefix, HISTNAMES[jj])
                #mur = xsection[0].mur(jjj)
                #muf = xsection[0].muf(jjj)
                #print hname, mur, muf
                for kk, xsect in enumerate(xsection):
                    cmd = 'xsect.%s(jjj, ii, which)' % name
                    y[kk] = eval(cmd)
                h = makeHist(hname, name, pt, y)
                hist.append(h)
        hfile.Write()
        hfile.Close()
#-----------------------------------------------------------------------------
def main():
    argv = sys.argv[1:]
    argc = len(argv)
    if argc < 1:
        print '''
    ./mkCIhists.py PDFset [PDFindex=all]
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
        PDFindexMin = 401
        PDFindexMax = 500

    print
    print "\t==> PDFset:            %s" % PDFset
    print "\t==> PDFset index(min): %d" % PDFindexMin
    print "\t==> PDFset index(max): %d" % PDFindexMax

    pt = getBins()
    hfile = []
    for index in xrange(PDFindexMin, PDFindexMax+1):
        dirname = '%s/%3.3d' % (pdfsetdir, index)        
        if index % 10 == 0:
            print dirname
            
        xsection= []
        for pT in pt[:-1]:
            filename= '%s/ci%4.4d.txt' % (dirname, pT)
            xsection.append( CIXsection(filename) )
            
        makeHistograms(dirname, 'bi',  6, pt, xsection)
        makeHistograms(dirname, 'aig', 6, pt, xsection)
        makeHistograms(dirname, 'ai',  6, pt, xsection)
        
        makeHistograms(dirname, 'bij', 9, pt, xsection)
        makeHistograms(dirname, 'aijg',9, pt, xsection)
        makeHistograms(dirname, 'aij', 9, pt, xsection)
        
        makeHistograms(dirname, 'bi4', 4, pt, xsection)
        makeHistograms(dirname, 'ai4g',4, pt, xsection)                      
        makeHistograms(dirname, 'ai4', 4, pt, xsection)                      
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
