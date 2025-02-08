#!/usr/bin/env python
#-----------------------------------------------------------------------------
# plotCISpectra.py
# read spectra and plot them
# created 12-Oct-2014 Harrison B. Prosper - a rewrite
#-----------------------------------------------------------------------------
import os, sys, re
from math import *
from histutil import *
from string import *
from glob import glob
from array import array
from time import sleep
from random import shuffle
from ROOT import *
#-----------------------------------------------------------------------------
LHAPATH = os.environ['LHAPDF_DATA_PATH']
YMIN    = 1.0e-10
YMAX    = 1.0e+8
#-----------------------------------------------------------------------------
gSystem.Load("libCI.so")
LUMI = 19340.0 # 1/pb
getname = re.compile('(?<=[0-9]/)[ab].+(?=[._])')

PDFS = {'ALL'  : 'CT10nlo, MSTW2008nlo, NNPDF23_nlo',
        'CT10' : 'CT10nlo',
        'MSTW' : 'MSTW2008nlo',
        'NNPDF': 'NNPDF23_nlo'}
#-----------------------------------------------------------------------------
def makePlot(hname, spectrum, pt, COLOR=kBlue):
    x = array('d')
    nbins = pt.size()-1
    for ii in xrange(nbins+1): x.append(pt[ii])
    h = TH1D(hname, '', nbins, x)
    h.GetXaxis().SetTitle('p_{T} (GeV)')
    h.GetYaxis().SetTitle('NP #otimes #d^{2}#sigma /dp_{T}dy (pb/GeV)')
    h.SetLineColor(COLOR)
    for ii in xrange(nbins):
        pT = (pt[ii+1]+pt[ii])/2
        h.SetBinContent(ii+1, spectrum(pT))
        h.SetBinError(ii+1, 0)
    return h
#-----------------------------------------------------------------------------
def main():
    print "\n\t<=== plotCISpectra.py ===>"
    setStyle()
    os.system('mkdir -p figs/CT10;'\
              'mkdir -p figs/MSTW;'\
              'mkdir -p figs/NNPDF;'\
              'mkdir -p figs/ALL')
    
    # --------------------------------------------------------
    argv = sys.argv[1:]
    argc = len(argv)
    if argc < 1:
        print '''
    ./plotCISpectra.py PDFset [true/reco] [number]
        '''
        sys.exit(0)

    PDFset = argv[0]
    if argc > 1:
        RECO = upper(argv[1])[0] == 'R'
    else:
        RECO = False

    if argc > 2:
        PDFsetmin = atoi(argv[2])
        PDFsetmax = PDFsetmin
    else:
        PDFsetmin =  0
        PDFsetmax =100

    # --------------------------------------------------------
    print "="*80
    print "\tPDFset:   %s" % PDFset
    print "\tPDFset(range): %d - %d" % (PDFsetmin, PDFsetmax)
    if RECO:
        print "\tRECO"
    else:
        print "\tTRUE"
    print "="*80
    # --------------------------------------------------------
    
    tdir = TDirectory('CI', 'C')
    hist  = []
    kolor = [kRed, kOrange, kYellow+1, kGreen+1, kBlue, kMagenta]
    jcolor=0       
    for index in xrange(PDFsetmin, PDFsetmax+1):
        dirname = '../%s/%3.3d' % (PDFset, index)
        print dirname
        
        rootfiles = glob('%s/*.root' % dirname)
        if RECO:
            rootfiles = filter(lambda x: find(x, '_smeared') > 0, rootfiles)
        else:
            rootfiles = filter(lambda x: find(x, '_smeared') < 0, rootfiles)

        # ----------------------------------------------------
        # sort files
        # ----------------------------------------------------
        filelist = []
        for kk, rootfile in enumerate(rootfiles):
            if not os.path.exists(rootfile):
                hutil.error("plotCISpectra.py",
                            "can't find rootfile %s" % rootfile)
            hfile = TFile(rootfile)
            if not hfile.IsOpen():
                hutil.error("plotCISpectra.py",
                            "can't open rootfile %s" % rootfile)
            name = getname.findall(rootfile)[0]
            filelist.append(('%4s' % name, rootfile))
        filelist.sort()
        
        # ----------------------------------------------------
        # plot
        # ----------------------------------------------------             
        figname = 'figs/%s/coeff_%3.3d' %  (PDFset, index)
        canvas = TCanvas(figname, figname, 10, 10, 1000, 800)
        canvas.Divide(6,10)

        rootfiles.sort()
        lastname = ''
        for kk, (name, rootfile) in enumerate(filelist):
            if lastname in ['ai5',  'ai43', 'aij8',
                            'aig5', 'ai4g3','aijg8',
                            'bi5',  'bi43', 'bij8']:
                jcolor += 1
                if jcolor >= len(kolor): jcolor = 0
            lastname = strip(name)
            color = kolor[jcolor]
            # get histograms from current file
            if not os.path.exists(rootfile):
                hutil.error("plotCISpectra.py",
                            "can't find rootfile %s" % rootfile)
            hfile = TFile(rootfile)
            if not hfile.IsOpen():
                hutil.error("plotCISpectra.py",
                            "can't open rootfile %s" % rootfile)
            print rootfile
            
            histnames = filter(lambda x: x[:3] == 'nlo',
                               hutil.histogramNames(hfile))
            canvas.cd(kk+1)
            option = 'c'
            for histname in histnames:
                print "\t%s" % histname
                h = hfile.Get(histname)
                tdir.cd()
                h = h.Clone('%s%d' % (histname, kk))
                h.SetLineColor(color)
                scribe = Scribe(0.5, 0.6, 0.24)
                h.Draw(option)
                scribe.write(name)
                hist.append(h)
                option = 'c same'
            canvas.Update()
            
            hfile.Close()
        canvas.SaveAs('.pdf')            
    gApplication.Run()        
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
