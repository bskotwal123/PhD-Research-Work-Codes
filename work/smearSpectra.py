#!/usr/bin/env python
#-----------------------------------------------------------------------------
# smearSpectra.py
# Apply the jet response function given by the CMS Inclusive Jet Team
# (see 6-Sep-2013 presentation by Sanmay Ganguly to the Exotica Group) to the
# specified true spectra 
# created 09-Oct-2014 Harrison B. Prosper - a rewrite
#-----------------------------------------------------------------------------
import os, sys, re, optparse
from math import *
from histutil import *
from string import *
from glob import glob
from array import array
from time import sleep
from ROOT import gSystem, gPad, TH1D, TFile, TCanvas, kFALSE, kTRUE
#-----------------------------------------------------------------------------
LHAPATH = os.environ['LHAPDF_DATA_PATH']
YMIN    = 1.0e-10
YMAX    = 1.0e+8
#-----------------------------------------------------------------------------
# default list of histograms to be smeared
HISTNAMES = '''
nlo_0.500_0.500
nlo_0.500_1.000
nlo_1.000_0.500
nlo_1.000_1.000
nlo_1.000_2.000
nlo_2.000_1.000
nlo_2.000_2.000
'''
HISTNAMES = split(strip(HISTNAMES))
DATAFILENAME  = '../data/data_8TeV_L19.71_V8.root'
JECUNFILENAME = '../corrections/Winter14_V5_DATA_UncertaintySources_AK7PF.txt'
EWKFILENAME   = '../corrections/ElectroWeakCorrection.root'
JECUNTOTAL    = 'TotalNoTime'
#-----------------------------------------------------------------------------
LUMI = 19710.0 # 1/pb
#-----------------------------------------------------------------------------
def decodeCommandLine():
    VERSION = '22-Nov-2014'
    USAGE = '''
    python smearSpectra.py [options] <pathname to unsmeared spectra>

    options
       -s<smearing>   jes+jer+pdf, jes+jer, or pdf [default=jes+jer+pdf]
    '''
    parser = optparse.OptionParser(usage=USAGE,
                                   version=VERSION)
    parser.add_option('-s', '--smearing',
                      action='store',
                      dest='smearing',
                      type='string',
                      default='jes+jer+pdf',
                      help='level of smearing to apply')
    options, args = parser.parse_args()
    if len(args) == 0:
        print USAGE
        sys.exit(0)
    PDFdir = args[0]
    if PDFdir[-1] == '/': PDFdir=PDFdir[:-1]
    return (upper(options.smearing), PDFdir)
#-----------------------------------------------------------------------------
def makePlot(hname, spectrum, pt, COLOR=kBlue):
    x = array('d')
    nbins = pt.size()-1
    for ii in xrange(nbins+1): x.append(pt[ii])
    h = TH1D(hname, '', nbins, x)
    h.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    h.GetYaxis().SetTitle('')
    h.SetLineColor(COLOR)
    h.SetLineWidth(1)
    for ii in xrange(nbins):
        # save cross section/bin
        h.SetBinContent(ii+1, spectrum(pt[ii], pt[ii+1]))
        h.SetBinError(ii+1, 0)
    return h
#-----------------------------------------------------------------------------
def main():
    print "\n\t<=== smearSpectra.py ===>"
    
    smearingLevel, PDFdir = decodeCommandLine()
    doJESJER = find(smearingLevel, 'J') > -1
    doPDF    = find(smearingLevel, 'P') > -1

    cmd = '(?<=%s/)[0-9]+' % PDFdir
    getPDFmember = re.compile(cmd)

    gSystem.Load("libCI.so")
    from ROOT import JetCorrectorParameters
    from ROOT import JetCorrectionUncertainty
    from ROOT import JetSpectrum, JetSpectrumSmeared
    from ROOT import hutil
          
    setStyle()
    
    # ------------------------------------------------------------------------
    # determine which files to use
    # if doJESJER and not doPDF, assume that we wish to
    # include JES and JER uncertainties only, in which
    # case we use member zero, with mur=muf=1
    # and sample JES, JER 500 times.
    # ------------------------------------------------------------------------
    doSmearing = True
    if   doJESJER and not doPDF:
        member   = '000'
        histnames= ['nlo_1.000_1.000']
        nsmears  = 500
        prefix   = 'JESJER'

    # ------------------------------------------------------------------------
    # if doPDF and not doJESJER, assume that we wish to
    # include PDF uncertainties only, in which case we
    # use all PDF members but do no JES JER smearing
    # ------------------------------------------------------------------------
    elif not doJESJER and doPDF:
        member   = '*'
        histnames= HISTNAMES
        nsmears  = 1
        prefix   = 'PDF'
        
    # ------------------------------------------------------------------------
    # if doJESJER and doPDF, assume that we wish to
    # include JES, JER, and PDF uncertainties, in
    # which case we use all PDF members and
    # for each sample JES and JER ONCE.
    # ------------------------------------------------------------------------
    elif doJESJER and doPDF:
        member   = '*'
        histnames= HISTNAMES
        nsmears  = 1
        prefix   = 'JESJERPDF'
    else:
        member   = '*'
        histnames= HISTNAMES
        nsmears  = 1
        prefix   = 'NONE'
        doSmearing = False
        
    # ------------------------------------------------------------------------
    # get input rootfiles and construct output file names
    # ------------------------------------------------------------------------
    fdirs = glob('%s/%s' % (PDFdir, member))
    fdirs.sort()
    rootfiles= []
    for fdir in fdirs:
        rootfiles += glob('%s/*.root' % fdir)
        
    outfiles = []
    dirnames = {}
    for rootfile in rootfiles:
        t = split(rootfile, '/')
        filename = t[-1]
        dirname  = '%s/%s' % (joinfields(t[:-1], '/'), prefix)
        outfile  = '%s/%s' % (dirname, filename)
        if not os.path.exists(dirname):
            os.system('mkdir -p %s' % dirname)
        outfiles.append(outfile)
        
    # which is it, fastCI or fastNLO?
    fastNLO = find(rootfiles[0], 'fastNLO') > 0
    
    print "="*80
    print "\tfirst input root file:   %s" % rootfiles[0]
    print "\tlast  input root file:   %s" % rootfiles[-1]
    print "\tnumber of output files:  %d" % len(outfiles)
    print "\tnumber of smearings:     %d" % nsmears
    print "\thistograms to smear:     %s" % histnames
    
    if doJESJER:
        print "\tinclude JES+JER uncertainties"
    if doPDF:
        print "\tinclude PDF uncertainties"
    if not doJESJER and not doPDF:
        print "\tno smearing"
        
    if fastNLO:
        print "\tfastNLO"
    else:
        print "\tfastCI"
    print "="*80

    # --------------------------------------------------------
    # get correction functions
    # --------------------------------------------------------
    if not os.path.exists(JECUNFILENAME):
        print "** can't open file %s" % JECUNFILENAME
        sys.exit(0)
  
    jetcor = JetCorrectorParameters(JECUNFILENAME, JECUNTOTAL)
    JESunc = JetCorrectionUncertainty(jetcor)
    JERunc = 0.1 # assume a flat 10% uncertainty in JER

    if not os.path.exists(EWKFILENAME):
        print "** can't open file %s" % EWKFILENAME
        sys.exit(0)            
    fEWK = TFile(EWKFILENAME)
    hEWK = fEWK.Get("IncJet_Ewk_Y0")
    
    doNPcor  = True  # non-perturbative corrections
    
    # --------------------------------------------------------
    # read normalvariates.
    # we do this so that we maintain the JES/JER correlation
    # between QCD and the CI cross sections. Any (PDF member,
    # mur, muf) triplet can go with any random choice of the
    # (JES, JER) doublets, but we want to make sure that the
    # same doublet is used for the same QCD and CI triplet.
    # we do this using the variate map (see below).
    # --------------------------------------------------------
    nxy = map(lambda x: map(atof, x),
              map(split, open("normalvariates.txt").readlines()))

    # --------------------------------------------------------      
    # read data
    # --------------------------------------------------------    
    hdfile = TFile(DATAFILENAME)
    if not hdfile.IsOpen():
        print "** can't open file %s" % DATAFILENAME
        sys.exit(0)
    hdata  = hdfile.Get('hdata')
    hdata.GetYaxis().SetTitle('#sigma / bin (pb)')
    hdata.Scale(1.0/LUMI)
    ptmin = 500.0
    ptmax =2500.0    
    pT = hutil.binlowedges(hdata)
    pT.push_back(ptmax)

    # --------------------------------------------------------
    # loop over files and smear selected histograms within
    # each file
    # --------------------------------------------------------              
    cspect = TCanvas('cspec', 'spectra', 10, 10, 500, 500)
    jxy = 0
    variates = {} # map between (PDF member, mur, muf) and (x, y)
    for index, rootfile in enumerate(rootfiles):
        if not os.path.exists(rootfile):
            print "\t*** can't find rootfile %s" % rootfile
            sys.exit(0)

        member = getPDFmember.findall(rootfile)
        if len(member) == 0:
            hutil.error('smearSpectra.py', "can't get PDF member from %s" % \
                        rootfile)
        member = member[0]

        if index % 100 == 0:
            print "%5d\t%s" % (index, rootfile)
         
        # open an output file for smeared histograms
        hfile = TFile(outfiles[index], 'recreate')
        hist = []
        for histname in histnames:
            spectrum = JetSpectrum(rootfile, histname, fastNLO, doNPcor, hEWK)
            
            for ii in xrange(nsmears):

                if doSmearing:
                    # get the correct mapping between
                    # (PDFmember, mur, muf) and (x, y)
                    key = '%s/%s/%3.3d' % (member, histname, ii)
                    if variates.has_key(key):
                        x, y = variates[key]
                    else:
                        if doJESJER:
                            x, y = nxy[jxy]
                            jxy += 1
                        else:
                            x, y = 0.0, 0.0
                        variates[key] = (x, y)

                    sspectrum = JetSpectrumSmeared(spectrum,
                                                   JESunc, JERunc,
                                                   x, y)
                else:
                    sspectrum = JetSpectrumSmeared(spectrum)

                hname = '%s_%3.3d' % (histname, ii)
               
                hfile.cd()
                h = makePlot(hname, sspectrum, pT, kBlue)
                hist.append(h)

                cspect.cd()
                if fastNLO:
                    gPad.SetLogy()
                    hdata.Draw('ep')
                    h.Draw('c same')
                else:
                    gPad.SetLogy(kFALSE)
                    h.Draw('c')
                cspect.Update()
            del spectrum
            del sspectrum
        hfile.Write()
        hfile.Close()
    sleep(5)
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
