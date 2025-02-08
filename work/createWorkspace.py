#!/usr/bin/env python
#-----------------------------------------------------------------------------
# loadSpectra.py
# read spectra and cache them in a simple server
# created 12-Oct-2014 Harrison B. Prosper - a rewrite
#         08-Feb-2015 HBP - add nominal QCD spectrum, that is, the
#                     spectrum computed with PDF member zero
#-----------------------------------------------------------------------------
import os, sys, re, optparse, histutil
from math import *
from string import *
from glob import glob
from array import array
from time import sleep
from random import shuffle, randint
from ROOT import gSystem, TFile, kFALSE, kTRUE, \
     RooWorkspace, RooMsgService, RooFit, RooDataSet, RooCmdArg
#-----------------------------------------------------------------------------
LHAPATH = os.environ['LHAPDF_DATA_PATH']
DATAFILE= '../data/data_8TeV_L19.71_V8.root'
#-----------------------------------------------------------------------------
def decodeCommandLine():
    VERSION = '05-Feb-2015'
    USAGE = '''
    python createWorkspace.py [options] <PDF sets>

    options
       -s<smearing>  jes+jer+pdf, jes+jer, pdf [def.=jes+jer+pdf]
       -o<root-filename>                       [def.=<PDF>_<s>_workspace.root]
    '''
    parser = optparse.OptionParser(usage=USAGE,
                                   version=VERSION)
    parser.add_option('-s', '--smearing',
                      action='store',
                      dest='smearing',
                      type='string',
                      default='jes+jer+pdf',
                      help='level of smearing to use')

    parser.add_option('-o', '--output',
                      action='store',
                      dest='filename',
                      type='string',
                      default='',
                      help='output file containing QCD and CI spectra')    
        
    options, args = parser.parse_args()
    if len(args) == 0:
        print USAGE
        sys.exit(0)
    PDFsets = args

    # create directory name containing smeared spectra
    directory = upper(replace(options.smearing, '+', ''))

    # create name of output file
    filename = options.filename
    if filename == '':
        if len(PDFsets) > 1:
            prefix = 'ALL'
        else:
            prefix = PDFsets[0]
        filename = '%s_%s_workspace.root' % (prefix, directory)
        
    return (directory, PDFsets, filename)
#-----------------------------------------------------------------------------
def main():
    print "\n\t\t\t=== createWorkspace.py ==="
    
    dirname, PDFsets, filename = decodeCommandLine()

    gSystem.Load("libCI")
    from ROOT import hutil, QCDSpectrum, CIXsection, CISpectrum, \
     RooInclusiveJetPdf
    
    # -------------------------------------
    # determine which files to use
    # -------------------------------------
    if   dirname == 'JESJER':
        member = '000'
    else:
        member = '*'

    # -------------------------------------                
    # get list of histograms
    # -------------------------------------            
    rootfile = '../fastNLO/CT10/000/%s/qcd.root' % dirname
    if not os.path.exists(rootfile):
        hutil.error("createWorkspace.py",
                    "can't find rootfile %s" % rootfile)
    hfile = TFile(rootfile)
    if not hfile.IsOpen():
        hutil.error("createWorkspace.py",
                    "can't open rootfile %s" % rootfile)
    histnames = hutil.histogramNames(hfile)
    hfile.Close()

    # -------------------------------------                
    # make a sample of histograms in which
    # the renormalization and factorization
    # scales are randomly selected.
    # Each histogram corresponds to a
    # randomly selected pair of scales,PDFs,
    # jet energy scale, and jet energy
    # resolution.
    # -------------------------------------
    # get input directories
    histcollection = []    
    for pdfset in PDFsets:
        QCDdirs = glob('../fastNLO/%s/%s/%s' % (pdfset, member, dirname))
        QCDdirs.sort()
        QCDdirs = QCDdirs[1:]
                           
        for qcddir in QCDdirs:
            # randomly pick renormalization and
            # factorization scales
            kk = randint(0, 6)
            histcollection.append((qcddir, histnames[kk]))

    # insert nominal QCD histogram at position 0
    qcddir = '../fastNLO/%s/%s/%s' % (pdfset, '000', dirname)
    fqcdnom = glob(qcddir)
    if len(fqcdnom) == 1:
        fqcdnom = fqcdnom[0]
    else:
        hutil.error("createWorkspace.py",
                    "can't find directory %s" % qcddir)      
    histname = 'nlo_1.000_1.000_000'
    histcollection.insert(0, (qcddir, histname))

    print "="*80
    print " number of spectra: %d" % len(histcollection)
    print " first spectrum:      ", histcollection[0]
    print " last  spectrum:      ", histcollection[-1]
    print "="*80
    
    # --------------------------------------
    # create RooFit workspace
    # --------------------------------------
    print "\t==> create workspace..."
    RooWorkspace.autoImportClassCode(kFALSE)
    RooWorkspace.addClassDeclImportDir('../CI')
    RooWorkspace.addClassImplImportDir('../CI')
    
    ws = RooWorkspace("CI")

    # Suppress info messages
    RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)    

    # -------------------------------------
    # create model parameters
    # Note: l = 1/Lambda^2
    # -------------------------------------    
    ws.factory('lambda[0, 0, 0.015]')
    l = ws.var("lambda")
    l.setBins(150)    
    # note upper case "S"; yes, annoying!
    l.SetTitle('#lambda = 1/#Lambda^{2} (TeV^{-2})')

    # create kapppa parameters
    record = []
    for ii in xrange(6):
        name = 'kappa%d' % ii
        ws.factory('%s[0,-1,1]' % name)
        record.append(name)
    ws.var('kappa0').setVal(-1) # default: LL model
    record = joinfields(record, ',')
    ws.defineSet('kappaset', record)
        
    # -------------------------------------            
    # create data set
    # -------------------------------------
    hdfile = TFile(DATAFILE)
    hdata  = hdfile.Get('hdata')
    hdata.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    hdata.GetYaxis().SetTitle('count / bin')
    nbins = hdata.GetNbinsX()
    getattr(ws, 'import')(hdata, 'hdata')
    
    record = [] 
    for ii in xrange(nbins):
        c = hdata.GetBinContent(ii+1)
        j = int(c + 10*sqrt(c))/10000
        cmax = (j+1)*10000           
        name = 'N%2.2d' % ii
        variable = '%s[%10.0f, 0,%10.0f]' % (name, c, cmax)
        ws.factory(variable)
        record.append (name)
        print "%5.0f\t%s" % (hdata.GetBinLowEdge(ii+1), variable)
        record.append (name)
    record = joinfields(record, ',')
    ws.defineSet('Nset', record)    

    # make a RooDataSet from the counts
    data = RooDataSet('data',
                      'counts per bin',
                      ws.set('Nset'))
    data.add(ws.set('Nset'))
    getattr(ws, 'import')(data, RooCmdArg())

    # -------------------------------------            
    # create model
    # -------------------------------------
    model = RooInclusiveJetPdf('model', 'p(D|#lambda, #kappa)',
                               ws.set('Nset'),
                               ws.var('lambda'),
                               ws.set('kappaset'))

    print "\t==> saving spectra..."
    hfile = TFile(filename, "recreate")
    qcdspectrum = []
    cispectrum  = []
    qcdrecords  = []
    cirecords   = []
    for index, (QCDdir, histname) in enumerate(histcollection):
        
        CIdir = replace(QCDdir, 'fastNLO', 'fastCI')
        if index % 50 == 0:
            print "%4d" % index, QCDdir, histname

        name = "QCD%3.3d" % index
        qcdrecords.append(name)
        qcdspectrum.append(QCDSpectrum(name, name, QCDdir, histname))

        name = "CI%3.3d" % index
        cirecords.append(name)
        cispectrum.append( CISpectrum(name, name, CIdir, histname) )
        
        model.add(qcdspectrum[-1], cispectrum[-1])

    getattr(ws, 'import')(model, RooCmdArg())
    # -------------------------------------                
    print "="*80    
    ws.Print()
    print "\t==> writing workspace to file: %s" % filename
    hfile.Close()
    ws.writeToFile(filename, kFALSE)
    print "\tdone!"    
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
