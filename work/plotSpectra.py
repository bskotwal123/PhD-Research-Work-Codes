#!/usr/bin/env python
#-----------------------------------------------------------------------------
# plotSpectra.py
# read spectra and plot them
# created 12-Oct-2014 Harrison B. Prosper - a rewrite
#-----------------------------------------------------------------------------
import os, sys, re, optparse
from math import *
from histutil import *
from string import *
from array import array
from time import sleep
from ROOT import *
#-----------------------------------------------------------------------------
LHAPATH = os.environ['LHAPDF_DATA_PATH']
YMIN    = 1.0e-10
YMAX    = 1.0e+8
#                        kappa
#                   1   2   3   4   5   6
KAPPA   = {'LL' : [-1,  0,  0,  0,  0,  0],
           'RR' : [ 0,  0,  0,  0, -1,  0],
           'VV' : [-1,  0, -2,  0, -1,  0],
           'AA' : [-1,  0, +2,  0, -1,  0],
           'V-A': [ 0,  0, -2,  0,  0,  0]}
MODEL   = ['LL', 'RR', 'VV', 'AA', 'V-A']

class Context:
    pass    
#-----------------------------------------------------------------------------
gSystem.Load("libCI.so")
LUMI = 19340.0 # 1/pb
getpdfdir   = re.compile('(?<=fastNLO/).+(?=/[0-9])|(?<=fastCI/).+(?=/[0-9])')
getsmearing = re.compile('(?<=/[0-9][0-9][0-9]/).+')
PDFS = {'ALL'  : 'CT10nlo, MSTW2008nlo, NNPDF23_nlo',
        'CT10' : 'CT10nlo',
        'MSTW' : 'MSTW2008nlo',
        'NNPDF': 'NNPDF23_nlo'}
#-----------------------------------------------------------------------------
def decodeCommandLine():
    VERSION = '09-Nov-2014'
    USAGE = '''
    python plotSpectra.py [options] <workspace-file>

    options
       -m<model>      LL, RR, VV, AA, V-A, QCD     [default=ALL]
       -L<Lambda>                                  [default=20 (TeV)]
    '''
    parser = optparse.OptionParser(usage=USAGE,
                                   version=VERSION)

    parser.add_option('-L', '--Lambda',
                      action='store',
                      dest='Lambda',
                      type='string',
                      default='20',
                      help='CI mass scale')

    parser.add_option('-m', '--model',
                      action='store',
                      dest='model',
                      type='string',
                      default='QCD LL RR VV AA V-A',
                      help='model')    
        
    options, args = parser.parse_args()
    if len(args) == 0:
        print USAGE
        sys.exit(0)
    filenames = args

    # Set up model
    Lambda = float(options.Lambda)
    kappa = vector('double')(6, 0)
    kappa[0] =-1
    model = upper(options.model)
    if model[0] == '-':
        model = model[1:]
        key = model
        sign=-1
        interf = 'd'
    else:
        key = model
        sign= 1
        interf = 'c'
    if KAPPA.has_key(key):
        for ii in xrange(kappa.size()):
            kappa[ii] = sign*KAPPA[key][ii]

    return (filenames, Lambda, kappa, model, interf)    
#-----------------------------------------------------------------------------
def makePlot(context):
    model    = context.model
    Lambda   = context.Lambda
    kappa    = context.kappa
    interf   = context.interf
    hdata    = context.hdata
    pT       = context.pT
    qcdspectrum = context.qcdspectrum
    cispectrum  = context.cispectrum

    lam   = 1.0/Lambda**2
    pTmin = pT[0]
    pTmax = pT[-1]
    nbins = hdata.GetNbinsX()
    # --------------------------------------------------------
    smearing = getsmearing.findall(qcdspectrum[0].dirname())[0]
    pdfdir   = getpdfdir.findall(qcdspectrum[0].dirname())[0]
    if   smearing == 'JESJER':
        label  = 'with JES+JER uncert.'
        postfix= '_JES_JER'         
    elif smearing == 'PDF':
        label  = 'with PDF uncert.'
        postfix= '_PDF'
    else:
        label  = 'with JES+JER+PDF uncert.'
        postfix= '_JES_JER_PDF'
    
    # --------------------------------------------------------
    print "="*80
    addCI = model != 'QCD'
    if addCI:
        postfix += '_%s_%3.1f_%s' % (model, Lambda, interf[0])
    else:
        postfix += '_QCD'
    
    pcQCD   = PercentileCurve(nbins)
    pcQCDCI = PercentileCurve(nbins)
    hQCD    = []
    hQCDCI  = []
    
    for index, QCD in enumerate(qcdspectrum):
        hqcd = QCD()
        hutil.divideByWidth(hqcd) # convert to a density
        pcQCD.add(hqcd)
        hQCD.append(hqcd)

        if addCI:
            CI  = cispectrum[index]
            hci = CI(lam, kappa)
            hutil.divideByWidth(hci) # convert to a density
            hci.Add(hqcd)
            pcQCDCI.add(hci)
            hQCDCI.append(hci)

        if index % 100 == 0:
            print "%4d" % index, hqcd.GetName()
    
    # --------------------------------------------------------
    # plot spectrum
    # --------------------------------------------------------
    name = 'figs/%s/%s_d2sigma_dpTdy%s' % (pdfdir,
                                           pdfdir, postfix)
    cspect = TCanvas(name, name, 10, 10, 500, 500)
    x = map(lambda i: (pT[i+1]+pT[i])/2, range(nbins))

    # decide what spectrum is to be plotted
    if addCI:
        pc   = pcQCDCI
        hist = hQCDCI
    else:
        pc   = pcQCD
        hist = hQCD
        
    # get median QCD curve. we shall use it as our reference curve
    medQCD = pcQCD(0.5)
    median= hdata.Clone('median')
    median.SetLineColor(kBlue)
    for ii, c in enumerate(medQCD):
        median.SetBinContent(ii+1, c)
        median.SetBinError(ii+1, 0)
        
    # compute percentile spectra
    curve = []
    for p in PERCENT: curve.append( pc(p) )
    p95 = mkpline(x, curve[0], curve[-1], hdata, color=kGreen)
    p68 = mkpline(x, curve[1], curve[-2], hdata, color=kYellow)
    p50 = mkgraph(x, curve[2],
                  "Jet p_{T} (GeV)",
                  "d^{2}#sigma /dp_{T}dy (pb/GeV)",
                  pTmin, pTmax,
                  ymin=YMIN,
                  ymax=YMAX,
                  color=kRed,
                  lwidth=1)
    
    cspect.cd()
    gPad.SetLogy()
    hdata.Draw("ep")
    p95.Draw('f same')
    p68.Draw('f same')
    p50.Draw('c same')
    median.Draw('c same')
    hdata.Draw("ep same")
    
    scribe = addTitle('CMS Preliminary  #surds=8TeV CI Search L=19.3/fb',
                      0.035)
    scribe.vspace()
    if addCI:
        scribe.write("%s(#Lambda=%3.1f TeV) %s" % (model,
                                                   Lambda,
                                                   interf), 0.04)    
    scribe.write(PDFS[pdfdir], 0.08)
    scribe.write(label, 0.12)

    xp = 0.55
    yp = 0.55
    xw = 0.16
    yw = 0.22
    lg = mklegend(xp, yp, xw, yw)
    lg.AddEntry(hdata, 'data', 'p')
    lg.AddEntry(median, 'QCD', 'l')
    lg.AddEntry(p50, 'median', 'l')
    lg.AddEntry(p68, '68%s' % '%', 'f')
    lg.AddEntry(p95, '95%s' % '%', 'f')
    lg.Draw('same')
    
    cspect.Update()
    cspect.SaveAs('.pdf')

    xdata = hdata.Integral()
    xsect = median.Integral()
    scale = xdata / xsect
    print "data/theory: %10.3f" % scale

    # --------------------------------------------------------
    # now plot the ratio
    # --------------------------------------------------------    
    pc  = PercentileCurve(nbins)
    for h in hist:
        h.Divide(median)
        pc.add(h)
        
    median.Scale(scale)
    hratio = hdata.Clone('hratio')
    hratio.Divide(median)
    hratio.SetMinimum(0)
    hratio.SetMaximum(3)
    if addCI:
        title = 'data / (QCD + CI)_{NLO}'
    else:
        title = 'data / QCD_{NLO}'

    hratio.GetYaxis().SetTitle(title)

    cratio = TCanvas('figs/%s/%s_data_over_theory%s' % (pdfdir,
                                                        pdfdir, postfix),
                     '%s - spectrum-ratio' % pdfdir,
                     515, 10, 500, 500)
   
    curve = []
    for p in PERCENT: curve.append( pc(p) )
    p95r = mkpline(x, curve[0], curve[-1], hratio, color=kGreen)
    p68r = mkpline(x, curve[1], curve[-2], hratio, color=kYellow)
    p50r = mkgraph(x, curve[2],
                   "Jet p_{T} (GeV)", title,
                   pTmin, pTmax,
                   ymin=0,
                   ymax=3,
                   color=kRed,
                   lwidth=1)

    xx = array('d'); xx.append(pTmin); xx.append(pTmax)
    yy = array('d'); yy.append(1); yy.append(1)
    qcdr = TGraph(2, xx, yy)
    qcdr.SetLineWidth(2)
    qcdr.SetLineColor(kBlue)
        
    cratio.cd()
    gPad.SetLogy(kFALSE)
    hratio.Draw("ep")
    p95r.Draw('f same')
    p68r.Draw('f same')
    p50r.Draw('c same')
    qcdr.Draw('c same')
    hratio.Draw("ep same")
    
    scriber = addTitle('CMS Preliminary  #surds=8TeV CI Search L=19.3/fb',
                       0.035)
    scriber.vspace()
    if addCI:
        scriber.write("%s(#Lambda=%3.1f TeV) %s" % (model,
                                                    Lambda,
                                                    interf), 0.04)
    scriber.write(PDFS[pdfdir], 0.04)
    scriber.write(label, 0.04)

    lg = mklegend(xp+0.1, yp, xw, yw)
    lg.AddEntry(hratio, 'data', 'p')
    lg.AddEntry(qcdr, 'QCD', 'l')
    lg.AddEntry(p50r, 'median', 'l')
    lg.AddEntry(p68r, '68%s' % '%', 'f')
    lg.AddEntry(p95r, '95%s' % '%', 'f')
    lg.Draw('same')
    
    cratio.Update()
    cratio.SaveAs('.pdf')    

    sleep(5)           
#-----------------------------------------------------------------------------
def main():
    print "\n\t<=== plotSpectra.py ===>"
    setStyle()
 
    filenames, Lambda, kappa, model, interf = decodeCommandLine()
    
    # --------------------------------------------------------      
    # read data
    # --------------------------------------------------------    
    hdfile = TFile('../data/data_8TeV_L19.34.root')
    hdata  = hdfile.Get('hdata').Clone('data')
    hdata.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    hdata.GetYaxis().SetTitle('d^{2}#sigma /dp_{T}dy (pb/GeV)')
    hutil.divideByWidth(hdata)
    hdata.Scale(1.0/LUMI)
    pT = hutil.binlowedges(hdata)
    pT.push_back(2500)

    # --------------------------------------------------------        
    print "\nloading spectra..."
    # --------------------------------------------------------    
    addCI = model != "QCD"
    qcdspectrum = []
    cispectrum  = []
    hfile = []
    for filename in filenames:
        hfile.append( TFile(filename) )
        nspectra = min(QCDSpectrum.count(hfile[-1]),
                       CISpectrum.count(hfile[-1]))
        for ii in xrange(nspectra):
            jj = ii + 1
            
            qcd = hfile[-1].Get("QCD%3.3d" % jj)
            if qcd == None: break
            qcdspectrum.append(qcd)

            if addCI:
                ci = hfile[-1].Get("CI%3.3d" % jj)
                if ci == None: break
                cispectrum.append(ci)
                            
            if jj % 100 == 0:
                print "%4d %s %s" % (jj, qcd.dirname(), qcd.histname())
                if addCI:
                    print "%5s %s %s" % ('', ci.dirname(), ci.histname())
    # --------------------------------------------------------
    # now plot
    # --------------------------------------------------------
                                       
    context = Context()
    context.Lambda= Lambda
    context.kappa = kappa
    context.hdata = hdata
    context.interf= ''    
    context.pT = pT
    context.qcdspectrum = qcdspectrum
    context.cispectrum  = cispectrum

    models = split(model)
    print "\nplotting..."
    for model in models:
        context.model = model
        print "\n\t<=== %s ===>" % context.model
         
        if context.model == 'QCD':
            makePlot(context)
        else:
            for ii in xrange(kappa.size()):
                kappa[ii] = KAPPA[model][ii]
            context.kappa = kappa
            context.interf= 'constructive'
            makePlot(context)

            for ii in xrange(kappa.size()):
                kappa[ii] =-KAPPA[model][ii]
            context.kappa = kappa
            context.interf= 'destructive'                
            makePlot(context)
    
    gApplication.Run()        
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
