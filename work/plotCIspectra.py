#!/usr/bin/env python
#-----------------------------------------------------------------------------
# plotCIspectra.py
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
YMAX    = 1.0e+08
#                        kappa
#                   1   2   3   4   5   6
KAPPA   = {'LL' : [-1,  0,  0,  0,  0,  0],
           'RR' : [ 0,  0,  0,  0, -1,  0],
           'VV' : [-1,  0, -2,  0, -1,  0],
           'AA' : [-1,  0, +2,  0, -1,  0],
           'V-A': [ 0,  0, -2,  0,  0,  0]}
MODEL   = ['LL', 'RR', 'VV', 'AA', 'V-A']

LAMBDA  = [10.0, 15.0, 20.0]
COLOR   = [kRed, kOrange+1, kGreen+1]
class Context:
    pass    
#-----------------------------------------------------------------------------
gSystem.Load("libCI.so")
LUMI = {'7':  5000.0,
        '8': 19710.0} # 1/pb

getpdfdir   = re.compile('(?<=fastNLO/).+(?=/[0-9])|'\
                         '(?<=fastCI/).+(?=/[0-9])')
getsmearing = re.compile('(?<=/[0-9][0-9][0-9]/).+')
PDFS = {'ALL'  : 'CT10nlo, MSTW2008nlo, NNPDF23_nlo',
        'CT10' : 'CT10nlo',
        'CTEQ6.6': 'CTEQ6.6',
        'MSTW' : 'MSTW2008nlo',
        'NNPDF': 'NNPDF23_nlo'}
#-----------------------------------------------------------------    
def check(o, message):
    if o == None:
        print "** %s **" % message
        sys.exit(0)
def waitForKey(message):
    print message
    sys.stdout.flush()
    raw_input("")
#-----------------------------------------------------------------
def decodeCommandLine():
    VERSION = '29-Nov-2014'
    USAGE = '''
    python plotCIspectra.py [options] <workspace-file>

    options
       -m<model>      LL, RR, VV, AA, V-A     [default=LL]
    '''
    parser = optparse.OptionParser(usage=USAGE,
                                   version=VERSION)

    parser.add_option('-m', '--model',
                      action='store',
                      dest='model',
                      type='string',
                      default='LL',
                      help='model')    
        
    options, args = parser.parse_args()
    if len(args) == 0:
        print USAGE
        sys.exit(0)
    filename = args[0]

    # Set up model
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

    return (filename, kappa, model, interf)    
#-----------------------------------------------------------------------------
def makePlots(context):
    model    = context.model
    Lambda   = context.Lambda
    kappa    = context.kappa
    interf   = context.interf
    hdata    = context.hdata
    pT       = context.pT
    is7TeV   = context.is7TeV
    lwidth   = context.lwidth
    color    = context.color
    
    qcdspectrum = context.qcdspectrum
    cispectrum  = context.cispectrum
    
    print "\n\t<=== %s ===> %s" % (context.model, context.interf)
    
    lam   = 1.0/Lambda**2
    pTmin = pT[0]
    pTmax = pT[-1]
    nbins = hdata.GetNbinsX()

    pcQCD   = PercentileCurve(nbins)
    pcQCDCI = PercentileCurve(nbins)
    hQCD    = []
    hQCDCI  = []
    
    for index, QCD in enumerate(qcdspectrum):
        hqcd = QCD()
        hutil.divideByWidth(hqcd) # convert to a density
        pcQCD.add(hqcd)
        hQCD.append(hqcd)

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

    x = map(lambda i: (pT[i+1]+pT[i])/2, range(nbins))
    pc   = pcQCDCI
    hist = hQCDCI
        
    # get median QCD curve. we shall use it as our reference curve
    medQCD = pcQCD(0.5)
    median = hdata.Clone('median')
    median.SetLineColor(kBlue)
    for ii, c in enumerate(medQCD):
        median.SetBinContent(ii+1, c)
        median.SetBinError(ii+1, 0)

    # hack to rescale 7TeV spectrum
    if is7TeV:
        hdata.Scale(median.Integral() / hdata.Integral())
                
    xdata = hdata.Integral()
    xsect = median.Integral()
    scale = xdata / xsect
    print "\n\t==> data/theory: %10.3f\n" % scale
    
    # compute median spectrum
    curve = pc(0.5)
    p50 = mkgraph(x, curve,
                  "Jet p_{T} (GeV)",
                  "d^{2}#sigma /dp_{T}dy (pb/GeV)",
                  pTmin, pTmax,
                  ymin=YMIN,
                  ymax=YMAX,
                  color=color,
                  lwidth=lwidth)

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
    hratio.SetMinimum(context.ymin)
    hratio.SetMaximum(context.ymax)
    if is7TeV:
        title = 'data / (QCD_{NLO} + CI_{LO})'
    else:
        title = 'data / (QCD + CI)_{NLO}'
    hratio.GetYaxis().SetTitle(title)
  
    curve = pc(0.5)
    p50r = mkgraph(x, curve,
                   "Jet p_{T} (GeV)", title,
                   pTmin, pTmax,
                   ymin=context.ymin,
                   ymax=context.ymax,
                   color=color,
                   lwidth=lwidth)

    xx = array('d'); xx.append(pTmin); xx.append(pTmax)
    yy = array('d'); yy.append(1); yy.append(1)
    qcdr = TGraph(2, xx, yy)
    qcdr.SetLineWidth(2)
    qcdr.SetLineColor(kBlue)
    return (hdata, median, p50, hratio, qcdr, p50r)
#-----------------------------------------------------------------------------
def drawPlots(context, plots):
    model    = context.model
    Lambda   = context.Lambda
    kappa    = context.kappa
    interf   = context.interf
    is7TeV   = context.is7TeV
    
    if is7TeV:
        smearing = ''
        pdfdir   = 'CTEQ6.6'
    else:
        smearing = getsmearing.findall(context.qcdspectrum[0].dirname())[0]
        pdfdir   = getpdfdir.findall(context.qcdspectrum[0].dirname())[0]

    if   smearing == 'NONE':
        label   = ''
        postfix = '_NONE' 
    elif smearing == 'JESJER':
        label   = 'JES+JER uncert.'
        postfix = '_JES_JER'         
    elif smearing == 'PDF':
        label   = 'PDF uncert.'
        postfix = '_PDF'
    else:
        label   = 'JES+JER+PDF uncert.'
        postfix = '_JES_JER_PDF'

    # --------------------------------------------------------
    postfix += '_%s_%s' % (model, interf[0])
    os.system('mkdir -p figs/%s' % pdfdir)
    
    context.specfile  = 'figs/%s/%s_d2sigma_dpTdy%s' % (pdfdir,
                                                        pdfdir,
                                                        postfix)
    context.ratiofile =  'figs/%s/%s_data_over_theory%s' % (pdfdir,
                                                            pdfdir,
                                                            postfix)
    print context.specfile

    # --------------------------------------------------------
         
    cspect = TCanvas(context.specfile,
                     context.specfile,
                     10, 10, 500, 500)

    hdata, median, p50, hratio, qcdr, p50r = plots[0]
    cspect.cd()
    gPad.SetLogy()
    hdata.Draw("ep")
    median.Draw('c same')

    xp = 0.23
    yp = 0.22
    xw = 0.22
    yw = 0.28
    lg = mklegend(xp, yp, xw, yw)
    lg.AddEntry(hdata,  'data', 'p')
    lg.AddEntry(median, 'QCD',  'l')

    p50 = []
    p50r= []
    for ii, (h, m, p, hr, mr, pr) in enumerate(plots):
        p.Draw('c same')
        p50.append( p )
        p50r.append( pr )
        lg.AddEntry(p50[-1], '#Lambda = %2.0f TeV' % LAMBDA[ii], 'l')
    hdata.Draw("ep same")

    energy = context.energy
    lumi   = context.lumi
    scribe = addTitle('CMS Preliminary  #surds=%sTeV CI Search L=%s/fb' % \
                      (energy, lumi),
                      0.035)
    scribe.vspace()
    scribe.write("%s %s" % (model, interf), 0.04)
    scribe.write(PDFS[pdfdir], 0.08)
    scribe.write(label, 0.12)

    lg.Draw('same')
    
    cspect.Update()
    cspect.SaveAs('.pdf')
    
    # --------------------------------------------------------
    cratio = TCanvas(context.ratiofile,
                     context.ratiofile,
                     515, 10, 500, 500)
    cratio.cd()
    gPad.SetLogy(kFALSE)
    hratio.Draw("ep")
    qcdr.Draw('c same')

    lgr = mklegend(xp, yp, xw, yw)
    lgr.AddEntry(hratio, 'data', 'p')
    lgr.AddEntry(qcdr, 'QCD', 'l')
    
    for ii, pr in enumerate(p50r):
        pr.Draw('c same')
        lgr.AddEntry(pr, '#Lambda = %2.0f TeV' % LAMBDA[ii], 'l')
    hratio.Draw("ep same")
    
    scriber = addTitle('CMS Preliminary  #surds=%sTeV CI Search L=%s/fb' % \
                       (energy, lumi),
                       0.035)
    scriber.vspace()
    scriber.write("%s %s" % (model, interf), 0.04)
    scriber.write(PDFS[pdfdir], 0.04)
    scriber.write(label, 0.04)

    lgr.Draw('same')
    
    cratio.Update()
    cratio.SaveAs('.pdf')
    waitForKey('\n==> Hit key to continue>> ')    
#-----------------------------------------------------------------------------
def main():
    print "\n\t<=== plotCIspectra.py ===>"
    setStyle()
 
    filename, kappa, model, interf = decodeCommandLine()

    context = Context()
    context.is7TeV = filename[:5] == 'CTEQ6'
    context.kappa  = kappa
    
    if context.is7TeV:
        context.energy = '7'
        context.lumi   = '5.0'
    else:
        context.energy = '8'
        context.lumi   = '19.7'
        
    # --------------------------------------------------------        
    print "\nloading workspace..."
    # --------------------------------------------------------    
    qcdspectrum = []
    cispectrum  = []
    hfile = TFile(filename)
    ws = hfile.Get('CI')
    check(ws, "can't get workspace CI")

    # --------------------------------------------------------      
    # get data
    # --------------------------------------------------------
    hdata = ws.obj('hdata')
    check(hdata, "can't get hdata")
    nbins   = hdata.GetNbinsX()
    hdata.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    hdata.GetYaxis().SetTitle('d^{2}#sigma /dp_{T}dy (pb/GeV)')
    hdata.SetMarkerSize(0.8)
    hdata.SetMarkerStyle(20)
    hutil.divideByWidth(hdata)
    hdata.Scale(1.0/LUMI[context.energy])
    pT = hutil.binlowedges(hdata)
    pT.push_back(pT[-1]+hdata.GetBinWidth(nbins))

    # --------------------------------------------------------      
    # get model
    # --------------------------------------------------------    
    pdf = ws.pdf('model')
    nspectra = pdf.size()
    print "number of spectra %d" % nspectra

    for ii in xrange(nspectra):
        qcd = pdf.QCD(ii)
        if qcd == None: break
        qcdspectrum.append(qcd)

        ci = pdf.CI(ii)
        if ci == None: break
        cispectrum.append(ci)

        if ii % 100 == 0:
            print "%4d %s %s" % (ii,
                                 qcd.dirname(),
                                 qcd.histname())
    
    # --------------------------------------------------------
    # now plot
    # --------------------------------------------------------
    context.hdata = hdata
    context.interf= ''    
    context.pT = pT
    context.qcdspectrum = qcdspectrum
    context.cispectrum  = cispectrum
    context.ymin = 0.8
    context.ymax = 1.2
    models = split(model)
    for model in models:
        context.model = model
                
        for interf, sign in [('constructive', 1),
                             ('destructive', -1)]:
            context.interf= interf
            for ii in xrange(kappa.size()):
                kappa[ii] = sign*KAPPA[model][ii]
            context.kappa = kappa
            
            plots = []
            for ii, Lambda in enumerate([10.0, 15.0, 20.0]):
                context.hdata  = hdata
                context.Lambda = Lambda
                context.color  = COLOR[ii]
                context.lwidth = 2
                plots.append( makePlots(context) )
            drawPlots(context, plots)    
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
