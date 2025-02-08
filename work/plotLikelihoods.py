#!/usr/bin/env python
#-----------------------------------------------------------------
# File:        plotLikelihoods.py
# Description: read CI 8 TeV workspace and plot likelihoods
# Created:     02-Mar-2014 HBP
#-----------------------------------------------------------------
import os,sys,re
from array import array
from time import sleep
from histutil import *
from string import *
from ROOT import *
#-----------------------------------------------------------------
WSPACE   = 'CI'
DEBUG    = 0
#                        kappa
#                   0   1   2   3   4   5
KAPPA   = {'LL' : [-1,  0,  0,  0,  0,  0],
           'RR' : [ 0,  0,  0,  0, -1,  0],
           'VV' : [-1,  0, -2,  0, -1,  0],
           'AA' : [-1,  0, +2,  0, -1,  0],
           'V-A': [ 0,  0, -2,  0,  0,  0]}
MODEL   = ['LL', 'RR', 'VV', 'AA', 'V-A']
#-----------------------------------------------------------------
def check(o, message):
    if o == None:
        print "** %s **" % message
        sys.exit(0)
def waitForKey():
    print "\t** hit any key to continue**"
    sys.stdout.flush()
    raw_input("")        
#-----------------------------------------------------------------        
def makePlot(ws, likelihood, x, hname, number, kappa,
             color=kBlack, lstyle=1, lwidth=1):

    likelihood.setNumber(number)
    xmin   = x.getMin()
    xmax   = x.getMax()
    xbins  = x.getBins()
    xstep  = (xmax-xmin)/xbins
    
    hl = TH1D(hname, "", xbins, xmin, xmax)
    hl.SetLineColor(color)
    hl.SetLineStyle(lstyle)
    hl.GetXaxis().SetTitle('#lambda (TeV^{-2})')
    hl.GetXaxis().SetNdivisions(505)

    hl.GetYaxis().SetTitle('p(#lambda | D)')
    hl.GetYaxis().SetNdivisions(505)
    hl.GetYaxis().SetTitleOffset(1.7)
    
    hl.SetLineWidth(lwidth)
    hl.SetMinimum(0.0)
    hl.SetMaximum(0.15)

        
    xtitle = 'Jet p_{T} (GeV)'
    ytitle = 'data / theory'

    hQCD = likelihood.QCD(0)()
    hqcd = likelihood.QCD(number)()

    hist = []
    print
    print hname
    for ii in xrange(xbins):
        l = xmin + ii*xstep
        x.setVal(l)
        p = likelihood.getVal()
        print "%10.5f\t%10.3e" % (l, p)
        
        hl.SetBinContent(ii+1, p)
        hl.SetBinError(ii+1, 0)        
    
        hci  = likelihood.CI(number)(l, kappa)
        h    = hqcd.Clone(hname+"_%3.3d" % ii)
        h.Add(hci)
        h.Divide(hQCD)
        h.SetLineColor(color+ii%5)
        h.SetLineStyle(lstyle)
        h.GetXaxis().SetTitle(xtitle)
        h.GetXaxis().SetNdivisions(505)

        h.GetYaxis().SetTitle(ytitle)
        h.GetYaxis().SetNdivisions(505)
        h.GetYaxis().SetTitleOffset(1.7)
        h.SetLineWidth(lwidth)

        h.SetMinimum(0.00)            
        h.SetMaximum(2.00)
        hist.append(h)

    return (hl, hist)
#-----------------------------------------------------------------
def setValues(ws, d, scale=1.0):
    if type(d) == type(""): d = ws.set(d)
    iterator = d.iterator()
    ii = 0
    while 1:
        o = iterator.Next()
        if  o == None: break
        name = o.GetName()
        ws.var(name).setVal(d[name].getVal()/scale)
        print "%5d\t%10.2f" % (ii, ws.var(name).getVal())
        ii += 1
#-----------------------------------------------------------------
def toVector(ws, setname):
    data = ws.set(setname)
    if data == None: return None
    iterator = data.iterator()
    vdata = vector('double')()
    while 1:
        o = iterator.Next()
        if  o == None: break
        name = o.GetName()
        vdata.push_back( ws.var(name).getVal() )
        #print '%10s %f' % (name, vdata.back())
    return vdata
#-----------------------------------------------------------------
def main():
    # --------------------------------------
    # load various codes needed for the
    # calculatioms
    # --------------------------------------
    gSystem.Load('libCI.so')

    # set up some standard graphics style
    setStyle()

    # --------------------------------------
    # load workspace into memory
    # --------------------------------------
    argv = sys.argv[1:]
    if len(argv) == 0:
        print '''
    Usage:
         ./plotLikelihoods.py <workspace-file> [models]
        '''
        sys.exit(0)
        
    filename = argv[0]
    wfile = TFile(filename)
    if not wfile.IsOpen():
        print "** can't open %s" % filename
        sys.exit(0)
    ws = wfile.Get(WSPACE)
    check(ws, "can't access workspace %s" % WSPACE)

    if len(argv) > 1:
        models = argv[1:]
    else:
        models = MODEL
    
    if filename[0:5] == 'CTEQ6':
        energy = '7'
        lumi = '5'
    else:
        energy = '8'
        lumi   = '19.7'
        
    prefix = nameonly(filename)
    dname  = split(prefix, '_')[0]
    
    # --------------------------------------
    # get model etc.
    # --------------------------------------
    model = ws.pdf('model')
    poi   = ws.var('lambda')
    poi.setRange(0, 0.015)
    poi.setBins(100)    
    kappa = vector('double')(6, 0)
        

    print "="*80
    print "input filename: %s" % filename
    print "workspace:      %s" % WSPACE
    print "models:         %s" % models
    print "="*80
      
    # --------------------------------------
    # create graphs of likelihoods
    # --------------------------------------
    xl = array('d'); xl.append(500); xl.append(2500)
    yl = array('d'); yl.append(1); yl.append(1)
    gline = TGraph(2, xl, yl)

    xmin   = poi.getMin()
    xmax   = poi.getMax()
    xbins  = poi.getBins()
    xstep  = (xmax-xmin)/xbins
    
    hl = TH1D("havg", "", xbins, xmin, xmax)
    hl.SetLineColor(kRed)
    hl.SetLineStyle(1)
    hl.GetXaxis().SetTitle('#lambda (TeV^{-2})')
    hl.GetXaxis().SetNdivisions(505)

    hl.GetYaxis().SetTitle('p(#lambda | D)')
    hl.GetYaxis().SetNdivisions(505)
    hl.GetYaxis().SetTitleOffset(1.7)
    
    hl.SetLineWidth(1)
    hl.SetMinimum(0.0)
    hl.SetMaximum(0.15)

    
    ll = 0.0
    L  = 19.71    # lumi
    sign = 1
    for key in models:
        if KAPPA.has_key(key):
            for ii in xrange(kappa.size()):
                kappa[ii] = sign*KAPPA[key][ii]
                vname = 'kappa%d' % ii
                ws.var(vname).setVal(kappa[ii])
                print 'kappa[%d] = %10.1f' % (ii, kappa[ii])

        hl.Reset()
        
        name  = key
        fname = 'figs/%s/%s_likelihood_%s' % (dname, prefix, name)
        clike = TCanvas(fname, fname, 10, 10, 500, 500)
        clike.cd()
        scribe1 = addTitle('CMS Preliminary  '\
                          '#surds=%sTeV CI Search L=%s/fb' % \
                          (energy, lumi),
                          0.035)

        
        fname = 'figs/%s/%s_spectra_%s' % (dname, prefix, name)
        cspec = TCanvas(fname, fname, 515, 10, 500, 500)
        cspec.cd()
        gPad.SetGridy()
        scribe2 = addTitle('CMS Preliminary  '\
                          '#surds=%sTeV CI Search L=%s/fb' % \
                          (energy, lumi),
                          0.035)

        # --------------------------------------
        # plot posterior density
        # --------------------------------------
        try:
            del hist
            del hlike
        except:
            pass

        model.setAsimov(True, L, ll)
        model.setBinRange(0)
        size  = model.size()
        color = [kRed, kOrange, kYellow+2, kGreen, kMagenta, kBlue]
        option= 'l'
        jj =-1
        for ii in xrange(1, size):
            jj += 1
            kk = jj % len(color)
            nn = ii+1
            hname = "h%3.3d" % nn
            
            hlike, hist = makePlot(ws, model, poi, hname, nn, kappa,
                                       color=color[kk], lstyle=1)
            hl.Add(hlike)

            if ii % 50 == 0:
                hlike.Scale(1.0/hlike.Integral())            
                clike.cd()
                hlike.Draw('l')
                clike.Update()
            
                cspec.cd()
                hist[0].Draw('l')
                gline.Draw('l same')
                for h in hist[1:]:
                    h.Draw('l same')
                cspec.Update()
                sleep(0.01)
        clike.cd()
        hl.Scale(1.0/hl.Integral())
        hl.SetMinimum(0.0)
        hl.SetMaximum(0.15)    
        hl.Draw('l')
        clike.Update()
        sleep(5)
    #gApplication.Run()
#----------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print
    print "ciao!"
