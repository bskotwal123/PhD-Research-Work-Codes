#!/usr/bin/env python
#-----------------------------------------------------------------
# File:        computeLimits.py
# Description: read CI 8 TeV workspace and compute limits
# Created:     02-Mar-2014 HBP
#              22-May-2014 HBP slighty modified version of program
#                          used to compute 7 TeV limits
#              09-Jun-2014 HBP add check function
#              15-Nov-2014 HBP use new pdf
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
#-----------------------------------------------------------------        
def makePlot(ws, likelihood, q,
             hname,
             color=kBlack,
             lstyle=1,
             lwidth=2,
             ymin=0.0,
             ymax=0.3):

    xmin   = q.getMin()
    xmax   = q.getMax()
    xbins  = q.getBins()
    xtitle = q.GetTitle()
    ytitle = likelihood.GetTitle()
    xstep  = (xmax-xmin)/xbins

    ## h = mkhist1(hname, xtitle, 'p(D | #lambda )',
    ##             xbins, xmin, xmax,
    ##             color=color, lstyle=lstyle)
    h = TH1D(hname, "", xbins, xmin, xmax)

    h.SetLineColor(color)
    h.SetLineStyle(lstyle)
    h.GetXaxis().SetTitle('#lambda (TeV^{-2})')
    h.GetXaxis().SetNdivisions(505)

    h.GetYaxis().SetTitle('p(#lambda | D)')
    h.GetYaxis().SetNdivisions(505)
    h.GetYaxis().SetTitleOffset(1.7)
    h.SetLineWidth(lwidth)

    for ii in xrange(xbins):
        x = xmin + (ii+0.5)*xstep
        q.setVal(x)
        y = likelihood.getVal()
        h.SetBinContent(ii+1, y)
        h.SetBinError(ii+1, 0)
    total = h.Integral()
    if total == 0:
        print "total is ZERO!"
        sys.exit(0)
    h.Scale(1.0/total)
    h.SetMinimum(ymin)
    h.SetMaximum(ymax)
    return h
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
         ./computeLimits.py <workspace-file>
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

    L  = 19.71    # lumi
    #L  = 5.00
    if filename[0:5] == 'CTEQ6':
        energy = '7'
        lumi = '%5.2f' %  5.0
    else:
        energy = '8'
        lumi   = '%5.2f' % L
        
    prefix = nameonly(filename)
    dname  = split(prefix, '_')[0]
    
    # --------------------------------------
    # get model etc.
    # --------------------------------------
    model = ws.pdf('model')
    Nset  = ws.set('Nset')
    poi   = ws.var('lambda')

    print "="*80
    print "computing Bayesian interval"
    print "input filename: %s" % filename
    print "workspace:      %s" % WSPACE
    print "models:         %s" % models
    print "="*80
      
    # --------------------------------------
    # create a graph of likelihood
    # --------------------------------------
    
    # create wrapper for model
    pdf   = PDFWrapper(model, Nset, poi)

    CL = 0.95
    model.setAsimov(True, L)
    model.setBinRange(0)
    ntrials = 5
    step    = 100
    for key in models:
        
        for sign in [1, -1]:
            if sign > 0:
                name = '%s_positive' % key
                fname = 'figs/%s/%s_likelihood_%s' % (dname, prefix, name)
                clike1 = TCanvas(fname, fname, 10, 10, 500, 500)
                clike = clike1
            else:
                name = '%s_negative' % key
                fname = 'figs/%s/%s_likelihood_%s' % (dname, prefix, name)
                clike2 = TCanvas(fname, fname, 515, 10, 500, 500)
                clike = clike2

            kappa = map(lambda x: sign*x, KAPPA[key])
            print "\n\tmodel: %5s\t%s" % (key, kappa)
            for ii in xrange(len(kappa)):
                vname = 'kappa%d' % ii
                if energy == '7':
                    ws.var(vname).setVal(-kappa[ii])
                else:
                    ws.var(vname).setVal(kappa[ii])
                
            # --------------------------------------
            # compute nominal limits
            # --------------------------------------
            try:
                del hnom
            except:
                pass
            
            model.setSize()     # use full sample of likelihoods           
            model.setNumber(0)  # use nominal cross section            
            hnom = makePlot(ws, model, poi, "hnom", color=kBlue,lstyle=2)
            swatch = TStopwatch()
            swatch.Start()
            model.setNumber(0)
            bayes = Bayes(pdf, poi.getMin(), poi.getMax())        
            limit = bayes.quantile(CL)
            if limit > 0:
                Limit1= 1.0/sqrt(limit)
            else:
                Limit1=-1.0
            print "\tLambda > %8.1f TeV @ %3.1f%s CL (stats. only)" % \
              (Limit1, 100*CL, '%')
            clike.cd()
            hnom.Draw('l')
            clike.Update()

            # --------------------------------------
            # compute nominal limits
            # --------------------------------------
            model.setNumber(-1) # include systematic uncertainties            
            Limit2 = []
            try:
                for h in havg:
                    del h
            except:
                pass
            havg = [None]*ntrials
            colors = [kBlue, kGreen, kOrange+1, kMagenta, kRed]
            for iii in xrange(ntrials):
                size = (iii+1)*step
                model.setSize(size)
                lwidth = 1
                if iii == ntrials-1: lwidth=3
                havg[iii]=makePlot(ws, model, poi, "havg%3.3d" % iii,
                                   color=colors[iii],
                                   lstyle=1,
                                   lwidth=lwidth)
                clike.cd()
                havg[iii].Draw('l same')            
                clike.Update()
                bayes.normalize()
                limit= bayes.quantile(CL)
                if limit > 0:
                    Limit2.append(1.0/sqrt(limit))
                    print "\tLambda > %8.1f TeV @ %3.1f%s CL (all uncert.)" % \
                    (Limit2[-1], 100*CL, '%'), size
                    print "\t\t\t==> real time: %8.3f s " % swatch.RealTime()
            # --------------------------------------
            # plot posterior density and limits
            # --------------------------------------
            clike.cd()
            scribe = addTitle('CMS Preliminary  '\
                              '#surds=%sTeV CI Search L=%s/fb' % \
                              (energy, lumi),
                              0.035)
            scribe.vspace()
            for iii in xrange(ntrials):
                size = (iii+1)*step
                scribe.write("%s(#kappa=%s) #Lambda > %3.1fTeV (%d)" % \
                            (key, kappa, Limit2[iii], size), 0.04)
                            
            clike.Update()
            clike.SaveAs('.pdf')
            
        sleep(5)
    #gApplication.Run()
#----------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print
    print "ciao!"
