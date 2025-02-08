#!/usr/bin/env python
#-------------------------------------------------------------------------------
import os, sys, re
from string import *
from glob import glob
from ROOT import *
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
def xsection(cixsec, l, kappa, which=1):
    bi  = vector('double')(6)
    ai  = vector('double')(6)
    aig = vector('double')(6)

    bij = vector('double')(9)
    aij = vector('double')(9)
    aijg= vector('double')(9)    

    bi4 = vector('double')(4)
    ai4 = vector('double')(4)
    ai4g= vector('double')(4)

    print "aij[3] = ", aij[3]
    
    xsec = [0]*9
    for ii in xrange(9):
        for jj in xrange(6):
            bi[jj] = cixsec.bi(ii, jj, which)
            aig[jj]= cixsec.aig(ii,jj, which)            
            ai[jj] = cixsec.ai(ii, jj, which)

        for jj in xrange(9):
            bij[jj] = cixsec.bij(ii, jj, which)
            aijg[jj]= cixsec.aijg(ii,jj, which)
            aij[jj] = cixsec.aij(ii, jj, which)


        for jj in xrange(4):
            bi4[jj] = cixsec.bi4(ii, jj, which)
            ai4g[jj]= cixsec.ai4g(ii,jj, which)            
            ai4[jj] = cixsec.ai4(ii, jj, which)

        xsec[ii] = CIXsection.evaluate(l, kappa,
                                       bi, aig, ai,
                                       bij,aijg,aij,
                                       bi4,ai4g,ai4)
    return xsec
#-------------------------------------------------------------------------------
def main():
    gSystem.Load('libCI.so')

    Lambda = 10.0
    l = 1.0/Lambda**2
    kappa = vector('double')(6, 0)
    kappa[0] =-1

    histdir = '../CT10/000'
    dirname = 'CT10_000'
    bins = [507, 1497]
    ibin = [20, 36]
    
    cifnames = []
    xsec = {}
    for b in bins:
        cifnames.append('%s/ci%4.4d.txt' % (histdir, b))
        
        fname = 'xsec%4.4d.txt' % b
        print "\n\t==> %s" % fname
        recs = map(split, open(fname).readlines()[6:])
        xsec[b] = []
        for index in MURMUF:
            xsec[b].append((recs[index][1], recs[index][0],
                            recs[index][2], recs[index][3]))

            print xsec[b][-1]            

    # --------------------------------------------------------------------
    cilo = []
    cinlo= []
    print '\n\t...loading LO CI histograms...'
    swatch = TStopwatch()
    for postfix in HISTNAMES:
        hname = 'lo%s' % postfix
        cilo.append( CISpectrum(histdir,  hname) )
        hname = 'nlo%s' % postfix
        cinlo.append( CISpectrum(histdir, hname) )
    print "time/constructor: %10.3fs" % (swatch.RealTime()/(2*len(HISTNAMES)))

    CI = cilo[0]
    h = CI(l, kappa)
    
    for ii in xrange(h.GetNbinsX()):
        print "%5d\t%10.0f\t%10.3e" % (ii+1,
                                       h.GetBinLowEdge(ii+1),
                                       h.GetBinContent(ii+1)*
                                       h.GetBinWidth(ii+1)/2)
    # --------------------------------------------------------------------

 
    out = open('testcixsec.txt', 'w')

    for index, filename in enumerate(cifnames):
        cixsec = CIXsection(filename)
        loxsec = cixsec(l, kappa, 0); loxsec = map(lambda i: loxsec[i], MURMUF)
        nloxsec= cixsec(l, kappa, 1); nloxsec= map(lambda i:nloxsec[i], MURMUF)
        xsection = xsec[bins[index]]
        
        #lo  = xsection(cixsec, l, kappa, 0)
        #nlo = xsection(cixsec, l, kappa, 1)
        
        record = "==> %s" % filename
        
        out.write('%s\n' % record)     
        print record
        
        record = " %6s %6s %10s %10s %10s %10s %10s %10s" % \
          ('mur', 'muf',
           'LO (pb)',
           'LO (CIX)',
           'LO (CIS)',
           'NLO (pb)',
           'NLO (CIX)',
           'NLO (CIS)')
        out.write('%s\n' % record)
        print record
        
        for ii in xrange(len(loxsec)):
            hlo  = cilo[ii](l, kappa)
            xlo = hlo.GetBinContent(ibin[index])*\
              hlo.GetBinWidth(ibin[index])/2
              
            hnlo = cinlo[ii](l, kappa)
            xnlo = hnlo.GetBinContent(ibin[index])*\
              hnlo.GetBinWidth(ibin[index])/2

            mur, muf, lo, nlo = xsection[ii]
            fmt = " %6s %6s" + 2*" %10s %10.3e %10.3e"
            record = fmt % \
              (mur, muf,
               lo,
               loxsec[ii],
               xlo,
               nlo,
               nloxsec[ii],
               xnlo)

            out.write('%s\n' % record)
            print record
    out.close()

    swatch = TStopwatch()
    N = 10000
    CI = cinlo[0]
    for ii in xrange(N):
        h = CI(l, kappa)
    t = 1000*swatch.RealTime()/N
    print "time/CI calculation: %10.2f ms" % t
#------------------------------------------------------------------------------ 
main()

