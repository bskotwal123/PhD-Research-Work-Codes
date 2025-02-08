#!/usr/bin/env python
from ROOT import *
gSystem.Load('libCI.so')

h = hutil()
h.hists.push_back(TH1D('h', '', 10, 0, 1))
print len(h.hists)

