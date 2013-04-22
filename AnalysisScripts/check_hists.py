#!/usr/bin/env python
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--infile",dest="infile",help="File to search")
parser.add_option("-g","--grep",dest="grep",default=[],action="append",help="Grep string (can take multiple arguments)")
(options,args)=parser.parse_args()

import ROOT as r

tf = r.TFile(options.infile)

for key in tf.GetListOfKeys():
  if 'th1f' in key.GetName():
    for grep in options.grep:
      if grep in key.GetName():
        th1f = tf.Get(key.GetName())
        print '%50s -- %10d -- %10.24f'%(th1f.GetName(),th1f.GetEntries(),th1f.Integral())
