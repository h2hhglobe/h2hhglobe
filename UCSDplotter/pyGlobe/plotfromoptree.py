import ROOT
import array

def parseInputfiles(filename):
    samples = dict()
    file = open(filename)
    lines = file.readlines()
    file.close()
    for l in lines:
        if ("#" not in l and "typ" in l):
            items = l.split()
            itype = -1
            name = ""
            for item in items:
                if (item.find("typ=") is not -1):
                    itype = int(item.split("=")[1])
                elif (item.find("Nam=") is not -1):
                    name = str(item.split("=")[1])
                #else:
                #    print "Item",item,"not parsed."
            samples[itype] = name

    return samples

def inputfileTree(mysamples):

    inputfiletree = ROOT.TTree("inputfiles", "globe inputfiles provenance information");

    nfiles = array.array('i', [1]) #junk
    nindfiles = array.array('i', [len(mysamples)])  #junk
    intlumi = array.array('f', [12.])  #junk
    itype = array.array('i', nindfiles[0]*[0])
    infoind = array.array('i', nfiles[0]*[1])  #junk
    histoindfromfiles = array.array('i', nfiles[0]*[1])  #junk
    inshortnames = ROOT.TClonesArray("TObjString", len(mysamples))
    infilenames = ROOT.TClonesArray("TObjString", len(mysamples)) #junk
        
    inputfiletree.Branch("nfiles", nfiles, "nfiles/I");
    inputfiletree.Branch("nindfiles", nindfiles, "nindfiles/I");
    inputfiletree.Branch("intlumi",  intlumi, "intlumi/F");
    inputfiletree.Branch("itype", itype, "itype[nindfiles]/I");
    inputfiletree.Branch("histoind", histoindfromfiles, "histoindfromfiles[nfiles]/I");
    inputfiletree.Branch("infoind", infoind, "infoind[nindfiles]/I");
    inputfiletree.Branch("inshortnames", "TClonesArray", ROOT.AddressOf(inshortnames), 32000, 0)
    inputfiletree.Branch("infilenames", "TClonesArray", ROOT.AddressOf(infilenames), 32000, 0)

    for i,s in enumerate(mysamples):
        itype[i] = s
        temp = ROOT.TObjString()
        temp.SetString(mysamples[s])
        inshortnames[i] = temp
    
    inputfiletree.Fill()

    return inputfiletree

def plotvariableTree(allHistos):

    plotvartree = ROOT.TTree("plotvariables","globe plotvariables provenance information");

    Nvar = array.array('i', [len(allHistos.name)])
    histoncat = array.array('i', Nvar[0]*[0])
    typplotall = array.array('i', [0])  #junk
    doplot = array.array('i', Nvar[0]*[0]) #junk
    h2d = array.array('i', Nvar[0]*[0]) #junk
    typplot = array.array('i', Nvar[0]*[0]) #junk
    histoncatindtonames = array.array('i', Nvar[0]*[0]) #junk
    nbinsx = array.array('i', Nvar[0]*[0])
    nbinsy = array.array('i', Nvar[0]*[0])
    lowlim = array.array('f', Nvar[0]*[0])
    highlim = array.array('f', Nvar[0]*[0])
    lowlim2 = array.array('f', Nvar[0]*[0])
    highlim2 = array.array('f', Nvar[0]*[0])
    xaxislabels = ROOT.TClonesArray("TObjString", Nvar[0])
    yaxislabels = ROOT.TClonesArray("TObjString", Nvar[0])
    plotvarnames = ROOT.TClonesArray("TObjString", Nvar[0])
        
    plotvartree.Branch("Nvar", Nvar, "Nvar/I");
    plotvartree.Branch("typplotall", typplotall, "typplotall/I");
    plotvartree.Branch("doplot", doplot, "doplot[Nvar]/I");
    plotvartree.Branch("h2d", h2d, "h2d[Nvar]/I");
    plotvartree.Branch("typplot", typplot, "typplot[Nvar]/I");
    plotvartree.Branch("histoncat", histoncat, "histoncat[Nvar]/I");
    plotvartree.Branch("histoncatindtonames", histoncatindtonames, "histoncatindtonames[Nvar]/I");
    plotvartree.Branch("nbinsx", nbinsx, "nbinsx[Nvar]/I");
    plotvartree.Branch("nbinsy", nbinsy, "nbinsy[Nvar]/I");
    plotvartree.Branch("lowlim", lowlim, "lowlim[Nvar]/F");
    plotvartree.Branch("highlim", highlim, "highlim[Nvar]/F");
    plotvartree.Branch("lowlim2", lowlim2, "lowlim2[Nvar]/F");
    plotvartree.Branch("highlim2", highlim2, "highlim2[Nvar]/F");
    plotvartree.Branch("xaxislabels", "TClonesArray", ROOT.AddressOf(xaxislabels), 32000, 0)
    plotvartree.Branch("yaxislabels", "TClonesArray", ROOT.AddressOf(yaxislabels), 32000, 0)
    plotvartree.Branch("plotvarnames", "TClonesArray", ROOT.AddressOf(plotvarnames), 32000, 0)
 
    for i in xrange(Nvar[0]):
        h2d[i] = 0
        doplot[i] = 1
        typplot[i] = 0
        histoncat[i] = allHistos.ncat[i]
        histoncatindtonames[i] = -1
        nbinsx[i]   = allHistos.xbins[i]
        nbinsy[i]   = allHistos.ybins[i]
        lowlim[i]   = allHistos.xmin[i]
        highlim[i]  = allHistos.xmax[i]
        lowlim2[i]  = allHistos.ymin[i]
        highlim2[i] = allHistos.ymax[i]

        temp = ROOT.TObjString()
        #print allHistos.xaxis[i]
        temp.SetString(allHistos.xaxis[i])
        #print temp
        xaxislabels[i] = temp
        #print xaxislabels[i]

        temp.SetString(allHistos.yaxis[i])
        yaxislabels[i] = temp

        temp.SetString(allHistos.name[i])
        plotvarnames[i] = temp

    plotvartree.Fill()

    return plotvartree
  
class histoContainer:
    def __init__(self):
        self.histo_type = []
        self.ncat = []
        self.xbins = []
        self.ybins = []
        self.xmin, self.xmax = [], []
        self.ymin, self.ymax = [], []
        self.name = []
        self.xaxis = []
        self.yaxis = []
        self.vars = []
        self.catTypes =[]
        self.histos = dict()
        
    def createHistos(self, samples):
        for s in samples:
            for i in xrange(self.ncat[-1]):
                plotkey = self.name[-1]+"_cat"+str(i)+"_"+samples[s]
                if (self.histo_type[-1] == 0): # TH1F
                    self.histos[plotkey] = ROOT.TH1F(plotkey, plotkey, self.xbins[-1], self.xmin[-1], self.xmax[-1])
                if (self.histo_type[-1] == 1): # TH1I
                    self.histos[plotkey] = ROOT.TH1I(plotkey, plotkey, self.xbins[-1], self.xmin[-1], self.xmax[-1])

    def __str__(self):
        s = "typ:%s" % (str(self.histo_type[0]))
        s = s + " ncat:%s" % (str(self.ncat[0]))
        s = s + " name:%s" % (self.name[0])
        return s

    def fillHisto(self, n, cat, sample, var, w=1.):
        key = n+"_cat"+str(cat)+"_"+sample
        self.histos[key].Fill(var, w)

def parsePlotvariables(filename, samples):
    h = histoContainer()

    file = open(filename)
    lines = file.readlines()
    file.close()

    for l in lines:
      if ("#" not in l and "htyp" in l):
        items = l.split()
        for item in items:
            if ("default" in item or item == ""):
                continue
            elif (item.find("htyp=") is not -1):
                h.histo_type.append(int(item.split("=")[1]))
            elif item.find("ncat=") is not -1:
                h.ncat.append(int(item.split("=")[1]))
            elif item.find("xbins=") is not -1:
                h.xbins.append(int(item.split("=")[1]))
            elif item.find("ybins=") is not -1:
                h.ybins.append(int(item.split("=")[1]))
            elif item.find("xmin=") is not -1:
                h.xmin.append(float(item.split("=")[1]))
            elif item.find("xmax=") is not -1:
                h.xmax.append(float(item.split("=")[1]))
            elif item.find("ymin=") is not -1:
                h.ymin.append(float(item.split("=")[1]))
            elif item.find("ymax=") is not -1:
                h.ymax.append(float(item.split("=")[1]))
            elif item.find("name=") is not -1:
                h.name.append(str(item.split("=")[1]))
            elif item.find("xaxis=") is not -1:
                dummy = str(item.split("=")[1])
                #print dummy
                dummy = dummy.replace("@"," ")
                dummy = dummy.replace("%","#")
                #print dummy
                h.xaxis.append(dummy)
            elif item.find("yaxis=") is not -1:
                dummy = str(item.split("=")[1])
                #print dummy
                dummy = dummy.replace("@"," ")
                dummy = dummy.replace("%","#")
                #print dummy
                h.yaxis.append(dummy)
            elif item.find("name=") is not -1:
                h.name.append(str(item.split("=")[1]))
            elif item.find("var=") is not -1:
                h.vars.append(str(item.split("=")[1]))
            elif item.find("catType=") is not -1:
                h.catTypes.append(str(item.split("=")[1]))
            #else:
            #    print "Item",item,"not parsed."
        h.createHistos(samples)
        #print h.name
    return h
