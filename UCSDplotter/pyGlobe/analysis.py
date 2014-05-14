import helpers
import ROOT
from helpers import getHighestSumPtPairs, tightElectronCharge, getAllPairs
from  makeWS import WSProducer
import plotfromoptree, table

class Analysis:
    def __init__(self, options):
        self.options = options
        self.counter = []

    def Run(self):
        self.initialize()
        self.loadTree()
        self.loop()
        self.finalize()

    def finalize(self):
        self.wsProducer.finalize()
        self.wsProducer.saveWS()

        self.counter.Print(self.samples)

        output = ROOT.TFile(self.options.output, "recreate")
        for h in self.allHistos.histos.values():
            h.Write()
        
        self.inputfiletree.Write()
        self.plotvariabletree.Write()
        output.Close()

    def initialize(self):
        print "Parsing inputfiles..."
        self.samples = plotfromoptree.parseInputfiles(self.options.inputfile)
        print "Parsing plotvariables..."
        self.allHistos = plotfromoptree.parsePlotvariables(self.options.plotvariables, self.samples)
        print "Writing inputfiles Tree..."
        self.inputfiletree = plotfromoptree.inputfileTree(self.samples)
        print "Writing plotvariables Tree..."
        self.plotvariabletree = plotfromoptree.plotvariableTree(self.allHistos)

        self.wsProducer = WSProducer()
        if (self.options.mode == "highestSumpt"):
            self.counter = table.table(2, 9)
            self.wsProducer.prepareDataSets(3)
        elif (self.options.mode == "allPairs"):
            self.counter = table.table(2, 1)
            self.wsProducer.prepareDataSets(3)

    def loadTree(self):
        self.file = ROOT.TFile(self.options.input)
        self.tree = self.file.Get("opttree")

        self.tree.SetBranchStatus("*",0)
        self.tree.SetBranchStatus("itype",1)
        self.tree.SetBranchStatus("weight",1)
        self.tree.SetBranchStatus("pairs", 1)
        self.tree.SetBranchStatus("mass",1)
        self.tree.SetBranchStatus("type",1)
        self.tree.SetBranchStatus("cat",1)
        self.tree.SetBranchStatus("sumpt", 1)
        self.tree.SetBranchStatus("ch1_1", 1)
        self.tree.SetBranchStatus("ch2_1", 1)
        self.tree.SetBranchStatus("ch3_1", 1)
        self.tree.SetBranchStatus("ch1_2", 1)
        self.tree.SetBranchStatus("ch2_2", 1)
        self.tree.SetBranchStatus("ch3_2", 1)
        self.tree.SetBranchStatus("id1", 1)
        self.tree.SetBranchStatus("id2", 1)
        self.tree.SetBranchStatus("iso1", 1)
        self.tree.SetBranchStatus("iso2", 1)
        self.tree.SetBranchStatus("met",1)
        self.tree.SetBranchStatus("njets30",1)

        if (self.options.mode == "highestSumpt"):
            self.tree.SetBranchStatus("njets20",1)
            self.tree.SetBranchStatus("btag", 1)
        
        if (self.options.mode == "allPairs"):
            self.tree.SetBranchStatus("eta1", 1)
            self.tree.SetBranchStatus("eta2", 1)
            self.tree.SetBranchStatus("phi1", 1)
            self.tree.SetBranchStatus("phi2", 1)

        self.entries = self.tree.GetEntries()

    def loop(self):
        for z in xrange(self.entries):
            self.tree.GetEntry(z)
    
            if (self.options.mode == "highestSumpt"):
                self.singlePairs()
            elif (self.options.mode == "allPairs"):
                self.doublePairs()

    def singlePairs(self):
        itype = self.tree.itype
        if (itype not in self.samples.keys()):
            return

        pairs = self.tree.pairs
        pairTypes = self.tree.type
        cats = self.tree.cat
        weight = self.tree.weight
        masses = self.tree.mass
        id1 = self.tree.id1
        id2 = self.tree.id2
        iso1 = self.tree.iso1
        iso2 = self.tree.iso2
        
        for p in getHighestSumPtPairs(pairs, pairTypes, id1, id2, iso1, iso2, self.tree.sumpt):
            if (p == -1):
                continue

            self.counter.Fill(0, itype, pairTypes[p]*3+cats[p], weight)
            nbtags = 0
            btag = self.tree.btag
            for j in xrange(self.tree.njets20):
                self.allHistos.fillHisto("jetet", pairTypes[p]*3+cats[p], self.samples[itype], self.tree.jetet[j], weight)
                self.allHistos.fillHisto("btag",  pairTypes[p]*3+cats[p], self.samples[itype], self.tree.btag[j], weight)
                if (btag[j] > 0.7):
                    nbtags += 1
            self.allHistos.fillHisto("nbtag", pairTypes[p]*3+cats[p], self.samples[itype], nbtags, weight)
                    
            if (not tightElectronCharge(self.tree.ch1_1[p], self.tree.ch2_1[p], self.tree.ch3_1[p], self.tree.ch1_2[p], self.tree.ch2_2[p], self.tree.ch3_2[p])):
                continue
            
            if (masses[p] > 12.):
                self.allHistos.fillHisto("mass_nocut", pairTypes[p]*3+cats[p], self.samples[itype], masses[p], weight)
                self.allHistos.fillHisto("met", pairTypes[p]*3+cats[p], self.samples[itype], self.tree.met, weight)
            
                self.allHistos.fillHisto("njet20", pairTypes[p]*3+cats[p], self.samples[itype], self.tree.njets20, weight)
                self.allHistos.fillHisto("njet30", pairTypes[p]*3+cats[p], self.samples[itype], self.tree.njets30, weight)

                if (self.tree.met < 40.):
                    continue

                if (self.tree.njets30 < 2):
                    continue

                self.counter.Fill(1, itype, pairTypes[p]*3+cats[p], weight)
                self.allHistos.fillHisto("mass", pairTypes[p]*3+cats[p], self.samples[itype], masses[p], weight)
                self.allHistos.fillHisto("id1",  pairTypes[p]*3+cats[p], self.samples[itype], id1[p],    weight)
                self.allHistos.fillHisto("id2",  pairTypes[p]*3+cats[p], self.samples[itype], id2[p],    weight)
                self.allHistos.fillHisto("iso1", pairTypes[p]*3+cats[p], self.samples[itype], iso1[p],   weight)
                self.allHistos.fillHisto("iso2", pairTypes[p]*3+cats[p], self.samples[itype], iso2[p],   weight)
        
                self.wsProducer.fillDataset(itype, pairTypes[p], cats[p], masses[p], weight)

    def doublePairs(self):
        itype = self.tree.itype
        if (itype not in self.samples.keys()):
            return

        pairs = self.tree.pairs
        if (pairs < 2):
            return
        pairTypes = self.tree.type
        #cats = self.tree.cat
        weight = self.tree.weight
        masses = self.tree.mass
        id1 = self.tree.id1
        id2 = self.tree.id2
        iso1 = self.tree.iso1
        iso2 = self.tree.iso2

        for p in getAllPairs(pairs, pairTypes, id1, id2, iso1, iso2, self.tree.eta1, self.tree.eta2, self.tree.phi1, self.tree.phi2, self.tree.sumpt):
            if (p == -1):
                continue
            #print list(self.tree.type)
            #print list(pairTypes)
                
            self.allHistos.fillHisto("type", 0, self.samples[itype], pairTypes[p],    1)
            self.allHistos.fillHisto("id1",  0, self.samples[itype], id1[p],    weight)
            self.allHistos.fillHisto("id2",  0, self.samples[itype], id2[p],    weight)
            self.allHistos.fillHisto("iso1", 0, self.samples[itype], iso1[p],   weight)
            self.allHistos.fillHisto("iso2", 0, self.samples[itype], iso2[p],   weight)

            self.counter.Fill(0, itype, 0, weight)
            if (not tightElectronCharge(self.tree.ch1_1[p], self.tree.ch2_1[p], self.tree.ch3_1[p], self.tree.ch1_2[p], self.tree.ch2_2[p], self.tree.ch3_2[p])):
                continue
            
            if (masses[p] > 8.):
                if (self.tree.met < 0.):
                    continue

                if (self.tree.njets30 < 0):
                    continue

                self.counter.Fill(1, itype, 0, weight)
                self.allHistos.fillHisto("mass", 0, self.samples[itype], masses[p], weight)
        
                self.wsProducer.fillDataset(itype, 1, 1, masses[p], weight)



