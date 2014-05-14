class table:
    tables = []
    nCheckPoints = 0
    nCats = 0

    def __init__(self, nCheckPoints, nCats):
        self.nCheckPoints = nCheckPoints
        self.nCats = nCats
        for i in xrange(nCheckPoints):
            temp = []
            for j in xrange(nCats):
                temp.append(dict())
            self.tables.append(temp)
    
    def Fill(self, checkPoint, itype, cat, weight):
        if (itype in self.tables[checkPoint][cat].keys()):
            self.tables[checkPoint][cat][itype] += weight
        else:
            self.tables[checkPoint][cat][itype] = weight

    def Print(self, samples, checkPoint = -1):
        if (checkPoint != -1):
            for c in xrange(self.nCats):
                nbkg = [0 for i in xrange(checkPoint)]
                print "CAT: ", c
                print "___________"
                keys = sorted(self.tables[checkPoint][c].keys())
                for k in keys:
                    print "%20s %.2f"%(samples[k], self.tables[checkPoint][c][k])
                    if (k > 0):
                        nbkg[checkPoint] += self.tables[checkPoint][c][k]
                print "%20s"%("TOT. Bkg"),
                for cp in  xrange(self.nCheckPoints):
                    print "%10.2f"%(nbkg[cp]),
                print ""
        else:
            for c in xrange(self.nCats):
                nbkg = [0 for i in xrange(self.nCheckPoints)]
                print "CAT: ", c
                print "___________"
                keys = sorted(self.tables[0][c].keys())
                for k in keys:
                    print "%20s"%(samples[k]),
                    for cp in xrange(self.nCheckPoints):
                        if (k in self.tables[cp][c].keys()):
                            print "%10.2f"%(self.tables[cp][c][k]),
                            if (k > 0):
                                nbkg[cp] += self.tables[cp][c][k]
                        else:
                            print "      0.00",
                    print ""
                print "%20s"%("TOT. Bkg"),
                for cp in  xrange(self.nCheckPoints):
                    print "%10.2f"%(nbkg[cp]),
                print ""
