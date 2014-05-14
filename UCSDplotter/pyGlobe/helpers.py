import ROOT
import sys

def getHighestSumPtPairs(nPairs, pairTypes, id1, id2, iso1, iso2, sumpts):
    bestSumpts = [0, 0, 0]
    bestPairs = [-1, -1, -1]
    for p in xrange(nPairs):
        goOn = False
        if (pairTypes[p] == 0):
            if (id1[p] ==11 and id2[p] == 11):
                goOn = True
        elif (pairTypes[p] == 1):
            if (id1[p] > 0.92 and id2[p] > 0.92 and
                iso1[p] < 0.1 and iso2[p] < 0.1):
                goOn = True
        else:
            if (id1[p] > 0.92 and iso1[p] < 0.1 and id2[p] == 11):
                goOn = True
        if (not goOn):
            continue
        if (sumpts[p] > bestSumpts[pairTypes[p]]):
            bestSumpts[pairTypes[p]] = sumpts[p]
            bestPairs[pairTypes[p]] = p

    return bestPairs

def getAllPairs(nPairs, pairTypes, id1, id2, iso1, iso2, eta1, eta2, phi1, phi2, sumpts):
    if (nPairs > 1):
        sortedPairs = sorted(zip(sumpts, range(nPairs)))
        sortedIndices = zip(*sortedPairs)[1]
        leptonList = []
        for p in sortedIndices:
            #goOn = False
            #if (pairTypes[p] == 0):
            #    if (id1[p] > 0 and id2[p] > 0):
            #        goOn = True
            #elif (pairTypes[p] == 1):
            #    if (id1[p] > 0.1 and id2[p] > 0.1 and
            #        iso1[p] < 0.5 and iso2[p] < 0.5):
            #        goOn = True
            #else:
            #    if (id1[p] > 0.1 and iso1[p] < 0.5 and id2[p] > 0):
            #        goOn = True
            #if (not goOn):
            #    continue
            leptonList.append((p,(eta1[p], phi1[p]),(eta2[p], phi2[p])))

        to_remove = []
        for i in xrange(len(leptonList)-1):
            for j in xrange(i+1, len(leptonList)):
                if (leptonList[i][1] == leptonList[j][1] or leptonList[i][1] == leptonList[j][2] 
                    or leptonList[i][2] == leptonList[j][1] or leptonList[i][2] == leptonList[j][2]):
                    to_remove.append(leptonList[j])
        for i in to_remove: 
            if (i in leptonList):
                leptonList.remove(i)
            
        finalList = [p[0] for p in leptonList]
        if (len(finalList) > 1):
            return finalList
        else:
            return []

def tightElectronCharge(ch11, ch21, ch31, ch12, ch22, ch32):
    ch1 = 19
    ch2 = 21
    if ((ch11 == ch21) and (ch11 == ch31)):
        ch1 = ch11
    if ((ch12 == ch22) and (ch12 == ch32)):
        ch2 = ch12

    return (ch1 == ch2)
