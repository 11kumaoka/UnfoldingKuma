import ROOT
from ROOT import TH1F

ROOT.gROOT.SetBatch(True)

def calcPtHardScalFactor(geneEventHistList, ptHardBin):
    nEventsAccSum = 0
    nEventsTotSum=0
    # nNEventsTot= GetNEvents(eventListName, bin, hNEventsTot, verbose, bAcceptedEventsOnly=False)
    # nEventsAcc = GetNEvents(eventListName, bin, hNEventsAcc, verbose, bAcceptedEventsOnly=True)
    # nEventsAccSum += nEventsAcc
    # nEventsTotSum += nNEventsTot
    # nEventsAccAvg = nEventsAccSum/PtHardBins
    # nNEventsTotAvg= nEventsTotSum/PtHardBins
    
    # hXsecPtHard      = geneEventHistList.FindObject("fHistXsection")
    hXsecPtHard      = geneEventHistList[1]
    # hTrialsPtHard    = geneEventHistList.FindObject("hHistTrials")
    hTrialsPtHard    = geneEventHistList[2]
    # Get NEvents histos from file, since we undo the entry after writing it
    # hNEventsTotal    = geneEventHistList.FindObject("fHistEventCount") 
    hNEventsTotal    = geneEventHistList[3]
    # hNEventsAccepted = f.Get("hNEventsAcc")

    # Compute: scale factor = xsec per event / trials per event
    nEventsTot = hNEventsTotal.GetEntries()
    nEventsAcc = hNEventsTotal.GetBinContent(ptHardBin+1)

    #entries in this case are number of files that were merged together
    xsec = hXsecPtHard.GetBinContent(ptHardBin+1) / hXsecPtHard.GetEntries() 
    #print("ooo Test entries: {0}".format(hXsecPtHard.GetEntries()))
    trials = 1.*hTrialsPtHard.GetBinContent(ptHardBin+1) / nEventsTot

    print(ptHardBin)
    scaleFactor = xsec/trials
    print('xSec = {0}, trial = {1}, nEvents = {2}, scaleFactor = {3}'\
        .format(xsec, trials, nEventsTot, scaleFactor))

    # if nEventsAcc>0:
        # eventScaleFactorAcc = nEventsAccAvg/nEventsAcc 
        # also scale to account that there are different number of events in each Pt-hard bin
    # else:
        # eventScaleFactorAcc = 0
    # if nEventsTot>0:
        # eventScaleFactorTot = nNEventsTotAvg/nEventsTot
    # else:
        # eventScaleFactorTot = 0
      # return scaleFactor * eventScaleFactorAcc
    
    
    return scaleFactor 