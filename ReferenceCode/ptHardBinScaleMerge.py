import ROOT
import argparse
import ctypes
import os
import gc

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

###################################################################################
# Main function
def ptHardBinScaleMerge():

    inFileDir = "/home/alidock/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/embedding/"
    inFileBaseName = inFileDir + "AnalysisResultsEmbFullCh18r21May06_"
    taskName = "AliAnalysisTaskEmcalEmbeddingHelper_histos"

    PtHardBins = 20
    ptHardL = [5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235]
    ptHardH = [7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, -1]
    
    # Create histogram of NEvents accepted and NEvents acc+rej, as a function of pT-hard bin
    hNEventsAcc = ROOT.TH1F("hNEventsAcc", "hNEventsAccepted", PtHardBins+1, 0, PtHardBins+1)
    hNEventsTot = ROOT.TH1F("hNEventsTot", "hNEventsTotal", PtHardBins+1, 0, PtHardBins+1)
    nEventsAccSum = 0
    nEventsTotSum = 0
    for bin in range(1, PtHardBins):
        hNEventsAcc.GetXaxis().SetBinLabel(bin, "%d-%d" % (ptHardL[bin-1],ptHardH[bin-1]))
        hNEventsTot.GetXaxis().SetBinLabel(bin, "%d-%d" % (ptHardL[bin-1],ptHardH[bin-1]))
    for bin in range(1, PtHardBins):
        nNEventsTot= GetNEvents(inFileBaseName, taskName, bin, hNEventsTot, bAccEvents=False)
        nEventsAcc = GetNEvents(inFileBaseName, taskName, bin, hNEventsAcc, bAccEvents=True)
        nEventsAccSum += nEventsAcc
        nEventsTotSum += nNEventsTot

    nEventsAccAvg = nEventsAccSum/PtHardBins
    nNEventsTotAvg= nEventsTotSum/PtHardBins
    
    lScaleFactor, lPtHardEventNum, hXSecPerEvent, hNTrialsPerEvent, hScaleFactor = \
        eachPtHardBinScaleFact(inFileBaseName, taskName, PtHardBins, ptHardL, ptHardH, \
            hNEventsAcc, hNEventsTot, nEventsAccAvg, nNEventsTotAvg)
    
    # outputFileName = 'hoge.root'
    # outputFile = ROOT.TFile(outputFileName, "RECREATE")
    # lOMainTree = ROOT.TList()
# 
    # lOMainTree.Add(hXSecPerEvent)
    # lOMainTree.Add(hNTrialsPerEvent)
    # lOMainTree.Add(hScaleFactor)
    # lOMainTree.Add(hNEventsAcc)
    # lOMainTree.Add(hNEventsTot)
# 
    # outputFile.cd()
    # lOMainTree.Write('mainTree', 1)
    # outputFile.Close()

    return lScaleFactor, lPtHardEventNum, hXSecPerEvent, hNTrialsPerEvent, hScaleFactor, hNEventsAcc, hNEventsTot


###################################################################################
# Given event list name eventListName, pT-hard bin number, and histogram hNEvents of appropriate form, 
# this function fills the number of events 
# (accepted events only if bAcceptedEventsOnly=True, otherwise all events)
def GetNEvents(inFileBaseName, taskName, bin, hNEvents, bAccEvents = True):
    
    inputFile = inFileBaseName + "{0}.root".format(bin)
    f = ROOT.TFile(inputFile, "READ")
    
    lHEventInfo = f.Get(taskName)
    f.Close()

    nEvents = 0
    if not lHEventInfo:
        print("ERROR no lHEventInfo found")

    # Look for the EventCutOutput from AliEventCuts, and if it doesn't exist, look for histo fHistEventCount
    # eventCutList = lHEventInfo.FindObject("EventCutOutput")
    eventCutList = lHEventInfo.FindObject("EventCuts")

    bCutHist = False
    if eventCutList: bCutHist = True
    
    histBin = 0
    if bCutHist: 
        hNEventsPtHard = eventCutList.FindObject("fCutStats")
        if bAccEvents: histBin = hNEventsPtHard.GetXaxis().FindBin("All cuts")
        else: histBin = hNEventsPtHard.GetXaxis().FindBin("No cuts")
    else: 
        hNEventsPtHard = lHEventInfo.FindObject("fHistEventCount")
        histBin = 1

    nEvents = hNEventsPtHard.GetBinContent(histBin) #Bin All cuts
    if (not bAccEvents) and (not bCutHist):
        nEvents += hNEventsPtHard.GetBinContent(2)

    hNEvents.Fill(bin-0.5, nEvents)
    

    return nEvents


def eachPtHardBinScaleFact(inFileBaseName, taskName, PtHardBins, ptHardL, ptHardH,\
    hNEventsAcc, hNEventsTot, nEventsAccAvg, nNEventsTotAvg):

    print("ooo Determine Scale Factors from single pT Hard Bins")
    lScaleFactor = list()
    lPtHardEventNum = list()
    hXSecPerEvent = ROOT.TH1F("hXSecPerEvent", "hXSecPerEvent", PtHardBins+1, 0, PtHardBins+1)
    hNTrialsPerEvent = ROOT.TH1F("hNTrialsPerEvent", "hNTrialsPerEvent", PtHardBins+1, 0, PtHardBins+1)
    hScaleFactor = ROOT.TH1F("hScaleFactor", "hScaleFactor", PtHardBins+1, 0, PtHardBins+1)
    for bin in range(1, PtHardBins):
        # Label histograms
        hXSecPerEvent.GetXaxis().SetBinLabel(bin, "%d-%d" % (ptHardL[bin-1],ptHardH[bin-1]))
        hNTrialsPerEvent.GetXaxis().SetBinLabel(bin, "%d-%d" % (ptHardL[bin-1],ptHardH[bin-1]))
        hScaleFactor.GetXaxis().SetBinLabel(bin, "%d-%d" % (ptHardL[bin-1],ptHardH[bin-1]))

    # Extract cross sections from pT Hard bins
    for bin in range(1, PtHardBins):
        
        # Open input file and get relevant lists
        inputFile = inFileBaseName + "{0}.root".format(bin)
        f = ROOT.TFile(inputFile, "READ")
        print("ooo Pt-hard bin %d" % (bin))
        print("  o File: {}".format(inputFile))
        
        lHEventInfo = f.Get(taskName)
        f.Close()
        
        print("ooo Computing scaling factors with list: " + lHEventInfo.GetName())

        # hXsecPtHard      = lHEventInfo.FindObject("hXsec")
        # hTrialsPtHard    = lHEventInfo.FindObject("hNtrials")
        hXsecPtHard      = lHEventInfo.FindObject("fHistXsection")
        hTrialsPtHard    = lHEventInfo.FindObject("fHistTrials")
        # hNEventsTotal    = f.Get("hNEventsTot") # Get NEvents histos from file, since we undo the entry after writing it
        # hNEventsAccepted = f.Get("hNEventsAcc")

        # Compute: scale factor = xsec per event / trials per event
        nEventsTot = hNEventsTot.GetBinContent(bin)
        nEventsAcc = hNEventsAcc.GetBinContent(bin)

        # xsec = hXsecPtHard.GetBinContent(1) / hXsecPtHard.GetEntries() 
        # trials = 1.*hTrialsPtHard.GetBinContent(1) / nEventsTot
        xsec = hXsecPtHard.GetBinContent(bin+1) / hXsecPtHard.GetEntries() 
        trials = 1.*hTrialsPtHard.GetBinContent(bin+1) / nEventsTot
        scaleFactor = xsec/trials
        # also scale to account that there are different number of events in each Pt-hard bin
        if nEventsAcc>0: eventScaleFactorAcc = nEventsAccAvg/nEventsAcc 
        else: eventScaleFactorAcc = 0
        if nEventsTot>0: eventScaleFactorTot = nNEventsTotAvg/nEventsTot
        else: eventScaleFactorTot = 0

        hXSecPerEvent.Fill(bin-0.5, xsec)
        hNTrialsPerEvent.Fill(bin-0.5, trials)
        hScaleFactor.Fill(bin-0.5, eventScaleFactorAcc*scaleFactor)

        print("xSec: {0}, trials: {1}, histVal: {2}, nEventsTot: {3}, nEventsTot: {4}, eventScalFactor: {5}, scaleFactor: {6}"\
            .format(xsec, trials, hTrialsPtHard.GetBinContent(1), nEventsTot, nEventsAcc, eventScaleFactorAcc, scaleFactor))

        lScaleFactor.append(eventScaleFactorAcc*scaleFactor)
        lPtHardEventNum.append(nEventsTot)
        f.Close()

    return lScaleFactor, lPtHardEventNum, hXSecPerEvent, hNTrialsPerEvent, hScaleFactor


if __name__ == "__main__":
    ptHardBinScaleMerge()