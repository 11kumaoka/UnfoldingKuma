import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TMath
import argparse
import ctypes
import os
import gc
from array import array
import numpy as np

import sys

import genePlotSets
import histSetting
import plotPerformanceHists
import PtRangeList

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

# centBinRangeList = [[1,2],[2,3],[3,5],[5,7],[7,9],[9, 11],[11, 13],[13,15],[15,17],[17,19]]

# numOfCentBin = 10
# centRangeList = [[0,5],[5,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80],[80,90]]
numOfCentBin = 5

lCentLabel = ['0-5','5-10','10-30','30-50','50-80']
rawJetCentBinLoopRange = [[0,1],[1,2],[2,4],[4,6],[6,9]]
embCentBinRange = [[0,1],[1,2],[2,3],[3,4],[4,5]]
centRangeList = [[0,5],[5,10],[10,30],[30,50],[50,90]]

ptHardL = [5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235]
ptHardH = [7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, -1]

lEPLabel = ['OutOfPlane', 'InPlane', 'Inclusive']

iniPtHardBin = 1 # default 1
PtHardBins = 21
# PtHardBins = 2

bRawJetExtract = 1

variation = 1
minPtGen = 10.
maxPtGen = 200.
minPtDet = 20. 
maxPtDet = 150.

# lrefPtHardScalFactor = [1.9044557347701444e-05, 9.653419162066905e-06] #??? test
# lrefPtHardNumOfEvent = [6635787.0, 6707788.0] #??? test

lrefPtHardScalFactor = [2.446091855262593e-07, 7.314975239763274e-08, 3.362231880594244e-08, 1.1895882634924634e-08, 4.040376437588744e-09, 8.080455194155069e-09, 4.750910913607525e-10, 1.566097040702298e-10, 3.4940423913705284e-10, 2.118136932704947e-11, 8.003818748583827e-12, 2.8357605981006887e-12, 1.3703850028545194e-12, 6.368826505468165e-13, 3.1137657459266875e-13, 2.221801809420425e-13, 8.647496708263527e-14, 2.7811870500499066e-14, 3.492944863521209e-14, 6.539761692063409e-14]

lrefPtHardNumOfEvent = [5883650.0, 6097914.0, 6526103.0, 6487913.0, 6542589.0, \
    6618306.0, 6659974.0, 6606850.0, 6600430.0, 6473709.0, 6788644.0, 6612948.0, \
        6467705.0, 6215084.0, 6635787.0, 6707788.0, 6710764.0, 5544633.0, 6355388.0, 6619566.]


# LHC18r
# lrefPtHardScalFactor = \
#     [2.456280125649383e-07, 7.61698249503861e-08, 3.22581122679987e-08, 1.1812252486320773e-08, \
#     4.049254578128142e-09, 1.6319721765403964e-09, 4.65763556935813e-10, 1.6123117645106004e-10,\
#         6.392822378921412e-11, 2.2751243702969846e-11, 8.132893034348747e-12,\
#         2.9160375006758812e-12, 1.483583932062433e-12, 7.274705162337206e-13,\
#             3.486126907680189e-13, 1.1722360279760496e-12, 9.228858393189566e-14, \
#                 4.5892010095811983e-14, 2.908187970843322e-14, 2.538214345770791e-13]
# lrefPtHardNumOfEvent = \
#     [298012544.0, 286666656.0, 300784576.0, 299847456.0, 298053152.0, 285056032.0, \
#     292419200.0, 288783872.0, 293465568.0, 284348704.0, 297642656.0, 296179264.0, \
#         290634720.0, 283121344.0, 285828544.0, 114284304.0, 290955040.0, 298294720.0,\
#             275080032.0, 104309320.0]


###################################################################################
# Main function (0/1/2/3/4, "LHC18q/LHC18r", 0/5/7, "98%: ''/94%: 'TrackEff094'"
def MergeScalePtHardPlots(centKind, LHCPeriod, leadingTrackPtCut,  diffSys):
# def main(centKind, LHCPeriod, leadingTrackPtCut, diffSys):
    
    inRawJetFile = '~/ALICE/cernbox/SWAN_projects/outputFiles/'+LHCPeriod+'/pass3/Ch/RawJet/'
    if (diffSys != 'V0C') and (diffSys != 'V0A') and (diffSys != 'BKGV2'): inRawJetFile += 'TrainAnaResTreeBKGWay.root'
    if (diffSys == 'BKGV2'): inRawJetFile += 'TrainAnaResTreeBKGWayV2.root'
    else:  inRawJetFile += 'TrainAnaResTreeBKGWayDiffV0.root'

    # inEmbFileDir = "./"
    # inEmbFileDir = "/Volumes/KumaSSD/TrainOutput/+LHCPeriod+/Embedding/8870/"
    # inEmbFileDir = '~/ALICE/cernbox/SWAN_projects/outputFiles/'+LHCPeriod+'/pass3/Ch/Embedding/Train/8870/'
    inEmbFileDir = '~/ALICE/cernbox/SWAN_projects/outputFiles/'+LHCPeriod+'/pass3/Ch/Embedding/Train/8922/'

    inEmbFileBaseName = inEmbFileDir + "AnalysisResultsPtHard"
    # outEmbFileDir = './'
    outEmbFileDir = '~/ALICE/cernbox/SWAN_projects/outputFiles/'+LHCPeriod+'/pass3/Ch/Embedding/'

    outEmbFile = outEmbFileDir + 'EmbedPtHardScaledResults'\
        +'_TrackPtCut'+str(leadingTrackPtCut)+'_'+diffSys+'_CentBin'+str(centKind)+'.root'
    # outEmbFile = outEmbFileDir + 'hoge.root'

    embHelperTaskName = "AliAnalysisTaskEmcalEmbeddingHelper_histos"
    # jetTaskName = 'AliAnalysisTaskEmbeddingJetWithEP_1_new_histos'
    jetTaskName = 'AliAnalysisTaskEmbeddingJetWithEP_R02PtCut' +str(leadingTrackPtCut)+ diffSys + '_histos'
    

    # Create histogram of NEvents accepted and NEvents acc+rej, as a function of pT-hard bin
    hNEventsAcc = ROOT.TH1F("hNEventsAcc", "hNEventsAccepted", PtHardBins+1, iniPtHardBin-1, PtHardBins+1)
    hNEventsTot = ROOT.TH1F("hNEventsTot", "hNEventsTotal", PtHardBins+1, iniPtHardBin-1, PtHardBins+1)
    nEventsAccSum = 0
    nEventsTotSum = 0
    for bin in range(iniPtHardBin, PtHardBins):
        hNEventsAcc.GetXaxis().SetBinLabel(bin, "%d-%d" % (ptHardL[bin-1],ptHardH[bin-1]))
        hNEventsTot.GetXaxis().SetBinLabel(bin, "%d-%d" % (ptHardL[bin-1],ptHardH[bin-1]))
    for bin in range(iniPtHardBin, PtHardBins):
        nNEventsTot = GetNEvents(inEmbFileBaseName, embHelperTaskName,\
            bin, hNEventsTot, bAccEvents=False)
        nEventsAcc  = GetNEvents(inEmbFileBaseName, embHelperTaskName,\
            bin, hNEventsAcc, bAccEvents=True)
        
        nEventsAccSum += nEventsAcc
        nEventsTotSum += nNEventsTot

    nEventsAccAvg = nEventsAccSum/(PtHardBins-iniPtHardBin)
    nNEventsTotAvg= nEventsTotSum/(PtHardBins-iniPtHardBin)
    
    lScaleFactor, lPtHardEventNum, hXSecPerEvent, hNTrialsPerEvent, hScaleFactor, hPtHardSpectOrig, hPtHardSpectScaled = \
        eachPtHardBinScaleFact(inEmbFileBaseName, embHelperTaskName, PtHardBins, ptHardL, ptHardH, \
            hNEventsAcc, hNEventsTot, nEventsAccAvg, nNEventsTotAvg)
    
    # return lScaleFactor, lPtHardEventNum, hXSecPerEvent, hNTrialsPerEvent, hScaleFactor, hNEventsAcc, hNEventsTot
    print('lScaleFactor')
    print(lScaleFactor)
    print('lPtHardEventNum')
    print(lPtHardEventNum)

    outRootFile = ROOT.TFile(outEmbFile, 'RECREATE')
    lMainTree = ROOT.TList()
    OFileStracture(lMainTree, numOfCentBin)
    
    lMainTree.Add(hXSecPerEvent)
    lMainTree.Add(hNTrialsPerEvent)
    lMainTree.Add(hScaleFactor)
    lMainTree.Add(hNEventsAcc)
    lMainTree.Add(hNEventsTot)
    lMainTree.Add(hPtHardSpectOrig)
    lMainTree.Add(hPtHardSpectScaled)

    print('== s == ExtractDataRawJet ===========')
    ExtractDataRawJet(inRawJetFile, leadingTrackPtCut, diffSys, centKind, lMainTree)
    print('== s == ScaleJetHists ===========')
    ScaleJetHists(lMainTree, lScaleFactor, PtHardBins, centKind, inEmbFileBaseName, jetTaskName, leadingTrackPtCut, outEmbFileDir)
    
    print('== s == PlotPerformanceHistEachCent ===========')
    for epBin in range(0, 2):
        PlotPerformanceHistEachCent(lMainTree[epBin], centKind, epBin)

        for centBin in range(embCentBinRange[centKind][0], embCentBinRange[centKind][1]):
            label = 'Cent'+str(centBin)
            histName = 'hRM_forUnfold_' + 'Cent{0}'.format(centBin)
            hOriginRM = lMainTree[epBin][centBin].FindObject(histName)

            EditRMForUF(lMainTree[epBin][centBin], hOriginRM, centBin, epBin, label)
    
    # print(len(lMainTree))
    # for centBin in range(embCentBinRange[centKind][0], embCentBinRange[centKind][1]):
    #     for histBin in range(0, len(lMainTree[0][centBin])):
    #         tempHist = lMainTree[0][centBin][histBin].Clone()
    #         # lMainTree[2][centBin].Add(lMainTree[0][centBin][histBin])
    #         lMainTree[2][centBin].Add(tempHist)
    #         lMainTree[2][centBin][histBin].Add(lMainTree[1][centBin][histBin])
    #         lMainTree[2][centBin].ls()
    # # print(len(lMainTree[0]))
    # for histBin in range(numOfCentBin, len(lMainTree[0])):
    #     tempHist = lMainTree[0][histBin].Clone()
    #     # lMainTree[2].Add(lMainTree[0][histBin])
    #     lMainTree[2].Add(tempHist)
    #     lMainTree[2][histBin].Add(lMainTree[1][histBin])

    outRootFile.cd()
    lMainTree.Write('mainTree', 1)
    outRootFile.Close()
    print('root ' + outEmbFile)


###################################################################################
def ExtractDataRawJet(inRawJetFile, leadingTrackPtCut, diffSys, centKind, lMainTree):
    print(diffSys)
    print(inRawJetFile)
    f = ROOT.TFile(inRawJetFile, "READ")
    
    diffSysName = ''
    if (diffSys == 'V0C') or (diffSys == 'V0A'): diffSysName = diffSys
    elif (diffSys == 'BKGNoFit'): diffSysName = 'kNoFit'
    elif (diffSys == 'BKGV2'): diffSysName = 'kV2'
    mainJetTaskName = 'AliAnalysisTaskRawJetWithEP_R02PtCut'+str(leadingTrackPtCut)+diffSysName+'_histos'
    lMainHJetDists = f.Get(mainJetTaskName)
    f.ls()
    print(mainJetTaskName)
    jetTaskName = 'Jet_AKTChargedR020_tracks_pT0150_pt_scheme'
    lHJetDists = lMainHJetDists.FindObject(jetTaskName)
    f.Close()
    
    if not lHJetDists:
        print("ERROR no "+jetTaskName+" found")

    fNormalisationHist = lMainHJetDists.FindObject('fNormalisationHist')
    lMainTree.Add(fNormalisationHist)
    hCentrality_raw = lMainHJetDists.FindObject('Centrality_raw')
    lMainTree.Add(hCentrality_raw)
    hNEventsCentrality = lMainHJetDists.FindObject('fHistCentrality')
    lMainTree.Add(hNEventsCentrality)

    lHJetDistsInc = lHJetDists.FindObject('Inclusive')
    if not lHJetDistsInc:
        print("ERROR no Inclusive found")
    
    jetHistName = 'hJetCorrPtVsEP2AngleVsCent'
    hRawJet3D = lHJetDistsInc.FindObject(jetHistName)
    if not hRawJet3D:
        print("ERROR no hJetCorrPtVsEP2AngleVsCent found")
    for epBin in range(0, 2):

        hRawJet = ROOT.TH1D()
        for centBin in range(rawJetCentBinLoopRange[centKind][0], rawJetCentBinLoopRange[centKind][1]):
            if bRawJetExtract == 1:
                lHJetDistsEP = lHJetDists.FindObject(lEPLabel[epBin])
                if not lHJetDistsEP: print("ERROR no Inclusive found")
                histName = "hJetCorrPtLocal_"+str(centBin)
                tempHRawJet = lHJetDistsEP.FindObject(histName)
                if centBin == rawJetCentBinLoopRange[centKind][0]: hRawJet = tempHRawJet
                else: hRawJet.Add(tempHRawJet)
                print('centBin'+ str(centBin))
                # lMainTree[epBin][centBin].Add(tempHRawJet)
                
            elif bRawJetExtract == 2:
                histName = 'hRawJet3D_Cent' + str(centBin) + '_EP' + str(epBin)
                hRawJet3DCP = hRawJet3D.Clone(histName)
                hRawJet3DCP.GetZaxis().SetRangeUser(centRangeList[centBin][0], centRangeList[centBin][1])
                hRawJet2D = hRawJet3DCP.Project3D("xy")
                if(epBin==0):
                    histName = 'hRawJet2DCP1_Cent' + str(centBin) + '_EP' + str(epBin)
                    hRawJet2DCP1 = hRawJet2D.Clone(histName)
                    hRawJet2DCP1.GetXaxis().SetRangeUser(0, TMath.Pi()/4)
                    hRawJet1 = hRawJet2DCP1.ProjectionY()
                    histName = 'hRawJet2DCP2_Cent' + str(centBin) + '_EP' + str(epBin)
                    hRawJet2DCP2 = hRawJet2D.Clone(histName)
                    hRawJet2DCP2.GetXaxis().SetRangeUser(TMath.Pi()*3/4, TMath.Pi()*5/4)
                    hRawJet2 = hRawJet2DCP2.ProjectionY()
                    histName = 'hRawJet2DCP3_Cent' + str(centBin) + '_EP' + str(epBin)
                    hRawJet2DCP3 = hRawJet2D.Clone(histName)
                    hRawJet2DCP3.GetXaxis().SetRangeUser(TMath.Pi()*7/4, TMath.Pi()*2)
                    hRawJet3 = hRawJet2DCP3.ProjectionY()

                    hRawJet1.Add(hRawJet2)
                    hRawJet1.Add(hRawJet3)

                    histName = 'hJetCorrPtLocal_' + str(centBin)
                    hRawJet1.SetName(histName)

                    if centBin == rawJetCentBinLoopRange[centKind][0]: hRawJet = hRawJet1
                    else: hRawJet.Add(hRawJet1)
                    # lMainTree[epBin][centBin].Add(hRawJet1)
                else:
                    histName = 'hRawJet2DCP1_Cent' + str(centBin) + '_EP' + str(epBin)
                    hRawJet2DCP1 = hRawJet2D.Clone(histName)
                    hRawJet1 = hRawJet2DCP1.ProjectionY()
                    hRawJet2DCP1.GetXaxis().SetRangeUser(TMath.Pi()/4, TMath.Pi()*3/4)
                    histName = 'hRawJet2DCP2_Cent' + str(centBin) + '_EP' + str(epBin)
                    hRawJet2DCP2 = hRawJet2D.Clone(histName)
                    hRawJet2DCP2.GetXaxis().SetRangeUser(TMath.Pi()*3/4, TMath.Pi()*7/4)
                    hRawJet2 = hRawJet2DCP2.ProjectionY()

                    histName = 'hJetCorrPtLocal_' + str(centBin)
                    hRawJet1.SetName(histName)
                    hRawJet1.Add(hRawJet2)
                    
                    if centBin == rawJetCentBinLoopRange[centKind][0]: hRawJet = hRawJet1
                    else: hRawJet.Add(hRawJet1)
                    # lMainTree[epBin][centBin].Add(hRawJet1)

        histName = "hJetCorrPtLocal_"+str(centKind)
        hRawJet.SetName(histName)
        lMainTree[epBin][centKind].Add(hRawJet)
        
###################################################################################
# Given event list name eventListName, pT-hard bin number, and histogram hNEvents of appropriate form, 
# this function fills the number of events 
# (accepted events only if bAcceptedEventsOnly=True, otherwise all events)
def GetNEvents(inEmbFileBaseName, taskName, bin, hNEvents, bAccEvents = True):
    
    inEmbFile = inEmbFileBaseName + "{0}.root".format(bin)
    f = ROOT.TFile(inEmbFile, "READ")
    
    lHEventInfo = f.Get(taskName)
    f.Close()
    
    nEvents = 0
    if not lHEventInfo:
        print("ERROR no lHEventInfo found")

    # Look for the EventCutOutput from AliEventCuts, and if it doesn't exist, look for histo fHistEventCount
    eventCutList = lHEventInfo.FindObject("EventCuts")

    bCutHist = False
    if eventCutList: bCutHist = True
    
    histBin = 0
    # if bCutHist: 
    #     hNEventsPtHard = eventCutList.FindObject("fCutStats")
    #     if bAccEvents: histBin = hNEventsPtHard.GetXaxis().FindBin("All cuts")
    #     else: histBin = hNEventsPtHard.GetXaxis().FindBin("No cuts")
    # else: 
    #     hNEventsPtHard = lHEventInfo.FindObject("fHistEventCount")
    #     histBin = 1
    # if (not bAccEvents) and (not bCutHist):
    #     nEvents += hNEventsPtHard.GetBinContent(2)

    hNEventsPtHard = lHEventInfo.FindObject("fHistEventCount")
    histBin = 1

    nEvents = hNEventsPtHard.GetBinContent(histBin) #Bin All cuts

    hNEvents.Fill(bin-0.5, nEvents)
    
    return nEvents


def eachPtHardBinScaleFact(inEmbFileBaseName, taskName, PtHardBins, ptHardL, ptHardH,\
    hNEventsAcc, hNEventsTot, nEventsAccAvg, nNEventsTotAvg):

    print("ooo Determine Scale Factors from single pT Hard Bins")
    lScaleFactor = list()
    lPtHardEventNum = list()
    hXSecPerEvent = ROOT.TH1F("hXSecPerEvent", "hXSecPerEvent", PtHardBins+1, iniPtHardBin-1, PtHardBins+1)
    hNTrialsPerEvent = ROOT.TH1F("hNTrialsPerEvent", "hNTrialsPerEvent", PtHardBins+1, iniPtHardBin-1, PtHardBins+1)
    hScaleFactor = ROOT.TH1F("hScaleFactor", "hScaleFactor", PtHardBins+1, iniPtHardBin-1, PtHardBins+1)
    hPtHardSpectScaled = ROOT.TH1F()
    hPtHardSpectOrig = ROOT.TH1F()

    for bin in range(iniPtHardBin, PtHardBins):
        # Label histograms
        hXSecPerEvent.GetXaxis().SetBinLabel(bin, "%d-%d" % (ptHardL[bin-1],ptHardH[bin-1]))
        hNTrialsPerEvent.GetXaxis().SetBinLabel(bin, "%d-%d" % (ptHardL[bin-1],ptHardH[bin-1]))
        hScaleFactor.GetXaxis().SetBinLabel(bin, "%d-%d" % (ptHardL[bin-1],ptHardH[bin-1]))

    # Extract cross sections from pT Hard bins
    for bin in range(iniPtHardBin, PtHardBins):
        
        # Open input file and get relevant lists
        inEmbFile = inEmbFileBaseName + "{0}.root".format(bin)
        f = ROOT.TFile(inEmbFile, "READ")
        print("ooo Pt-hard bin %d" % (bin))
        print("  o File: {}".format(inEmbFile))
        
        lHEventInfo = f.Get(taskName)
        f.Close()
        
        print("ooo Computing scaling factors with list: " + lHEventInfo.GetName())

        # hXsecPtHard      = lHEventInfo.FindObject("hXsec")
        # hTrialsPtHard    = lHEventInfo.FindObject("hNtrials")
        hXsecPtHard      = lHEventInfo.FindObject("fHistXsection")
        hTrialsPtHard    = lHEventInfo.FindObject("fHistTrials")

        # Compute: scale factor = xsec per event / trials per event
        nEventsTot = hNEventsTot.GetBinContent(bin)
        nEventsAcc = hNEventsAcc.GetBinContent(bin)

        xsec = hXsecPtHard.GetBinContent(bin+1) / hXsecPtHard.GetEntries() 
        trials = 1.*hTrialsPtHard.GetBinContent(bin+1) / nEventsTot
        # xsec = hXsecPtHard.GetBinContent(bin+1) * hXsecPtHard.GetEntries() 
        # trials = hTrialsPtHard.GetBinContent(bin+1)

        print('xsec: '+str(xsec))
        print('trials: '+str(trials))
        
        scaleFactor = xsec/trials # Check ????

        # scaleFactor = xsec/(trials*nEventsAcc) # Check ????
        # scaleFactor = lrefPtHardScalFactor[bin-1] / lrefPtHardNumOfEvent[bin-1]# Check ????

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
        # lScaleFactor.append(nEventsAccAvg*scaleFactor)
        # lScaleFactor.append(scaleFactor)
        lPtHardEventNum.append(nEventsTot)

        origHPtHardSpectScaled = lHEventInfo.FindObject("fHistPtHard")

        cpHPtHardSpectOrig = origHPtHardSpectScaled.Clone()
        cpHPtHardSpectOrig.SetName("hPtHardSpectOrig".format(bin))
        if bin == 1: hPtHardSpectOrig = cpHPtHardSpectOrig
        else : hPtHardSpectOrig.Add(cpHPtHardSpectOrig)

        cpHPtHardSpectScaled = origHPtHardSpectScaled.Clone()
        cpHPtHardSpectScaled.SetName("hPtHardSpectScaled".format(bin))
        cpHPtHardSpectScaled.Scale(scaleFactor)
        if bin == 1: hPtHardSpectScaled = cpHPtHardSpectScaled
        else : hPtHardSpectScaled.Add(cpHPtHardSpectScaled)

        f.Close()

    return lScaleFactor, lPtHardEventNum, hXSecPerEvent, hNTrialsPerEvent, hScaleFactor, hPtHardSpectOrig, hPtHardSpectScaled


def ScaleJetHists(lMainTree, lScaleFactor, PtHardBins, centKind, inEmbFileBaseName, jetTaskName, leadingTrackPtCut, outputDir):
    
    lAllScaledGenJetPtDist = list()
    lAllScaledDetJetPtDist = list()
    for epBin in range(0, 2):
        lAllScaledGenJetPtDistForEP = list()
        lAllScaledDetJetPtDistForEP = list()
        lAllScaledGenJetPtDist.append(lAllScaledGenJetPtDistForEP)
        lAllScaledDetJetPtDist.append(lAllScaledDetJetPtDistForEP)
    
    for ptHardBin in range(iniPtHardBin, PtHardBins):
        inEmbFile = inEmbFileBaseName + "{0}.root".format(ptHardBin)
        print('  READ ' + inEmbFile)
        f = ROOT.TFile(inEmbFile, "READ")
        lHEventInfo = f.Get(jetTaskName)
        f.Close()
        
        if not lHEventInfo: print('ERROR no '+jetTaskName+' found')

        # Look for the EventCutOutput from AliEventCuts, and if it doesn't exist, look for histo fHistEventCount
        lHybRawJet = lHEventInfo.FindObject("hybridRawJet")
        lParRawJet = lHEventInfo.FindObject("particleRawJet")
        
        ptHardScaleFactor = lScaleFactor[ptHardBin-1]
        for epBin in range(0, 2):

            epLabel = lEPLabel[epBin]
            lHybRawJetEP = lHybRawJet.FindObject(epLabel)
            lParRawJetEP = lParRawJet.FindObject(epLabel)
            
            for centBin in range(embCentBinRange[centKind][0], embCentBinRange[centKind][1]):
                print('    '+lEPLabel[epBin]+' Scale Each Jet distributions Centbin'+str(centBin))
                hTempHybRawJetEP = lHybRawJetEP.FindObject('hJetCorrPtLocal_'+str(centBin))
                hHybRawJetEP = hTempHybRawJetEP.Clone()
                histName = 'hHybJetCorrPtLocal_'+str(centBin)
                hHybRawJetEP.SetName(histName)
                hHybRawJetEP.GetXaxis().SetRangeUser(0, 250)
                reBin = 1 # 5
                histSetting.histLabelSetting(hHybRawJetEP, histName, \
                    '#it{p}_{T, det}^{jet} [GeV/#it{c}]', 'count', centBin, 1, reBin)
                addHistIntoList(lMainTree[epBin][centBin], hHybRawJetEP, histName, \
                    ptHardScaleFactor, ptHardBin)
                
                hTempParRawJetEP = lParRawJetEP.FindObject('hJetCorrPtLocal_'+str(centBin))
                hParRawJetEP = hTempParRawJetEP.Clone()
                histName = 'hParJetCorrPtLocal_'+str(centBin)
                hParRawJetEP.SetName(histName)
                hParRawJetEP.GetXaxis().SetRangeUser(0, 250)
                reBin = 1 # 5
                histSetting.histLabelSetting(hParRawJetEP, histName, \
                    '#it{p}_{T, gen}^{jet} [GeV/#it{c}]', 'count', centBin, 1, reBin)
                addHistIntoList(lMainTree[epBin][centBin], hParRawJetEP, histName, \
                    ptHardScaleFactor, ptHardBin)
            
            listName = 'MatchedJetHisto_' + epLabel
            lMatchingJet = lHEventInfo.FindObject(listName)
            
            print('    '+lEPLabel[epBin]+' Sum hResponseMatrixDiff Plots')
            histName = 'hResponseMatrixDiff'
            hTemp = lMatchingJet.FindObject(histName)
            hTemp.GetAxis(5).SetRangeUser(centRangeList[centKind][0], centRangeList[centKind][1])
            hTemp.GetAxis(1).SetRangeUser(10., ptHardH[ptHardBin-1]*4) # ????
            histName = 'hRM_' + 'Cent{0}'.format(centKind)
            hEachCentRM = hTemp.Projection(0, 1, "")
            hEachCentRM.SetName(histName)
            addHistIntoList(lMainTree[epBin][centKind], hEachCentRM, histName, ptHardScaleFactor, ptHardBin)

            histName = 'hPartJetPt_'+str(centBin)+'pTHardBin'+str(ptHardBin)
            tempHParRawJetEPRM = hEachCentRM.Clone()
            tempHParRawJetEP = hTemp.Projection(1, "")
            tempHParRawJetEP.SetName(histName)
            tempHParRawJetEP.Scale(ptHardScaleFactor)
            lAllScaledGenJetPtDist[epBin].append(tempHParRawJetEP)
            histName = 'hHybJetPt_'+str(centBin)+'pTHardBin'+str(ptHardBin)
            tempHHybRawJetEPRM = hEachCentRM.Clone()
            tempHHybRawJetEP = hTemp.Projection(0, "")
            tempHHybRawJetEP.SetName(histName)
            tempHHybRawJetEP.Scale(ptHardScaleFactor)
            lAllScaledDetJetPtDist[epBin].append(tempHHybRawJetEP)

            print('    '+lEPLabel[epBin]+' Sum Each Jet Performace Plots')
            histName = 'hJESshift'
            hTemp = lMatchingJet.FindObject(histName)
            addHistIntoList(lMainTree[epBin], hTemp, histName, 1, ptHardBin)
            histName = 'hEmbDeltaPt'
            hTemp = lMatchingJet.FindObject(histName)
            addHistIntoList(lMainTree[epBin], hTemp, histName, 1, ptHardBin)
            histName = 'hJESshiftHybDet'
            hTemp = lMatchingJet.FindObject(histName)
            addHistIntoList(lMainTree[epBin], hTemp, histName, 1, ptHardBin)
            histName = 'hEmbDeltaPtHybDet'
            hTemp = lMatchingJet.FindObject(histName)
            addHistIntoList(lMainTree[epBin], hTemp, histName, 1, ptHardBin)
            histName = 'hJESshiftDetPar'
            hTemp = lMatchingJet.FindObject(histName)
            addHistIntoList(lMainTree[epBin], hTemp, histName, 1, ptHardBin)
            histName = 'hEmbDeltaPtDetPar'
            hTemp = lMatchingJet.FindObject(histName)
            addHistIntoList(lMainTree[epBin], hTemp, histName, 1, ptHardBin)

            histName = 'hNEFVsPt'
            hTemp = lMatchingJet.FindObject(histName)
            addHistIntoList(lMainTree[epBin], hTemp, histName, ptHardScaleFactor, ptHardBin)
            histName = 'hZLeadingVsPt'
            hTemp = lMatchingJet.FindObject(histName)
            addHistIntoList(lMainTree[epBin], hTemp, histName, ptHardScaleFactor, ptHardBin)
            histName = 'hMatchingDistance'
            hTemp = lMatchingJet.FindObject(histName)
            addHistIntoList(lMainTree[epBin], hTemp, histName, ptHardScaleFactor, ptHardBin)
            histName = 'fHistJetMatchingQA'
            # hTemp = lMatchingJet.FindObject(histName)
            # addHistIntoList(lMainTree[epBin], hTemp, histName, 1, ptHardBin)
        
    for epBin in range(0, 2):
        histName = 'hParJetCorrPtLocal_'+str(centKind)
        tempFullGenPt = lMainTree[epBin][centKind].FindObject(histName)
        lAllScaledGenJetPtDist[epBin].append(tempFullGenPt)
        histName = 'hHybJetCorrPtLocal_'+str(centKind)
        tempFullHybPt = lMainTree[epBin][centKind].FindObject(histName)
        lAllScaledDetJetPtDist[epBin].append(tempFullHybPt)

        tempGenDistList = list()
        tempDetDistList = list()
        
        for l in range(0,len(lAllScaledGenJetPtDist[epBin])): tempGenDistList.append(lAllScaledGenJetPtDist[epBin][l])
        for l in range(0,len(lAllScaledDetJetPtDist[epBin])): tempDetDistList.append(lAllScaledDetJetPtDist[epBin][l])
        #[kRed, kBlue, kGreen+3, kOrange+8, kAzure+10, kSpring+4, kViolet+1, kPink+9, kSpring-9, kOrange+8, kAzure+1, kTeal+4, kRed-7, kBlue-7, kGreen-2, kOrange-2, kViolet+7, kSpring+4, kTeal-5,kBrack]
        lTempFillColor = [1, 632, 600, 416+3, 800+8, 860+10, 900+10, 820+4, 880+1, 900+9, 820-9, 800+1, 860+1, 840+4, 632-7, 600-7, 416-2, 800-2, 880+7, 820+4, 840-5]
        plotStyleList = list()
        for l in range(0, len(tempGenDistList)): plotStyleList.append(0)
        colorList = [1, 600, 416+3, 900+9, 800+8, 880+1]
        hBlankForPtBinCompare = ROOT.TObject()
        legALICE = ROOT.TLegend(0.3,0.7,0.88,0.93,'')
        genePlotSets.addALICELegend(hBlankForPtBinCompare,legALICE,lCentLabel[centBin],'', 0.03)
        legJet = ROOT.TLegend(0.3,0.7,0.88,0.93,'')
        genePlotSets.addJetLegend(hBlankForPtBinCompare, legJet, 0.2, leadingTrackPtCut, 0.03)
        lNamePtBin = list()
        for ptHBin in range(1, len(tempGenDistList)): lNamePtBin.append('pT bin '+str(ptHBin))
        lNamePtBin.append('Total')
        
        tempGenDistList.reverse()
        tempDetDistList.reverse()
        lNamePtBin.reverse()
        
        legPtBin = ROOT.TLegend(0.3,0.7,0.88,0.93,'pT bins')
        genePlotSets.addSomeHistsLegend(legPtBin, tempGenDistList, lNamePtBin, 0.03)
        label = 'ParLvPtBinCompare'+ '_TrackPtCut'+str(leadingTrackPtCut)+'_Cent'+str(centBin)+lEPLabel[epBin]
        genePlotSets.overwrightSomePlots(hBlankForPtBinCompare, tempGenDistList, \
            0, legALICE,legJet,legPtBin, lTempFillColor, plotStyleList, label, '', outputDir)

        legPtBin = ROOT.TLegend(0.3,0.7,0.88,0.93,'pT bins')
        genePlotSets.addSomeHistsLegend(legPtBin, tempDetDistList, lNamePtBin, 0.03)
        label = 'HybLvPtBinCompare'+ '_TrackPtCut'+str(leadingTrackPtCut)+'_Cent'+str(centBin)+lEPLabel[epBin]
        genePlotSets.overwrightSomePlots(hBlankForPtBinCompare, tempDetDistList, \
            0, legALICE,legJet,legPtBin, lTempFillColor, plotStyleList, label, 'p ', outputDir)


#22  devHistForEachCent   222222222222222222222222222222222222222222222222222222
def PlotPerformanceHistEachCent(lMainTree, centKind, epBin):
    print('=== s === Add Performance Hists  #############################')
    
    # baseMatchedJetHists = lMainTree.FindObject('hResponseMatrixDiff')
    histName = 'hRM_' + 'Cent{0}'.format(centKind)
    baseRM = lMainTree[centKind].FindObject(histName)
    baseHJESshift = lMainTree.FindObject('hJESshift')
    baseHEmbDeltaPt = lMainTree.FindObject('hEmbDeltaPt')
    baseHJESshiftHybDet = lMainTree.FindObject('hJESshiftHybDet')
    baseHEmbDeltaPtHybDet = lMainTree.FindObject('hEmbDeltaPtHybDet')
    baseHJESshiftDetPar = lMainTree.FindObject('hJESshiftDetPar')
    baseHEmbDeltaPtDetPar = lMainTree.FindObject('hEmbDeltaPtDetPar')
    baseHNEFVsPt = lMainTree.FindObject('hNEFVsPt')
    baseHZLeadingVsPt = lMainTree.FindObject('hZLeadingVsPt')
    baseHMatchingDistance = lMainTree.FindObject('hMatchingDistance')

    for centBin in range(embCentBinRange[centKind][0], embCentBinRange[centKind][1]):
        ptRangeDict, ptBinArrayDict = PtRangeList.eachJetPtBinDef(0,centBin)
        ########################################################################
        # print('1. Set Up RM to Each Centrality for all histograms        #####')
        # tempOriginMatchedJetHists = baseMatchedJetHists.Clone()
        # histName = 'RMHistsCent'+str(centBin)
        # tempOriginMatchedJetHists.SetName(histName)
        # tempOriginMatchedJetHists.GetAxis(5).SetRangeUser(\
        #     centRangeList[centBin][0], centRangeList[centBin][1])
        #######################################################################

        ###  Add RM hist  ######################################################
        # print('1. Add Original RM                                        #####')
        # histName = 'hRM_' + 'Cent{0}'.format(centBin)
        # hTempEachCentRM = tempOriginMatchedJetHists.Clone()
        # hEachCentRM = hTempEachCentRM.Projection(0, 1, "")
        # hEachCentRM.SetName(histName)
        # lMainTree[centBin].Add(hEachCentRM)
        ########################################################################

        ###  Add RM hist  ######################################################
        print('2. Add rebin and normalize RM                             #####')
        histName = 'hRM_forUnfold_' + 'Cent{0}'.format(centBin)
        # hTempEachCentRM_forUnfold = tempOriginMatchedJetHists.Clone()
        # hEachCentRM_forUnfold = hTempEachCentRM_forUnfold.Projection(0, 1, "")
        hEachCentRM_forUnfold = baseRM.Clone()
        hEachCentRM_forUnfold.SetName(histName)
        lMainTree[centBin].Add(hEachCentRM_forUnfold)
        ########################################################################

        ###  Add Matched Gen Jet hist  #########################################
        print('3. Add Matched Particle Level Jet hist Centrality         #####')
        histName = 'hMatchGenJetPt_' + 'Cent{0}'.format(centBin)
        # hTempMatchedGenJet = tempOriginMatchedJetHists.Clone()
        # hMatchGenJetPt = hTempMatchedGenJet.Projection(0, histName)
        hTempMatchedGenJet = baseRM.Clone()
        hMatchGenJetPt = hTempMatchedGenJet.ProjectionY()
        hMatchGenJetPt.SetName(histName)
        reBin = 1 # 10
        histSetting.histLabelSetting(hMatchGenJetPt, histName, \
            '#it{p}_{T, gen}^{Matched jet} [GeV/#it{c}]', 'count', centBin, 1, reBin)
        hMatchGenJetPt.Draw('AP')
        lMainTree[centBin].Add(hMatchGenJetPt)
        ########################################################################

        ###  Add Matched Det Jet hist  #########################################
        print('4. Add Matched Detector Level Jet hist Centrality         #####')
        histName = 'hMatchDetJetPt_' + 'Cent{0}'.format(centBin)
        # hTempMatchedDetJet = tempOriginMatchedJetHists.Clone()
        # hMatchDetJetPt = hTempMatchedDetJet.Projection(1, histName)
        hTempMatchedDetJet = baseRM.Clone()
        hMatchDetJetPt = hTempMatchedDetJet.ProjectionX()
        hMatchDetJetPt.SetName(histName)
        reBin = 1 # 10
        histSetting.histLabelSetting(hMatchDetJetPt, histName, \
            '#it{p}_{T, det}^{Matched jet} [GeV/#it{c}]', 'count', centBin, 1, reBin)
        hMatchDetJetPt.Draw('AP')
        lMainTree[centBin].Add(hMatchDetJetPt)
        ##############################################################################

        ########################################################################
        # print('5. Add Matched Jets Delta R hist Centrality               #####')
        # histName = 'hMatchJetDeltaR_' + 'Cent{0}'.format(centBin)
        # hTempMatchJetDeltaR = tempOriginMatchedJetHists.Clone()
        # hMatchJetDeltaR = hTempMatchJetDeltaR.Projection(2, histName)
        # hMatchJetDeltaR.SetName(histName)
        # histSetting.histLabelSetting(hMatchJetDeltaR, histName, '#Delta R', 'count', centBin, 1, 1)
        # hMatchJetDeltaR.Draw('AP')
        # lMainTree[centBin].Add(hMatchJetDeltaR)
        ########################################################################

        ########################################################################
        # print('6. Add Matched Jets Delta Angle hist Centrality           #####')
        # histName = 'hMatchJetDAngle_' + 'Cent{0}'.format(centBin)
        # hTemphMatchJetDAngle = tempOriginMatchedJetHists.Clone()
        # hMatchJetDAngle = hTemphMatchJetDAngle.Projection(4, histName)
        # hMatchJetDAngle.SetName(histName)
        # histSetting.histLabelSetting(hMatchJetDAngle, histName, 'angle', 'count', centBin, 1, 1)
        # hMatchJetDAngle.Draw('AP')
        # lMainTree[centBin].Add(hMatchJetDAngle)
        ########################################################################

        ###  Add hJESshiftEMCal hist   #########################################
        print('7. Add Jet Energy Scale Shift hists                       #####')
        lJESshiftPtRange = [[20,30], [50,70], [100,120]]
        lJESshiftDistHists = ROOT.TList()
        lJESshiftDistHists.SetName('hJESshiftList')
        
        histName = 'hJESshift_' + 'Cent{0}'.format(centBin)
        hJESshift = baseHJESshift.Clone(histName)
        hJESshift.GetXaxis().SetRangeUser(centRangeList[centBin][0], centRangeList[centBin][1])
        
        histName = 'hJESshiftForProj_' + 'Cent{0}'.format(centBin) 
        hJESshift_clone = hJESshift.Clone()
        hJESshift_clone.SetName(histName)
        hJESshift_proj = hJESshift_clone.Project3D("zyeo")
        hJESshift_prof = hJESshift_proj.ProfileX()
        
        xTitle = '#it{p}_{T, particle}^{jet}'
        yTitle = "#frac{#it{p}_{T}^{det} - #it{p}_{T}^{particle}}{#it{p}_{T}^{particle}}"
        histSetting.histLabelSetting(hJESshift_prof, histName, xTitle, yTitle, centBin, 1, reBin)
        lMainTree[centBin].Add(hJESshift_prof)

        print('8.1 Add Jet Energy Scale Shift projected hists for three jet pt range   #####')
        lPtRangeColor = [1, 632, 600]
        for ptRangeKind in range(0, 2):
            histName = 'hJESshiftDist_' + 'Cent{0}'.format(centBin)\
                + '_ptRange{}'.format(ptRangeKind)
            hJESshif_clone = hJESshift.Clone()
            hJESshif_clone.SetName(histName)
            hJESshif_clone.GetYaxis().SetRangeUser(\
                lJESshiftPtRange[ptRangeKind][0], lJESshiftPtRange[ptRangeKind][1])
            hJESshif_clone_proj = hJESshif_clone.Project3D('ze')

            xTitle = "#frac{#it{p}_{T}^{det} - #it{p}_{T}^{gen}}{#it{p}_{T}^{gen}}"
            yTitle = 'probability'
            histSetting.histLabelSetting(hJESshif_clone_proj, histName, xTitle, yTitle, \
                ptRangeKind, 1, reBin)
            hJESshif_clone_proj.SetLineColor(lPtRangeColor[ptRangeKind])
            hJESshif_clone_proj.SetMarkerColor(lPtRangeColor[ptRangeKind])
            lMainTree[centBin].Add(hJESshif_clone_proj)
        ########################################################################

        ###  Add JER1 hist  ####################################################
        print('9. Add Jet Energy Resolution hists 1 way                  #####')
        histName = 'hJESForJER3D_' + 'Cent{0}'.format(centBin)
        hJESshift_clone = hJESshift.Clone(histName)
        histName = 'hJESForJER2D_' + 'Cent{0}'.format(centBin)
        hJESshift_proj2D = hJESshift_clone.Project3D("zyeo")
        hJESshift_proj2DCP = hJESshift_proj2D.Clone(histName)
        
        histName = 'JER1_cent{}'.format(centBin)+'_'+lEPLabel[epBin]
        hEachCentJER = plotPerformanceHists.getJER1(hJESshift_proj2DCP, ptBinArrayDict['mcGen'], histName)

        xTitle = '#it{p}_{T, gen}^{jet} [GeV/#it{c}]'
        yTitle = "#sigma(#frac{#it{p}_{T}^{det} - #it{p}_{T}^{particle}}{#it{p}_{T}^{particle}})"
        histSetting.histLabelSetting(hEachCentJER, histName, xTitle, yTitle, centBin, 1, reBin)
        lMainTree[centBin].Add(hEachCentJER)
        ########################################################################

        ###  Add JER2 hist  ####################################################
        # print('10. Add Jet Energy Resolution hists 2 way                 #####')
        # histName = 'hRM_forJER2_' + 'Cent{0}'.format(centBin)+'_'+lEPLabel[epBin]
        # # hTempEachCentRM_forJER = tempOriginMatchedJetHists.Clone()
        # # hEachCentRM_forJER = hTempEachCentRM_forJER.Projection(0, 1, "")
        # hEachCentRM_forJER = baseRM.Clone()
        # hEachCentRM_forJER.SetName(histName)
        # histName = 'JER2_cent{}'.format(centBin)+'_'+lEPLabel[epBin]
        # hEachCentJER = plotPerformanceHists.getJER2(hEachCentRM_forJER,histName)

        # xTitle = '#it{p}_{T, gen}^{jet} [GeV/#it{c}]'
        # yTitle = "#frac{#sigma(#it{p}_{T}^{gen})}{#it{p}_{T}^{gen}}"
        # histSetting.histLabelSetting(hEachCentJER, histName, xTitle, yTitle, centBin, 1, reBin)
        # lMainTree[centBin].Add(hEachCentJER)
        ########################################################################

        ###  Add hNEFVsPt hist   ###############################################
        # print('10. Add NEFVsPt                                           #####')
        # histName = 'hNEFVsPt_' + 'Cent{0}'.format(centBin)
        # hNEFVsPt = baseHNEFVsPt.Clone()
        # hNEFVsPt.SetName(histName)
        # hNEFVsPt.GetXaxis().SetRangeUser(centRangeList[centBin][0], centRangeList[centBin][1])
        # hNEFVsPt_proj = hNEFVsPt.Project3D("zyeo")
        # hNEFVsPt_profX = hNEFVsPt_proj.ProfileX()

        # xTitle = '#it{p}_{T, gen}^{jet}'
        # yTitle = "NEF"
        # histSetting.histLabelSetting(hNEFVsPt_profX, histName, xTitle, yTitle, centBin, 1, reBin)
        # lMainTree[centBin].Add(hNEFVsPt_profX)
        ##############################################################################

        ###  Add hZLeadingVsPt hist   ################################################
        # histName = 'hZLeadingVsPt_' + 'Cent{0}'.format(centBin)
        # hZLeadingVsPt = baseHZLeadingVsPt.Clone()
        # hZLeadingVsPt.SetName(histName)
        # hZLeadingVsPt.GetXaxis().SetRangeUser(centRangeList[centBin][0], centRangeList[centBin][1])
        # hZLeadingVsPt_proj = hZLeadingVsPt.Project3D("zyeo")
        # hZLeadingVsPt_profX = hZLeadingVsPt_proj.ProfileX()

        # xTitle = '#it{p}_{T, gen}^{jet} [GeV/#it{c}]'
        # yTitle = "zLeading #it{p}_{T}"
        # histSetting.histLabelSetting(hZLeadingVsPt_profX, histName, xTitle, yTitle, centBin, 1, reBin)
        # lMainTree[centBin].Add(hZLeadingVsPt_profX)
        ##############################################################################


def EditRMForUF(outList, hOriginRM, centBin, epBin, label):
    ptRangeDict, ptBinArrayDict = PtRangeList.eachJetPtBinDef(0,centBin)
    hRM_rebin = rebinRM(outList, hOriginRM, epBin, ptBinArrayDict)
    hRM_norm = normalizeRM(outList,hRM_rebin, epBin, outList, ptBinArrayDict)
    
    return hRM_norm

################################################################################
# Rebin the response matrix to have variable binning                          ##
################################################################################
def rebinRM(outList, hOriginRM, epBin, ptBinArrayDict):
    genBinArray = ptBinArrayDict['mcGen']
    nGenBins = len(genBinArray )-1
    detBinArray = ptBinArrayDict['mcDet']
    nDetBins = len(detBinArray )-1
    histname = "{}_Rebinned".format(hOriginRM.GetName())+'_'+lEPLabel[epBin]
    title = histname + ";#it{p}_{T,corr}^{det} (GeV/#it{c});#it{p}_{T}^{truth} (GeV/#it{c})"
    hRM_rebin = ROOT.TH2D(histname, title, nDetBins, detBinArray, nGenBins, genBinArray)

    # Loop over all bins in fine-binned response matrix, 
    # and fill appropriate bin in new response matrix
    # Assume that the bin edges overlap appropriately
    for detBin in range(1, hOriginRM.GetNbinsX() + 1):
        for genBin in range(1, hOriginRM.GetNbinsY() + 1):
            oldContent = hOriginRM.GetBinContent(detBin, genBin)

            # Find the bin that should be filled in the new histogram, and fill it
            # Need to get (x,y) location from bins (ibin, jbin)
            x = hOriginRM.GetXaxis().GetBinCenter(detBin)
            y = hOriginRM.GetYaxis().GetBinCenter(genBin)
            hRM_rebin.Fill(x, y, oldContent)
    
    # Assume 0 errors on response matrix
    for bin in range(1, hRM_rebin.GetNcells() + 1):
        hRM_rebin.SetBinError(bin, 0)

    outList.Add(hRM_rebin)

    return hRM_rebin

################################################################################
# Normalize response matrix                                                   ##
# Normalize the pT-truth projection to 1                                      ##
################################################################################
def normalizeRM(outList,hRM_rebin, epBin, outputDir, ptBinArrayDict):
    genBinArray = ptBinArrayDict['mcGen']
    nGenBins = len(genBinArray)-1
    detBinArray = ptBinArrayDict['mcDet']
    nDetBins = len(detBinArray)-1
    # Make projection onto pT-true axis (y-axis), and scale appropriately
    # Do exclude under and overflow bins
    hGenProjBefore = hRM_rebin.ProjectionY("_py",1,hRM_rebin.GetNbinsX()) 
    hGenProjBefore.SetName("hGenProjectionBefore")

    histname = "{}Normed".format(hRM_rebin.GetName())+'_'+lEPLabel[epBin]
    title = histname + ";#it{p}_{T,corr}^{det} (GeV/#it{c});#it{p}_{T}^{truth} (GeV/#it{c})"
    hRM_norm = ROOT.TH2D(histname, title, nDetBins, detBinArray, nGenBins, genBinArray)

    # Loop through truth-level bins, and apply normalization factor to all bins.
    nBinsY = hRM_rebin.GetNbinsY() # pT-gen
    nBinsX = hRM_rebin.GetNbinsX() # pT-det
    for genBin in range(1,nBinsY+1):
        normFactor = hGenProjBefore.GetBinContent(genBin)
        if normFactor > 0:
            genBinCenter = hGenProjBefore.GetXaxis().GetBinCenter(genBin)

            for detBin in range(1,nBinsX+1):
                binContent = hRM_rebin.GetBinContent(detBin, genBin)
                # hRM_rebin.SetBinContent(detBin, genBin, binContent/normFactor)
                hRM_norm.SetBinContent(detBin, genBin, binContent/normFactor)
        
    # Plot response matrix
    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    c.cd().SetLeftMargin(0.15)

    # hRM_rebin.Draw("colz")
    hRM_norm.Draw("colz")
    line = ROOT.TLine(minPtDet,0,minPtDet,250)
    line.SetLineColor(0)
    line.SetLineStyle(2)
    line.Draw("same")
    line2 = ROOT.TLine(maxPtDet,0,maxPtDet,250)
    line2.SetLineColor(0)
    line2.SetLineStyle(2)
    line2.Draw("same")
    line3 = ROOT.TLine(0,minPtGen,100,minPtGen)
    line3.SetLineColor(0)
    line3.SetLineStyle(2)
    line3.Draw("same")
    line4 = ROOT.TLine(0,maxPtGen,100,maxPtGen)
    line4.SetLineColor(0)
    line4.SetLineStyle(2)
    line4.Draw("same")
    
    # c.SaveAs('hoge.root')
    c.Close()

    outList.Add(hRM_norm)

    return hRM_norm


def addHistIntoList(lTree, inputHist, histName, scaleFactor, ptHardBin):
    inputHist.Sumw2()
    inputHist.Scale(scaleFactor)

    if ptHardBin == 1: 
        lTree.Add(inputHist)
    else :
        tempHist = lTree.FindObject(histName)
        # tempHist = lTree.Get(histName)
        tempHist.Add(inputHist)
    

################################################################################
def OFileStracture(lMainTree, numOfCentBin):
    for epBin in range(0, 3):
        listName = lEPLabel[epBin]
        lEPList = ROOT.TList()
        lEPList.SetName(listName)
        for centBin in range(0, numOfCentBin):
            lEachCentHists = ROOT.TList()
            lEachCentHists.SetName('lCent{0}'.format(centBin))
            lEPList.Add(lEachCentHists)
            
        lMainTree.Add(lEPList)
################################################################################

if __name__ == "__main__":
    args = sys.argv
    # centKind, LHCPeriod, leadingTrackPtCut, diffSys = sys.argv
    # MergeScalePtHardPlots()
    # MergeScalePtHardPlots(centKind, LHCPeriod, leadingTrackPtCut, diffSys)
    print(args[1], args[2], args[3],args[4])
    MergeScalePtHardPlots(int(args[1]), args[2], int(args[3]),args[4])