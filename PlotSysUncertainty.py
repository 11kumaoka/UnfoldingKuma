import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TMath, TGraphErrors, TGraphMultiErrors
import argparse
import ctypes
import os
import gc
from array import array
import numpy as np

import PtRangeList
import genePlotSets
import histSetting
import plotPerformanceHists

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

numOfCentBin = 4
lCentLabel = ['0-5','5-10','10-30','30-50','50-80']

leadingTrackPtCut = 5 # 0, 5, 7 GeV/c

lValKindName = ['lChJetV2', 'lIncChPbPbJetYeild', 'lChJetRAA', 'lOEPChPbPbJetYeild', 'lIEPChPbPbJetYeild']


# lSysType = ['Nominal','iteKM1', 'iteKP1','DetLvLcut', 'DetLvUcut', 'DiffPrior', 'bkgWay', 'V0Diff']
lSysType = ['Nominal','iteKM1', 'iteKP1','DetLvLcut', 'DetLvUcut', 'DiffPrior',\
    'BKGNoFit', 'BKGV2', 'V0C', 'V0A']

lEPLabel = ['OutOfPlane', 'InPlane', 'Inclusive']
lUnfoldType = ['Bayes', 'SVD']
unfoldType = 'Bayes'

uncAbs = 1 # Ratio: 0, Abs: 1

# lPsi2Reso = [0.4513, 0.6388, 0.6801, 0.6463, 0.5609, 0.4455, 0.2860, 0.1400, 0.0477]
# lPsi2Reso = [0.617, 0.798, 0.815, 0.803, 0.736, 0.623, 0.431, 0.238, 0.082]
# 0-5, 5-10, 10-30, 30-50, 50-80
lPsi2Reso = [0.4513, 0.6388, 0.6632, 0.5032,  0.1579]
# lPsi2Reso = [0.617, 0.798, 0.809, 0.68, 0.250]


# inUnfoldedFileDir = './'
inUnfoldedFileDir = '/Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18qr/pass3/Ch/Unfolded/'

# outFinalFileDir = '~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/FinalPlots/'
outFinalFileDir = '/Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18qr/pass3/Ch/FinalPlots/'
# outFinalFileDir = './FinalPlots/'


lFillMerker = [20, 21, 22, 23, 34, 29, 33, 31, 24, 24, 26, 32, 20, 21, 22, 23, 34, 29]
lNoFillMerker = [24, 25, 26, 32, 28, 30, 27, 24, 25, 26]

#[kBrack, kRed, kBlue, kGreen+3, kOrange+8, kAzure+10, kSpring+4, kViolet+1, kPink+9, kSpring-9, kOrange+8, kAzure+1, kTeal+4, kRed-7, kBlue-7, kGreen-2, kOrange-2, kViolet+7, kSpring+4, kTeal-5]
lFillColor = [1, 632, 600, 416+3, 800+8, 860+10, 900+10, 820+4, 880+1, 900+9, 820-9, 800+1, 860+1, 840+4, 632-7, 600-7, 416-2, 800-2, 880+7, 820+4, 840-5]

lBkgSysUncErr = [[1.7754592434432792, 2.8046837280277077, 3.412777197542915, 2.1635127562552214, \
        13.210083385018262, 0.9124633802012357, 3.2648648648648657], \
    [2.2590301304044056, 3.345589270985158, 2.7947686217378713, 2.6134988719298167, \
        0.6349590956078255, 1.7493155319004983, 1.0710999281092741], \
    [2.1786806311149878, 2.200098969310484, 1.5525326719450536, 3.53059313347085, \
        1.5380989787902595, 0.8879668049792533, 0.7414965986394557], \
    [1.0180839303016933, 0.8720927017184233, 0.5818151665579498, 0.6840665579027901, \
        0.8788701792504074, 0.7788553869790669, 1.0339489885664028]]

lV0DiffSysUncErr = [[0.5850256766345653, 0.4441944711573965, 0.4886183863596415, \
        0.567446184760034, 19.048366586274938, 0.20437710437710443, 0.5798491151378411], \
    [0.4353385105741184, 0.4443557629511435, 0.48484707736804766, 0.3111777748793532, \
        0.10662784152113035, 0.3949969416862457, 0.4627170582226761], \
    [0.3272861619207525, 0.29838228128600036, 0.24542852533268028, 0.7909063489669733, \
        0.023008801105955518, 0.1916564185034527, 0.3077955185831211], \
    [0.15321061229985566, 0.2386214750195241, 0.12569887604335997, 0.11536357604104279, \
        0.2846680304063439, 0.7015631880228973, 0.1986779287550497]]


###################################################################################
# Main function
def PlotSysUncertainty():
    inUnfoldedFileName = inUnfoldedFileDir \
        + 'UnfoldedPtDists'+ '_TrackPtCut'+str(leadingTrackPtCut)+'_.root'
    inputFile = ROOT.TFile(inUnfoldedFileName, "READ")
    iMainTree = inputFile.Get('mainTree')
    
    oEpTree = iMainTree.FindObject(lEPLabel[0])
    iEpTree = iMainTree.FindObject(lEPLabel[1])
    incTree = iMainTree.FindObject(lEPLabel[2])

    inUnfoldedFileNameEff94 = inUnfoldedFileDir \
        + 'UnfoldedPtDists' + '_TrackPtCut'+str(leadingTrackPtCut)+'_TrackEff094.root'
    inputFileEff94 = ROOT.TFile(inUnfoldedFileNameEff94 , "READ")
    iMainTreeEff94 = inputFileEff94.Get('mainTree')
    oEpTreeEff94 = iMainTreeEff94.FindObject(lEPLabel[0])
    iEpTreeEff94 = iMainTreeEff94.FindObject(lEPLabel[1])
    incTreeEff94 = iMainTreeEff94.FindObject(lEPLabel[2])

    # === BKG No Fit ========= 
    inUnfoldedFileNameBKGNoFit = inUnfoldedFileDir \
        + 'UnfoldedPtDists' + '_TrackPtCut'+str(leadingTrackPtCut)+'_BKGNoFit.root'
    inputFileBKGNoFit = ROOT.TFile(inUnfoldedFileNameBKGNoFit, "READ")
    iMainTreeBKGNoFit = inputFileBKGNoFit.Get('mainTree')
    oEpTreeBKGNoFit = iMainTreeBKGNoFit.FindObject(lEPLabel[0])
    iEpTreeBKGNoFit = iMainTreeBKGNoFit.FindObject(lEPLabel[1])
    incTreeBKGNoFit = iMainTreeBKGNoFit.FindObject(lEPLabel[2])

    # === BKG V2 Fit ========= 
    inUnfoldedFileNameBKGV2 = inUnfoldedFileDir \
        + 'UnfoldedPtDists' + '_TrackPtCut'+str(leadingTrackPtCut)+'_BKGV2.root'
    inputFileBKGV2 = ROOT.TFile(inUnfoldedFileNameBKGV2, "READ")
    iMainTreeBKGV2 = inputFileBKGV2.Get('mainTree')
    oEpTreeBKGV2 = iMainTreeBKGV2.FindObject(lEPLabel[0])
    iEpTreeBKGV2 = iMainTreeBKGV2.FindObject(lEPLabel[1])
    incTreeBKGV2 = iMainTreeBKGV2.FindObject(lEPLabel[2])
    
    # === V0C ========= 
    inUnfoldedFileNameV0C = inUnfoldedFileDir \
        + 'UnfoldedPtDists' + '_TrackPtCut'+str(leadingTrackPtCut)+'_V0C.root'
    inputFileV0C = ROOT.TFile(inUnfoldedFileNameV0C, "READ")
    iMainTreeV0C = inputFileV0C.Get('mainTree')
    oEpTreeV0C = iMainTreeV0C.FindObject(lEPLabel[0])
    iEpTreeV0C = iMainTreeV0C.FindObject(lEPLabel[1])
    incTreeV0C = iMainTreeV0C.FindObject(lEPLabel[2])
    
    # === V0A ========= 
    inUnfoldedFileNameV0A = inUnfoldedFileDir \
        + 'UnfoldedPtDists' + '_TrackPtCut'+str(leadingTrackPtCut)+'_V0A.root'
    inputFileV0A = ROOT.TFile(inUnfoldedFileNameV0A, "READ")
    iMainTreeV0A = inputFileV0A.Get('mainTree')
    oEpTreeV0A = iMainTreeV0A.FindObject(lEPLabel[0])
    iEpTreeV0A = iMainTreeV0A.FindObject(lEPLabel[1])
    incTreeV0A = iMainTreeV0A.FindObject(lEPLabel[2])


    oFileName = outFinalFileDir + 'FinalResultsWithSys.root'
    oFile = ROOT.TFile(oFileName, "RECREATE")
    lOMainTree = ROOT.TList()
    OFileStracture(lOMainTree, numOfCentBin)
    oFile.cd()
    
    # 0: ChJetV2, PbPbYield, ChJetRAA,  PbPbYieldOEP, PbPbYieldIEP
    lMainAllCentObs = list()
    for valKindBin in range(0, 5):
        lAllCentObs = list()
        for centBin in range(0, numOfCentBin):
            lEachCentObs = list()
            lAllCentObs.append(lEachCentObs)
        
        lMainAllCentObs.append(lAllCentObs)


    for centBin in range(0, numOfCentBin):
        
        lDTargedObs, lRawJetObs = ExtractTargeObservable(\
            oEpTree, iEpTree, incTree, \
            oEpTreeEff94, iEpTreeEff94, incTreeEff94,\
            oEpTreeBKGNoFit, iEpTreeBKGNoFit, incTreeBKGNoFit, \
            oEpTreeBKGV2, iEpTreeBKGV2, incTreeBKGV2,\
            oEpTreeV0C, iEpTreeV0C, incTreeV0C,\
            oEpTreeV0A, iEpTreeV0A, incTreeV0A,\
            centBin)
        
        
        ptRangeDict, ptBinArrayDict = PtRangeList.eachJetPtBinDef(0, centBin)
        pTBinArray = ptBinArrayDict['reported']
        dhRelaSysUncOEP, dlRelaSysUncErrOEP = GetSysUncErrObj(lMainAllCentObs[3],lDTargedObs[3],\
            3, pTBinArray, centBin,lOMainTree,oFile)
        dhRelaSysUncIEP, dlRelaSysUncErrIPE = GetSysUncErrObj(lMainAllCentObs[4],lDTargedObs[4],\
            4, pTBinArray, centBin,lOMainTree,oFile)
        # print(dlRelaSysUncErrOEP)
        GetSysUncErrObjV2(lMainAllCentObs[0],lDTargedObs[0],\
            lDTargedObs[3]['Nominal'], lDTargedObs[4]['Nominal'],\
            dlRelaSysUncErrOEP, dlRelaSysUncErrIPE, \
                0, pTBinArray, centBin,lOMainTree,oFile)

        ptRangeDict, ptBinArrayDict = PtRangeList.eachJetPtBinDef(1, centBin)
        pTBinArray = ptBinArrayDict['reported']
        GetSysUncErrObj(lMainAllCentObs[1],lDTargedObs[1], 1, pTBinArray, centBin,lOMainTree,oFile)

        lOMainTree[0][centBin].Add(lRawJetObs[0])
        lOMainTree[1][centBin].Add(lRawJetObs[1])

    CalcChRAA(lMainAllCentObs, lOMainTree[2])
    
    PlotAllCentTargetObs(lMainAllCentObs[0], 0) # ChJetV2
    PlotAllCentTargetObs(lMainAllCentObs[2], 2) # ChJetRAA

    CompareWithAnotherCentral(lMainAllCentObs[0][0])
    CompareWithAnotherSemiCentral(lMainAllCentObs[0][3])

    oFile.cd()
    lOMainTree.Write('mainTree', 1)
    oFile.Close()

    print('root ' + outFinalFileDir)

def ExtractTargeObservable(oEpTree, iEpTree, incTree, \
    oEpTreeEff94, iEpTreeEff94, incTreeEff94,\
    oEpTreeBKGNoFit, iEpTreeBKGNoFit, incTreeBKGNoFit, \
    oEpTreeBKGV2, iEpTreeBKGV2, incTreeBKGV2,\
    oEpTreeV0C, iEpTreeV0C, incTreeV0C,\
    oEpTreeV0A, iEpTreeV0A, incTreeV0A,\
    centBin):
    dChJetV2 = {}
    dJetYield = {}
    dJetYieldOEP = {}
    dJetYieldIEP = {}

    treeName = 'lCent' + str(centBin)
    oEpCentTree = oEpTree.FindObject(treeName)
    iEpCentTree = iEpTree.FindObject(treeName)
    incCentTree = incTree.FindObject(treeName)
    treeName = 'lUnfoldedJets'
    oEpCentUFJetTree = oEpCentTree.FindObject(treeName)
    iEpCentUFJetTree = iEpCentTree.FindObject(treeName)
    incCentUFJetTree = incCentTree.FindObject(treeName)

    treeName = 'lJetDistQA'
    oEpCentRawJetTree = oEpCentTree.FindObject(treeName)
    iEpCentRawJetTree = iEpCentTree.FindObject(treeName)
    incCentRawJetTree = incCentTree.FindObject(treeName)
    RawJetName = 'hMeasuredJetYieldPerGeV'
    oEpRawJetPtYield = oEpCentRawJetTree.FindObject(RawJetName)
    iEpRawJetPtYield = iEpCentRawJetTree.FindObject(RawJetName)
    incRawJetPtYield = incCentRawJetTree.FindObject(RawJetName)
    
    ptRangeDict, ptBinArrayDict = PtRangeList.eachJetPtBinDef(1, centBin)
    pTBinArray = ptBinArrayDict['reported']
    nBins = len(pTBinArray) - 1
    incRawJetPtYield = incRawJetPtYield.Rebin(nBins, RawJetName, pTBinArray)

    hRawChJetV2 = CalcChJetV2(oEpRawJetPtYield, iEpRawJetPtYield,centBin,'Raw')
    ptRangeDict, ptBinArrayDict = PtRangeList.eachJetPtBinDef(0, centBin)
    pTBinArray = ptBinArrayDict['reported']
    nBins = len(pTBinArray) - 1
    RawJetV2Name = 'hRawChargedJetV2Cent'+str(centBin)
    hRawChJetV2 = hRawChJetV2.Rebin(nBins, RawJetV2Name, pTBinArray)

    baseUFJetName = 'hUnfoldedJetPt'
    for sysTypeBin in range(0, len(lSysType)-4): 
        UFJetName = baseUFJetName +'_'+ unfoldType +'_'+ lSysType[sysTypeBin]
        oEpUFJetPtYield = oEpCentUFJetTree.FindObject(UFJetName)
        iEpUFJetPtYield = iEpCentUFJetTree.FindObject(UFJetName)
        incUFJetPtYield = incCentUFJetTree.FindObject(UFJetName)
        if not oEpUFJetPtYield: print("ERROR no "+UFJetName+" found")
        histName = 'hOEpncPbPbJetYield'+'_'+ unfoldType +'_'+ lSysType[sysTypeBin]
        hOEpUFJetPtYield= oEpUFJetPtYield.Clone(histName)
        dJetYieldOEP[lSysType[sysTypeBin]] = hOEpUFJetPtYield
        histName = 'hIEpIncPbPbJetYield'+'_'+ unfoldType +'_'+ lSysType[sysTypeBin]
        hIEpIncPbPbJetYield = iEpUFJetPtYield.Clone(histName)
        dJetYieldIEP[lSysType[sysTypeBin]] = hIEpIncPbPbJetYield
        
        hChJetV2 = CalcChJetV2(oEpUFJetPtYield, iEpUFJetPtYield,centBin, lSysType[sysTypeBin])
        dChJetV2[lSysType[sysTypeBin]] = hChJetV2
        
        histName = 'hIncPbPbJetYield'+'_'+ unfoldType +'_'+ lSysType[sysTypeBin]
        hIncPbPbJetYield = incUFJetPtYield.Clone(histName)
        dJetYield[lSysType[sysTypeBin]] = hIncPbPbJetYield

    
    # == s == track Efficiency 94%  ===================================
    treeName = 'lCent' + str(centBin)
    oEpCentTreeEff94 = oEpTreeEff94.FindObject(treeName)
    iEpCentTreeEff94 = iEpTreeEff94.FindObject(treeName)
    incCentTreeEff94 = incTreeEff94.FindObject(treeName)
    treeName = 'lUnfoldedJets'
    oEpCentUFJetTreeEff94 = oEpCentTreeEff94.FindObject(treeName)
    iEpCentUFJetTreeEff94 = iEpCentTreeEff94.FindObject(treeName)
    incCentUFJetTreeEff94 = incCentTreeEff94.FindObject(treeName)

    UFJetName = baseUFJetName +'_'+ unfoldType +'_'+ lSysType[0]
    oEpUFJetPtYieldEff94 = oEpCentUFJetTreeEff94.FindObject(UFJetName)
    iEpUFJetPtYieldEff94 = iEpCentUFJetTreeEff94.FindObject(UFJetName)
    histName = 'hOEpncPbPbJetYield'+'_'+ unfoldType +'_'+ 'trackEff94'
    hOEpUFJetPtYieldEff94 = oEpUFJetPtYieldEff94.Clone(histName)
    dJetYieldOEP['trackEff94'] = hOEpUFJetPtYieldEff94
    histName = 'hIEpIncPbPbJetYield'+'_'+ unfoldType +'_'+ 'trackEff94'
    hIEpIncPbPbJetYieldEff94 = iEpUFJetPtYieldEff94.Clone(histName)
    dJetYieldIEP['trackEff94'] = hIEpIncPbPbJetYieldEff94
    
    hChJetV2Eff94 = CalcChJetV2(oEpUFJetPtYieldEff94,iEpUFJetPtYieldEff94,centBin,'trackEff94')
    dChJetV2['trackEff94'] = hChJetV2Eff94
    
    incUFJetPtYieldEff94 = incCentUFJetTree.FindObject(UFJetName)
    histName = 'hIncPbPbJetYield'+'_'+ unfoldType +'_'+ lSysType[0]
    hIncPbPbJetYieldEff94 = incUFJetPtYieldEff94.Clone(histName)
    dJetYield['trackEff94'] = hIncPbPbJetYieldEff94
    # == e == track Efficiency 94%  ===================================
    
    # == s == BKGNoFit  ===================================
    treeName = 'lCent' + str(centBin)
    oEpCentTreeBKGNoFit = oEpTreeBKGNoFit.FindObject(treeName)
    iEpCentTreeBKGNoFit = iEpTreeBKGNoFit.FindObject(treeName)
    incCentTreeBKGNoFit = incTreeBKGNoFit.FindObject(treeName)
    treeName = 'lUnfoldedJets'
    oEpCentUFJetTreeBKGNoFit = oEpCentTreeBKGNoFit.FindObject(treeName)
    iEpCentUFJetTreeBKGNoFit = iEpCentTreeBKGNoFit.FindObject(treeName)
    incCentUFJetTreeBKGNoFit = incCentTreeBKGNoFit.FindObject(treeName)

    UFJetName = baseUFJetName +'_'+ unfoldType +'_'+ lSysType[0]
    oEpUFJetPtYieldBKGNoFit = oEpCentUFJetTreeBKGNoFit.FindObject(UFJetName)
    iEpUFJetPtYieldBKGNoFit = iEpCentUFJetTreeBKGNoFit.FindObject(UFJetName)
    
    histName = 'hOEpncPbPbJetYield'+'_'+ unfoldType +'_'+ 'BKGNoFit'
    hOEpUFJetPtYieldBKGNoFit = oEpUFJetPtYieldBKGNoFit.Clone(histName)
    dJetYieldOEP['BKGNoFit'] = hOEpUFJetPtYieldBKGNoFit
    histName = 'hIEpIncPbPbJetYield'+'_'+ unfoldType +'_'+ 'BKGNoFit'
    hIEpIncPbPbJetYieldBKGNoFit = iEpUFJetPtYieldBKGNoFit.Clone(histName)
    dJetYieldIEP['BKGNoFit'] = hIEpIncPbPbJetYieldBKGNoFit
    
    hChJetV2BKGNoFit = CalcChJetV2(oEpUFJetPtYieldBKGNoFit,iEpUFJetPtYieldBKGNoFit,centBin,'BKGNoFit')
    dChJetV2['BKGNoFit'] = hChJetV2BKGNoFit
    
    incUFJetPtYieldBKGNoFit = incCentUFJetTree.FindObject(UFJetName)
    histName = 'hIncPbPbJetYield'+'_'+ unfoldType +'_'+ lSysType[0]
    hIncPbPbJetYieldBKGNoFit = incUFJetPtYieldBKGNoFit.Clone(histName)
    dJetYield['BKGNoFit'] = hIncPbPbJetYieldBKGNoFit
    # == e == BKGNoFit  ===================================
    
    # == s == BKGV2Fit  ===================================
    treeName = 'lCent' + str(centBin)
    oEpCentTreeBKGV2 = oEpTreeBKGV2.FindObject(treeName)
    iEpCentTreeBKGV2 = iEpTreeBKGV2.FindObject(treeName)
    incCentTreeBKGV2 = incTreeBKGV2.FindObject(treeName)
    treeName = 'lUnfoldedJets'
    oEpCentUFJetTreeBKGV2 = oEpCentTreeBKGV2.FindObject(treeName)
    iEpCentUFJetTreeBKGV2 = iEpCentTreeBKGV2.FindObject(treeName)
    incCentUFJetTreeBKGV2 = incCentTreeBKGV2.FindObject(treeName)

    UFJetName = baseUFJetName +'_'+ unfoldType +'_'+ lSysType[0]
    oEpUFJetPtYieldBKGV2 = oEpCentUFJetTreeBKGV2.FindObject(UFJetName)
    iEpUFJetPtYieldBKGV2 = iEpCentUFJetTreeBKGV2.FindObject(UFJetName)
    
    histName = 'hOEpncPbPbJetYield'+'_'+ unfoldType +'_'+ 'BKGV2'
    hOEpUFJetPtYieldBKGV2 = oEpUFJetPtYieldBKGV2.Clone(histName)
    dJetYieldOEP['BKGV2'] = hOEpUFJetPtYieldBKGV2
    histName = 'hIEpIncPbPbJetYield'+'_'+ unfoldType +'_'+ 'trackBKGV2'
    hIEpIncPbPbJetYieldBKGV2 = iEpUFJetPtYieldBKGV2.Clone(histName)
    dJetYieldIEP['BKGV2'] = hIEpIncPbPbJetYieldBKGV2
    
    hChJetV2BKGV2 = CalcChJetV2(oEpUFJetPtYieldBKGV2,iEpUFJetPtYieldBKGV2,centBin,'BKGV2')
    dChJetV2['BKGV2'] = hChJetV2BKGV2
    

    incUFJetPtYieldBKGV2 = incCentUFJetTree.FindObject(UFJetName)
    histName = 'hIncPbPbJetYield'+'_'+ unfoldType +'_'+ lSysType[0]
    hIncPbPbJetYieldBKGV2 = incUFJetPtYieldBKGV2.Clone(histName)
    dJetYield['BKGV2'] = hIncPbPbJetYieldBKGV2
    # == e == BKGV2Fit  ===================================

    # == s == V0CFit  ===================================
    treeName = 'lCent' + str(centBin)
    oEpCentTreeV0C = oEpTreeV0C.FindObject(treeName)
    iEpCentTreeV0C = iEpTreeV0C.FindObject(treeName)
    incCentTreeV0C = incTreeV0C.FindObject(treeName)
    treeName = 'lUnfoldedJets'
    oEpCentUFJetTreeV0C = oEpCentTreeV0C.FindObject(treeName)
    iEpCentUFJetTreeV0C = iEpCentTreeV0C.FindObject(treeName)
    incCentUFJetTreeV0C = incCentTreeV0C.FindObject(treeName)

    UFJetName = baseUFJetName +'_'+ unfoldType +'_'+ lSysType[0]
    oEpUFJetPtYieldV0C = oEpCentUFJetTreeV0C.FindObject(UFJetName)
    iEpUFJetPtYieldV0C = iEpCentUFJetTreeV0C.FindObject(UFJetName)

    histName = 'hOEpncPbPbJetYield'+'_'+ unfoldType +'_'+ 'V0C'
    hOEpUFJetPtYieldV0C = oEpUFJetPtYieldV0C.Clone(histName)
    dJetYieldOEP['V0C'] = hOEpUFJetPtYieldV0C
    histName = 'hIEpIncPbPbJetYield'+'_'+ unfoldType +'_'+ 'V0C'
    hIEpIncPbPbJetYieldV0C = iEpUFJetPtYieldV0C.Clone(histName)
    dJetYieldIEP['V0C'] = hIEpIncPbPbJetYieldV0C
    
    hChJetV2V0C = CalcChJetV2(oEpUFJetPtYieldV0C,iEpUFJetPtYieldV0C,centBin,'V0C')
    dChJetV2['V0C'] = hChJetV2V0C

    incUFJetPtYieldV0C = incCentUFJetTree.FindObject(UFJetName)
    histName = 'hIncPbPbJetYield'+'_'+ unfoldType +'_'+ lSysType[0]
    hIncPbPbJetYieldV0C = incUFJetPtYieldV0C.Clone(histName)
    dJetYield['V0C'] = hIncPbPbJetYieldV0C
    # == e == V0CFit  ===================================

    # == s == V0AFit  ===================================
    treeName = 'lCent' + str(centBin)
    oEpCentTreeV0A = oEpTreeV0A.FindObject(treeName)
    iEpCentTreeV0A = iEpTreeV0A.FindObject(treeName)
    incCentTreeV0A = incTreeV0A.FindObject(treeName)
    treeName = 'lUnfoldedJets'
    oEpCentUFJetTreeV0A = oEpCentTreeV0A.FindObject(treeName)
    iEpCentUFJetTreeV0A = iEpCentTreeV0A.FindObject(treeName)
    incCentUFJetTreeV0A = incCentTreeV0A.FindObject(treeName)

    UFJetName = baseUFJetName +'_'+ unfoldType +'_'+ lSysType[0]
    oEpUFJetPtYieldV0A = oEpCentUFJetTreeV0A.FindObject(UFJetName)
    iEpUFJetPtYieldV0A = iEpCentUFJetTreeV0A.FindObject(UFJetName)
    
    histName = 'hOEpncPbPbJetYield'+'_'+ unfoldType +'_'+ 'V0A'
    hOEpUFJetPtYieldV0A = oEpUFJetPtYieldV0A.Clone(histName)
    dJetYieldOEP['V0A'] = hOEpUFJetPtYieldV0A
    histName = 'hIEpIncPbPbJetYield'+'_'+ unfoldType +'_'+ 'V0A'
    hIEpIncPbPbJetYieldV0A = iEpUFJetPtYieldV0A.Clone(histName)
    dJetYieldIEP['V0A'] = hIEpIncPbPbJetYieldV0A

    hChJetV2V0A = CalcChJetV2(oEpUFJetPtYieldV0A,iEpUFJetPtYieldV0A,centBin,'V0A')
    dChJetV2['V0A'] = hChJetV2V0A

    incUFJetPtYieldV0A = incCentUFJetTree.FindObject(UFJetName)
    histName = 'hIncPbPbJetYield'+'_'+ unfoldType +'_'+ lSysType[0]
    hIncPbPbJetYieldV0A = incUFJetPtYieldV0A.Clone(histName)
    dJetYield['V0A'] = hIncPbPbJetYieldV0A
    # == e == V0AFit  ===================================


    lDTargedObs = [dChJetV2, dJetYield, dJetYield, dJetYieldOEP, dJetYieldIEP]
    lRawJetObs = [hRawChJetV2, incRawJetPtYield]
    
    return lDTargedObs, lRawJetObs

def GetSysUncErrObj(lAllValWithErr, dTargedObs, valKindBin, pTBinArray, centBin, lOMainTree, oFile):

    # == s == Calculate systematic uncertainty  ===========================
    valDire = ['ChJetV2', 'ChJetRAA', 'ChJetInc', 'ChJetOEP', 'ChJetIEP']

    dhAbsSysUnc = {}
    dlAbsSysUncErr  = {}
    dhRelaSysUnc = {}
    dlRelaSysUncErr = {}
    if uncAbs == 1:
        dhAbsSysUnc['prior'], dlAbsSysUncErr['prior'] = AbsSysUncErrCalc1(dTargedObs['Nominal'],\
            dTargedObs['DiffPrior'], 'prior', 1, pTBinArray, centBin)
        dhAbsSysUnc['ptBinRange'], dlAbsSysUncErr['ptBinRange']  = AbsSysUncErrCalc2(dTargedObs['Nominal'], \
            dTargedObs['DetLvLcut'], dTargedObs['DetLvUcut'], 'ptBinRange', 2, pTBinArray, centBin)
        dhAbsSysUnc['iteration'], dlAbsSysUncErr['iteration'] = AbsSysUncErrCalc2(dTargedObs['Nominal'], \
            dTargedObs['iteKM1'], dTargedObs['iteKP1'], 'iteration', 3, pTBinArray,centBin)

        dhAbsSysUnc['trackEff'], dlAbsSysUncErr['trackEff'] = AbsSysUncErrCalc1(dTargedObs['Nominal'], \
            dTargedObs['trackEff94'], 'trackEff', 4, pTBinArray, centBin)
        
        dhAbsSysUnc['bkgWay'], dlAbsSysUncErr['bkgWay'] \
            = AbsSysUncErrCalc2(dTargedObs['Nominal'], dTargedObs['BKGNoFit'], dTargedObs['BKGV2'],\
                'bkgWay', 5, pTBinArray,centBin)

        print(dTargedObs)

        dhAbsSysUnc['V0Diff'], dlAbsSysUncErr['V0Diff'] \
            = AbsSysUncErrCalc2(dTargedObs['Nominal'], dTargedObs['V0C'], dTargedObs['V0A'],\
                'V0Diff', 6, pTBinArray,centBin)

        # $$$$ Temp Raw Uncertain Error $$$$
        # lBkgSysUncErrEachCent = array('d', lBkgSysUncErr[centBin])
        # dhAbsSysUnc['bkgWay'], dlAbsSysUncErr['bkgWay'] = AbsSysUncErrCalcTempRaw(dTargedObs['Nominal'], \
        #     lBkgSysUncErrEachCent, 'bkgWay', 5, pTBinArray, centBin)
        # lV0DiffSysUncErrEachCent = array('d', lV0DiffSysUncErr[centBin])
        # dhAbsSysUnc['V0Diff'], dlAbsSysUncErr['V0Diff'] = AbsSysUncErrCalcTempRaw(dTargedObs['Nominal'], \
        #    lV0DiffSysUncErrEachCent, 'bV0Diff', 6, pTBinArray, centBin)

        
        SumUpAllAbsSysErr(lAllValWithErr, dTargedObs['Nominal'], dhAbsSysUnc, dlAbsSysUncErr ,pTBinArray,centBin)
        
        hNameSysErr = 'h'+valDire[valKindBin]+'SysErr' + '_Cent' + str(centBin)
        hTargedObsSysUncErr = SetupTGraphErrors(dTargedObs['Nominal'],dlAbsSysUncErr['total'],\
            hNameSysErr, pTBinArray, centBin)

    elif uncAbs == 0:
        dhRelaSysUnc['prior'], dlRelaSysUncErr['prior'] = RelativeSysUncErrCalc1(dTargedObs['Nominal'],\
            dTargedObs['DiffPrior'], 'prior', 1, pTBinArray, centBin)
        dhRelaSysUnc['ptBinRange'], dlRelaSysUncErr['ptBinRange']  = RelativeSysUncErrCalc2(dTargedObs['Nominal'], \
            dTargedObs['DetLvLcut'], dTargedObs['DetLvUcut'], 'ptBinRange', 2, pTBinArray, centBin)
        dhRelaSysUnc['iteration'], dlRelaSysUncErr['iteration'] = RelativeSysUncErrCalc2(dTargedObs['Nominal'], \
            dTargedObs['iteKM1'], dTargedObs['iteKP1'], 'iteration', 3, pTBinArray,centBin)

        dhRelaSysUnc['trackEff'], dlRelaSysUncErr['trackEff'] = RelativeSysUncErrCalc1(dTargedObs['Nominal'], \
            dTargedObs['trackEff94'], 'trackEff', 4, pTBinArray, centBin)
        
        print('CheeeeeeeeeeeeeeeeecKuma002')
        SumUpAllRelaSysErr(lAllValWithErr, dTargedObs['Nominal'], dhRelaSysUnc, dlRelaSysUncErr ,pTBinArray,centBin)
        
        hNameSysErr = 'h'+valDire[valKindBin]+'SysErr' + '_Cent' + str(centBin)
        hTargedObsSysUncErr = SetupTGraphErrors(dTargedObs['Nominal'], dlRelaSysUncErr ['total'],\
            hNameSysErr, pTBinArray, centBin)
    
    # == s == Calculate systematic uncertainty  ===========================

    # == s == plot systematic error #######################################
    lhSysUncError = list()
    if uncAbs == 1: 
        for key in dhAbsSysUnc.values(): lhSysUncError.append(key)
    elif uncAbs == 0: 
        for key in dhRelaSysUnc.values(): lhSysUncError.append(key)
    lhSysUncError.reverse()
    lNameSysUncError = ['Reweighting Prior', 'Truncation', 'Iterations', \
        'Tracking Efficiency', 'Background Estiamtion Way', 'V0 Detectors','Total']
    lNameSysUncError.reverse()
    nBin = (pTBinArray[-1] - pTBinArray[0])/10
    lRange = [[pTBinArray[0], pTBinArray[-1]],[0.0, 2.]]
    lTitle = ['hBlankForSysError','#it{p}_{T, ch jet} [GeV/#it{c}]','systematic error ratio']
    hBlankEachCentSysErr = genePlotSets.makeBlank1DHist(nBin, lRange, lTitle)
    legALICE = ROOT.TLegend(0.15,0.8,0.25,0.9,'')
    legJet = ROOT.TLegend(0.15,0.65,0.25,0.75,'')
    legSysUncError = ROOT.TLegend(0.6,0.7,0.8,0.88,'systematic error ratio')
    genePlotSets.addALICELegend(hBlankEachCentSysErr, legALICE, lCentLabel[centBin], \
        'Work in progress', 0.03)
    genePlotSets.addJetLegend(hBlankEachCentSysErr, legJet, 0.2, leadingTrackPtCut, 0.03)        
    genePlotSets.addSomeHistsLegend(legSysUncError, lhSysUncError, lNameSysUncError, 0.03)
    
    label = valDire[valKindBin]+'/SysUncError'+'_TrackPtCut'+str(leadingTrackPtCut)+'_Cent'+str(centBin)
    plotStyleList = list()
    # [kViolet+1, kBlue, kAzure+1, kGreen+3, kOrange+8, kPink+9]
    colorList = [880+1, 600, 860+1, 416+3, 800-3, 800+8, 900+9]
    for l in range(0, len(lhSysUncError)): plotStyleList.append(0)
    genePlotSets.overwrightSomePlots(hBlankEachCentSysErr, lhSysUncError, \
        0, legALICE,legJet,legSysUncError, lFillColor, plotStyleList, label, '', outFinalFileDir)
    # == e == plot systematic error #######################################

    oFile.cd()
    for key in dTargedObs.values(): lOMainTree[valKindBin][centBin][0].Add(key)
    for key in dhAbsSysUnc.values(): lOMainTree[valKindBin][centBin][1].Add(key)
    # for key in dhRelaSysUnc.values(): lOMainTree[valKindBin][centBin][1].Add(key)
    lOMainTree[valKindBin][centBin][2].Add(hTargedObsSysUncErr)
    
    return dhAbsSysUnc, dlAbsSysUncErr
    # return dhRelaSysUnc, dlRelaSysUncErr 

def AbsSysUncErrCalc1(hNominal, hDiff, diffName, diffNum, pTBinArray ,centBin):
    histName = "hNominalCP_" + diffName + '_Cent'+ str(centBin)
    hNominalCP = hNominal.Clone(histName)
    
    histname = "hDiff_" + diffName + '_Cent'+ str(centBin)
    title = histname + ";#it{p}_{T,ch jet} (GeV/#it{c}); Ratio Systematic Uncertainty [%]"
    
    nBins = len(pTBinArray)-1
    hSysUncPer = ROOT.TH1D(histname, title, nBins, pTBinArray)
    
    lSysUncErr = list()
    for iBin in range(1, hNominalCP.GetNbinsX() + 1):
        nomiVal = hNominalCP.GetBinContent(iBin)
        deltaVal = nomiVal - hDiff.GetBinContent(iBin)
        SysUnc = np.sqrt((deltaVal/nomiVal)*(deltaVal/nomiVal))
        # SysUnc = np.sqrt((deltaVal)*(deltaVal))*100
        fillVal = hNominalCP.GetXaxis().GetBinCenter(iBin)

        if (fillVal > pTBinArray[0]) and (fillVal < pTBinArray[-1]):
            hSysUncPer.Fill(fillVal, SysUnc)

            SysUncErr = np.sqrt(deltaVal*deltaVal)
            lSysUncErr.append(SysUncErr)
            # lSysUncErr.append(SysUnc)

    for iBin in range(1, hSysUncPer.GetNbinsX() + 1): hSysUncPer.SetBinError(iBin, 0.)
    hSysUncPer.GetXaxis().SetTitle("#it{p}_{T, ch jet} [GeV/#it{c}]")
    hSysUncPer.GetYaxis().SetTitle("Ratio Systematic Uncertainty")
    # hSysUncPer.SetMarkerStyle(21)
    hSysUncPer.SetMarkerColor(lFillColor[diffNum])
    hSysUncPer.SetLineColor(lFillColor[diffNum])
    hSysUncPer.SetLineWidth(4)
    
    lSysUncErr = array('d',lSysUncErr)

    return hSysUncPer, lSysUncErr

def AbsSysUncErrCalc2(hNominal, hDiff1, hDiff2, diffName, diffNum, pTBinArray, centBin):
    histName = "hNominalCP_" + diffName + '_Cent'+ str(centBin)
    hNominalCP = hNominal.Clone(histName)
    
    histname = "hDiff_" + diffName + '_Cent'+ str(centBin)
    title = histname + ";#it{p}_{T,ch jet} (GeV/#it{c}) Ratio Systematic Uncertainty [%]"
    
    nBins = len(pTBinArray)-1
    hSysUncPer = ROOT.TH1D(histname, title, nBins, pTBinArray)
    
    lSysUncErr = list()
    for iBin in range(1, hNominalCP.GetNbinsX() + 1):
        nomiVal = hNominalCP.GetBinContent(iBin)
        deltaVal1 = nomiVal - hDiff1.GetBinContent(iBin)
        deltaVal2 = nomiVal - hDiff2.GetBinContent(iBin)
        SysUnc1 = (deltaVal1/nomiVal)*(deltaVal1/nomiVal)
        SysUnc2 = (deltaVal2/nomiVal)*(deltaVal2/nomiVal)
        SysUncTot = np.sqrt(SysUnc1 + SysUnc2)
        # SysUnc = np.sqrt((deltaVal)*(deltaVal))*100
        fillVal = hNominalCP.GetXaxis().GetBinCenter(iBin)

        if (fillVal > pTBinArray[0]) and (fillVal < pTBinArray[-1]):
            hSysUncPer.Fill(fillVal, SysUncTot)

            SysUncErr = np.sqrt(deltaVal1*deltaVal1+deltaVal2*deltaVal2)
            lSysUncErr.append(SysUncErr)
            # lSysUncErr.append(SysUncTot)

    for iBin in range(1, hSysUncPer.GetNbinsX() + 1): hSysUncPer.SetBinError(iBin, 0.)
    hSysUncPer.GetXaxis().SetTitle("#it{p}_{T, ch jet} [GeV/#it{c}]")
    hSysUncPer.GetYaxis().SetTitle("Ratio Systematic Uncertainty")
    # hSysUncPer.SetMarkerStyle(21)
    hSysUncPer.SetMarkerColor(lFillColor[diffNum])
    hSysUncPer.SetLineColor(lFillColor[diffNum])
    hSysUncPer.SetLineWidth(4)

    lSysUncErr = array('d',lSysUncErr)

    return hSysUncPer, lSysUncErr

def AbsSysUncErrCalcTempRaw(hNominal, lTemRawSysErr, diffName, diffNum, pTBinArray, centBin):
    histName = "hNominalCP_" + diffName + '_Cent'+ str(lBkgSysUncErr)
    hNominalCP = hNominal.Clone(histName)
    
    histname = "hDiff_" + diffName + '_Cent'+ str(centBin)
    title = histname + ";#it{p}_{T,ch jet} (GeV/#it{c}); Ratio Systematic Uncertainty [%]"
    
    nBins = len(pTBinArray)-1
    hSysUncPer = ROOT.TH1D(histname, title, nBins, pTBinArray)
    
    lSysUncErr = list()
    print('centBin: '+str(centBin))
    print('hNominalCP.GetNbinsX():'+ str(hNominalCP.GetNbinsX()))
    print(pTBinArray)
    fillPtBin = 0
    for iBin in range(1, hNominalCP.GetNbinsX() + 1):
        fillVal = hNominalCP.GetXaxis().GetBinCenter(iBin)

        if (fillVal > pTBinArray[0]) and (fillVal < pTBinArray[-1]):
            nomiVal = hNominalCP.GetBinContent(iBin)
            print('pT = '+str(fillVal))
            print(lTemRawSysErr[fillPtBin])
            deltaVal = nomiVal*(lTemRawSysErr[fillPtBin])
            SysUnc = np.sqrt((deltaVal/nomiVal)*(deltaVal/nomiVal))
            # SysUnc = np.sqrt((deltaVal)*(deltaVal))*100

            hSysUncPer.Fill(fillVal, SysUnc)

            SysUncErr = np.sqrt(deltaVal*deltaVal)
            lSysUncErr.append(SysUncErr)
            # lSysUncErr.append(SysUnc)

            fillPtBin += 1

    for iBin in range(1, hSysUncPer.GetNbinsX() + 1): hSysUncPer.SetBinError(iBin, 0.)
    hSysUncPer.GetXaxis().SetTitle("#it{p}_{T, ch jet} [GeV/#it{c}]")
    hSysUncPer.GetYaxis().SetTitle("Ratio Systematic Uncertainty")
    # hSysUncPer.SetMarkerStyle(21)
    hSysUncPer.SetMarkerColor(lFillColor[diffNum])
    hSysUncPer.SetLineColor(lFillColor[diffNum])
    hSysUncPer.SetLineWidth(4)
    
    lSysUncErr = array('d',lSysUncErr)

    return hSysUncPer, lSysUncErr


def RelativeSysUncErrCalc1(hNominal, hDiff, diffName, diffNum, pTBinArray ,centBin):
    histName = "hNominalRelativeCP_" + diffName + '_Cent'+ str(centBin)
    hNominalCP = hNominal.Clone(histName)
    
    histname = "hRelativeDiff_" + diffName + '_Cent'+ str(centBin)
    title = histname + ";#it{p}_{T,ch jet} (GeV/#it{c}); Ratio Systematic Uncertainty [%]"
    
    nBins = len(pTBinArray)-1
    hSysUncPer = ROOT.TH1D(histname, title, nBins, pTBinArray)
    
    lSysUncErr = list()
    for iBin in range(1, hNominalCP.GetNbinsX() + 1):
        nomiVal = hNominalCP.GetBinContent(iBin)
        deltaVal = hDiff.GetBinContent(iBin) - nomiVal
        SysUnc = deltaVal/nomiVal
        fillVal = hNominalCP.GetXaxis().GetBinCenter(iBin)

        if (fillVal > pTBinArray[0]) and (fillVal < pTBinArray[-1]):
            hSysUncPer.Fill(fillVal, SysUnc)

            SysUncErr = np.sqrt(deltaVal*deltaVal)
            lSysUncErr.append(SysUncErr)
            # lSysUncErr.append(SysUnc)

    for iBin in range(1, hSysUncPer.GetNbinsX() + 1): hSysUncPer.SetBinError(iBin, 0.)
    hSysUncPer.GetXaxis().SetTitle("#it{p}_{T, ch jet} [GeV/#it{c}]")
    hSysUncPer.GetYaxis().SetTitle("Ratio Systematic Uncertainty")
    # hSysUncPer.SetMarkerStyle(21)
    hSysUncPer.SetMarkerColor(lFillColor[diffNum])
    hSysUncPer.SetLineColor(lFillColor[diffNum])
    hSysUncPer.SetLineWidth(4)
    
    lSysUncErr = array('d',lSysUncErr)

    return hSysUncPer, lSysUncErr

def RelativeSysUncErrCalc2(hNominal, hDiff1, hDiff2, diffName, diffNum, pTBinArray, centBin):
    histName = "hNominalRelativeCP_" + diffName + '_Cent'+ str(centBin)
    hNominalCP = hNominal.Clone(histName)
    
    histname = "hRelativeDiff_" + diffName + '_Cent'+ str(centBin)
    title = histname + ";#it{p}_{T,ch jet} (GeV/#it{c}); Ratio Systematic Uncertainty [%]"
    
    nBins = len(pTBinArray)-1
    hSysUncPer = ROOT.TH1D(histname, title, nBins, pTBinArray)
    
    lSysUncErr = list()
    for iBin in range(1, hNominalCP.GetNbinsX() + 1):
        nomiVal = hNominalCP.GetBinContent(iBin)
        deltaVal1 = nomiVal - hDiff1.GetBinContent(iBin)
        deltaVal2 = nomiVal - hDiff2.GetBinContent(iBin)
        SysUnc1 = deltaVal1/nomiVal
        SysUnc2 = deltaVal2/nomiVal
        SysUncTot = SysUnc1 + SysUnc2
        fillVal = hNominalCP.GetXaxis().GetBinCenter(iBin)

        if (fillVal > pTBinArray[0]) and (fillVal < pTBinArray[-1]):
            hSysUncPer.Fill(fillVal, SysUncTot)

            SysUncErr = deltaVal1+deltaVal2
            lSysUncErr.append(SysUncErr)
            # lSysUncErr.append(SysUncTot)

    for iBin in range(1, hSysUncPer.GetNbinsX() + 1): hSysUncPer.SetBinError(iBin, 0.)
    hSysUncPer.GetXaxis().SetTitle("#it{p}_{T, ch jet} [GeV/#it{c}]")
    hSysUncPer.GetYaxis().SetTitle("Ratio Systematic Uncertainty")
    # hSysUncPer.SetMarkerStyle(21)
    hSysUncPer.SetMarkerColor(lFillColor[diffNum])
    hSysUncPer.SetLineColor(lFillColor[diffNum])
    hSysUncPer.SetLineWidth(4)

    lSysUncErr = array('d',lSysUncErr)

    return hSysUncPer, lSysUncErr

def CalcChJetV2(hOEPUFJet, hIEPUFJet, centBin, label):
    hChJetV2Nume = hIEPUFJet - hOEPUFJet
    hChJetV2Deno = hIEPUFJet + hOEPUFJet
    
    # == s == stat error proper gation Calculation      ========================
    lv2Error = list()
    for iBin in range(1, hChJetV2Nume.GetNbinsX() + 1):
        # numeErr = hChJetV2Nume.GetBinError(iBin)
        # denoErr = hChJetV2Deno.GetBinError(iBin)
        numeContent = hChJetV2Nume.GetBinContent(iBin)
        denoContent = hChJetV2Deno.GetBinContent(iBin)
        # multErr = np.sqrt((numeErr/numeContent)*(numeErr/numeContent)\
        #     + (denoErr/denoContent)*(denoErr/denoContent))
        
        IEPUFJetContent = hIEPUFJet.GetBinContent(iBin)
        OEPUFJetContent = hOEPUFJet.GetBinContent(iBin)
        IEPUFJetErr = hIEPUFJet.GetBinError(iBin)
        OEPUFJetErr = hOEPUFJet.GetBinError(iBin)
        
        # multErr = np.sqrt((IEPUFJetContent+OEPUFJetContent)*numeContent+(IEPUFJetErr+OEPUFJetErr)*denoContent)

        DenoErr = (IEPUFJetContent+OEPUFJetContent)*(IEPUFJetContent+OEPUFJetContent)\
            *(IEPUFJetContent+OEPUFJetContent)*(IEPUFJetContent+OEPUFJetContent)
        ONumeErr = 4*(OEPUFJetContent*OEPUFJetContent)*(IEPUFJetErr*IEPUFJetErr)
        INumeErr = 4*(IEPUFJetContent*IEPUFJetContent)*(OEPUFJetErr*OEPUFJetErr) 
        CoEffError =  TMath.Pi() / (4 * lPsi2Reso[centBin])
        multErr = np.sqrt(ONumeErr/DenoErr + INumeErr/DenoErr)

        lv2Error.append(multErr)
        # print('Err multErr:'+str(multErr)+', nume: '+str(numeErr)+', deno: '+str(denoErr)+' numeRatio:'+str(numeErr/numeContent)+', denoRatio:'+str(denoErr/denoContent))
    # == e == stat error proper gation Calculation      ========================

    hChJetV2Nume.Divide(hChJetV2Deno)
    jetV2Coeff =  TMath.Pi() / (4 * lPsi2Reso[centBin])
    hChJetV2Nume.Scale(jetV2Coeff)
    
    histName  = 'hChJetV2'+'_Cent'+str(centBin)+'_'+label
    histTitle = 'hChJetV2'+'_Cent'+str(centBin)+'_'+label
    hChJetV2Nume.SetName(histName)
    hChJetV2Nume.SetTitle(histTitle)
    hChJetV2Nume.GetXaxis().SetTitle('#it{p}_{T, ch jet} [GeV/#it{c}]')
    hChJetV2Nume.GetYaxis().SetTitle('#it{v}_{2}^{jet}')

    for iBin in range(1, hChJetV2Nume.GetNbinsX() + 1):
        # error = lv2Error[iBin-1]/hChJetV2Nume.GetBinContent(iBin)
        error = lv2Error[iBin-1]
        hChJetV2Nume.SetBinError(iBin, error)

    return hChJetV2Nume

def CalcChRAA(lAllValWithErr, lOMainTree):
    for centBin in range(0, numOfCentBin):
        lAllIncPbPbWithErr = lAllValWithErr[1][centBin]
        
        # HEP https://www.hepdata.net/record/ins1733689  bins: 18
        # lPPCrossVal = [0.222851, 0.102591, 0.053269, 0.030328, 0.0185155, 0.010319, \
        #     0.00503213, 0.00279515, 0.00164066, 0.00102357, 0.000529236, 0.000219935, \
        #         7.91E-05, 2.36E-05, 8.71E-06, 3.75E-06, 1.54E-06, 5.99E-07]
        # lPPBinArray = [5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100]
        # lPPBinCent = [5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17, 19, 22.5,\
        #     27.5, 35, 45, 55, 65, 77.5, 92.5]
        # lPPBinWidth = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 5, 5, 10, 10, 10, 10, 15, 15]
        # lPPStatErr = [0.000176827, 0.000116457, 8.18E-05, 5.71E-05, 3.82E-05, 2.27E-05,\
        #     1.69E-05, 1.83E-05, 1.94E-05, 1.96E-05, 1.51E-05, 8.77E-06, 4.18E-06, \
        #         1.58E-06, 7.07E-07, 3.57E-07, 1.66E-07, 7.11E-08]
        # lPPSysErr = [0.00824549, 0.00379588, 0.00197095, 0.00112214, 0.000685073, \
        #     0.000381804, 0.000191221, 0.000111806, 6.73E-05, 4.40E-05, 2.49E-05, \
        #         1.12E-05, 4.75E-06, 1.72E-06, 7.49E-07, 3.71E-07, 1.79E-07, 7.84E-08]

        lPPCrossVal = [7.91E-05, 2.36E-05, 8.71E-06, 3.75E-06, 1.54E-06, 5.99E-07]
        lPPBinArray = [30, 40, 50, 60, 70, 85, 100]
        lPPBinCent = [35, 45, 55, 65, 77.5, 92.5]
        lPPBinWidth = [5, 5, 5, 5, 7.5, 7.5]
        lPPStatErr = [4.18E-06, 1.58E-06, 7.07E-07, 3.57E-07, 1.66E-07, 7.11E-08]
        lPPSysErr = [4.75E-06, 1.72E-06, 7.49E-07, 3.71E-07, 1.79E-07, 7.84E-08]

        lPPCrossVal = array('d', lPPCrossVal)
        lPPBinArray = array('d', lPPBinArray)
        lPPBinCent = array('d', lPPBinCent)
        lPPBinWidth = array('d', lPPBinWidth)
        lPPStatErr = array('d', lPPStatErr) 
        lPPSysErr = array('d', lPPSysErr)

        lChJetRAA = list()
        lChJetRAAStat = list()
        lChJetRAASys = list()
        nPPBins = len(lPPCrossVal)
        
        nPbPbBins = len(lAllIncPbPbWithErr[0])
        diffNumOfBins = nPPBins - nPbPbBins
        
        lChJetRAABinCent = list()
        lChJetRAABinWidth = list()
        lChJetRAABinArray = list()
        iPbPbBin = 0
        for iBin in range(0, nPPBins):
            if iBin < diffNumOfBins: continue
            
            chJetRAA = lAllIncPbPbWithErr[0][iPbPbBin]/lPPCrossVal[iBin]
            
            incPbPbYieldStatErrRatio = lAllIncPbPbWithErr[3][iPbPbBin]/lAllIncPbPbWithErr[0][iPbPbBin]
            incPPYieldStatErrRatio = lPPStatErr[iBin]/lPPCrossVal[iBin]
            chJetRAAStatErr = np.sqrt(incPbPbYieldStatErrRatio *incPbPbYieldStatErrRatio \
                + incPPYieldStatErrRatio*incPPYieldStatErrRatio)
            
            incPbPbYieldSysErrRatio = lAllIncPbPbWithErr[5][iPbPbBin]/lAllIncPbPbWithErr[0][iPbPbBin]
            incPPYieldSysErrRatio = lPPSysErr[iBin]/lPPCrossVal[iBin]
            chJetRAASysErr = np.sqrt(incPbPbYieldSysErrRatio*incPbPbYieldSysErrRatio \
                + incPPYieldSysErrRatio*incPPYieldSysErrRatio)
            
            lChJetRAA.append(chJetRAA)
            lChJetRAABinCent.append(lPPBinCent[iBin])
            lChJetRAABinWidth.append(lPPBinWidth[iBin])
            lChJetRAAStat.append(chJetRAAStatErr)
            lChJetRAABinArray.append(lPPBinArray[iBin])
            lChJetRAASys.append(chJetRAASysErr)
            
            iPbPbBin += 1

        lChJetRAABinArray.append(lPPBinArray[-1])

        lChJetRAA = array('d',lChJetRAA)
        lChJetRAABinCent = array('d', lChJetRAABinCent)
        lChJetRAABinWidth = array('d', lChJetRAABinWidth)
        lChJetRAAStat = array('d',lChJetRAAStat)
        lChJetRAABinArray = array('d', lChJetRAABinArray)
        lChJetRAASys = array('d',lChJetRAASys)

        lAllValWithErr[2][centBin].append(lChJetRAA)
        lAllValWithErr[2][centBin].append(lChJetRAABinCent)
        lAllValWithErr[2][centBin].append(lChJetRAABinWidth)
        lAllValWithErr[2][centBin].append(lChJetRAAStat)
        lAllValWithErr[2][centBin].append(lChJetRAABinArray)
        lAllValWithErr[2][centBin].append(lChJetRAASys)
        
        histName  = 'hChJetRAA'+'_Cent'+str(centBin)
        hChJetRAA = ROOT.TH1D(histName, histName, nPbPbBins-1, lChJetRAABinArray)
        hChJetRAA.GetXaxis().SetTitle('#it{p}_{T, ch jet} [GeV/#it{c}]')
        hChJetRAA.GetYaxis().SetTitle('#it{R}_{AA}^{ ch jet}')
        for ptBin in range(0, len(lChJetRAA)):
            hChJetRAA.SetBinContent(ptBin+1, lChJetRAA[ptBin])
            hChJetRAA.SetBinError(ptBin+1, lChJetRAAStat[ptBin])
        lOMainTree[centBin][0].Add(hChJetRAA)

        hChRAAErr= ROOT.TGraphErrors(nPbPbBins, lChJetRAABinCent, lChJetRAA, \
            lChJetRAABinWidth, lChJetRAASys)
        lOMainTree[centBin][2].Add(hChRAAErr)

def SumUpAllAbsSysErr(lAllValWithErr, hNominal, dhAbsSysUnc, dlAbsSysUncErr , pTBinArray, centBin):
    histName = "hNominalCP_" +'ToTal' + '_Cent'+ str(centBin)
    hNominalCP = hNominal.Clone(histName)

    lNominalVal = list()
    lNominalBinCent = list()
    lNominalBinWidth = list()
    lNominalValStat = list()

    lPriorSysUncErr = dlAbsSysUncErr ['prior']
    lTotalSysUncErr = list()
    
    for ptBin in range(0, len(lPriorSysUncErr)):
        totalSysUncErr = 0.
        for lSysUncErr in dlAbsSysUncErr .values():
            # print('EachError'+str(lSysUncErr[ptBin]))
            totalSysUncErr += lSysUncErr[ptBin]*lSysUncErr[ptBin]
        totalSysUncErr = np.sqrt(totalSysUncErr)
        # print('Total Error: ' +str(totalSysUncErr))
        lTotalSysUncErr.append(totalSysUncErr)

    lTotalSysUncErr = array('d',lTotalSysUncErr)

    histname = "hTotalSysUnc" + '_CentBin' +str(centBin)
    title = histname + ";#it{p}_{T,ch jet} (GeV/#it{c}) Relative Systematic Uncertainty [%]"
    
    nBins = len(pTBinArray)-1
    hTotSysUncPer = ROOT.TH1D(histname, title, nBins, pTBinArray)
    tempPtBin = 0
    for iBin in range(1, hNominalCP.GetNbinsX() + 1):
        nomiVal = hNominalCP.GetBinContent(iBin)
        fillVal = hNominalCP.GetXaxis().GetBinCenter(iBin)

        if (fillVal > pTBinArray[0]) and (fillVal < pTBinArray[-1]):
            SysUnc = lTotalSysUncErr[tempPtBin]
            SysUnc /= np.sqrt(nomiVal*nomiVal)
            hTotSysUncPer.Fill(fillVal, SysUnc)
            tempPtBin += 1
            
            lNominalVal.append(nomiVal)
            lNominalBinCent.append(fillVal)
            binWidth = hNominalCP.GetBinWidth(iBin)*0.5
            lNominalBinWidth.append(binWidth)
            nomiValErr = hNominalCP.GetBinError(iBin)
            lNominalValStat.append(nomiValErr)
            
    for iBin in range(1, hTotSysUncPer.GetNbinsX() + 1): hTotSysUncPer.SetBinError(iBin, 0.)
    hTotSysUncPer.GetXaxis().SetTitle("#it{p}_{T, ch jet} [GeV/#it{c}]")
    hTotSysUncPer.GetYaxis().SetTitle("Relative Systematic Uncertainty")
    # hSysUncPer.SetMarkerStyle(21)
    hTotSysUncPer.SetMarkerColor(1)
    hTotSysUncPer.SetLineColor(1)
    hTotSysUncPer.SetLineWidth(4)

    dhAbsSysUnc['total'] = hTotSysUncPer
    dlAbsSysUncErr ['total'] = lTotalSysUncErr
    
    lNominalVal      = array('d', lNominalVal)
    lNominalBinCent  = array('d', lNominalBinCent)
    lNominalBinWidth = array('d', lNominalBinWidth)
    lNominalValStat  = array('d', lNominalValStat)

    lAllValWithErr[centBin].append(lNominalVal)
    lAllValWithErr[centBin].append(lNominalBinCent)
    lAllValWithErr[centBin].append(lNominalBinWidth)
    lAllValWithErr[centBin].append(lNominalValStat)
    lAllValWithErr[centBin].append(pTBinArray)
    lAllValWithErr[centBin].append(lTotalSysUncErr)

    print('Centrality Bin ' + str(centBin))
    print(dlAbsSysUncErr)

def SumUpAllRelaSysErr(lAllValWithErr, hNominal, dhRelaSysUnc, dlRelaSysUncErr , pTBinArray, centBin):
    histName = "hNominalRelativeCP_" +'ToTal' + '_Cent'+ str(centBin)
    hNominalCP = hNominal.Clone(histName)

    lNominalVal = list()
    lNominalBinCent = list()
    lNominalBinWidth = list()
    lNominalValStat = list()

    lPriorSysUncErr = dlRelaSysUncErr ['prior']
    lTotalSysUncErr = list()
    
    for ptBin in range(0, len(lPriorSysUncErr)):
        totalSysUncErr = 0.
        for lSysUncErr in dlRelaSysUncErr .values():
            # print('EachError'+str(lSysUncErr[ptBin]))
            totalSysUncErr += lSysUncErr[ptBin]
        # print('Total Error: ' +str(totalSysUncErr))
        lTotalSysUncErr.append(totalSysUncErr)

    lTotalSysUncErr = array('d',lTotalSysUncErr)

    histname = "hTotalSysUnc" + '_CentBin' +str(centBin)
    title = histname + ";#it{p}_{T,ch jet} (GeV/#it{c}); Ratio Systematic Uncertainty [%]"
    
    nBins = len(pTBinArray)-1
    hTotSysUncPer = ROOT.TH1D(histname, title, nBins, pTBinArray)
    tempPtBin = 0
    for iBin in range(1, hNominalCP.GetNbinsX() + 1):
        nomiVal = hNominalCP.GetBinContent(iBin)
        fillVal = hNominalCP.GetXaxis().GetBinCenter(iBin)

        if (fillVal > pTBinArray[0]) and (fillVal < pTBinArray[-1]):
            SysUnc = lTotalSysUncErr[tempPtBin]
            SysUnc /= nomiVal
            hTotSysUncPer.Fill(fillVal, SysUnc)
            tempPtBin += 1
            
            lNominalVal.append(nomiVal)
            lNominalBinCent.append(fillVal)
            # print('fillVal '+str(fillVal))
            binWidth = hNominalCP.GetBinWidth(iBin)*0.5
            lNominalBinWidth.append(binWidth)
            nomiValErr = hNominalCP.GetBinError(iBin)
            lNominalValStat.append(nomiValErr)
            
    for iBin in range(1, hTotSysUncPer.GetNbinsX() + 1): hTotSysUncPer.SetBinError(iBin, 0.)
    hTotSysUncPer.GetXaxis().SetTitle("#it{p}_{T, ch jet} [GeV/#it{c}]")
    hTotSysUncPer.GetYaxis().SetTitle("Ratio Systematic Uncertainty")
    # hSysUncPer.SetMarkerStyle(21)
    hTotSysUncPer.SetMarkerColor(1)
    hTotSysUncPer.SetLineColor(1)
    hTotSysUncPer.SetLineWidth(4)

    dhRelaSysUnc['total'] = hTotSysUncPer
    dlRelaSysUncErr ['total'] = lTotalSysUncErr
    
    lNominalVal      = array('d', lNominalVal)
    lNominalBinCent  = array('d', lNominalBinCent)
    lNominalBinWidth = array('d', lNominalBinWidth)
    lNominalValStat  = array('d', lNominalValStat)

    lAllValWithErr[centBin].append(lNominalVal)
    lAllValWithErr[centBin].append(lNominalBinCent)
    lAllValWithErr[centBin].append(lNominalBinWidth)
    lAllValWithErr[centBin].append(lNominalValStat)
    lAllValWithErr[centBin].append(pTBinArray)
    lAllValWithErr[centBin].append(lTotalSysUncErr)

    print('Centrality Bin ' + str(centBin))
    print(dlRelaSysUncErr)

def SetupTGraphErrors(hNominal, lSysErr, errHistName, pTBinArray, centBin):
    histName = 'hNominalCP' + '_Cent'+ str(centBin)
    hNominalCP = hNominal.Clone(histName)
    lPtBinCenter = list()
    lBinWidth = list()
    lTargetVal = list()
    lStatErr = list()
    tempPtBin = 0

    for iBin in range(1, hNominalCP.GetNbinsX() + 1):
        ptVal = hNominalCP.GetXaxis().GetBinCenter(iBin)

        if (ptVal > pTBinArray[0]) and (ptVal < pTBinArray[-1]):
            lPtBinCenter.append(ptVal)
            binWidth = hNominalCP.GetBinWidth(iBin)*0.5
            lBinWidth.append(binWidth)
            targetVal = hNominalCP.GetBinContent(iBin)
            lTargetVal.append(targetVal)
            statErr = hNominalCP.GetBinError(iBin)
            lStatErr.append(statErr)
            tempPtBin += 0
    
    nBins = len(pTBinArray)-1
    lPtBinCenter = array('d',lPtBinCenter)
    lBinWidth = array('d',lBinWidth)
    lTargetVal = array('d', lTargetVal)
    lStatErr = array('d',lStatErr)
    
    # errHistName, errHistName, \
    hErr= ROOT.TGraphErrors(nBins, 
        lPtBinCenter, lTargetVal, lBinWidth, lSysErr)

    return hErr

def MergeDataAndError(hNominal, lSysErr, finalHistName, pTBinArray, centBin):
    histName = 'hNominalCP' + '_Cent'+ str(centBin)
    hNominalCP = hNominal.Clone(histName)
    lPtBinCenter = list()
    lBinWidth = list()
    lTargetVal = list()
    lStatErr = list()
    tempPtBin = 0
    for iBin in range(1, hNominalCP.GetNbinsX() + 1):
        ptVal = hNominalCP.GetXaxis().GetBinCenter(iBin)

        if (ptVal > pTBinArray[0]) and (ptVal < pTBinArray[-1]):
            lPtBinCenter.append(ptVal)
            binWidth = hNominalCP.GetBinWidth(iBin)*0.5
            lBinWidth.append(binWidth)
            targetVal = hNominalCP.GetBinContent(iBin)
            lTargetVal.append(targetVal)
            statErr = hNominalCP.GetBinError(iBin)
            lStatErr.append(statErr)
            tempPtBin += 0
    
    nBins = len(pTBinArray)-1
    lPtBinCenter = array('d',lPtBinCenter)
    lBinWidth = array('d',lBinWidth)
    lTargetVal = array('d', lTargetVal)
    lStatErr = array('d',lStatErr)
    
    hFinalPlotWithErr = ROOT.TGraphMultiErrors(finalHistName, finalHistName, \
        nBins, lPtBinCenter, lTargetVal, lBinWidth, lBinWidth, lStatErr, lStatErr)
    hFinalPlotWithErr.AddYError(nBins, lSysErr, lSysErr)

    return hFinalPlotWithErr

def PlotFinalFigure1(hFinalVal, color, label):
    
    hFinalVal.SetMarkerStyle(20)
    hFinalVal.SetMarkerColor(color)
    hFinalVal.GetAttLine(0).SetLineColor(color)
    hFinalVal.GetAttLine(0).SetLineWidth(4)
    
    hFinalVal.GetAttLine(1).SetLineColor(color)
    # hFinalVal.GetAttFill(1).SetFillStyle(3001)
    # hFinalChJetV2.GetAttFill(1).SetFillColorAlpha(858, 0.8)
    
    hFinalVal.GetXaxis().SetTitle('#it{p}_{T, ch jet} [GeV/#it{c}]')
    hFinalVal.GetYaxis().SetTitle('#it{v}_{2}^{jet}')

    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    c.cd().SetLeftMargin(0.15)
    c.cd().SetRightMargin(0.15)

    # Set pad and histo arrangement
    myPad = ROOT.TPad("myPad", "The pad",0,0,1,1)
    myPad.SetLeftMargin(0.21)
    myPad.SetTopMargin(0.04)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.15)
    myPad.Draw()
    myPad.cd()

    hFinalVal.Draw("a p s ; ; 5 s=0.5")
    # hFinalChJetV2.Draw("epz same")

    canvName = outFinalFileDir + 'c' + label + '.root'
    c.SaveAs(canvName)
    c.Close()

    return hFinalVal
    # return c

def PlotFinalFigure2(hFinalVal, hSysErr, color, label):
    
    hFinalVal.SetMarkerStyle(20)
    hFinalVal.SetMarkerColor(color)
    hFinalVal.SetLineColor(color)
    hFinalVal.SetLineWidth(4)
    
    hSysErr.SetLineColor(color)
    hSysErr.SetLineWidth(2)
    hSysErr.SetFillStyle(3001)
    hSysErr.SetFillColorAlpha(color, 5.0)
    hSysErr.SetFillColor(color)
    
    hFinalVal.GetXaxis().SetTitle('#it{p}_{T, ch jet} [GeV/#it{c}]')
    hFinalVal.GetYaxis().SetTitle('#it{v}_{2}^{jet}')

    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    c.cd().SetLeftMargin(0.15)
    c.cd().SetRightMargin(0.15)

    # Set pad and histo arrangement
    myPad = ROOT.TPad("myPad", "The pad",0,0,1,1)
    myPad.SetLeftMargin(0.21)
    myPad.SetTopMargin(0.04)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.15)
    myPad.Draw()
    myPad.cd()

    hFinalVal.Draw("EP")
    hSysErr.Draw("5 same")
    
    
    canvName = outFinalFileDir + 'c' + label + '.root'
    c.SaveAs(canvName)
    c.Close()

    return hFinalVal
    # return c

def PlotAllCentTargetObs(lAllValWithErr, valKindBin):
    lhTargetObsAllCent = list()
    lhTargetObsSysErrAllCent = list()
    for centBin in range(0, numOfCentBin):
        ptRangeDict, ptBinArrayDict = PtRangeList.eachJetPtBinDef(valKindBin, centBin)
        pTBinArray = ptBinArrayDict['reported']
        
        print(lAllValWithErr[centBin])
        print('lAllValWithErr[centBin][4]')
        print(lAllValWithErr[centBin][4])
        print('lAllValWithErr[centBin][0]')
        print(lAllValWithErr[centBin][0])
            
        nBins = len(lAllValWithErr[centBin][4])
        histName = ''
        if valKindBin==0: histName = 'hChJetV2_Cent'+str(centBin)+'_ForPlot'
        elif valKindBin==2: histName = 'hChJetRAA_Cent'+str(centBin)+'_ForPlot'
        hTargetObs = ROOT.TH1D(histName, histName, nBins-1, lAllValWithErr[centBin][4])
        hTargetObs.GetXaxis().SetTitle('#it{p}_{T, ch jet} [GeV/#it{c}]')
        yTitle = ''
        if valKindBin==0: yTitle = '#it{v}_{2}^{ ch jet}'
        elif valKindBin==2: yTitle = '#it{R}_{AA}^{ ch jet}'
        hTargetObs.GetYaxis().SetTitle(yTitle)
        for ptBin in range(0, nBins-1):
            hTargetObs.SetBinContent(ptBin+1, lAllValWithErr[centBin][0][ptBin])
            hTargetObs.SetBinError(ptBin+1, lAllValWithErr[centBin][3][ptBin])

        lhTargetObsAllCent.append(hTargetObs)

        hTargetObsSysErr= ROOT.TGraphErrors(nBins, lAllValWithErr[centBin][1], \
            lAllValWithErr[centBin][0], lAllValWithErr[centBin][2], lAllValWithErr[centBin][5])
        # [kViolet+1, kAzure+1, kOrange+8, kGreen+3,  kPink+9]
        colorList = [880+1, 860+7, 800+7, 416+2, 800-3,  900+9]
        color = colorList[centBin]
        hTargetObsSysErr.SetLineColor(color)
        hTargetObsSysErr.SetLineWidth(2)
        hTargetObsSysErr.SetFillStyle(3001)
        hTargetObsSysErr.SetFillColorAlpha(color, 5.0)
        hTargetObsSysErr.SetFillColor(color)
        lhTargetObsSysErrAllCent.append(hTargetObsSysErr)
    
    nBin = 120
    lRange = [[30, 150],[-0.01, 0.2]]
    if valKindBin==1:
        lRange[1][0] = -0.1
        lRange[1][1] = 1.2

    lTitle = list()
    if valKindBin==0: 
        histName = 'hBlankForChJetV2AllCent'
        lTitle = [histName ,'#it{p}_{T, ch jet} [GeV/#it{c}]','#it{v}_{2}^{ ch jet}']
    elif valKindBin==2:
        histName = 'hBlankForChJetRAAAllCent'
        lTitle = [histName ,'#it{p}_{T, ch jet} [GeV/#it{c}]','#it{R}_{AA}^{ ch jet}']
    hBlankTargetObsAllCent = genePlotSets.makeBlank1DHist(nBin, lRange, lTitle)
    legALICE = ROOT.TLegend(0.15,0.8,0.25,0.9,'')
    legJet = ROOT.TLegend(0.15,0.65,0.25,0.75,'')
    legTargetObsAllCent = ROOT.TLegend(0.6,0.7,0.8,0.88,'Centrality[%]')

    genePlotSets.addALICELegend(hBlankTargetObsAllCent, legALICE, lCentLabel[centBin], 'Work in progress', 0.03)
    genePlotSets.addJetLegend(hBlankTargetObsAllCent, legJet, 0.2, leadingTrackPtCut, 0.03)
    genePlotSets.addSomeHistsLegend(legTargetObsAllCent, lhTargetObsAllCent, lCentLabel, 0.03)
    label = ''
    if valKindBin==0: label = 'ChJetV2/ChJetV2AllCent'+ '_TrackPtCut'+str(leadingTrackPtCut)
    elif valKindBin==2: label = 'ChJetRAA/ChJetRAAAllCent'+ '_TrackPtCut'+str(leadingTrackPtCut)
    # [kViolet+1, kAzure+1, kOrange+8, kGreen+3,  kPink+9]
    colorList = [880+1, 860+7, 800+7, 416+2, 800-3,  900+9]
    styleList = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20]
    genePlotSets.overwrightSomePlotsWithErr(hBlankTargetObsAllCent, \
        lhTargetObsAllCent, lhTargetObsSysErrAllCent,\
            legALICE,legJet,legTargetObsAllCent, colorList, styleList, \
                label, 'EP X0', outFinalFileDir)


def CompareWithAnotherCentral(lAllCentChJetV2Val):
    # == s == prepare this measurement (5.02 TeV ALICE Pb-Pb)  =================
    ptBinArray5TeV = lAllCentChJetV2Val[4]
    nBins = len(ptBinArray5TeV) - 1

    hChJetV5TeV = ROOT.TH1D("hChJetV5TeV", "hChJetV5TeV", nBins, ptBinArray5TeV)
    for ptBin in range(0, nBins):
        hChJetV5TeV.Fill(ptBinArray5TeV[ptBin]+0.01, lAllCentChJetV2Val[0][ptBin])
        hChJetV5TeV.SetBinError(ptBin+1, lAllCentChJetV2Val[3][ptBin])

    sysErr5TeV = ROOT.TGraphErrors(nBins, lAllCentChJetV2Val[1], \
        lAllCentChJetV2Val[0], lAllCentChJetV2Val[2], lAllCentChJetV2Val[5])
    color = 880+1
    genePlotSets.setHistLooksWSysE(hChJetV5TeV, sysErr5TeV, 20, color, 0)
    # == e == prepare this measurement (5.02 TeV ALICE Pb-Pb)  =================


    # == s == prepare paper figure (2.7 TeV ALICE Pb-Pb)  =================
    # https://arxiv.org/pdf/1509.07334.pdf
    # https://www.hepdata.net/record/ins1394678
    ChJetV23TeV = [0.0388, 0.0413, 0.0528, 0.0648, 0.0738, 0.0799, 0.0843]
    ptBin3TeV = [35., 45., 55., 65., 75., 85., 95.]
    binWidth3TeV = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    statErr3TeV = [0.00819, 0.011, 0.0148, 0.0203, 0.027, 0.0345, 0.0421]
    binArray3TeV = [30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100]
    shapeSysErr3TeV = [0.0322, 0.0247, 0.0181, 0.0141, 0.0193, 0.0287, 0.03]

    ChJetV23TeV = array('d',ChJetV23TeV)
    shapeSysErr3TeV = array('d',shapeSysErr3TeV)
    ptBin3TeV = array('d',ptBin3TeV)
    binWidth3TeV = array('d',binWidth3TeV)
    statErr3TeV = array('d',statErr3TeV)
    binArray3TeV = array('d',binArray3TeV)

    nBins = len(binArray3TeV) - 1
    hChJetV23TeV = ROOT.TH1D("hChJetV23TeV", "hChJetV23TeV", nBins, binArray3TeV)
    for ptBin in range(0, nBins):
        hChJetV23TeV.Fill(binArray3TeV[ptBin]+0.01, ChJetV23TeV[ptBin])
        hChJetV23TeV.SetBinError(ptBin+1, statErr3TeV[ptBin])

    color = 600+2
    sysErr3TeV = ROOT.TGraphErrors(7, ptBin3TeV, ChJetV23TeV, binWidth3TeV, shapeSysErr3TeV)
    genePlotSets.setHistLooksWSysE(hChJetV23TeV, sysErr3TeV, 20, color, 0)
    # == e == prepare paper figure (2.7 TeV ALICE Pb-Pb)  =================


    # == s == prepare paper figure (5.02 TeV ATLAS Pb-Pb)  =================
    # https://www.hepdata.net/record/ins1967021
    # https://arxiv.org/pdf/2111.06606.pdf
    ChJetV2ATLAS5TeV = [0.00528595, -0.00455564, -0.00148288, 0.00398145, 0.00114189, -0.00458343, 0.00171849, -0.0108833, -0.00454052]
    ptBinATLAS5TeV = [75, 84, 94.5, 113, 142, 179, 225.5, 283.5, 357]
    binWidthATLAS5TeV = [4.0, 5.0, 5.5, 13.0, 16.0, 21.0, 25.5, 32., 41]
    statErrATLAS5TeV = [0.00294283, 0.00322651, 0.00298647, 0.00280327, 0.00264858, 0.00433627, 0.00739941, 0.0130912, 0.0236564]
    binArrayATLAS5TeV = [71., 79., 89., 100., 126., 158., 200., 251., 316., 398]
    shapeSysErrATLAS5TeV = [0.00532905, 0.00500799, 0.0033307, 0.00184235, 0.00328851, 0.00246496, 0.00134662, 0.00191156, 0.0026559]

    ChJetV2ATLAS5TeV = array('d',ChJetV2ATLAS5TeV)
    shapeSysErrATLAS5TeV = array('d',shapeSysErrATLAS5TeV)
    ptBinATLAS5TeV = array('d',ptBinATLAS5TeV)
    binWidthATLAS5TeV = array('d',binWidthATLAS5TeV)
    statErrATLAS5TeV = array('d',statErrATLAS5TeV)
    binArrayATLAS5TeV = array('d',binArrayATLAS5TeV)

    nBins = len(binArrayATLAS5TeV) - 1
    hChJetV2ATLAS5TeV = ROOT.TH1D("hChJetV2ATLAS5TeV", "hChJetV2ATLAS5TeV", nBins, binArrayATLAS5TeV)
    for ptBin in range(0, nBins):
        hChJetV2ATLAS5TeV.Fill(binArrayATLAS5TeV[ptBin]+0.01, ChJetV2ATLAS5TeV[ptBin])
        hChJetV2ATLAS5TeV.SetBinError(ptBin+1, statErrATLAS5TeV[ptBin])

    sysErrATLAS5TeV = ROOT.TGraphErrors(7, ptBinATLAS5TeV, ChJetV2ATLAS5TeV, binWidthATLAS5TeV, shapeSysErrATLAS5TeV)
    color = 632+2
    genePlotSets.setHistLooksWSysE(hChJetV2ATLAS5TeV, sysErrATLAS5TeV, 20, color, 0)
    # == e == prepare paper figure (5.02 TeV ATLAS Pb-Pb)  =================


    lhChJetV2Compare = list()
    lhChJetV2Compare.append(hChJetV23TeV)
    lhChJetV2Compare.append(hChJetV2ATLAS5TeV)
    lhChJetV2Compare.append(hChJetV5TeV)
    lhChJetV2SysErrCompare = list()
    lhChJetV2SysErrCompare.append(sysErr3TeV)
    lhChJetV2SysErrCompare.append(sysErrATLAS5TeV)
    lhChJetV2SysErrCompare.append(sysErr5TeV)

    # nBin = 120
    # lRange = [[30, 150],[-0.01, 0.2]]
    nBin = 350
    lRange = [[30, 380],[-0.01, 0.2]]
    lTitle = ['hBlanckForChJetV2Compare','#it{p}_{T, ch jet} [GeV/#it{c}]','#it{v}_{2}^{ ch jet}']
    hBlankChJetV2Compare = genePlotSets.makeBlank1DHist(nBin, lRange, lTitle)
    legALICE = ROOT.TLegend(0.15,0.8,0.25,0.9,'')
    legJet = ROOT.TLegend(0.15,0.65,0.25,0.75,'')
    legChJetV2Compare = ROOT.TLegend(0.6,0.7,0.8,0.88,'')

    genePlotSets.addALICELegend(hBlankChJetV2Compare, legALICE, lCentLabel[0], 'Work in progress', 0.03)
    genePlotSets.addJetLegend(hBlankChJetV2Compare, legJet, 0.2, leadingTrackPtCut, 0.03)
    lHistLabel = ['#sqrt{#it{s}_{NN}} = 2.76 TeV (ALICE Run1)', '#sqrt{#it{s}_{NN}} = 5.02 TeV (ATLAS)', '#sqrt{#it{s}_{NN}} = 5.02 TeV (This Measurement)']
    genePlotSets.addSomeHistsLegend(legChJetV2Compare, lhChJetV2Compare, lHistLabel, 0.03)
    label = 'CentralChJetV2Compare'+ '_TrackPtCut'+str(leadingTrackPtCut)

    colorList = [600+2, 632+2, 880+1]
    styleList = [20, 20, 20]

    genePlotSets.overwrightSomePlotsWithErr(hBlankChJetV2Compare, \
        lhChJetV2Compare, lhChJetV2SysErrCompare,\
            legALICE,legJet,legChJetV2Compare, colorList, styleList, \
                label, 'EP X0', outFinalFileDir)


    # ChJetV2ATLAS5TeV = [0.0253191, 0.0223062, 0.01972, 0.0184535, 0.0180286, 0.00906937, 0.0125028, 0.00406064, 0.0144245]
    # ptBinATLAS5TeV = [75, 84, 94.5, 113, 142, 179, 225.5, 283.5, 357]
    # binWidthATLAS5TeV = [4.0, 5.0, 5.5, 13.0, 16.0, 21.0, 25.5, 32., 41]
    # statErrATLAS5TeV = [0.00240366, 0.0026685, 0.00240415, 0.00230443, 0.0021526, 0.00350096, 0.00595128, 0.0104246, 0.0191402]
    # binArrayATLAS5TeV = [71., 79., 89., 100., 126., 158., 200., 251., 316., 398]
    # shapeSysErrATLAS5TeV = [0.00513033, 0.00485812, 0.00418428, 0.00378757, 0.0038047, 0.00400346, 0.00479492, 0.00528187, 0.00924739]


def CompareWithAnotherSemiCentral(lAllCentChJetV2Val):
    ptBinArray5TeV = lAllCentChJetV2Val[4]
    nBins = len(ptBinArray5TeV) - 1

    hChJetV5TeV = ROOT.TH1D("hChJetV5TeV", "hChJetV5TeV", nBins, ptBinArray5TeV)
    for ptBin in range(0, nBins):
        hChJetV5TeV.Fill(ptBinArray5TeV[ptBin]+0.01, lAllCentChJetV2Val[0][ptBin])
        hChJetV5TeV.SetBinError(ptBin+1, lAllCentChJetV2Val[3][ptBin])

    sysErr5TeV = ROOT.TGraphErrors(nBins, lAllCentChJetV2Val[1], \
        lAllCentChJetV2Val[0], lAllCentChJetV2Val[2], lAllCentChJetV2Val[5])
    color = 416+2
    genePlotSets.setHistLooksWSysE(hChJetV5TeV, sysErr5TeV, 20, color, 0)

    # == s == prepare paper figure (2.7 TeV ALICE Pb-Pb)  =================
    # https://arxiv.org/pdf/1509.07334.pdf
    # https://www.hepdata.net/record/ins1394678
    ChJetV23TeV = [0.0811, 0.0914, 0.0756, 0.0625, 0.0575, 0.0588, 0.0611]
    ptBin3TeV = [25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0]
    binWidth3TeV = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    statErr3TeV = [0.00753, 0.0114, 0.0165, 0.0213, 0.0286, 0.0383, 0.0489]
    binArray3TeV = [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]
    shapeSysErr3TeV = [0.0392, 0.0291, 0.0179, 0.0188, 0.0111, 0.0154, 0.0148]

    ChJetV23TeV = array('d',ChJetV23TeV)
    shapeSysErr3TeV = array('d',shapeSysErr3TeV)
    ptBin3TeV = array('d',ptBin3TeV)
    binWidth3TeV = array('d',binWidth3TeV)
    statErr3TeV = array('d',statErr3TeV)
    binArray3TeV = array('d',binArray3TeV)

    nBins = len(binArray3TeV) - 1
    hChJetV23TeV = ROOT.TH1D("hChJetV23TeV", "hChJetV23TeV", nBins, binArray3TeV)
    for ptBin in range(0, nBins):
        hChJetV23TeV.Fill(binArray3TeV[ptBin]+0.01, ChJetV23TeV[ptBin])
        hChJetV23TeV.SetBinError(ptBin+1, statErr3TeV[ptBin])

    color = 600+2
    sysErr3TeV = ROOT.TGraphErrors(7, ptBin3TeV, ChJetV23TeV, binWidth3TeV, shapeSysErr3TeV)
    genePlotSets.setHistLooksWSysE(hChJetV23TeV, sysErr3TeV, 20, color, 0)
    # == e == prepare paper figure (2.7 TeV ALICE Pb-Pb)  =================

    # == s == prepare paper figure (5.02 TeV ATLAS Pb-Pb)  =================
    # https://www.hepdata.net/record/ins1967021
    # https://arxiv.org/pdf/2111.06606.pdf
    ChJetV2ATLAS5TeV = [0.0467919, 0.0384123, 0.0357318, 0.0323667, 0.0275154, 0.0254492, 0.016756, 0.00858585, 0.0176957]
    ptBinATLAS5TeV = [75., 84., 94.5, 113., 142., 179., 225.5, 283.5, 357]
    binWidthATLAS5TeV = [4.0, 5.0, 5.5, 13.0, 16.0, 21.0, 25.5, 32., 41]
    statErrATLAS5TeV = [0.00168948, 0.00182129, 0.00162882, 0.00157211, 0.00149234, 0.00244447, 0.00423172, 0.00752307, 0.0135238]
    binArrayATLAS5TeV = [71., 79., 89., 100., 126., 158., 200., 251., 316., 398]
    shapeSysErrATLAS5TeV = [0.00536816, 0.00459791, 0.00245665, 0.00194695, 0.00214304, 0.00216564, 0.00192541, 0.00173636, 0.00174023]

    ChJetV2ATLAS5TeV = array('d',ChJetV2ATLAS5TeV)
    shapeSysErrATLAS5TeV = array('d',shapeSysErrATLAS5TeV)
    ptBinATLAS5TeV = array('d',ptBinATLAS5TeV)
    binWidthATLAS5TeV = array('d',binWidthATLAS5TeV)
    statErrATLAS5TeV = array('d',statErrATLAS5TeV)
    binArrayATLAS5TeV = array('d',binArrayATLAS5TeV)

    nBins = len(binArrayATLAS5TeV) - 1
    hChJetV2ATLAS5TeV = ROOT.TH1D("hChJetV2ATLAS5TeV", "hChJetV2ATLAS5TeV", nBins, binArrayATLAS5TeV)
    for ptBin in range(0, nBins):
        hChJetV2ATLAS5TeV.Fill(binArrayATLAS5TeV[ptBin]+0.01, ChJetV2ATLAS5TeV[ptBin])
        hChJetV2ATLAS5TeV.SetBinError(ptBin+1, statErrATLAS5TeV[ptBin])

    color = 632+2
    sysErrATLAS5TeV = ROOT.TGraphErrors(7, ptBinATLAS5TeV, ChJetV2ATLAS5TeV, binWidthATLAS5TeV, shapeSysErrATLAS5TeV)
    genePlotSets.setHistLooksWSysE(hChJetV2ATLAS5TeV, sysErrATLAS5TeV, 20, color, 0)
    # == e == prepare paper figure (5.02 TeV ATLAS Pb-Pb)  =================
    
    #   //get Caitie's results
    color = 900+10
    file1 = ROOT.TFile('v2check_R02_Apr20.root')
    R2v2Caitie = file1.Get('myv2')
    R2v2Caitie.SetLineColor(color)
    R2v2Caitie.SetLineWidth(2)
    R2v2Caitie.SetMarkerColor(color)
    R2v2Caitie.SetMarkerSize(0.8)
    R2v2Caitie.SetMarkerStyle(20)
    tempErrList = [0,0,0,0,0,0]
    tempErrList = [0,0,0,0,0,0]
    tempErrWidth = [5,10,10,20,20,20]
    tempErrList  = array('d', tempErrList)
    sysErrCaitie = ROOT.TGraphErrors(6, tempErrList, \
        tempErrList, tempErrList, tempErrList)

    lhChJetV2Compare = list()
    lhChJetV2Compare.append(hChJetV23TeV)
    lhChJetV2Compare.append(hChJetV2ATLAS5TeV)
    lhChJetV2Compare.append(R2v2Caitie)
    lhChJetV2Compare.append(hChJetV5TeV)
    lhChJetV2SysErrCompare = list()
    lhChJetV2SysErrCompare.append(sysErr3TeV)
    lhChJetV2SysErrCompare.append(sysErrATLAS5TeV)
    lhChJetV2SysErrCompare.append(sysErrCaitie)
    lhChJetV2SysErrCompare.append(sysErr5TeV)

    # nBin = 120
    # lRange = [[30, 150],[-0.01, 0.2]]
    nBin = 350
    lRange = [[30, 380],[-0.01, 0.2]]
    lTitle = ['hBlanckForChJetV2Compare','#it{p}_{T, ch jet} [GeV/#it{c}]','#it{v}_{2}^{ ch jet}']
    hBlankChJetV2Compare = genePlotSets.makeBlank1DHist(nBin, lRange, lTitle)
    legALICE = ROOT.TLegend(0.15,0.8,0.25,0.9,'')
    legJet = ROOT.TLegend(0.15,0.65,0.25,0.75,'')
    legChJetV2Compare = ROOT.TLegend(0.6,0.7,0.8,0.88,'')

    genePlotSets.addALICELegend(hBlankChJetV2Compare, legALICE, lCentLabel[3], 'Work in progress', 0.03)
    genePlotSets.addJetLegend(hBlankChJetV2Compare, legJet, 0.2, leadingTrackPtCut, 0.03)
    lHistLabel = ['#sqrt{#it{s}_{NN}} = 2.76 TeV (ALICE Run1)', '#sqrt{#it{s}_{NN}} = 5.02 TeV (ATLAS Centrality 20-40%)', 'Caitie Estimation ', '#sqrt{#it{s}_{NN}} = 5.02 TeV (This Measurement)']
    genePlotSets.addSomeHistsLegend(legChJetV2Compare, lhChJetV2Compare, lHistLabel, 0.03)
    label = 'ChJetV2SemiCentralCompare'+ '_TrackPtCut'+str(leadingTrackPtCut)

    colorList = [600+2, 632+2, 900+10, 416+2]
    styleList = [20, 20, 20, 20]

    genePlotSets.overwrightSomePlotsWithErr(hBlankChJetV2Compare, \
        lhChJetV2Compare, lhChJetV2SysErrCompare,\
            legALICE,legJet,legChJetV2Compare, colorList, styleList, \
                label, 'EP X0', outFinalFileDir)

def V2StatErrorPropergate(hOEPUFJet, hIEPUFJet):
    lv2Error = list()
    for iBin in range(1, hChJetV2Nume.GetNbinsX() + 1):
        IEPUFJetContent = hIEPUFJet.GetBinContent(iBin)
        OEPUFJetContent = hOEPUFJet.GetBinContent(iBin)
        IEPUFJetErr = hIEPUFJet.GetBinError(iBin)
        OEPUFJetErr = hOEPUFJet.GetBinError(iBin)
        
        DenoErr = (IEPUFJetContent+OEPUFJetContent)*(IEPUFJetContent+OEPUFJetContent)\
            *(IEPUFJetContent+OEPUFJetContent)*(IEPUFJetContent+OEPUFJetContent)
        ONumeErr = 4*(OEPUFJetContent*OEPUFJetContent)*(IEPUFJetErr*IEPUFJetErr)
        INumeErr = 4*(IEPUFJetContent*IEPUFJetContent)*(OEPUFJetErr*OEPUFJetErr) 
        CoEffError =  TMath.Pi() / (4 * lPsi2Reso[centBin])
        multErr = np.sqrt(ONumeErr/DenoErr + INumeErr/DenoErr)

        lv2Error.append(multErr)

def V2SysErrorPropergate(hV2, hOEPUFJet, hIEPUFJet, lRelaSysUncErrOEP, lRelaSysUncErrIEP, centBin, pTBinArray, diffNum, diffName):
    histname = "hDiff_" + diffName + '_Cent'+ str(centBin)
    title = histname + ";#it{p}_{T,ch jet} (GeV/#it{c}); Ratio Systematic Uncertainty [%]"
    
    nBins = len(pTBinArray)-1
    hSysUncPer = ROOT.TH1D(histname, title, nBins, pTBinArray)
    
    lSysUncErr  = list()
    errBin = 0
    for iBin in range(1, hV2.GetNbinsX() + 1):

        fillVal = hV2.GetXaxis().GetBinCenter(iBin)
        if (fillVal > pTBinArray[0]) and (fillVal < pTBinArray[-1]):
            IEPUFJetContent = hIEPUFJet.GetBinContent(iBin)
            OEPUFJetContent = hOEPUFJet.GetBinContent(iBin)
            
            IEPUFJetErr = lRelaSysUncErrIEP[errBin]
            OEPUFJetErr = lRelaSysUncErrOEP[errBin]
            DenoErr = (IEPUFJetContent+OEPUFJetContent)*(IEPUFJetContent+OEPUFJetContent)\
                *(IEPUFJetContent+OEPUFJetContent)*(IEPUFJetContent+OEPUFJetContent)
            ONumeErr = 4*(OEPUFJetContent*OEPUFJetContent)*(IEPUFJetErr*IEPUFJetErr)
            INumeErr = 4*(IEPUFJetContent*IEPUFJetContent)*(OEPUFJetErr*OEPUFJetErr) 
            # CoEffError =  TMath.Pi() / (4 * lPsi2Reso[centBin])
            multErr = np.sqrt(ONumeErr/DenoErr + INumeErr/DenoErr)

            ratio = multErr/hV2.GetBinContent(iBin)
            hSysUncPer.Fill(fillVal, ratio)
            # hSysUncPer.Fill(fillVal, multErr)
            # SysUncErr = np.sqrt(deltaVal*deltaVal)
            SysUncErr = multErr*hV2.GetBinContent(iBin)
            # SysUncErr = multErr
            lSysUncErr.append(SysUncErr)

            errBin+=1

    for iBin in range(1, hSysUncPer.GetNbinsX() + 1): hSysUncPer.SetBinError(iBin, 0.)
    hSysUncPer.GetXaxis().SetTitle("#it{p}_{T, ch jet} [GeV/#it{c}]")
    hSysUncPer.GetYaxis().SetTitle("Ratio Systematic Uncertainty")
    # hSysUncPer.SetMarkerStyle(21)
    hSysUncPer.SetMarkerColor(lFillColor[diffNum])
    hSysUncPer.SetLineColor(lFillColor[diffNum])
    hSysUncPer.SetLineWidth(4)
    
    lSysUncErr = array('d',lSysUncErr)

    return hSysUncPer, lSysUncErr



def GetSysUncErrObjV2(lAllValWithErr, dTargedObs, \
    hJetYeildOEP, hJetYeildIEP, dlRelaSysUncErrOEP, dlRelaSysUncErrIEP,\
        valKindBin, pTBinArray, centBin, lOMainTree, oFile):

    # == s == Calculate systematic uncertainty  ===========================
    valDire = ['ChJetV2', 'ChJetRAA', 'ChJetInc', 'ChJetOEP', 'ChJetIEP']

    dhAbsSysUnc = {}
    dlAbsSysUncErr  = {}
    dhRelaSysUnc = {}
    dlRelaSysUncErr = {}
    
    dhRelaSysUnc['prior'], dlRelaSysUncErr['prior'] = V2SysErrorPropergate(\
        dTargedObs['Nominal'], hJetYeildOEP, hJetYeildIEP, \
        dlRelaSysUncErrOEP['prior'], dlRelaSysUncErrIEP['prior'], centBin, pTBinArray, 1, 'prior')

    dhRelaSysUnc['ptBinRange'], dlRelaSysUncErr['ptBinRange'] = V2SysErrorPropergate(dTargedObs['Nominal'], hJetYeildOEP, hJetYeildIEP, \
        dlRelaSysUncErrOEP['ptBinRange'], dlRelaSysUncErrIEP['ptBinRange'], centBin, pTBinArray, 2, 'ptBinRange')

    dhRelaSysUnc['iteration'], dlRelaSysUncErr['iteration'] = V2SysErrorPropergate(dTargedObs['Nominal'], hJetYeildOEP, hJetYeildIEP, \
        dlRelaSysUncErrOEP['iteration'], dlRelaSysUncErrIEP['iteration'], centBin, pTBinArray, 3, 'iteration')

    dhRelaSysUnc['trackEff'], dlRelaSysUncErr['trackEff'] = V2SysErrorPropergate(dTargedObs['Nominal'], hJetYeildOEP, hJetYeildIEP, \
        dlRelaSysUncErrOEP['trackEff'], dlRelaSysUncErrIEP['trackEff'], centBin, pTBinArray, 4,'trackEff')

    dhRelaSysUnc['bkgWay'], dlRelaSysUncErr['bkgWay'] = V2SysErrorPropergate(dTargedObs['Nominal'], hJetYeildOEP, hJetYeildIEP, \
        dlRelaSysUncErrOEP['bkgWay'], dlRelaSysUncErrIEP['bkgWay'], centBin, pTBinArray, 5, 'bkgWay')

    dhRelaSysUnc['V0Diff'], dlRelaSysUncErr['V0Diff'] = V2SysErrorPropergate(dTargedObs['Nominal'], hJetYeildOEP, hJetYeildIEP, \
        dlRelaSysUncErrOEP['V0Diff'], dlRelaSysUncErrIEP['V0Diff'], centBin, pTBinArray, 6,'V0Diff')
    
    SumUpAllRelaSysErr(lAllValWithErr, dTargedObs['Nominal'], dhRelaSysUnc, \
        dlRelaSysUncErr ,pTBinArray,centBin)

    print('CheeeeeeeeeeeeeeeeecKumaV2')
    print(dlRelaSysUncErr)
    
    hNameSysErr = 'h'+valDire[valKindBin]+'SysErr' + '_Cent' + str(centBin)
    hTargedObsSysUncErr = SetupTGraphErrors(dTargedObs['Nominal'], dlRelaSysUncErr ['total'],\
        hNameSysErr, pTBinArray, centBin)

    # == s == Calculate systematic uncertainty  ===========================

    # == s == plot systematic error #######################################
    print('CheeeeeeeeeeeeeeeeecKuma111')
    print(dhRelaSysUnc)
    lhSysUncError = list()
    for key in dhRelaSysUnc.values(): lhSysUncError.append(key)
    lhSysUncError.reverse()
    lNameSysUncError = ['Reweighting Prior', 'Truncation', 'Iterations', \
        'Tracking Efficiency', 'Background Estiamtion Way', 'V0 Detectors','Total']
    lNameSysUncError.reverse()
    nBin = (pTBinArray[-1] - pTBinArray[0])/10
    lRange = [[pTBinArray[0], pTBinArray[-1]],[0.0, 2.]]
    lTitle = ['hBlankForSysError','#it{p}_{T, ch jet} [GeV/#it{c}]','systematic error ratio']
    hBlankEachCentSysErr = genePlotSets.makeBlank1DHist(nBin, lRange, lTitle)
    legALICE = ROOT.TLegend(0.15,0.8,0.25,0.9,'')
    legJet = ROOT.TLegend(0.15,0.65,0.25,0.75,'')
    legSysUncError = ROOT.TLegend(0.6,0.7,0.8,0.88,'systematic error ratio')
    genePlotSets.addALICELegend(hBlankEachCentSysErr, legALICE, lCentLabel[centBin], \
        'Work in progress', 0.03)
    genePlotSets.addJetLegend(hBlankEachCentSysErr, legJet, 0.2, leadingTrackPtCut, 0.03)        
    genePlotSets.addSomeHistsLegend(legSysUncError, lhSysUncError, lNameSysUncError, 0.03)
    
    label = valDire[valKindBin]+'/SysUncError'+'_TrackPtCut'+str(leadingTrackPtCut)+'_Cent'+str(centBin)
    plotStyleList = list()
    # [kViolet+1, kBlue, kAzure+1, kGreen+3, kOrange+8, kPink+9]
    colorList = [880+1, 600, 860+1, 416+3, 800-3, 800+8, 900+9]
    for l in range(0, len(lhSysUncError)): plotStyleList.append(0)
    genePlotSets.overwrightSomePlots(hBlankEachCentSysErr, lhSysUncError, \
        0, legALICE,legJet,legSysUncError, lFillColor, plotStyleList, label, '', outFinalFileDir)
    # == e == plot systematic error #######################################

    oFile.cd()
    for key in dTargedObs.values(): lOMainTree[valKindBin][centBin][0].Add(key)
    # for key in dhAbsSysUnc.values(): lOMainTree[valKindBin][centBin][1].Add(key)
    for key in dhRelaSysUnc.values(): lOMainTree[valKindBin][centBin][1].Add(key)
    lOMainTree[valKindBin][centBin][2].Add(hTargedObsSysUncErr)
    
    return dhRelaSysUnc, dlRelaSysUncErr 



################################################################################
def OFileStracture(lMainTree, numOfCentBin):
    for valKindBin in range(0, 5):
        lValTree = ROOT.TList()
        lValTree.SetName(lValKindName[valKindBin])
        for centBin in range(0, numOfCentBin):
            lEachCentHists = ROOT.TList()
            lEachCentHists.SetName('lCent{0}'.format(centBin))

            lTargetVal = ROOT.TList()
            lTargetVal.SetName('lTargetVal')
            lEachCentHists.Add(lTargetVal)

            lSysUncErr = ROOT.TList()
            lSysUncErr.SetName('lSysUncErr')
            lEachCentHists.Add(lSysUncErr)

            lSysUncErrGrough = ROOT.TList()
            lSysUncErrGrough.SetName('lSysUncErrGrough')
            lEachCentHists.Add(lSysUncErrGrough)

            lValTree.Add(lEachCentHists)

        lMainTree.Add(lValTree)
################################################################################

if __name__ == "__main__":
    PlotSysUncertainty()


