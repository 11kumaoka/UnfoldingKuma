#! /usr/bin/env python
import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F
import argparse
import ctypes
import os
import gc
from array import array
import numpy as np

import histSetting
import outputFileStracture
import calcPtHardScaleFactor
import ptHardBinScaleMerge
import makeRMProcess
import plotPerformanceHists
import countEachPtHardNEvents

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)


inputBasePassName = '/home/alidock/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/embedding/'
inputEmbFilePassName = inputBasePassName + 'AnalysisResultsEmbFullCh18r21May06'

inputDataBaseDir = '/home/alidock/OutPlot/LHC18qr/'
inputDataFileName = inputDataBaseDir + 'AnalysisResultsVariousCut.root'
FullCh = 'Ch' # 'Ch'/'Full'

ptHardBinList = [5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150,\
  169, 190, 212, 235, -1]

lrefPtHardScalFactor = [13.256237105043024, 3.994951508183162, 1.901382604892938, 0.6997084477251105,\
  0.23859307106048133, 0.08838464195186756, 0.026609353549318118, 0.009002896389014694,\
    0.0036894610775540953, 0.0012366270093495048, 0.0004840940273109418, 0.0001721894632700725, \
      8.419234025525398e-05, 3.914934018299818e-05, 1.9044557347701444e-05, 9.653419162066905e-06, \
        5.223266284853783e-06, 2.7229473408774812e-06, 1.4609950982940283e-06, 1.8795100140338064e-06]

lrefPtHardNumOfEvent = [5883650.0, 6097914.0, 6526103.0, 6487913.0, 6542589.0, \
  6618306.0, 6659974.0, 6606850.0, 6600430.0, 6473709.0, 6788644.0, 6612948.0, \
    6467705.0, 6215084.0, 6635787.0, 6707788.0, 6710764.0, 5544633.0, 6355388.0, 6619566.]

initPtHardBin = 1 # ordinaly 0

# resoPara = [2, 3, 4, 5]
resoPara = [2]
centBinRangeList = [[1,3], [3,7], [7,11], [11,15], [15, 21]]
centRangeList = [[0,10], [10,30], [30,50], [50,70], [70, 100]]

#00  Main function   00000000000000000000000000000000000000000000000000000000000
# Get a list of all the output list names, in order to scale all of them
def createPlotEmb():
  
  outRootFile = ROOT.TFile('1hJetPerfomancePlots' +FullCh+ '.root', 'RECREATE')
  lMainTree = ROOT.TList()

  for rBin in range(0, 1):
    lEachRJetTree = ROOT.TList()
    lEachRJetTree.SetName('jetTree_R{}'.format(resoPara[rBin]))
    lMainTree.Add(lEachRJetTree)

    outputFileStracture.outputFileStracture(lMainTree[rBin], resoPara[rBin])
  
  taskNameEmbHelper = 'AliAnalysisTaskEmcalEmbeddingHelper_histos'
  hNEventsAcc, hNEventsTot = \
    countEachPtHardNEvents.extractNEvetns(inputEmbFilePassName, taskNameEmbHelper, ptHardBinList)
  lMainTree.Add(hNEventsAcc)
  lMainTree.Add(hNEventsTot)
  outRootFile.cd()

  taskBaseName = 'AliAnalysisTaskEmcalJetSample'\
      + '_tracks' + '_caloClusters' + '_emcalCells' + '_leadingTrackCut0_histos'
  # taskBaseName = 'AliAnalysisTaskEmcalJetPerformance'\
      # + '_tracks' + '_caloClustersCombined' + '_02' + '_5GeV_histos'
  fCutStats, Centrality_raw, fNormalisationHist \
    = countEachPtHardNEvents.extractDataNEvetns(inputDataFileName, taskBaseName)
    # = countEachPtHardNEvents.extractDataNEvetnsV2(inputDataFileName, taskBaseName)
    
  lMainTree.Add(fCutStats)
  lMainTree.Add(Centrality_raw)
  lMainTree.Add(fNormalisationHist)
  outRootFile.cd()

  # lScaleFactor, lPtHardEventNum, hXSecPerEvent, hNTrialsPerEvent, hScaleFactor, hNEventsAcc, hNEventsTot\
    # = ptHardBinScaleMerge.ptHardBinScaleMerge() # scale factor ver 1
  # lScaleFactor = lrefPtHardScalFactor# scale factor ver 2

  for rBin in range(0, 1):
    lhCheckEachJetPtBinPart = [[],[],[],[],[]]
    lhCheckEachJetPtBinHyb = [[],[],[],[],[]]
    for ptHardBin in range(1, 19): #0,20
      # if((ptHardBin != 0) and (ptHardBin != 1)): continue
      # scaleFactor = lrefPtHardScalFactor[ptHardBin-1]
      scaleFactor = lrefPtHardScalFactor[ptHardBin-1] / lrefPtHardNumOfEvent[ptHardBin-1]
      histExtraction(inputEmbFilePassName, outRootFile, lMainTree[rBin], \
        lhCheckEachJetPtBinPart, lhCheckEachJetPtBinHyb, ptHardBin, rBin, scaleFactor)
  
    checkEachPtHardDist(lMainTree[rBin], lhCheckEachJetPtBinPart, lhCheckEachJetPtBinHyb)

  lMainTree.Write('mainTree', 1)
  outRootFile.Close()
#0000000000000000000000000000000000000000000000000000000000000000000000000000000


#11  histExtraction   1111111111111111111111111111111111111111111111111111111111
def histExtraction(inputEmbFilePassName, outRootFile, lMainTree, \
  lhCheckEachJetPtBinPart, lhCheckEachJetPtBinHyb, ptHardBin, rBin, ptHardScaleFactor):
  # f = ROOT.TFile("AnalysisResults.root", "READ")

  f = ROOT.TFile(inputEmbFilePassName + "_{}.root".format(ptHardBin), "READ")
  # f = ROOT.TFile("TrainAnalysisResults.root", "READ")
  performanceListKeys = f.GetListOfKeys()
  performanceListNames = []

  for key in performanceListKeys:
    name = key.GetName()
    # print(name)
  
  # resoPara = [2, 3, 4, 5]
  eventList = 'AliAnalysisTaskEmcalJetPerformance_'
  eventList += 'tracks_caloClustersCombined_JetPerformance_'
  eventList +=  FullCh + 'JetR0' + str(resoPara[rBin]) + '_histos'

  performanceList = f.Get(eventList)
  #print(performanceList)
  if not performanceList:
    print("ERROR no list found")

  eventFilterTask = "AliAnalysisTaskEmcalEmbeddingHelper_histos"
  eventFilterTaskList = f.Get(eventFilterTask)

  f.Close()

  # performanceList.Print()
  outRootFile.cd()
  
  # ptHardScaleFactor = calcPtHardScaleFactor.calcPtHardScalFactor(eventFilterTaskList, ptHardBin)
  # ptHardScaleFactor = 1

  ptHardDistTemp = eventFilterTaskList[3]
  ptHardDistTemp.Sumw2()
  ptHardDistTemp.Scale(ptHardScaleFactor)

  if ptHardBin == initPtHardBin: 
    lMainTree.Add(ptHardDistTemp)
  else :
    ptHardBinDist = lMainTree.FindObject('fHistPtHard')
    ptHardBinDist.Add(ptHardDistTemp)

  #Extract generation level pT distribution
  # performanceList[4].Print() #truthLevelJets  #old:? new: 4
  genJetBranch = performanceList[23] # new: 23
  genJetHistos = genJetBranch[0]
  hGenJetAreaPtCent = genJetHistos[3] # use generation level pT distribution
  
  #Extract detector level pT distribution
  # performanceList[5].Print() #detLevelJets  #old:? new: 5
  detJetBranch = performanceList[22] # old: 5 new: 22
  detJetHistos = detJetBranch[0]
  hDetJetAreaPtCent = detJetHistos[3] # use detector level pT distribution

  #Extract a branch including Matching histgrams 
  # performanceList[11].Print() #MatchedJeistograms  #old:11 new: 30
  if FullCh is 'Full':tempHistNum = 30
  elif FullCh is 'Ch':tempHistNum = 28
  else: tempHistNum = 9
  matchedJetBranch = performanceList[tempHistNum]
  # matchedJetBranch = performanceList.FindObject('MatchedJetHistograms')
  # matchedJetBranch.Print()
  fHistJetMatchingQA = matchedJetBranch[5] #use mathing hists list
  # hResponseMatrixEMCalDiff = matchedJetBranch[0]

  if ptHardBin is initPtHardBin: fPtRangeBinDev = open('ptRangeBinDev.txt', 'w')
  for centBin in range(0, 5):
    if ptHardBin is initPtHardBin: dataRawJetExtranct(lMainTree[centBin], FullCh, centBin, rBin, fPtRangeBinDev)
    # if ptHardBin is initPtHardBin: dataRawJetExtranctV2(lMainTree[centBin], centBin, rBin)
    outRootFile.cd()
    devHistForEachCent(lMainTree[centBin], hGenJetAreaPtCent, hDetJetAreaPtCent, \
      matchedJetBranch, lhCheckEachJetPtBinPart, lhCheckEachJetPtBinHyb, centBin, ptHardBin, ptHardScaleFactor)
    
  if ptHardBin is initPtHardBin: fPtRangeBinDev.close()
#1111111111111111111111111111111111111111111111111111111111111111111111111111111


#11  histExtraction   1111111111111111111111111111111111111111111111111111111111
def dataRawJetExtranct(lMainTree, FullCh, centBin, rBin, fPtRangeBinDev):
  TAADiv10 = [23.26, 14.4, 8.767,  5.086, 2.747, 1.352,  0.5992,  0.2385, 0.08383, 0.02527]
  
  f = ROOT.TFile(inputDataFileName, "READ")

  # taskName = 'AliAnalysisTaskEmcalJetSample' + '_tracks' + '_caloClusters' + '_emcalCells' + '_leadingTrackCut0_histos'
  taskName = 'AliAnalysisTaskEmcalJetSample' + '_tracks' + '_caloClusters' + '_emcalCells' + '_leadingPtCut5_histos'

  task = f.Get(taskName)
  if not task : print("ERROR no task found")

  lkeys = f.GetListOfKeys()
  lkeysName = []
  # print(" the following lists:")
  for key in lkeys:
    name = key.GetName()
    # print(name)
  f.Close()

  if FullCh is 'Full': jetTreeName = 'Jet_AKTFullR0{}0'.format(resoPara[rBin])\
    + '_tracks_pT0150' + '_caloClusters_E0300_pt_scheme'
  elif FullCh is 'Ch': jetTreeName = 'Jet_AKTChargedR0{}0'.format(resoPara[rBin])\
    + '_tracks_pT0150' + '_pt_scheme'
    
  print(jetTreeName)
  jetList = task.FindObject(jetTreeName)
  if not jetList: print("ERROR no list found")

  lCentMergeBin = [[0,1],[1,3],[3,5],[5,7],[7,9]]

  initCentBin = lCentMergeBin[centBin][0]
  endCentBin = lCentMergeBin[centBin][1]
  
  hEachCentRawJetCorrPt_rebin, lPtBins = ptRangeDef(jetList, centBin, initCentBin, endCentBin)
  lMainTree.Add(hEachCentRawJetCorrPt_rebin)

  fPtRangeBinDev.write('cent bin {} : ['.format(centBin))
  lPtBins_str = map(str, lPtBins)
  lPtBins_write = ','.join(lPtBins_str) + ']\n'
  fPtRangeBinDev.write(lPtBins_write)
  
  for cent10DivBin in range(initCentBin, endCentBin):
    hTempMeasJetPt = jetList.FindObject('histJetPt_{}'.format(initCentBin))
    hTempMeasJetCorrPt = jetList.FindObject('histJetCorrPt_{}'.format(initCentBin))
    hEachCentJetPt = hTempMeasJetPt.Clone()
    hEachCentJetCorrPt = hTempMeasJetCorrPt.Clone()

    hEachCentJetPt.Scale(1/TAADiv10[cent10DivBin])
    hEachCentJetCorrPt.Scale(1/TAADiv10[cent10DivBin])

    if cent10DivBin == initCentBin :
      hEachMergeCentJetPt = hEachCentJetPt
      hEachMergeCentJetPt.SetName('hMeasJetPt_Cent{}'.format(centBin))
      hEachMergeCentJetCorrPt = hEachCentJetCorrPt
      hEachMergeCentJetCorrPt.SetName('hMeasJetCorrPt_Cent{}'.format(centBin))
    else:
      hEachMergeCentJetPt.Add(hEachCentJetPt)
      hEachMergeCentJetCorrPt.Add(hEachCentJetCorrPt)

    print('centBin:{0}, initCentBin:{1}, endCentBin:{2}, cent10DivBin:{3}'\
      .format(centBin, initCentBin, endCentBin, cent10DivBin))

  lMainTree.Add(hEachMergeCentJetPt)
  lMainTree.Add(hEachMergeCentJetCorrPt)

#1111111111111111111111111111111111111111111111111111111111111111111111111111111


#22  devHistForEachCent   222222222222222222222222222222222222222222222222222222
def devHistForEachCent(lMainTree, baseGenJetHists, baseDetJetHists, \
  baseMatchedJetHists, lhCheckEachJetPtBinPart, lhCheckEachJetPtBinHyb, \
    centBin, ptHardBin, ptHardScaleFactor):

  histNum = 3

  ###  Add generation level jet pT hist  #######################################
  histName = 'hGenJetPt_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  tempOriginGenJetHists = baseGenJetHists.Clone()
  tempOriginGenJetHists.SetName(histName)
  tempOriginGenJetHists.GetXaxis().SetRangeUser(\
    centRangeList[centBin][0], centRangeList[centBin][1])
  tempOriginGenJetHists.GetYaxis().SetRangeUser(0, 250)

  hGenJetPtDist = tempOriginGenJetHists.Project3D("ye")
  reBin = 1 # 5
  histSetting.histLabelSetting(hGenJetPtDist, histName, \
    '#it{p}_{T, gen}^{jet}', 'count', centBin, 1, reBin)

  addHistIntoList(lMainTree, hGenJetPtDist, histNum, histName, \
    ptHardScaleFactor, ptHardBin)

  temphGenJetPtDist = hGenJetPtDist.Clone()
  tempName = histName+'_notScale'
  temphGenJetPtDist.SetName(tempName)
  lhCheckEachJetPtBinPart[centBin].append(temphGenJetPtDist)
  
  histNum += 1
  ##############################################################################

  ###  Add detector level jet pT hist  #######################################
  histName = 'hDetJetPt_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  tempOriginDetJetHists = baseDetJetHists.Clone()
  tempOriginDetJetHists.SetName(histName)
  tempOriginDetJetHists.GetXaxis().SetRangeUser(\
    centRangeList[centBin][0], centRangeList[centBin][1])
  tempOriginDetJetHists.GetYaxis().SetRangeUser(0, 250)

  # hDetJetPtDist = tempOriginDetJetHists.ProjectionY("det"+str(centBin),0,-1,0,-1,"e")
  hDetJetPtDist = tempOriginDetJetHists.Project3D("ye")
  reBin = 1 # 5
  histSetting.histLabelSetting(hDetJetPtDist, histName, \
    '#it{p}_{T, det}^{jet}', 'count', centBin, 1, reBin)

  addHistIntoList(lMainTree, hDetJetPtDist, histNum, histName, \
    ptHardScaleFactor, ptHardBin)

  temphDetJetPtDist = hDetJetPtDist.Clone()
  tempName = histName+'_notScale'
  temphDetJetPtDist.SetName(tempName)
  lhCheckEachJetPtBinHyb[centBin].append(temphDetJetPtDist)
  
  histNum += 1
  ##############################################################################

  ##############################################################################
  tempOriginMatchedJetHists = baseMatchedJetHists[0].Clone()
  histName = 'RMHistsCent'+str(centBin)
  tempOriginMatchedJetHists.SetName(histName)
  tempOriginMatchedJetHists.GetAxis(5).SetRangeUser(\
    centRangeList[centBin][0], centRangeList[centBin][1])
  ##############################################################################

  ###  Add Matched Gen Jet hist  ###############################################
  histName = 'hMatchGenJetPt_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  hTempMatchedGenJet = tempOriginMatchedJetHists.Clone()
  hMatchGenJetPt = hTempMatchedGenJet.Projection(0, histName)
  hMatchGenJetPt.SetName(histName)
  reBin = 1 # 10
  histSetting.histLabelSetting(hMatchGenJetPt, histName, \
    '#it{p}_{T, gen}^{Matched jet}', 'count', centBin, 0, reBin)
  hMatchGenJetPt.Draw('AP')

  addHistIntoList(lMainTree, hMatchGenJetPt, histNum, histName, \
    ptHardScaleFactor, ptHardBin)
  histNum += 1
  ##############################################################################

  ###  Add Matched Det Jet hist  ###############################################
  histName = 'hMatchDetJetPt_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  hTempMatchedDetJet = tempOriginMatchedJetHists.Clone()
  hMatchDetJetPt = hTempMatchedDetJet.Projection(1, histName)
  hMatchDetJetPt.SetName(histName)
  reBin = 1 # 10
  histSetting.histLabelSetting(hMatchDetJetPt, histName, \
    '#it{p}_{T, det}^{Matched jet}', 'count', centBin, 0, reBin)
  hMatchDetJetPt.Draw('AP')

  addHistIntoList(lMainTree, hMatchDetJetPt, histNum, histName, \
    ptHardScaleFactor, ptHardBin)
  histNum += 1
  ##############################################################################

  ##############################################################################
  histName = 'hMatchJetDeltaR_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  hTempMatchJetDeltaR = tempOriginMatchedJetHists.Clone()
  hMatchJetDeltaR = hTempMatchJetDeltaR.Projection(2, histName)
  hMatchJetDeltaR.SetName(histName)
  histSetting.histLabelSetting(hMatchJetDeltaR, histName, '#Delta R', 'count', centBin, 0, 1)
  hMatchJetDeltaR.Draw('AP')

  addHistIntoList(lMainTree, hMatchJetDeltaR, histNum, histName, \
    ptHardScaleFactor, ptHardBin)
  histNum += 1
  ##############################################################################

  ##############################################################################
  histName = 'hMatchJetDAngle_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  hTemphMatchJetDAngle = tempOriginMatchedJetHists.Clone()
  hMatchJetDAngle = hTemphMatchJetDAngle.Projection(4, histName)
  hMatchJetDAngle.SetName(histName)
  histSetting.histLabelSetting(hMatchJetDAngle, histName, 'angle', 'count', centBin, 0, 1)
  hMatchJetDAngle.Draw('AP')

  addHistIntoList(lMainTree, hMatchJetDAngle, histNum, histName, \
    ptHardScaleFactor, ptHardBin)
  histNum += 1
  ##############################################################################

  ###  Add RM hist  ############################################################
  histName = 'RM_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  hTempEachCentRM = tempOriginMatchedJetHists.Clone()
  hEachCentRM = hTempEachCentRM.Projection(0, 1, "")
  hEachCentRM.SetName(histName)

  addHistIntoList(lMainTree, hEachCentRM, histNum, histName, \
    ptHardScaleFactor, ptHardBin)
  histNum += 1
  ##############################################################################

  ###  Add RM hist  ############################################################
  histName = 'hRM_forUnfold_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  hTempEachCentRM_forUnfold = baseMatchedJetHists[0].Clone()
  hEachCentRM_forUnfold = hTempEachCentRM_forUnfold.Projection(0, 1, "")
  hEachCentRM_forUnfold.SetName(histName)

  addHistIntoList(lMainTree, hEachCentRM_forUnfold, histNum, histName, \
    ptHardScaleFactor, ptHardBin)
  histNum += 1
  ##############################################################################

  ###  Add JER hist  ###########################################################
  histName = 'hRM_forJER_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  hTempEachCentRM_forJER = tempOriginMatchedJetHists.Clone()
  hEachCentRM_forJER = hTempEachCentRM_forJER.Projection(0, 1, "")
  hEachCentRM_forJER.SetName(histName)
  # hEachCentRM_forJER.Sumw2()
  # hEachCentRM_forJER.Scale(ptHardScaleFactor)
  hEachCentJER = plotPerformanceHists.getJER(hEachCentRM_forJER, 'JER_cent{}'.format(centBin))

  xTitle = '#it{p}_{T, gen}^{jet}'
  yTitle = "#frac{#sigma(#it{p}_{T}^{gen})}{#it{p}_{T}^{gen}}"
  histSetting.histLabelSetting(hEachCentJER, histName, xTitle, yTitle, centBin, 1, reBin)
  if ptHardBin == initPtHardBin: lMainTree.Add(hEachCentJER)
  else :
    tempHist = lMainTree[histNum]
    # tempHist = lMainTree.FindObject('JER_cent{}'.format(centBin))
    tempHist.Add(hEachCentJER)
  histNum += 1
  ##############################################################################

  ###  Add Det Gen Efficiency hist  ############################################
  histName = 'hDetGenJetEff_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  hDetGenJetEff = hMatchDetJetPt.Clone()
  hDetGenJetEff.SetName(histName)
  hDetGenJetEff.Divide(hMatchDetJetPt, hMatchGenJetPt, 1., 1., 'B')

  xTitle = '#it{p}_{T, gen}^{jet}'
  yTitle = "efficiency"
  histSetting.histLabelSetting(hDetGenJetEff, histName, xTitle, yTitle, centBin, 1, reBin)
  addHistIntoList(lMainTree, hDetGenJetEff, histNum, histName, 1, ptHardBin)
  histNum += 1
  ##############################################################################

  ###  Add Centrality hist  ####################################################
  histName = 'hCentHistName_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  hTempCentHist = tempOriginMatchedJetHists.Clone()
  hTempCentHist.SetName(histName)
  hCentHist = hTempCentHist.Projection(5, histName)
  
  addHistIntoList(lMainTree, hCentHist, histNum, histName, ptHardScaleFactor, ptHardBin)
  histNum += 1
  ##############################################################################

  ###  Add hJESshiftEMCal hist   ###############################################
  lJESshiftPtRange = [[20,30], [50,70], [100,120]]
  if ptHardBin == initPtHardBin:
    lJESshiftDistHists = ROOT.TList()
    lJESshiftDistHists.SetName('hJESshiftList')
    lMainTree.Add(lJESshiftDistHists)

  histName = 'hJESshiftEMCal_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  originMatchedJetHists_clone = baseMatchedJetHists.Clone()
  originMatchedJetHists_clone.SetName(histName)
  hJESshift = originMatchedJetHists_clone[1]
  hJESshift.GetXaxis().SetRangeUser(centRangeList[centBin][0], centRangeList[centBin][1])

  histName = 'hJESshift_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  hJESshift_clone = hJESshift.Clone()
  hJESshift_clone.SetName(histName)
  hJESshift_proj = hJESshift_clone.Project3D("zyeo")
  hJESshift_prof = hJESshift_proj.ProfileX()
  
  xTitle = '#it{p}_{T, gen}^{jet}'
  yTitle = "#frac{#it{p}_{T}^{det} - #it{p}_{T}^{gen}}{#it{p}_{T}^{gen}}"
  histSetting.histLabelSetting(hJESshift_prof, histName, xTitle, yTitle, centBin, 1, reBin)
  if ptHardBin == initPtHardBin: 
    lMainTree[histNum].Add(hJESshift_prof)
  else :
    tempHist = lMainTree[histNum][0]
    # tempHist = lTree.FindObject(histName)
    tempHist.Add(hJESshift_prof)

  lPtRangeColor = [1, 632, 600]
  for ptRangeKind in range(0, 3):
    histName = 'hJESshiftDist_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)\
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
      
    if ptHardBin == initPtHardBin: 
      lMainTree[histNum].Add(hJESshif_clone_proj)
    else :
      tempHist = lMainTree[histNum][ptRangeKind+1]
      # tempHist = lTree.FindObject(histName)
      tempHist.Add(hJESshif_clone_proj)
  
  histNum += 1
  ##############################################################################

  ###  Add hNEFVsPt hist   #####################################################
  histName = 'hNEFVsPt_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  originMatchedJetHists_clone = baseMatchedJetHists.Clone()
  originMatchedJetHists_clone.SetName(histName)
  hNEFVsPt = originMatchedJetHists_clone[2]
  hNEFVsPt.GetXaxis().SetRangeUser(centRangeList[centBin][0], centRangeList[centBin][1])
  hNEFVsPt_proj = hNEFVsPt.Project3D("zyeo")
  hNEFVsPt_profX = hNEFVsPt_proj.ProfileX()

  xTitle = '#it{p}_{T, gen}^{jet}'
  yTitle = "NEF"
  histSetting.histLabelSetting(hNEFVsPt_profX, histName, xTitle, yTitle, centBin, 1, reBin)
  addHistIntoList(lMainTree, hNEFVsPt_profX, histNum, histName, ptHardScaleFactor, ptHardBin)
  histNum += 1
  ##############################################################################

  ###  Add hZLeadingVsPt hist   ################################################
  histName = 'hZLeadingVsPt_' + 'Cent{0}ptHard{1}'.format(centBin, ptHardBin)
  originMatchedJetHists_clone = baseMatchedJetHists.Clone()
  originMatchedJetHists_clone.SetName(histName)
  hZLeadingVsPt = originMatchedJetHists_clone[2]
  hZLeadingVsPt.GetXaxis().SetRangeUser(centRangeList[centBin][0], centRangeList[centBin][1])
  hZLeadingVsPt_proj = hZLeadingVsPt.Project3D("zyeo")
  hZLeadingVsPt_profX = hZLeadingVsPt_proj.ProfileX()

  xTitle = '#it{p}_{T, gen}^{jet}'
  yTitle = "zLeading #it{p}_{T}"
  histSetting.histLabelSetting(hZLeadingVsPt_profX, histName, xTitle, yTitle, centBin, 1, reBin)
  addHistIntoList(lMainTree, hZLeadingVsPt_profX, histNum, histName, \
    ptHardScaleFactor, ptHardBin)
  histNum += 1
  ##############################################################################


#333333333333333333333333333333333333333333333333333333333333333333333333333333333
def checkEachPtHardDist(lMainTree, lhCheckEachJetPtBinPart, lhCheckEachJetPtBinHyb):
  lFillColor = [1, 632, 600, 416+3, 800+8, 860+10, 900+10, 820+4, 880+1, 900+9, \
    820-9, 800+1, 860+1, 840+4, 632-7, 600-7, 416-2, 800-2, 880+7, 820+4, 840-5]

  for centBin in range(0, 5):
    canvasName = "canv_ptHardBinParL_cent{}".format(centBin)
    cPartJet = ROOT.TCanvas(canvasName, canvasName,800,850)
    cPartJet.Draw()
    cPartJet.cd()

    padName = 'pad_ptHardBinParL_cent{}'.format(centBin)
    padPartJet = ROOT.TPad(padName, padName,0,0,1,1)
    padPartJet.Draw()
    padPartJet.SetLogy()
    padPartJet.cd()
  
    legName = 'ptHardBinParL_cent{}'.format(centBin)
    legPartJet = ROOT.TLegend(0.7,0.6,0.88,0.93,legName)
    legPartJet.SetFillColor(10)
    legPartJet.SetBorderSize(0)
    legPartJet.SetFillStyle(0)
    legPartJet.SetTextSize(0.04)

    lMainTree[centBin][3].SetLineColor(lFillColor[0])
    lMainTree[centBin][3].SetMarkerColor(lFillColor[0])
    lMainTree[centBin][3].SetMarkerStyle(8)
    lMainTree[centBin][3].SetMarkerSize(0.8)
    lMainTree[centBin][3].Draw()
    for ptBin in range(1, 18):
      lhCheckEachJetPtBinPart[centBin][ptBin].SetLineColor(lFillColor[ptBin+1])
      lhCheckEachJetPtBinPart[centBin][ptBin].SetMarkerColor(lFillColor[ptBin+1])
      lhCheckEachJetPtBinPart[centBin][ptBin].Draw('same')
    
    cPartJet.SaveAs('ptDistCheck/' + canvasName + '.root')
    cPartJet.Close()

  for centBin in range(0, 5):
    canvasName = "canv_ptHardBinHybL_cent{}".format(centBin)
    cHybJet = ROOT.TCanvas(canvasName, canvasName,800,850)
    cHybJet.Draw()
    cHybJet.cd()

    padName = 'pad_ptHardBinHybL_cent{}'.format(centBin)
    padHybJet = ROOT.TPad(padName, padName,0,0,1,1)
    padHybJet.Draw()
    padHybJet.SetLogy()
    padHybJet.cd()
  
    legName = 'ptHardBinHybL_cent{}'.format(centBin)
    legHybJet = ROOT.TLegend(0.7,0.6,0.88,0.93,legName)
    legHybJet.SetFillColor(10)
    legHybJet.SetBorderSize(0)
    legHybJet.SetFillStyle(0)
    legHybJet.SetTextSize(0.04)

    lMainTree[centBin][4].SetLineColor(lFillColor[0])
    lMainTree[centBin][4].SetMarkerColor(lFillColor[0])
    lMainTree[centBin][4].SetMarkerStyle(8)
    lMainTree[centBin][4].SetMarkerSize(0.8)
    lMainTree[centBin][4].Draw()
    for ptBin in range(0, 18):
      lhCheckEachJetPtBinHyb[centBin][ptBin].SetLineColor(lFillColor[ptBin+1])
      lhCheckEachJetPtBinHyb[centBin][ptBin].SetMarkerColor(lFillColor[ptBin+1])
      lhCheckEachJetPtBinHyb[centBin][ptBin].Draw('same')
    
    cHybJet.SaveAs('ptDistCheck/' + canvasName + '.root')
    cHybJet.Close()
#33333333333333333333333333333333333333333333333333333333333333333333333333333333


#333  addHistIntoList   33333333333333333333333333333333333333333333333333333333
def addHistIntoList(lTree, inputHist, histNum, histName, ptHardScaleFactor, ptHardBin):
  inputHist.Sumw2()
  inputHist.Scale(ptHardScaleFactor)

  if ptHardBin == initPtHardBin: 
    lTree.Add(inputHist)
  else :
    tempHist = lTree[histNum]
    # tempHist = lTree.FindObject(histName)
    tempHist.Add(inputHist)
#3333333333333333333333333333333333333333333333333333333333333333333333333333333



# check pt bin width start  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def ptRangeDef(jetList, centBin, initCentBin, endCentBin):
  for cent10DivBin in range(initCentBin, endCentBin):
    hTempMeasRawJetCorrPt = jetList.FindObject('histJetCorrPt_{}'.format(initCentBin))
    hEachCentRawJetCorrPt = hTempMeasRawJetCorrPt.Clone()

    if cent10DivBin == initCentBin :
      hEachMergeCentRawJetCorrPt = hEachCentRawJetCorrPt
      hEachMergeCentRawJetCorrPt.SetName('hMeasRawJetCorrPt_Cent{}'.format(centBin))
    else :
      hEachMergeCentRawJetCorrPt.Add(hEachCentRawJetCorrPt)

  # for ptBin in range(0, 250):
      # binNum = hEachMergeCentRawJetCorrPt.FindBin(ptBin)
      # print('pt: {0}, bin: {1}'.format(ptBin, binNum))

  binWidth5 = hEachMergeCentRawJetCorrPt.FindBin(5)\
    - hEachMergeCentRawJetCorrPt.FindBin(0)
  # print('binWidth5: {}'.format(binWidth5))

  nBins = hEachMergeCentRawJetCorrPt.GetXaxis().GetNbins()
  binContents = 0
  finalBin = False
  lPtBin = []
  lCountBin = []
  # print('nBins: {}'.format(nBins))
  lPtBin.append(250.)
  for bin in range(nBins, 500, -binWidth5):
    for bin5range in range(0,binWidth5):
      if(bin > 500):
          binContents += hEachMergeCentRawJetCorrPt.GetBinContent(bin)
          bin -= 1
          # print('bin : {0}, binContents: {1}'.format(bin, binContents))
      else: 
          binContents += hEachMergeCentRawJetCorrPt.GetBinContent(bin)
          finalBin = True

    if (binContents > 50) or finalBin:
      lCountBin.append(binContents)
      inputPt = hEachMergeCentRawJetCorrPt.GetBinLowEdge(bin+1)
      lPtBin.append(inputPt)
      # print('bin : {0}, pT : {1}, binContents : {2}, finalBin : {3}'\
        # .format(bin+1, inputPt, binContents, finalBin))
      binContents = 0
      
  lPtBin.reverse()
  lPtBin_d = array('d', lPtBin)
  lCountBin.reverse()
  lCountBin_d = array('d', lCountBin)
  print(lPtBin_d)
  print(lCountBin_d)
  histName = 'hMeasRawJetCorrPt_Cent{}_rebin'.format(centBin)
  hEachCentRawJetCorrPt_rebin = ROOT.TH1D(histName, histName, len(lPtBin_d)-1, lPtBin_d)
  for ptBin in range(0, len(lPtBin_d)-1):
    hEachCentRawJetCorrPt_rebin.Fill(lPtBin_d[ptBin], lCountBin_d[ptBin])
    # print('pt : {0}, count: {1}'.format(lPtBin_d[ptBin], lCountBin_d[ptBin]))

  return hEachCentRawJetCorrPt_rebin, lPtBin
# check pt bin width end      $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#11  histExtraction   1111111111111111111111111111111111111111111111111111111111
def dataRawJetExtranctV2(lMainTree, FullCh, centBin, rBin):
  TAADiv10 = [23.26, 14.4, 8.767,  5.086, 2.747, 1.352,  0.5992,  0.2385, 0.08383, 0.02527]
  # resoPara = [2, 3, 4, 5]

  f = ROOT.TFile(inputDataFileName, "READ")

  cutPt = 5
  taskName = 'AliAnalysisTaskEmcalJetPerformance'\
    + '_tracks' + '_caloClustersCombined' \
      + '_0{}'.format(resoPara[rBin]) + '_{}'.format(cutPt) + 'GeV_histos'

  task = f.Get(taskName)
  if not task : print("ERROR no task found")

  lkeys = f.GetListOfKeys()
  lkeysName = []
  # print(" the following lists:")
  for key in lkeys:
    name = key.GetName()
    print(name)
  f.Close()

  if FullCh is 'Full': jetTreeName = 'Jet_AKTFullR0{}0'.format(resoPara[rBin])\
    + '_tracks_pT0150' + '_caloClusters_E0300_pt_scheme'
  elif FullCh is 'Ch': jetTreeName = 'Jet_AKTChargedR0{}0'.format(resoPara[rBin])\
    + '_tracks_pT0150' + '_pt_scheme'
  # jetList = task.FindObject(jetTreeName)
  jetList = task[6]
  if not jetList: print("ERROR no list found")

  # lJetHistgrams = jetList[2]

  lCentMergeBin = [[0,1],[2,3],[3,5],[5,7],[7,9]]
  initCentBin = lCentMergeBin[centBin][0]
  endCentBin = lCentMergeBin[centBin][1]

  # hTempMeasJetCorrPt = lJetHistgrams[3]
  hTempMeasJetCorrPt = jetList[2]

  for cent10DivBin in range(initCentBin, endCentBin):
    hEachCentJetCorrPt3D = hTempMeasJetCorrPt.Clone()
    hEachCentJetCorrPt3D.GetXaxis().SetRangeUser(initCentBin, endCentBin)
    hEachCentJetCorrPt = hEachCentJetCorrPt3D.Project3D("y")
    hEachCentJetCorrPt.Scale(1/TAADiv10[cent10DivBin])

    if cent10DivBin == initCentBin:
      hEachMergeCentJetCorrPt = hEachCentJetCorrPt.Clone()
      hEachMergeCentJetCorrPt.SetName('hMeasJetCorrPt_Cent{}'.format(centBin))
    else :
      hEachMergeCentJetCorrPt.Add(hEachCentJetCorrPt)
    
    print('centBin:{0}, initCentBin:{1}, endCentBin:{2}, cent10DivBin:{3}'\
      .format(centBin, initCentBin, endCentBin, cent10DivBin))

  hEachMergeCentJetPt = hEachMergeCentJetCorrPt.Clone()
  hEachMergeCentJetPt.SetName('hEachMergeCentJetPt')
  lMainTree.Add(hEachMergeCentJetPt)
  lMainTree.Add(hEachMergeCentJetCorrPt)
#1111111111111111111111111111111111111111111111111111111111111111111111111111111


if __name__ == "__main__":
  createPlotEmb()


