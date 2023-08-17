#! /usr/bin/env python

# General
import os
import sys
import argparse
import itertools
import math
from array import *
import numpy

import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F

# Load the RooUnfold library
# ROOT.gSystem.Load("$ALICE_WORK_DIR/osx_x86-64/RooUnfold/latest/lib/libRooUnfold.dylib")
# ROOT.gSystem.Load("$ALIBUILD_WORK_DIR/slc7_x86-64/RooUnfold/latest/lib/libRooUnfold.so")
ROOT.gSystem.Load("/persist/sw/slc7_x86-64/RooUnfold/latest/lib/libRooUnfold.so")

from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
from ROOT import RooUnfoldSvd
from ROOT import TCanvas
from ROOT import gROOT

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

import UnfoldPlots
import PtRangeList

isClosureTest = False
isClosure  = False

#[kBrack, kRed, kBlue, kGreen+3, kOrange+8, kAzure+10, kSpring+4, kViolet+1, kPink+9, kSpring-9, kOrange+8, kAzure+1, kTeal+4, kRed-7, kBlue-7, kGreen-2, kOrange-2, kViolet+7, kSpring+4, kTeal-5]
lFillColor = [1, 632, 600, 416+3, 800+8, 860+10, 900+10, 820+4, 880+1, 900+9, 820-9, 800+1, 860+1, 840+4, 632-7, 600-7, 416-2, 800-2, 880+7, 820+4, 840-5]


lEPLabel = ['OutOfPlane', 'InPlane', 'Inclusive']
leadingTrackPtCut = 5 # 0, 5, 7 GeV/c
diffSys = 'V0A' # 98%: '', 94%: 'TrackEff094', 'BKGNoFit', 'BKGV2', 'V0C', 'V0A'
lDiffSys = ['TrackEff094', 'BKGNoFit', 'BKGV2', 'V0C', 'V0A', '']

JetR = 2 # [2, 3, 4, 5]
lPtCut = [0, 3, 5, 7, 10]
# lCentDivKind = [[0, 1], [1, 2, 3], [3, 4, 5], [5, 6, 7], [7, 8, 9, 10]]
lCentDivKind = [[0, 5], [5, 10], [10, 30], [30, 50], [50, 80]]

# pT range 5 sigma [32., 26., 20., 10., 5.5]

numOfCentBin = 5

# mb^-1 (value taken from http://cds.cern.ch/record/2636623)
# TAA = [23.26, 14.4, 8.767,  5.086, 2.747, 1.352,  0.5992,  0.2385, 0.08383, 0.02527]
TAA = [26.08, 20.44, 11.58, 3.92, 0.89] # [0-5, 5-10, 10-30, 30-50, 50-80]

regPara = 6 # The number of unfold iteration 

inputDir = '/Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18qr/pass3/Ch/Embedding/'
# inputDir = './'

outputDir = '/Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18qr/pass3/Ch/Unfolded/'
# outputDir = './'
oPerformanceDir = outputDir + 'UnfoldPerformancePlot/'


### Main START  ################################################################
def UnfoldJetPtDist():
    for diffBin in range(0, 6):
        diffSys = lDiffSys[diffBin]

        inputFileName = 'EmbedPtHardScaledResults'\
            +'_TrackPtCut'+str(leadingTrackPtCut)+'_'+diffSys+'_Ver1.root'
        outputFileName =  outputDir + 'UnfoldedPtDists'\
            +'_TrackPtCut'+str(leadingTrackPtCut)+'_'+diffSys+'.root'

        ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;"); #..to supress a lot of standard output
        baseLabel = ''

        #Things to set for
        smearErrors    = True
        imposePrior    = False
        useRMflexibleR = False  
        matchingDistance = 0.6   #in % of jet radius
        
        inputFile = ROOT.TFile(inputDir + inputFileName, "READ")
        mainTree = inputFile.Get('mainTree')

        outputFile = ROOT.TFile(outputFileName, "RECREATE")
        lOMainTree = ROOT.TList()
        OFileStracture(lOMainTree, numOfCentBin)
        outputFile.cd()

        for epBin in range(0, 3):
            obsKind = 0
            if epBin == 2: obsKind = 1
            for centBin in range(0, numOfCentBin):
                ptRangeDict, ptBinArrayDict = PtRangeList.eachJetPtBinDef(obsKind, centBin)

                label = ''
                label = baseLabel
                # label = label + '_cent{}'.format(centBin) + '_' +lEPLabel[epBin]

                nEventsData, nEventsResponse, vertexEfficiency = eventInfoExtract(mainTree, centBin)
                scalePraDict = {'nEventsData':nEventsData, 'nEventsResponse':nEventsResponse,\
                    'vertexEfficiency':vertexEfficiency}
                print(scalePraDict)
                
                # mainTree[epBin][centBin].ls()
                hRMDict = getResponseMatrix(mainTree[epBin][centBin], centBin, ptRangeDict, ptBinArrayDict, label)

                jetSpectDict = getJetSpectrums(mainTree[epBin][centBin], hRMDict, \
                    ptRangeDict, ptBinArrayDict, scalePraDict, centBin, label)
                
                responsDict = getResponce(mainTree[epBin][centBin], hRMDict, jetSpectDict, label)
                
                label =  lEPLabel[epBin]+'/Cent'+str(centBin)+'/'
                hKinEff = UnfoldPlots.plotKinematicEfficiency(hRMDict['main'], hRMDict['uncut'], \
                    jetSpectDict, ptRangeDict, ptBinArrayDict, label, oPerformanceDir)
                lOMainTree[epBin][centBin].Add(hKinEff)

                unfoldType = 'Bayes'
                hUnfoldedDistDictBayes = UnfoldJetSpectrum(hRMDict, responsDict, jetSpectDict, hKinEff,\
                    ptRangeDict, ptBinArrayDict, unfoldType, label, epBin, centBin)
                unfoldType = 'SVD'
                # hUnfoldedDistDictSVD = UnfoldJetSpectrum(hRMDict, responsDict, jetSpectDict, hKinEff, \
                #     ptRangeDict, ptBinArrayDict, unfoldType, label, epBin, centBin)
                
                outputFile.cd()
                for key in hUnfoldedDistDictBayes.values():
                    lOMainTree[epBin][centBin][0].Add(key)
                # for key in hUnfoldedDistDictSVD.values():
                #     lOMainTree[epBin][centBin][0].Add(key)
                for key in hRMDict.values():
                    lOMainTree[epBin][centBin][1].Add(key)
                for key in jetSpectDict.values():
                    lOMainTree[epBin][centBin][2].Add(key)
                for key in responsDict.values():
                    lOMainTree[epBin][centBin][3].Add(key)
                    # print(key)
                inputFile.cd()

        UnfoldPlots.plotRMRatio(lOMainTree, numOfCentBin)

        inputFile.cd()

        outputFile.cd()
        lOMainTree.Write('mainTree', 1)
        outputFile.Close()

        print('root ' + outputDir)
### Main END  ##################################################################


################################################################################
###      get Response Matrix                                                 ###
###############################################################################
def getResponseMatrix(mainCentTree, centBin, ptRangeDict, ptBinArrayDict, label):
    # hOrigRMName = 'hRM_forUnfold' + '_Cent' + str(centBin) + '_RebinnedNormed'
    hOrigRMName = 'hRM' + '_Cent' + str(centBin)
    hRM = mainCentTree.FindObject(hOrigRMName) # No cut
    # histName = 'hRM' + label
    # hRM.SetName(histName)

    hRMUncut = rebinResponseMatrix(hRM, ptBinArrayDict, 'mcFullDet')
    histName = 'hRM_uncut' + '_Cent' + str(centBin)+ label
    hRMUncut.SetName(histName)

    hRMMain = rebinResponseMatrix(hRM, ptBinArrayDict, 'mcDet')
    histName = 'hRM_Main' + '_Cent' + str(centBin)+ label
    hRMMain.SetName(histName)
    
    histName = 'hRM_reBin' + '_Cent' + str(centBin)+ label
    hRMReBin = hRMMain.Clone(histName)
    hRMNorm = normalizeResponseMatrix(hRMReBin)
    histName = 'hRM_nomali' + '_Cent' + str(centBin)+ label
    hRMNorm.SetName(histName)

    hRMDetLvLcut = rebinResponseMatrix(hRM, ptBinArrayDict, 'mcDetL')
    histName = 'hRM_DetLvLcut' + '_Cent' + str(centBin)+ label
    hRMDetLvLcut.SetName(histName)
    # hRMDetLvLcut = normalizeResponseMatrix(hRMDetLvLcut)  #?????

    hRMDetLvUcut = rebinResponseMatrix(hRM, ptBinArrayDict, 'mcDetU')
    histName = 'hRM_DetLvUcut' + '_Cent' + str(centBin)+ label
    hRMDetLvUcut.SetName(histName)
    # hRMDetLvUcut = normalizeResponseMatrix(hRMDetLvUcut)  #?????

    # Get response matrix from response file (Measured, True) to be used for the unfolding, 
    # with pT-det range cut to desired range, and re-bin. 
    # Also get response matrices for refolding/closure tests.
    hRMReFold = hRMMain.Clone()
    histName = 'hRM_for_Refolding' + '_Cent' + str(centBin)+ label
    hRMReFold.SetName(histName)
    
    hRMReFoldComp = hRMMain.Clone()
    histName = 'hRM_for_RefoldComp' + '_Cent' + str(centBin)+ label
    hRMReFoldComp.SetName(histName)
    
    hRMDict = {'Orig':hRM, 'rebin':hRMReBin, 'main':hRMMain, 'DetLvLcut':hRMDetLvLcut, 'DetLvUcut':hRMDetLvUcut, \
        'nomali':hRMNorm, 'uncut':hRMUncut, 'reFold':hRMReFold, 'reFoldComp':hRMReFoldComp}
    
    return hRMDict
################################################################################


################################################################################
###      get JetSpectrums                                                    ###
################################################################################
def getJetSpectrums(mainCentTree, hRMDict,ptRangeDict,ptBinArrayDict,scalePraDict,centBin, label):
    dEta = 2 * (0.9 - (0.1* JetR))
    hOrigRawJetName = 'hJetCorrPtLocal_' + str(centBin)
    # lCentChange = [0,1,2,4,6]
    # hOrigRawJetName = 'hJetCorrPtLocal_' + str(lCentChange[centBin])
    hOrigMeasuredJetPerBin = mainCentTree.FindObject(hOrigRawJetName)

    hMeasuredJetPerBin = hOrigMeasuredJetPerBin.Clone() # w/ background subtraction
    histName = 'hMeasuredJet' + label
    hMeasuredJetPerBin.SetName('hMeasuredJet')
    nBins = len(ptBinArrayDict['mcDet']) - 1
    binArray = ptBinArrayDict['mcDet']
    hMeasuredJetPerBin.Sumw2()

    #rebin the measured input spectrum to the final binning 
    # and save it for statistical error comparision
    hMeasuredJet_rebin = hMeasuredJetPerBin.Clone()
    histName = 'hMeasuredJet_rebin' + label
    hMeasuredJet_rebin.SetName(histName)
    nBins = len(ptBinArrayDict['mcDet']) - 1
    binArray = ptBinArrayDict['mcDet']
    hMeasuredJet_rebin.Sumw2()
    hMeasuredJet_rebin = hMeasuredJet_rebin.Rebin(nBins, histName, binArray)

    hMeasuredJetYieldPerBin = hMeasuredJet_rebin.Clone()
    histName = 'hMeasuredJetYield' + label
    hMeasuredJetYieldPerBin.SetName(histName)
    hMeasuredJetYieldPerBin.Sumw2()
    hMeasuredJetYieldPerBin.Scale(1./scalePraDict['nEventsData'])
    hMeasuredJetYieldPerBin.Scale(1./TAA[centBin])
    # hMeasuredJetYieldPerBin.Scale(1./dEta)
    hMeasuredJetYieldPerBin.Scale(dEta)

    #rebin the measured input spectrum to the final binning 
    # and save it for statistical error comparision
    hMeasuredJetYield_LCut = hMeasuredJetPerBin.Clone()
    histName = 'hMeasuredJetYield_LCut' + label
    hMeasuredJetYield_LCut.SetName(histName)
    nBins = len(ptBinArrayDict['mcDetL']) - 1
    binArray = ptBinArrayDict['mcDetL']
    hMeasuredJetYield_LCut.Sumw2()
    hMeasuredJetYield_LCut = hMeasuredJetYield_LCut.Rebin(nBins, histName, binArray)
    hMeasuredJetYield_LCut.Scale(1./scalePraDict['nEventsData'])
    hMeasuredJetYield_LCut.Scale(1./TAA[centBin])
    # hMeasuredJetYield_LCut.Scale(1./dEta)
    hMeasuredJetYield_LCut.Scale(dEta)

    hMeasuredJetYield_UCut = hMeasuredJetPerBin.Clone()
    histName = 'hMeasuredJetYield_UCut' + label
    hMeasuredJetYield_UCut.SetName(histName)
    nBins = len(ptBinArrayDict['mcDetU']) - 1
    binArray = ptBinArrayDict['mcDetU']
    hMeasuredJetYield_UCut.Sumw2()
    hMeasuredJetYield_UCut = hMeasuredJetYield_UCut.Rebin(nBins, histName, binArray)
    hMeasuredJetYield_UCut.Scale(1./scalePraDict['nEventsData'])
    hMeasuredJetYield_UCut.Scale(1./TAA[centBin])
    # hMeasuredJetYield_UCut.Scale(1./dEta)
    hMeasuredJetYield_UCut.Scale(dEta)

    # Get the truth-level jet spectrum (matched) from response matrix (already re-binned)
    # Do exclude under and overflow bins
    hMCGenJetPerBin = hRMDict['main'].ProjectionY()
    hMCGenJetPerBin.SetName('hMCGenJetPerBin' + label)

    hMCGenJetUncutPerBin = hRMDict['uncut'].ProjectionY()
    hMCGenJetUncutPerBin.SetName('hMCGenJetUncutPerBin' + label)
    hMCGenJetUncutPerBin.Sumw2()
    
    hMCGenJetYieldPerBin = hRMDict['main'].ProjectionY()
    histName = 'hMCGenJetYieldPerBin' + label
    hMCGenJetYieldPerBin.SetName(histName)
    hMCGenJetYieldPerBin.Sumw2()
    # nBins = len(ptBinArrayDict['mcDet']) - 1
    # binArray = ptBinArrayDict['mcDet']
    nBins = len(ptBinArrayDict['mcGen']) - 1
    binArray = ptBinArrayDict['mcGen']
    hMCGenJetYieldPerBin = hMCGenJetYieldPerBin.Rebin(nBins, histName, binArray)
    # hMCGenJetYieldPerBin.Scale(1./scalePraDict['nEventsResponse'])

    hDiffPriorPerBin, hApplyChagePriorRatio \
        = changePrior(hRMDict, hOrigMeasuredJetPerBin, ptBinArrayDict, centBin)
    hApplyChagePriorRatio.SetName('hApplyChagePriorRatio'+'CentBin'+str(centBin))

    # Get the det-level jet spectrum (matched) from response matrix projection, 
    # after cutting the pT-det range (already re-binned)
    # Note that this is potentially sensitive to the low-pT cutoff of the RM 
    # -- but it is only used for "det-level" Raa plot as a sanity check
    hMCDetJetYieldPerBin = hRMDict['main'].ProjectionX()
    hMCDetJetYieldPerBin.SetName('hMCDetJetYieldPerBin' + label)
    hMCDetJetYieldPerBin.Sumw2()
    # hMCDetJetYieldPerBin.Scale(1./scalePraDict['nEventsResponse'])

    jetSpectDict = {'measYieldPerBin':hMeasuredJetYieldPerBin, 
        'measYieldPerBin_LCut':hMeasuredJetYield_LCut, 'measYieldPerBin_UCut':hMeasuredJetYield_UCut,\
            'mcGenPerBin':hMCGenJetPerBin, 'mcGenJetUncutPerBin':hMCGenJetUncutPerBin, \
                'mcGenJetDiffPriorPerBin':hDiffPriorPerBin,
                    'mcGenYieldPerBin':hMCGenJetYieldPerBin, 'mcDetYieldPerBin':hMCDetJetYieldPerBin,
                    'priorChangeRatio': hApplyChagePriorRatio}
    
    # For Measured spectra, divide by bin width to produce per GeV spectra
    hMeasuredJetYieldPerGeV = jetSpectDict['measYieldPerBin'].Clone()
    histName = "hMeasuredJetYieldPerGeV" + label
    hMeasuredJetYieldPerGeV.SetName(histName)
    hMeasuredJetYieldPerGeV.Scale(1., "width")
    jetSpectDict['measYieldPerGeV'] =  hMeasuredJetYieldPerGeV

    hMCGenJetYieldPerGeV = jetSpectDict['mcGenYieldPerBin'].Clone()
    hMCGenJetYieldPerGeV.SetName('hMCGenJetYieldPerGeV' + label)
    hMCGenJetYieldPerGeV.Scale(1., "width")
    jetSpectDict['mcGenYieldPerGeV'] =  hMCGenJetYieldPerGeV
    
    hMCDetJetYieldPerGeV = jetSpectDict['mcDetYieldPerBin'].Clone()
    hMCDetJetYieldPerGeV.SetName('hMCDetJetYieldPerGeV' + label)
    hMCDetJetYieldPerGeV.Scale(1., "width")
    jetSpectDict['mcDetYieldPerGeV'] =  hMCDetJetYieldPerGeV
    
    return jetSpectDict 
################################################################################


################################################################################
## retrun Response matrix                                                   ####
################################################################################
def getResponce(mainTree, hRMDict, jetSpectDict, label):
    hJetMCGenUncutPerBin = jetSpectDict['mcGenJetUncutPerBin']
    hJetMCGenDiffPriorPerBin = jetSpectDict['mcGenJetDiffPriorPerBin']

    histName = hRMDict['main'].GetName() + '_Response' + label
    response = RooUnfoldResponse(0, hJetMCGenUncutPerBin, hRMDict['main'], histName, histName)
    response.UseOverflow(False) # ??????????????????????????????????????????????

    histName = hRMDict['main'].GetName() + '_ResponseDiffPrior' + label
    responseDiffPrior = RooUnfoldResponse(0, hJetMCGenDiffPriorPerBin, hRMDict['changePrior'], histName, histName)
    responseDiffPrior.UseOverflow(False) # ??????????????????????????????????????????????
    
    histName = hRMDict['DetLvLcut'].GetName() + '_Response' + label
    responseDetLvLcut = RooUnfoldResponse(0, hJetMCGenUncutPerBin, hRMDict['DetLvLcut'], histName, histName)
    responseDetLvLcut.UseOverflow(False) # ??????????????????????????????????????????????
    histName = hRMDict['DetLvUcut'].GetName() + '_Response' + label
    responseDetLvUcut = RooUnfoldResponse(0, hJetMCGenUncutPerBin, hRMDict['DetLvUcut'], histName, histName)
    responseDetLvUcut.UseOverflow(False) # ??????????????????????????????????????????????

    histName = 'hRMUncut' + label
    respUncut = RooUnfoldResponse(0, 0, hRMDict['uncut'], histName, histName)
    histName = 'hRMReFold' + label
    respReFold = RooUnfoldResponse(0, 0, hRMDict['reFold'], histName, histName)
    histName = 'hRMReFoldComp' + label
    respReFoldComp = RooUnfoldResponse(0, 0, hRMDict['reFoldComp'], histName, histName)

    responsDict = {'main':response, 'diffPrior':responseDiffPrior, \
        'DetLvLcut':responseDetLvLcut,'DetLvUcut':responseDetLvUcut,\
        'uncut':respUncut, 'reFold':respReFold, 'reFoldComp':respReFoldComp}

    return responsDict
################################################################################


################################################################################
###  prepareToUnfold                                                        ####
################################################################################
def setPrior(hRM, hJetSpectTruUncutPerBin, label):
    # projection to the intended prior distribution
    # Keep response matrix as per-bin probabilities (i.e. don't scale by bin width)
    # Scale also hJetSpectrumTrueUncutPerBin accordingly, 
    # soas to preserve the kinematic efficiency 
    # (which we will use to create the RooUnfoldResponse object)

    # Loop through truth-level bins (of RM and uncut truth spectrum), 
    # and normalize each bin (of RM and uncut truth spectrum) to the prior
    nBinsY = hRM.GetNbinsY() # pT-gen
    nBinsX = hRM.GetNbinsX() # pT-det

    # Determine the scale factor for each bin and then scale the response and the truth spectrum
    
    for truthBin in range(1,nBinsY+1):
        #this is a cut projection
        # truthBinCenter  = hTruthProjectionBefore.GetXaxis().GetBinCenter(truthBin)
        #this is a cut projection
        # truthProjVal    = hTruthProjectionBefore.GetBinContent(truthBin)            
        truthBinCenter  = hJetSpectTruUncutPerBin.GetXaxis().GetBinCenter(truthBin)
        truthProjVal    = hJetSpectTruUncutPerBin.GetBinContent(truthBin)
        priorScalingVal = 1
        # the unfolded final result as a prior
        # priorScalingVal = math.pow(truthBinCenter, powerLawOffset)

        # Set uncut truth spectrum to prior (to preserve the kinematic efficiency)
        uncutContent = hJetSpectTruUncutPerBin.GetBinContent(truthBin)
        # hJetSpectTruUncutPerBin.SetBinContent(truthBin, uncutContent * priorScalingVal)
        hJetSpectTruUncutPerBin.SetBinContent(truthBin, uncutContent)
        # Set RM to prior
        for detBin in range(1,nBinsX+1):
            binContent = hRM.GetBinContent(detBin, truthBin)
            hRM.SetBinContent(detBin, truthBin, binContent * priorScalingVal)

    # Set up the RooUnfoldResponse object
    # response = RooUnfoldResponse(0, 0, hResponseMatrix, "hResponseMatrixMain", "hResponseMatrixMain")
    # One can supply the measured and truth distributions is to incorporate fakes and inefficiency
    # For the truth-level, 
    # we pass in the projection before the RM was cut to a specific pT-det range,
    # in order to account for the kinematic efficiency
################################################################################

def changePrior(hRMOrig, hMeasuredJetSpectPerBinOrig, ptBinArrayDict, centBin):
    nBins = len(ptBinArrayDict['mcDet']) - 1
    binArray = ptBinArrayDict['mcDet']

    hRM = hRMOrig['uncut'].Clone('hRMForChangePrior')

    hMCDet = hRM.ProjectionX()
    hMCDet.Sumw2()
    hMCDet = hMCDet.Rebin(nBins, 'hMCDet_px', binArray)
    hMCDet.Scale(1., "width")
    hMCDet.Scale(1./hMCDet.Integral())
    # hMCDet.SaveAs('hogeMCDet.root')
    
    hApplyRatio = hMeasuredJetSpectPerBinOrig.Clone('hMeasForChangePrior')
    hApplyRatio.Sumw2()
    hApplyRatio = hApplyRatio.Rebin(nBins, 'hApplyRatio', binArray)
    if centBin==3:hApplyRatio.SetName('hApplyRatio')
    hApplyRatio.Scale(1./hApplyRatio.Integral())
    hApplyRatio.Scale(1., "width")
    # if centBin==3:hApplyRatio.SaveAs('hogeMeas.root')
    hApplyRatio.Divide(hMCDet)
    # if centBin==3:hApplyRatio.SaveAs('hogeRatio.root')
    
    fRatioFit = ROOT.TF1('fRatioFit',"[0]*x*x + [1]*x + [2]",binArray[0],binArray[-1])
    hApplyRatio.Fit('fRatioFit')
    fRatioFitP0 = fRatioFit.GetParameter(0) 
    fRatioFitP1 = fRatioFit.GetParameter(1)
    fRatioFitP2 = fRatioFit.GetParameter(2)

    hMCGenReweightPrior = hRM.ProjectionY()
    hMCGenReweightPrior.Sumw2()
    nBins = len(ptBinArrayDict['mcGen']) - 1
    binArray = ptBinArrayDict['mcGen']
    hMCGenReweightPrior = hMCGenReweightPrior.Rebin(nBins, 'ReweightPrior_py', binArray)
    hMCGenReweightPrior.SetName('ReweightPrior_py')

    nXBins = hRM.GetNbinsX()
    nYBins = hMCGenReweightPrior.GetNbinsX()
    for ptTruBin in range(1,nYBins+1):
        genPtVal = hMCGenReweightPrior.GetBinCenter(ptTruBin)
        genCount = hMCGenReweightPrior.GetBinContent(ptTruBin)
        genPtError = hMCGenReweightPrior.GetBinError(ptTruBin)
        scaleVal = fRatioFitP0*genPtVal*genPtVal + fRatioFitP1*genPtVal +fRatioFitP2
        genCountScaled = genCount * scaleVal
        # print('scaleVal = '+str(scaleVal))
        # print('genCountScaled = '+str(genCountScaled))
        
        hMCGenReweightPrior.SetBinContent(ptTruBin, genCountScaled)
        # hMCGenReweightPrior.SetBinError(ptTruBin, genPtError)
        hMCGenReweightPrior.SetBinError(ptTruBin, 0.)

        for ptDetBin in range(1,nXBins+1):
            hRMContent = hRM.GetBinContent(ptDetBin,ptTruBin)
            hRMErr = hRM.GetBinError(ptDetBin,ptTruBin)
            
            hRM.SetBinContent(ptDetBin, ptTruBin, scaleVal*hRMContent)
            hRM.SetBinError(ptDetBin,ptTruBin, hRMErr)

    hRMReWeight = rebinResponseMatrix(hRM, ptBinArrayDict, 'mcDet')
    # hRMReWeight = normalizeResponseMatrix(hRMReWeight)  #?????
    hRMOrig['changePrior'] = hRMReWeight
    
    return hMCGenReweightPrior, hApplyRatio

################################################################################
################################################################################
def eventInfoExtract(mainTree, centBin):
    # Get the number of accepted events in data in the specified centrality range
    # hCentEvents = mainTree.FindObject("Centrality_rawData")
    # hCentEvents = mainTree.FindObject("Centrality_raw")
    hCentEvents = mainTree.FindObject("fHistCentrality")
    hCentEventsClone = hCentEvents.Clone()
    # hCentEventsClone.SetName("Centrality_raw"+'_Cent{}'.format(centBin))
    hCentEventsClone.SetName("fHistCentrality"+'_Cent{}'.format(centBin))
    hCentEventsClone.GetXaxis().SetRangeUser(\
        lCentDivKind[centBin][0], lCentDivKind[centBin][1])
    nEventsData = hCentEventsClone.Integral()
    print("N data events (in {}-{}% centrality): {}"\
        .format(lCentDivKind[centBin][0], lCentDivKind[centBin][1], nEventsData))

    # Get the average number of accepted events per pT-hard bin (in the embedded data)
    # This is used to scale the response spectra, due to convention in scaling the pT hard bins
    histNEventResponse = mainTree.FindObject("hNEventsAcc")
    if not histNEventResponse:
        print("couldn't find NEventResponse :(")

    PtHardBins = 21
    nEventsAccSum = 0.
    for bin in range(0,PtHardBins):
        nEventsAccSum += histNEventResponse.GetBinContent(bin+1)

    nEventsResponse = nEventsAccSum/PtHardBins

    # hNormalizationHist = mainTree.FindObject("fNormalisationHistData")
    hNormalizationHist = mainTree.FindObject("fNormalisationHist")

    #take only events that fulfill general selection criteria mostly pile-up 
    # as a reference for the vertex eff. calculation
    allEvents          = hNormalizationHist.GetBinContent(\
        hNormalizationHist.GetXaxis().FindBin("Event selection"))
    EvtsWithGoodVertex = hNormalizationHist.GetBinContent(\
        hNormalizationHist.GetXaxis().FindBin("Vertex reconstruction and quality"))
    # do not take more than 3 numbers after decimal
    vertexEfficiency   = (round(EvtsWithGoodVertex/allEvents,3))

    print("N response events (avg per pT-hard bin): %d" % nEventsResponse)

    return nEventsData, nEventsResponse, vertexEfficiency
################################################################################


################################################################################
# Rebin the response matrix to have variable binning                          ##
################################################################################
def rebinResponseMatrix(hRM, ptBinArrayDict, detRangeType):
    histname = "{}NewRebinned_{}".format(hRM.GetName(), detRangeType)
    title = histname + ";#it{p}_{T,corr}^{det} (GeV/#it{c});#it{p}_{T}^{truth} (GeV/#it{c})"
    genNBins = len(ptBinArrayDict['mcGen']) - 1
    genBinArray = ptBinArrayDict['mcGen']
    detNBins = len(ptBinArrayDict[detRangeType]) - 1
    detBinArray = ptBinArrayDict[detRangeType]
    hRM_rebin = ROOT.TH2D(histname, title, detNBins, detBinArray, genNBins, genBinArray)
    
    # Loop over all bins in fine-binned response matrix, 
    # and fill appropriate bin in new response matrix
    # Assume that the bin edges overlap appropriately
    for ibin in range(1, hRM.GetNbinsX() + 1):
        for jbin in range(1, hRM.GetNbinsY() + 1):
            oldContent = hRM.GetBinContent(ibin, jbin)
            
            # Find the bin that should be filled in the new histogram, and fill it
            # Need to get (x,y) location from bins (ibin, jbin)
            x = hRM.GetXaxis().GetBinCenter(ibin)
            y = hRM.GetYaxis().GetBinCenter(jbin)
            if ((x > ptBinArrayDict[detRangeType][0]) and (x < ptBinArrayDict[detRangeType][-1])) and ((y > ptBinArrayDict['mcGen'][0])and(y < ptBinArrayDict['mcGen'][-1])):
                hRM_rebin.Fill(x, y, oldContent)

            #print "Adding value {} from bin ({},{}) = ({},{})".format(oldContent, ibin, jbin, x, y)
            #newBin = hResponseMatrixNew.FindBin(x,y)
            #print "New bin content: {}".format(hResponseMatrixNew.GetBinContent(newBin))
    
    # Assume 0 errors on response matrix
    for bin in range(1, hRM_rebin.GetNcells() + 1):
        hRM_rebin.SetBinError(bin, 0)

    return hRM_rebin
################################################################################


################################################################################
# Normalize response matrix                                                   ##
# Normalize the pT-truth projection to 1                                      ##
################################################################################
def normalizeResponseMatrix(hResponseMatrix):
    # Make projection onto pT-true axis (y-axis), and scale appropriately
    # Do exclude under and overflow bins
    hGenProjBefore = hResponseMatrix.ProjectionY("_py", 1, hResponseMatrix.GetNbinsX()) 
    hGenProjBefore.SetName("hGenProjectionBeforeOf" + hResponseMatrix.GetName())

    # Loop through truth-level bins, and apply normalization factor to all bins.
    nBinsY = hResponseMatrix.GetNbinsY() # pT-gen
    nBinsX = hResponseMatrix.GetNbinsX() # pT-det
    for genBin in range(1,nBinsY+1):
        normFactor = hGenProjBefore.GetBinContent(genBin)
        if normFactor > 0:
            genBinCenter = hGenProjBefore.GetXaxis().GetBinCenter(genBin)

            for detBin in range(1,nBinsX+1):
                binContent = hResponseMatrix.GetBinContent(detBin, genBin)
                hResponseMatrix.SetBinContent(detBin, genBin, binContent/normFactor)

    return hResponseMatrix
################################################################################

################################################################################
# Unfold jet Spect, Bayes or SVD (specify type)                              ###
################################################################################
def UnfoldJetSpectrum(hRMDict, responsDict, jetSpectDict, hKinEff, \
        ptRangeDict, ptBinArrayDict, unfoldType, label, epBin, centBin):
    
    lUnfoldedDistsDict = {}

    response = responsDict['main']
    hMeasJetYieldPerBin = jetSpectDict['measYieldPerBin']
    # hMeasJetYieldPerBin = jetSpectDict['mcDetYieldPerBin'] # trivial test
    
    # Set some info based on specified unfolding type
    if "Bayes" in unfoldType: regParaName = "nIter"
    elif "SVD" in unfoldType: regParaName = "k"
    
    label =  lEPLabel[epBin]+'/Cent'+str(centBin)+'/Nomial/'
    hUFJetNomial, hCovNomial = UnfoldSpectrum1D(responsDict['main'], hMeasJetYieldPerBin, \
            regPara, unfoldType, 'Nominal')
    histName = 'hUnfoldedJetPt_'+unfoldType+'_Nominal'
    hUFJetNomial.SetName(histName)
    applyKinEff(hUFJetNomial, hKinEff)
    hUFJetNomial.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])
    hUFJetNomial.Scale(1., "width")
    lUnfoldedDistsDict['Nominal'] = hUFJetNomial
    
    
    # == s == regularization parameter difference Tests 222222222222222222222222
    hUnfKDevelopDict = {}
    hUnfoldRegParaKindDict = {}
    #errorType = RooUnfold.kErrors     #default in RooUnfold
    #errorType = RooUnfold.kCovariance #used for preliminary
    #recomended by Leticia, see hadron jet note 
    # https://alice-notes.web.cern.ch/system/files/notes/analysis/251/2017-Aug-11-analysis_note-HadronJet_analysis_note.pdf Section 5
    errorType = RooUnfold.kCovToy
    lHUFJetKDev = list()
    lHRFJetKDev = list()
    for i in range(1,regPara+3):
        label = 'iteK' + str(i)
        hUFJet, hCov = UnfoldSpectrum1D(responsDict['main'], hMeasJetYieldPerBin, \
            i, unfoldType, label)
        lHUFJetKDevForReFold = hUFJet.Clone(histName)
        applyKinEff(hUFJet, hKinEff)
        hUFJet.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])
        hUFJet.Scale(1., "width")
        lHUFJetKDev.append(hUFJet)
        histName = 'hReFold'+lEPLabel[epBin]+'_Ite'+str(i)+'_Cent'+str(centBin)
        
        hRFJet = responsDict['main'].ApplyToTruth(lHUFJetKDevForReFold)
        hRFJet.Scale(1., "width")
        lHRFJetKDev.append(hRFJet)
        
        # == s == Plots Unfolding QA figures 111111111111111111111111111111111111111
        label =  lEPLabel[epBin]+'/Cent'+str(centBin)+'/'+ unfoldType+'/Iteration'+str(i)+'/'
        if ("SVD" in unfoldType): UnfoldPlots.svdDVectorPlot(unfold, i, label)
        # if ("SVD" in unfoldType) and (i == regPara): 
            # UnfoldPlots.svdDVectorPlot(unfold, regPara, label)
        
        # Apply RM to unfolded result
        # UnfoldPlots.plotKinematicEfficiency(hRMDict['main'], hRMDict['uncut'], jetSpectDict, \
        #         ptRangeDict, ptBinArrayDict['mcGen'], label,oPerformanceDir)
        
        UnfoldPlots.plotCorrelationCoefficients(hCov, i, label,oPerformanceDir)

        print('Plot Fold Test ======')
        UnfoldPlots.plotFoldedTest(responsDict['main'], hUFJet, jetSpectDict['measYieldPerGeV'], \
            i, regParaName, unfoldType, ptRangeDict, ptBinArrayDict, label,oPerformanceDir)
        
        # Refolding test -- unfold Measured spectrum with response1, 
        # then apply response2 to unfolded result, and compare to the Measured spectrum.
        print('Plot ReFold Test ======')
        UnfoldPlots.plotRefoldingTest(responsDict['reFold'], responsDict['reFoldComp'], \
            hMeasJetYieldPerBin, jetSpectDict['measYieldPerGeV'], hKinEff,\
                i, regParaName, unfoldType, ptRangeDict, ptBinArrayDict, label,oPerformanceDir)
        
        # Unfolding test -- unfold the smeared Measured-level result with response, \
        # and compare to truth-level MC.
        print('Plot Closure Test ======')
        UnfoldPlots.plotClosureTest(responsDict['main'], \
            jetSpectDict['mcDetYieldPerBin'], jetSpectDict['mcGenJetUncutPerBin'], hKinEff,\
                i, regParaName, unfoldType, ptRangeDict, ptBinArrayDict, label,oPerformanceDir)
        
        # Plot unfolded results and Measured spectrum, to see the effect of unfolding
        # Divide by bin width to create per GeV spectrum. Cloning doesn't seem to work here.
        print('Plot Effect Of Unfolding Test ======')
        UnfoldPlots.plotEffectOfUnfolding(hUFJet, jetSpectDict['measYieldPerGeV'], \
            i, regParaName, unfoldType, ptRangeDict, ptBinArrayDict, label,oPerformanceDir)
        # == s == Plots Unfolding QA figures 111111111111111111111111111111111111111

        # Store the spectra near the optimal value of k, for later plotting the ratio
        hUnfoldRegParaKindNames = ["iteKM1", "Nominal", "iteKP1"]
        hUnfoldRegParaKindLineColor = [2, 4, 1]

        if (abs(i - regPara) == 1):
            histKindNum = i - regPara + 1
            tempHist = hUFJet.Clone()
            histName = 'hUnfoldedJetPt_'+unfoldType+'_' + hUnfoldRegParaKindNames[histKindNum]
            tempHist.SetName(histName)
            tempHist.SetLineColor(hUnfoldRegParaKindLineColor[histKindNum])
            hUnfoldRegParaKindDict[hUnfoldRegParaKindNames[histKindNum]] = tempHist
            lUnfoldedDistsDict[hUnfoldRegParaKindNames[histKindNum]] = tempHist
    
    label =  lEPLabel[epBin]+'/Cent'+str(centBin)+'/'+ unfoldType

    # kViolet+4 kViolet+1 kAzure-3 kAzure+10 kGreen+3 kSpring-8 kSpring+9 kOrange kOrange+7 kOrange+10
    lColor = [880+4, 880+1, 860-3, 860+10, 416+3,  820-8, 820+9, 800+7, 800+10]
    UnfoldPlots.plotUnfoldedKDevelopedSpectra(0, lHUFJetKDev, regPara, \
        ptRangeDict,ptBinArrayDict,label,lColor,oPerformanceDir)
    UnfoldPlots.plotUnfoldedKDevelopedSpectra(1, lHRFJetKDev, regPara, \
        ptRangeDict,ptBinArrayDict,label,lColor,oPerformanceDir)
    # == e == regularization parameter difference Tests 222222222222222222222222

    # == s == detector level pT range difference Tests 3333333333333333333333333
    hUFJetLcut, hCovLcut = UnfoldSpectrum1D(responsDict['DetLvLcut'], jetSpectDict['measYieldPerBin_LCut'], regPara, "Bayes", "DetLvLcut")
    histName = 'hUnfoldedJetPt_'+unfoldType+'_DetLvLcut'
    applyKinEff(hUFJetLcut, hKinEff)
    hUFJetLcut.SetName(histName)
    hUFJetUcut, hCovUcut = UnfoldSpectrum1D(responsDict['DetLvUcut'], jetSpectDict['measYieldPerBin_UCut'], regPara, "Bayes", "DetLvUcut")
    histName = 'hUnfoldedJetPt_'+unfoldType+'_DetLvUcut'
    hUFJetUcut.SetName(histName)
    applyKinEff(hUFJetUcut, hKinEff)
    # hUFJetLcut.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
    hUFJetLcut.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])
    hUFJetLcut.Scale(1., "width")
    lUnfoldedDistsDict['DetLvLcut'] = hUFJetLcut
    # hUFJetUcut.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
    hUFJetUcut.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])
    hUFJetUcut.Scale(1., "width")
    lUnfoldedDistsDict['DetLvUcut'] = hUFJetUcut
    # == s == detector level pT range difference Tests 3333333333333333333333333
    
    # == s == different prior Unfolding  44444444444444444444444444444444444444
    hUFJetDiffPrior, hCovDiffPrior = UnfoldSpectrum1D(responsDict['diffPrior'], \
        hMeasJetYieldPerBin, regPara, "Bayes", "DiffPrior")
    # hUFJetDiffPrior.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
    hUFJetDiffPrior.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])
    histName = 'hUnfoldedJetPt_'+unfoldType+'_DiffPrior'
    hUFJetDiffPrior.SetName(histName)
    applyKinEff(hUFJetDiffPrior, hKinEff)
    hUFJetDiffPrior.Scale(1., "width")
    lUnfoldedDistsDict['diffPrior'] = hUFJetDiffPrior
    # == e == different prior Unfolding  44444444444444444444444444444444444444

    return lUnfoldedDistsDict
################################################################################


################################################################################
# Unfold jet Spect, Bayes or SVD (specify type)                              ###
################################################################################
def UnfoldSpectrum1D(respons, jetSpect, iteNum, unfoldType, label):
    
    # Set some info based on specified unfolding type
    if "Bayes" in unfoldType: regParaName = "nIter"
    elif "SVD" in unfoldType: regParaName = "k"
    
    # Loop over values of regularization parameter
    # Bayes: number of iterations
    # SVD: k
    # errorType = RooUnfold.kErrors     #default in RooUnfold
    # errorType = RooUnfold.kCovariance #used for preliminary
    # errorType = RooUnfold.kCovToy     
    #recomended by Leticia, see hadron jet note 
    # https://alice-notes.web.cern.ch/system/files/notes/analysis/251/2017-Aug-11-analysis_note-HadronJet_analysis_note.pdf Section 5
    errorType = RooUnfold.kCovariance
    
    # Set up the SVD/Bayesian unfolding object
    if "Bayes" in unfoldType:
        unfold = RooUnfoldBayes(respons, jetSpect, iteNum)
        #unfold.SetNToys(1000)
    if "SVD" in unfoldType:
        unfold = RooUnfoldSvd(respons, jetSpect, iteNum)
        unfold.SetNToys(1000)
        # unfold.SetNToys(10)
        
        
    # Perform the unfolding
    # Produces the truth distribution, with errors, 
    # PerBin (will scale by bin width below, after refolding checks)
    hUnfoldedJetSpect = unfold.Hreco(errorType) 

    # 1 -- (default) sqrt of cov matrix diagonals 
    # (for SVD, uses toy MCs to account for response matrix errors)
    # 3 -- sqrt of cov matrix from toy MC tests
    hUnfoldedJetSpect.SetName("hUnfoldedJetSpect{}".format(label)+unfoldType)

    # Plot Pearson correlation coeffs for each k/nIter, 
    # to get a measure of the correlation between the bins
    covMatrix = unfold.Ereco(errorType) # Get the covariance matrix
    
    # Plot unfolded results and Measured spectrum, to see the effect of unfolding
    # Divide by bin width to create per GeV spectrum. Cloning doesn't seem to work here.
    # hUnfoldedJetSpect.Scale(1., "width")
    
    return hUnfoldedJetSpect, covMatrix
################################################################################


################################################################################
##+++  reweightPythia                                                        ###
################################################################################
def reweightPythia(jetSpectDict):
    hJetSpectTruPerBin.Scale(0.5) # To approximate RAA
    for bin in range(1, hJetSpectTruPerBin.GetNbinsX()):
        x = hJetSpectTruPerBin.GetBinCenter(bin)
        if x < 40:
            sf = 1/(1 + 0.01*x+0.5)
            hJetSpectTruPerBin.SetBinContent(bin, hJetSpectTruPerBin.GetBinContent(bin)*sf)
    hFoldedPythiaTruth = response.ApplyToTruth(hJetSpectTruPerBin)
    unfoldPlots.reweightTruthSpect(hFoldedPythiaTruth, hJetSpectMeasPerBin)
################################################################################

#######################################################################################
# Get errors from measured spectrum, rebinned to errors on the truth-level bining   ###
#######################################################################################
def getMeasuredErrors(hJetSpectMeasuredPerBin):
    dictErrors = {}
    for bin in range(1, hJetSpectMeasuredPerBin.GetNbinsX()+1):
        binVal = hJetSpectMeasuredPerBin.GetBinCenter(bin)
        content = hJetSpectMeasuredPerBin.GetBinContent(bin)
        error = hJetSpectMeasuredPerBin.GetBinError(bin)

        if content != 0: #????????kuma add ????????????
            dictErrors[binVal] = error/content

    return dictErrors
################################################################################

################################################################################
# Smear historgam for error testing                                         　##
################################################################################
def smearPoints(hist,fRandom):
    hnew = hist.Clone("hnew")

    for bin in range(1, hnew.GetNbinsX()):
        contgr  = hist.GetBinContent(bin)
        errg    = hist.GetBinError(bin)
        content = fRandom.Gaus(contgr, errg)
        hnew.SetBinContent(bin,content)
        #hnew.SetBinError(bin,0.)  #For Bayes OK
        hnew.SetBinError(bin,errg) #For SVd errors are needed

    return hnew
################################################################################

################################################################################
# Smear spectrum according to the error bars on the measured spectrum        ###
#################################################################################
def smearSpectrum(h, measuredErrors):
    h_smear = h.Clone()
    histName = h.GetName() + '_smear'
    h_smear.SetName(histName)

    # Zero errors in all bins
    for bin in range(1, h_smear.GetNbinsX() + 1):
        h_smear.SetBinError(bin, 0)

    # Loop through relevant bins and smear content according to errors in measured data, 
    # and set new errors
    r = ROOT.TRandom3(0)
    for binValue, error in measuredErrors.items():
        hBin = h_smear.FindBin(binValue)
        content = h_smear.GetBinContent(hBin)
        errorNew = error * content
        contentNew = content + r.Gaus(0, errorNew)
        h_smear.SetBinContent(hBin, contentNew)
        h_smear.SetBinError(hBin, errorNew)

    # If error was not set (because we ran out of statistics in data), 
    # set it to something manually
    enableSetManually = False
    if enableSetManually:
        for bin in range(1, h_smear.GetNbinsX() + 1):
            error = h_smear.GetBinError(bin)
            content = h_smear.GetBinContent(bin)
            if bin > 10 and error < 1e-15:
                h_smear.SetBinError(bin, 0.3*content)

    return h_smear
################################################################################

def applyKinEff(hUFJet, hKinEff):
    nBins = hUFJet.GetNbinsX()
    for ptBin in range(1,nBins+1):
        xPtVal = hUFJet.GetXaxis().GetBinCenter(ptBin)
        kinEffBin = hKinEff.FindBin(xPtVal)
        kinEffVal = hKinEff.GetBinContent(kinEffBin)
        ufYVal = hUFJet.GetBinContent(ptBin)
        corrUFYVal = ufYVal
        corrUFYVal = ufYVal*kinEffVal
        # corrUFYVal = ufYVal/kinEffVal
        
        hUFJet.SetBinContent(ptBin, corrUFYVal)
        ufYErr = hUFJet.GetBinError(ptBin)
        hUFJet.SetBinError(ptBin, ufYErr)

################################################################################
def OFileStracture(lMainTree, numOfCentBin):
    for epBin in range(0, 3):
        listName = lEPLabel[epBin]
        lEPList = ROOT.TList()
        lEPList.SetName(listName)
        for centBin in range(0, numOfCentBin):
            lEachCentHists = ROOT.TList()
            lEachCentHists.SetName('lCent{0}'.format(centBin))

            lUnfoldedJets = ROOT.TList()
            lUnfoldedJets.SetName('lUnfoldedJets')
            lEachCentHists.Add(lUnfoldedJets)

            lRMHists = ROOT.TList()
            lRMHists.SetName('lRMQA')
            lEachCentHists.Add(lRMHists)

            lJetDist = ROOT.TList()
            lJetDist.SetName('lJetDistQA')
            lEachCentHists.Add(lJetDist)

            lResponseHists = ROOT.TList()
            lResponseHists.SetName('lResponseQA')
            lEachCentHists.Add(lResponseHists)

            lEPList.Add(lEachCentHists)
            
        
        lMainTree.Add(lEPList)
################################################################################


if __name__ == "__main__":
    UnfoldJetPtDist()


