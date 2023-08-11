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

import unfoldPlots

ispp = False
isClosureTest = False
isClosure  = False

lJetR = [2, 3, 4, 5]
lPtCut = [0, 3, 5, 7, 10]
lCentDivKind = [[0, 1], [1, 2, 3], [3, 4, 5], [5, 6, 7], [7, 8, 9, 10]]

TAADiv10 = [23.26, 14.4, 8.767,  5.086, 2.747, 1.352,  0.5992,  0.2385, 0.08383, 0.02527]
# 0: ptMin, 1:ptMax
# ptRangeDict = {'mcGen':[0,250], 'mcDet':[0,250], 'meas':[0,250], 'reported':[0,250]}
ptRangeDict = {'mcGen':[5,250], 'mcDet':[20,120], 'meas':[0,250], 'reported':[0,250]}

### Main START  ################################################################
def unfoldingProcessForJ():
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;"); #..to supress a lot of standard output

    baseLabel = 'RMresuts'

    #Things to set for
    smearErrors    = True
    imposePrior    = False
    useRMflexibleR = False  
    matchingDistance = 0.6   #in % of jet radius

    regParam = 4
    powerLawOffset = 0

    # print('ptRangeMax = {}'.format(ptRangeDict['mcGen'][1]))
    ptBinArrayDict = eachJetPtBinDef(4, ptRangeDict)
    
    # inputBaseDir = '/home/alidock/OutPlot/LHC18r/pass3/'
    inputDir = '/home/alidock/JetAnalysisKuma/embedding/'
    inputFileName = '1hJetPerfomancePlotsFull_J.root'
    inputFile = ROOT.TFile(inputDir + inputFileName, "READ")

    mainTree = inputFile.Get('mainTree')

    # outputDir = '/home/alidock/cernbox/SWAN_projects/outputFiles/LHC15/pass3/UnfoldProcess/'
    # outputDir = '/home/alidock/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/UnfoldProcess/'
    outputDir = '/home/alidock/cernbox/SWAN_projects/outputFiles/LHC18qr/UnfoldProcess/'
    outputFileName = 'EachSpectrams.root'
    outputFile = ROOT.TFile(outputDir + outputFileName, "RECREATE")
    lOMainTree = ROOT.TList()

    for rBin in range(0, 1):
        label = ''
        label = baseLabel + '_R{}'.format(rBin)
        lOEachRJetTree = ROOT.TList()
        lOEachRJetTree.SetName('jetTree_R{}'.format(rBin))
        lOMainTree.Add(lOEachRJetTree)
        for centBin in range(0, 1):
            label = label + '_cent{}'.format(centBin)
            lOEachCentJetTree = ROOT.TList()
            lOEachCentJetTree.SetName('jetTree_R{}_Cent{}'.format(rBin, centBin))
            lOMainTree[rBin].Add(lOEachCentJetTree)

            nEventsData, nEventsResponse, vertexEfficiency = eventInfoExtract(mainTree, centBin)
            scalPraDict = {'nEventsData':nEventsData, 'nEventsResponse':nEventsResponse,\
                'vertexEfficiency':vertexEfficiency}

            hRMDict = getResponseMatrix(mainTree[rBin][centBin], \
                ptRangeDict, ptBinArrayDict, label)

            jetSpectDict = getJetSpectrums(mainTree[rBin][centBin], \
                hRMDict, ptRangeDict, ptBinArrayDict, scalPraDict, centBin, label)
            
            unfoldPlots.plotRmKinEff(hRMDict['main'], jetSpectDict['mcGenJetUncutPerBin'], \
                mainTree[rBin][centBin], label, 'bef')
            if imposePrior:
                setPrior(hRMDict['main'], jetSpectDict['mcGenJetUncutPerBin'], label)
            unfoldPlots.plotRmKinEff(hRMDict['main'], jetSpectDict['mcGenJetUncutPerBin'], \
                mainTree[rBin][centBin], label, 'aft')

            responsDict = getResponce(mainTree[rBin][centBin], \
                hRMDict, jetSpectDict['mcGenJetUncutPerBin'], label)

            unfoldPlots.plotBeforeProcess(responsDict['main'], jetSpectDict, label)
            unfoldPlots.plotKinematicEfficiency(jetSpectDict, ptRangeDict, label)

            unfoldType = 'Bayes'
            unfoldJetSpectrum(hRMDict, responsDict, jetSpectDict, ptRangeDict, ptBinArrayDict, \
                powerLawOffset, unfoldType, label, lOMainTree[rBin][centBin])
            unfoldType = 'SVD'
            # unfoldJetSpectrum(hRMDict, responsDict, jetSpectDict, ptRangeDict, ptBinArrayDict, \
                # powerLawOffset, unfoldType, label, lOMainTree[rBin][centBin])

            outputFile.cd()
            for key in hRMDict.values():
                lOMainTree[rBin][centBin].Add(key)
                print(key)
            for key in jetSpectDict.values():
                lOMainTree[rBin][centBin].Add(key)
                # print(key)
            for key in responsDict.values():
                lOMainTree[rBin][centBin].Add(key)
                # print(key)
            inputFile.cd()

    outputFile.cd()
    lOMainTree.Write('mainTree', 1)
    outputFile.Close()

    print(outputDir)
### Main END  ##################################################################


################################################################################
###      get Response Matrix                                                 ###
################################################################################
def getResponseMatrix(mainTree, ptRangeDict, ptBinArrayDict, label):
    hRM = mainTree[3] # No cut
    histName = 'RM' + label
    hRM.SetName(histName)

    hRMFineBin = hRM.Clone()
    histName = 'RM_uncut_finebinning' + label
    hRMFineBin.SetName(histName)
    hRMFineBin.GetXaxis().SetRangeUser(ptRangeDict['mcGen'][0], ptRangeDict['mcGen'][1])
    hRMFineBin.GetYaxis().SetRangeUser(0, ptRangeDict['mcGen'][1])
    hRMFineBin = normalizeResponseMatrix(hRMFineBin)

    hRMUncut = hRM.Clone()
    histName = 'RM_uncut' + label
    hRMUncut.SetName(histName)
    hRMUncut.GetXaxis().SetRangeUser(ptRangeDict['mcGen'][0], ptRangeDict['mcGen'][1])
    hRMUncut.GetYaxis().SetRangeUser(0, ptRangeDict['mcGen'][1])

    hRMForKinEff = hRM.Clone()
    histName = 'For_KinEff' + label
    hRMForKinEff.SetName(histName)
    hRMForKinEff.GetXaxis().SetRangeUser(ptRangeDict['mcDet'][0], ptRangeDict['mcDet'][1])
    hRMForKinEff.GetYaxis().SetRangeUser(0, ptRangeDict['mcGen'][1])

    hRMForKinEffUnCut = hRM.Clone()
    histName = 'For_KinEff_Uncut' + label
    hRMForKinEffUnCut.SetName(histName)
    hRMForKinEffUnCut.GetXaxis().SetRangeUser(-100, ptRangeDict['mcGen'][1])
    hRMForKinEffUnCut.GetYaxis().SetRangeUser(0, ptRangeDict['mcGen'][1])


    hRMMain = rebinResponseMatrix(hRM, ptBinArrayDict)
    histName = 'RM_Main' + label
    hRMMain.SetName(histName)
    hRMMain.GetXaxis().SetRangeUser(ptRangeDict['mcDet'][0], ptRangeDict['mcDet'][1])
    hRMMain.GetYaxis().SetRangeUser(ptRangeDict['mcGen'][0], ptRangeDict['mcGen'][1])
    # hRMMain = normalizeResponseMatrix(hRMMain) # Kuma ????????????????????????

    hRMReBin = hRMMain.Clone()
    histName = 'RM_ReBin' + label
    hRMReBin.SetName(histName)
    hRMReBin = normalizeResponseMatrix(hRMReBin)

    hRMNoKinEff = hRMMain.Clone()
    histName = 'RM_W/O_Kinematic_Efficiency' + label
    hRMNoKinEff.SetName(histName)
    
    # Get response matrix from response file (Measured, True) to be used for the unfolding, 
    # with pT-det range cut to desired range, and re-bin. 
    # Also get response matrices for refolding/closure tests.
    hRMReFold = hRMMain.Clone()
    histName = 'RM_for_Refolding' + label
    hRMReFold.SetName(histName)
    
    hRMReFoldComp = hRMMain.Clone()
    histName = 'RM_for_RefoldComp' + label
    hRMReFoldComp.SetName(histName)
    
    hRMDict = {'main':hRMMain, 'fineBin':hRMFineBin, 'reBin':hRMReBin, 'uncut':hRMUncut, \
        'forKinEff':hRMForKinEff, 'forKinEffUnCut':hRMForKinEffUnCut, 'noKinEff':hRMNoKinEff, \
            'reFold':hRMReFold, 'reFoldComp':hRMReFoldComp}
    
    return hRMDict
################################################################################


################################################################################
###      get JetSpectrums                                                    ###
################################################################################
def getJetSpectrums(mainTree, hRMDict, ptRangeDict, ptBinArrayDict, scalPraDict, centBin, label):
    hMeaseardJetPerBin = mainTree[2].Clone() # w/ background subtraction
    histName = 'hMeasuredJet' + label
    hMeaseardJetPerBin.SetName('hMeasuredJet')
    nBins = len(ptBinArrayDict['mcDet']) - 1
    binArray = ptBinArrayDict['mcDet']
    hMeaseardJetPerBin = hMeaseardJetPerBin.Rebin(nBins, histName, binArray)

    # hMeaseardJetPerBin.SetAxisRange(ptRangeDict['mcDet'][0], ptRangeDict['mcDet'][1],'X')

    #rebin the measured input spectrum to the final binning 
    # and save it for statistical error comparision
    hMeaseardJet_rebin = hMeaseardJetPerBin.Clone()
    histName = 'hMeasuredJet_rebin' + label
    hMeaseardJet_rebin.SetName(histName)
    nBins = len(ptBinArrayDict['mcGen']) - 1
    binArray = ptBinArrayDict['mcGen']
    hMeaseardJet_rebin = hMeaseardJet_rebin.Rebin(nBins, histName, binArray)

    # hMCGenJetPerBin = mainTree[4] # w/ matching
    # hMCDetJetPerBin = mainTree[5] # w/ matching

    # Get the truth-level jet spectrum (matched) from response matrix (already re-binned)
    # Do exclude under and overflow bins
    hMCGenJetPerBin = hRMDict['main'].ProjectionY()
    hMCGenJetPerBin.SetName('hMCGenJetPerBin' + label)

    hMCGenJetPerFineBin = hRMDict['main'].ProjectionY("_py",1,hRMDict['fineBin'].GetNbinsX())
    hMCGenJetPerFineBin.SetName('hMCGenJetPerFineBin' + label)

    hMCGenJetUncutPerBin = hRMDict['uncut'].ProjectionY()
    hMCGenJetUncutPerBin = hMCGenJetUncutPerBin.Rebin(\
        len(ptBinArrayDict['mcGen'])-1, "{}_NewBinning".format(hMCGenJetUncutPerBin.GetName()), \
            ptBinArrayDict['mcGen'])
    hMCGenJetUncutPerBin.SetName("hMCGenJetUncutPerBin" + label)

    hMCGenJetCutPerBinFoKinEff = hRMDict['forKinEff'].ProjectionY()
    hMCGenJetCutPerBinFoKinEff.SetName("hJetSpectrumTrueCutPerBin")
    hMCGenJetCutPerBinFoKinEff.Rebin(5)
    hMCGenJetCutPerBinFoKinEff.Scale(1., "with")

    hMCGenJetUnCutPerBinFoKinEff = hRMDict['forKinEffUnCut'].ProjectionY()
    hMCGenJetUnCutPerBinFoKinEff.SetName("hJetSpectrumTrueUnCutPerBin")
    hMCGenJetUnCutPerBinFoKinEff.Rebin(5)
    hMCGenJetUnCutPerBinFoKinEff.Scale(1., "with")

    # Get the det-level jet spectrum (matched) from response matrix projection, 
    # after cutting the pT-det range (already re-binned)
    # Note that this is potentially sensitive to the low-pT cutoff of the RM 
    # -- but it is only used for "det-level" Raa plot as a sanity check
    hMCDetJetPerBin = hRMDict['main'].ProjectionX()
    hMCDetJetPerBin.SetName('hMCDetJetPerBin' + label)


    jetSpectDict = {'measPerBin':hMeaseardJetPerBin, 'measPerBin_rebin':hMeaseardJet_rebin,\
            'mcGenPerBin':hMCGenJetPerBin, \
                'mcGenCutPerBinForKin':hMCGenJetCutPerBinFoKinEff, \
                    'mcGenUnCutPerBinForKin':hMCGenJetUnCutPerBinFoKinEff,\
                        'mcGenJetUncutPerBin':hMCGenJetUncutPerBin, \
                            'mcGenJetPerFineBin':hMCGenJetPerFineBin, \
                                'mcDetPerBin':hMCDetJetPerBin}

    hMCGenJetPerFineGeV = hRMDict['fineBin'].ProjectionY("_py",1,hRMDict['fineBin'].GetNbinsX())
    hMCGenJetPerFineGeV.SetName('hMCGenJetPerGeV' + label)
    hMCGenJetPerFineGeV.Scale(1., "width")
    jetSpectDict['mcGenJetPerFineGeV'] =  hMCGenJetPerFineGeV

    prepareScaleSpectra(mainTree, jetSpectDict, scalPraDict, centBin)

    # For Measured spectra, divide by bin width to produce per GeV spectra
    hMeaseardJetPerGeV = jetSpectDict['measPerBin'].Clone()
    histName = "hMeaseardJetPerGeV" + label
    hMeaseardJetPerGeV.SetName(histName)
    hMeaseardJetPerGeV.Scale(1., "width")
    jetSpectDict['measPerGeV'] =  hMeaseardJetPerGeV

    hMeaseardJetPerGeV_rebin = hMeaseardJetPerGeV.Clone()
    histName = "hMeaseardJetPerGeV_rebin" + label
    nBins = len(ptBinArrayDict['mcGen']) - 1
    binArray = ptBinArrayDict['mcGen']
    hMeaseardJetPerGeV_rebin = hMeaseardJetPerGeV_rebin.Rebin(nBins, histName, binArray)
    hMeaseardJetPerGeV_rebin.Scale(1., "width")
    
    jetSpectDict['measPerGeV_rebin'] =  hMeaseardJetPerGeV_rebin

    hMCGenJetPerGeV = jetSpectDict['mcGenPerBin'].Clone()
    hMCGenJetPerGeV.SetName('hMCGenJetPerGeV' + label)
    hMCGenJetPerGeV.Scale(1., "width")
    jetSpectDict['mcGenPerGeV'] =  hMCGenJetPerGeV

    hMCDetJetPerGeV = hMCDetJetPerBin.Clone()
    hMCDetJetPerGeV.SetName('hMCDetJetPerGeV' + label)
    hMCDetJetPerGeV.Scale(1., "width")
    jetSpectDict['mcDetPerGeV'] =  hMCDetJetPerGeV
    
    return jetSpectDict
################################################################################


################################################################################
## retrun Response matrix                                                   ####
################################################################################
def getResponce(mainTree, hRMDict, hJetMCGenUncutPerBin, label):

    histName = hRMDict['main'].GetName() + '_Response' + label
    response = RooUnfoldResponse(0, hJetMCGenUncutPerBin, hRMDict['main'], histName, histName)
    response.UseOverflow(False) # ??????????????????????????????????????????????

    histName = 'hRMUncut' + label
    respUncut = RooUnfoldResponse(0, 0, hRMDict['uncut'], histName, histName)
    histName = 'hRMNoKinEff' + label
    respNoKinEff = RooUnfoldResponse(0, 0, hRMDict['noKinEff'], histName, histName)
    histName = 'hRMReFold' + label
    respReFold = RooUnfoldResponse(0, 0, hRMDict['reFold'], histName, histName)
    histName = 'hRMReFoldComp' + label
    respReFoldComp = RooUnfoldResponse(0, 0, hRMDict['reFoldComp'], histName, histName)

    responsDict = {'main':response, 'uncut':respUncut,'noKinEff':respNoKinEff, \
        'reFold':respReFold, 'reFoldComp':respReFoldComp}

    return responsDict
################################################################################

################################################################################
# Prepare to unfold jet spectrum                                              ##
# Returns RooUnfoldResponse object.                                           ##
################################################################################
def prepareScaleSpectra(mainTree, jetSpectDict, scalPraDict, centBin):

    # Normalization (after re-binning has been done)
    if jetSpectDict['measPerBin']:
        if ispp and not isClosure: # Normalize the data spectrum by the visible MB cross-section
            #(mb) V0AND cross section (value taken from https://cds.cern.ch/record/2648933)
            visibleMBCrossSection = 50.87 
            jetSpectDict['measPerBin'].Scale(visibleMBCrossSection)
            jetSpectDict['measPerBin'].Scale(scalPraDict['vertexEff'])
            print("Scale the pp spectrum by visibleMBCrossSection ({}) and vertex efficiency ({})"\
                .format(visibleMBCrossSection,scalPraDict['vertexEff']))
        elif isClosure:
            # Normalize by avg Nevents per pT hard bin, below (due to scaling script convention)
            pass 
        else: 
            pass
            # Normalize the data spectrum by Nevents, and by TAA
            # Taa = 23.07        # mb^-1 (value taken from http://cds.cern.ch/record/2636623)
            # jetSpectDict['measPerBin'].Scale(1./Taa)
            # print("Scale the PbPb spectrum by TAA ({})".format(Taa))

        jetSpectDict['measPerBin'].Scale(1./scalPraDict['nEventsData'])
        print("Scale the Data spectrum by nEvents ({})".format(scalPraDict['nEventsData']))

    # Normalize the truth spectrum by avg Nevents per pT hard bin (in the embedded data), 
    # to form the pp cross-section
    # jetSpectDict['mcGenPerBin'].Scale(1./scalPraDict['nEventsResponse'])
    print("Scale MC Particle Level spectrum by nEventsResponse ({})"\
        .format(scalPraDict['nEventsResponse']))
    # jetSpectDict['mcGenJetPerFineBin'].Scale(1./scalPraDict['nEventsResponse'])
    # jetSpectDict['mcGenJetPerFineGeV'].Scale(1./scalPraDict['nEventsResponse'])

    # Normalize MC det spectrum by avg Nevents per pT hard bin (in the embedded data)
    # jetSpectDict['mcDetPerBin'].Scale(1./scalPraDict['nEventsResponse'])
    print("Scale MC Detector Level spectrum by nEventsResponse ({})"\
        .format(scalPraDict['nEventsResponse']))


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
        #truthBinCenter  = hTruthProjectionBefore.GetXaxis().GetBinCenter(truthBin)
        #this is a cut projection
        #truthProjVal    = hTruthProjectionBefore.GetBinContent(truthBin)            
        truthBinCenter  = hJetSpectTruUncutPerBin.GetXaxis().GetBinCenter(truthBin)
        truthProjVal    = hJetSpectTruUncutPerBin.GetBinContent(truthBin)
        priorScalingVal = 1
        #ispp=False #if you want to use the smearing for pp
        # For pp take as an alternative 1) the POWEG spectrum as prior 2) 
        # the unfolded final result as a prior
        if ispp:
            if hAlterPriorSpect:
                alterProjVal = hAlterPriorSpectSF.GetBinContent(\
                    hAlterPriorSpectSF.FindBin(truthBinCenter))
                if truthProjVal>0:
                    priorScalingVal  = alterProjVal / truthProjVal
                else:
                    priorScalingVal  = 0
            else:
                print("Error: in this setting there should be a file provided \
                    which contains an alternative prior distribution")
        # For PbPb, scale the power law exponent
        else:
            priorScalingVal = math.pow(truthBinCenter, powerLawOffset)

        # Set uncut truth spectrum to prior (to preserve the kinematic efficiency)
        uncutContent = hJetSpectTruUncutPerBin.GetBinContent(truthBin)
        hJetSpectTruUncutPerBin.SetBinContent(truthBin, uncutContent * priorScalingVal)
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



################################################################################
################################################################################
def eventInfoExtract(mainTree, centBin):
    # Get the number of accepted events in data in the specified centrality range
    # evtCutList = mainTree.FindObject("EventCutOutput")
    if ispp:
        hCutStats  = mainTree.FindObject("fCutStatsData")
        nEventsData= hCutStats.GetBinContent(hCutStats.GetXaxis().FindBin("All cuts")) #Bin AllCuts
        print("N data events accepted: {}"\
            .format(lCentDivKind[centBin][0], lCentDivKind[centBin][-1], nEventsData))
    else:
        # hRhoVsCent = jetHistogramList.FindObject("hRhoVsCent") #James
        # hRhoVsCent.GetXaxis().SetRangeUser(minCent, maxCent)
        # histNEventData = hRhoVsCent.ProjectionX()
        # nEventsData = histNEventData.Integral()
        hCentEvents = mainTree.FindObject("Centrality_rawData")
        hCentEventsClone = hCentEvents.Clone()
        hCentEventsClone.SetName("Centrality_raw"+'_Cent{}'.format(centBin))
        hCentEventsClone.GetXaxis().SetRangeUser(\
            lCentDivKind[centBin][0]*10, lCentDivKind[centBin][-1]*10)
        nEventsData = hCentEventsClone.Integral()
        print("N data events (in {}0-{}0% centrality): {}"\
            .format(lCentDivKind[centBin][0], lCentDivKind[centBin][-1], nEventsData))

    # Get the average number of accepted events per pT-hard bin (in the embedded data)
    # This is used to scale the response spectra, due to convention in scaling the pT hard bins
    histNEventResponse = mainTree.FindObject("hNEventsAcc")
    if not histNEventResponse:
        print("couldn't find NEventResponse :(")
    PtHardBins = 20
    nEventsAccSum = 0.
    for bin in range(0,PtHardBins):
        nEventsAccSum += histNEventResponse.GetBinContent(bin+1)

    nEventsResponse = nEventsAccSum/PtHardBins
    if isClosureTest:
        nEventsData = nEventsResponse

    hNormalizationHist = mainTree.FindObject("fNormalisationHistData")
    #allEvents          = hNormalizationHist.GetBinContent(\
    # hNormalizationHist.GetXaxis().FindBin("No cuts"))
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
# Rebin the response matrix to have variable binning                          ##
################################################################################
def rebinResponseMatrix(hResponseMatrix, ptBinArrayDict):

    histname = "{}_Rebin".format(hResponseMatrix.GetName())
    title = histname + ";#it{p}_{T,corr}^{det} (GeV/#it{c});#it{p}_{T}^{truth} (GeV/#it{c})"
    nBinGen = len(ptBinArrayDict['mcGen']) - 1
    binArrayGen = ptBinArrayDict['mcGen']
    nBinDet = len(ptBinArrayDict['mcDet']) - 1
    binArrayDet = ptBinArrayDict['mcDet']
    hResponseMatrix_rebin = ROOT.TH2D(histname, title, \
        nBinsDet, binArrayDet, nBinsGen, binArrayGen)
    
    # Loop over all bins in fine-binned response matrix, 
    # and fill appropriate bin in new response matrix
    # Assume that the bin edges overlap appropriately
    for ibin in range(1, hResponseMatrix.GetNbinsX() + 1):
        for jbin in range(1, hResponseMatrix.GetNbinsY() + 1):

            oldContent = hResponseMatrix.GetBinContent(ibin, jbin)

            # Find the bin that should be filled in the new histogram, and fill it
            # Need to get (x,y) location from bins (ibin, jbin)
            x = hResponseMatrix.GetXaxis().GetBinCenter(ibin)
            y = hResponseMatrix.GetYaxis().GetBinCenter(jbin)
            hResponseMatrix_rebin.Fill(x, y, oldContent)
    
    # Assume 0 errors on response matrix
    for bin in range(1, hResponseMatrix_rebin.GetNcells() + 1):
        hResponseMatrix_rebin.SetBinError(bin, 0)

    return hResponseMatrix_rebin
################################################################################

################################################################################
# Unfold jet Spect, Bayes or SVD (specify type)                              ###
################################################################################
def unfoldJetSpectrum(hRMDict, responsDict, jetSpectDict, ptRangeDict, ptBinArrayDict, \
                powerLawOffset, unfoldType, label, lOTree):

    response = responsDict['main']
    regPara = 3
    # Set some info based on specified unfolding type
    if "Bayes" in unfoldType: regParaName = "nIter"
    elif "SVD" in unfoldType: regParaName = "k"


    # Smear MC truth spectrum according to the error bars on the Measured spectrum
    MeasuredErrors = getMeasuredErrors(jetSpectDict['measPerBin'])

    # Store errors on Measured spectrum, 
    # for closure test (stored as dictionary {bin:error} in same binning)
    hJetSpectMCDetPerBinSmear = smearSpectrum(jetSpectDict['mcDetPerBin'], MeasuredErrors)

    # reweightPythia(jetSpectDict)

    #- - - - - - - - - - - - - - - - - - - - - - - - - -
    # change the Measured spectrum within the statistical errors
    # Do that 1000 times and unfold
    RegulForStatTest= 3
    errorType       = RooUnfold.kCovariance #used for prel.

    fRandom = ROOT.TRandom3(0) #James
    for i in range(1, 10):  #1000
        # fRandom = ROOT.TRandom3(i)
        hJetSpectMeasPerBin_T = smearPoints(jetSpectDict['measPerBin'],fRandom)
        histName = "hJetSpectMeasuredPerBin{}Round{}".format(label,i) + unfoldType
        hJetSpectMeasPerBin_T.SetName(histName)
    
        # Set up the unfolding object
        if "Bayes" in unfoldType: unfold = RooUnfoldBayes(responsDict['main'], \
            hJetSpectMeasPerBin_T, RegulForStatTest)
        elif "SVD" in unfoldType: unfold = RooUnfoldSvd(responsDict['main'], \
            hJetSpectMeasPerBin_T, RegulForStatTest)

        # Perform the unfolding
        hunf = unfold.Hreco(errorType) # Produces the truth distribution, with errors, PerBin
        histName = "hJetSpectUnfoldedPerGeV{}Round{}".format(label,i) + unfoldType
        hJetSpectUnfoldedPerGeV_T = hunf.Clone(histName)

        # Plot unfolded results and Measured spectrum, to see the effect of unfolding
        # Divide by bin width to create per GeV spectrum.
        hJetSpectUnfoldedPerGeV_T.Scale(1., "width") 

    
    # Loop over values of regularization parameter
    # Bayes: number of iterations
    # SVD: k
    hUnfKDevelopDict = {}
    hUnfoldRegParaKindDict = {}
    for i in range(1,regPara+3):
        #errorType = RooUnfold.kErrors     #default in RooUnfold
        #errorType = RooUnfold.kCovariance #used for preliminary
        #recomended by Leticia, see hadron jet note 
        # https://alice-notes.web.cern.ch/system/files/notes/analysis/251/2017-Aug-11-analysis_note-HadronJet_analysis_note.pdf Section 5
        errorType = RooUnfold.kCovToy     

        # Set up the SVD/Bayesian unfolding object
        if "Bayes" in unfoldType:
            unfold = RooUnfoldBayes(responsDict['main'], jetSpectDict['measPerBin'], i)
            #unfold.SetNToys(1000)
        if "SVD" in unfoldType:
            unfold = RooUnfoldSvd(responsDict['main'], jetSpectDict['measPerBin'], i)
            # unfold.SetNToys(1000)
            unfold.SetNToys(10)

        # Perform the unfolding
        # Produces the truth distribution, with errors, 
        # PerBin (will scale by bin width below, after refolding checks)
        hJetSpectUnfoldedPerGeV = unfold.Hreco(errorType) 

        # 1 -- (default) sqrt of cov matrix diagonals 
        # (for SVD, uses toy MCs to account for response matrix errors)
        # 3 -- sqrt of cov matrix from toy MC tests
        hJetSpectUnfoldedPerGeV.SetName("hJetSpectUnfoldedPerGeV{}".format(label)+unfoldType)
        jetSpectDict['unfoldPerBin'] = hJetSpectUnfoldedPerGeV

        # Plot Pearson correlation coeffs for each k/nIter, 
        # to get a measure of the correlation between the bins
        covMatrix = unfold.Ereco(errorType) # Get the covariance matrix
        unfoldPlots.plotCorrelationCoefficients(covMatrix, i, label)
        
        if ("SVD" in unfoldType) and (i is 1): 
            unfoldPlots.svdDVectorPlot(unfold, regPara, label)
        
        # Apply RM to unfolded result, 
        # and check that I obtain Measured spectrum (simple technical check)
        unfoldPlots.plotResultFolded(ispp, responsDict['main'], \
            hJetSpectUnfoldedPerGeV, jetSpectDict['measPerGeV'], \
            i, regParaName, unfoldType, label)
        
        
        # Refolding test -- unfold Measured spectrum with response1, 
        # then apply response2 to unfolded result, and compare to the Measured spectrum.
        unfoldPlots.plotRefoldingTest(responsDict['reFold'], responsDict['reFoldComp'], \
            jetSpectDict['measPerBin'], jetSpectDict['measPerGeV'], \
                i, regParaName, unfoldType, label)
        """
        # Unfolding test -- unfold the smeared Measured-level result with response, \
        # and compare to truth-level MC.
        unfoldPlots.plotClosureTest(responsDict['noKinEff'], \
            jetSpectDict['mcDetPerBin'], jetSpectDict['mcGenPerBin'], \
            i, regParaName, unfoldType ,label)
        """
        # Plot unfolded results and Measured spectrum, to see the effect of unfolding
        # Divide by bin width to create per GeV spectrum. Cloning doesn't seem to work here.
        hJetSpectUnfoldedPerGeV.Scale(1., "width") 

        unfoldPlots.plotEffectOfUnfolding(\
            hJetSpectUnfoldedPerGeV, jetSpectDict['measPerGeV_rebin'], \
            i, regParaName, unfoldType, label)
        
        #Save all the different spectra for unfolding QA
        hUnfKDevelopDict['hUnfoldedSpectra_k{}'.format(i)]   \
            = hJetSpectUnfoldedPerGeV.Clone("hUnfoldedSpectra_k{}".format(i))

        # Store the spectra near the optimal value of k, for later plotting the ratio
        hUnfoldRegParaKindNames = ["hHigherkResult", "hMainResult", "hLowerkResult"]
        hUnfoldRegParaKindLineColor = [2, 4, 1]
        
        
        if (abs(i - regPara) < 2) or (regPara==2):
            histKindNum = i - regPara + 1
            print('histKindNum = {}, histName = '.format(histKindNum) + hUnfoldRegParaKindNames[histKindNum])
            # This is for the special case of having a very low regularization parameter
            # in this case take the two higher k values as variations
            if (regPara==2) and (i is regPara+2):
                histKindNum = 1
            tempHist = hJetSpectUnfoldedPerGeV.Clone()
            tempHist.SetName(hUnfoldRegParaKindNames[histKindNum])
            tempHist.SetLineColor(hUnfoldRegParaKindLineColor[histKindNum])
            tempHist.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
            hUnfoldRegParaKindDict[hUnfoldRegParaKindNames[histKindNum]] = tempHist

    unfoldPlots.plotUnfoldedKDevelopedSpectra(hUnfKDevelopDict, hUnfoldRegParaKindDict, \
        regPara, ptRangeDict, unfoldType, label, lOTree)
    

################################################################################



################################################################################
# Rebin the response matrix to have variable binning                          ##
################################################################################
def rebinResponseMatrix(hRM, ptBinArrayDict):
    histname = "{}NewRebinned".format(hRM.GetName())
    title = histname + ";#it{p}_{T,corr}^{det} (GeV/#it{c});#it{p}_{T}^{truth} (GeV/#it{c})"
    genNBins = len(ptBinArrayDict['mcGen']) - 1
    genBinArray = ptBinArrayDict['mcGen']
    detNBins = len(ptBinArrayDict['mcDet']) - 1
    detBinArray = ptBinArrayDict['mcDet']
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
        #print "Setting val {} to have relative error {}".format(binValue, errorNew/contentNew)

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


################################################################################
# Smear spectrum according to the error bars on the measured spectrum        ###
#################################################################################
def returnEventStatus():
    # Get the number of accepted events in data in the specified centrality range
    evtCutList = outputListData.FindObject("EventCutOutput")
    if ispp:
        hCutStats  = evtCutList.FindObject("fCutStats")
        nEventsData= hCutStats.GetBinContent(hCutStats.GetXaxis().FindBin("All cuts")) #Bin AllCuts
        print("N data events accepted: {}".format(minCent, maxCent, nEventsData))
    else:
        hRhoVsCent = jetHistogramList.FindObject("hRhoVsCent")
        hRhoVsCent.GetXaxis().SetRangeUser(minCent, maxCent)
        histNEventData = hRhoVsCent.ProjectionX()
        nEventsData = histNEventData.Integral()
        print("N data events (in {}-{}% centrality): {}".format(minCent, maxCent, nEventsData))

    # Get the average number of accepted events per pT-hard bin (in the embedded data)
    # This is used to scale the response spectra, due to convention in scaling the pT hard bins
    histNEventResponse = fResponse.Get("hNEventsAcc")
    if not histNEventResponse:
        print("couldn't find NEventResponse :(")
    PtHardBins = 20
    nEventsAccSum = 0.
    for bin in range(0,PtHardBins):
        nEventsAcc = histNEventResponse.GetBinContent(bin+1)
        nEventsAccSum += nEventsAcc
    nEventsAccAvg = nEventsAccSum/PtHardBins
    nEventsResponse = nEventsAccAvg
    if isClosureTest:
        nEventsData = nEventsResponse

    hNormalizationHist = evtCutList.FindObject("fNormalisationHist")
    #allEvents          = hNormalizationHist.GetBinContent(hNormalizationHist.GetXaxis().FindBin("No cuts"))
    #take only events that fulfill general selection criteria mostly pile-up 
    # as a reference for the vertex eff. calculation
    allEvents          = hNormalizationHist.GetBinContent(\
        hNormalizationHist.GetXaxis().FindBin("Event selection"))
    EvtsWithGoodVertex = hNormalizationHist.GetBinContent(\
        hNormalizationHist.GetXaxis().FindBin("Vertex reconstruction and quality"))
    vertexEfficiency   = (round(EvtsWithGoodVertex/allEvents,3)) # do not take more than 3 numbers after decimal

    print("N response events (avg per pT-hard bin): %d" % nEventsResponse)

    return nEventsResponse, allEvents, EvtsWithGoodVertex, vertexEfficiency
################################################################################


################################################################################
def eachJetPtBinDef(kind, ptRangeDict):
    binArrayGen = None
    binArrayDet = None

    if kind is 1:
        binArrayGen = ([ptRangeDict['mcGen'][0], 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75,\
            80, 85, 90, 95, 100, 105, 110, ptRangeDict['mcGen'][1]])
        binArrayDet = ([ptRangeDict['mcDet'][0], 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75,\
            80, 85, 90, 95, 100, 105, 110, ptRangeDict['mcDet'][1]])
    elif kind is 2:
        binArrayGen = ([ptRangeDict['mcGen'][0], 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75,\
            80, 85, 90, 95, 100, 105, 110, ptRangeDict['mcGen'][1]])
        binArrayDet = ([ptRangeDict['mcDet'][0], 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75,\
            80, 85, 90, 95, 100, 105, 110, ptRangeDict['mcDet'][1]])
    elif kind is 3:
        binArrayGen  = ([5, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 140, 190, 250])
        binArrayDet = ([0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,\
            60.0,65.0,70.0,75.0,80.0,85.0,95.0,105.0,125.0,250.0])
    elif kind is 4:
        binArrayGen = ([5, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 140, 190, 250])
        binArrayDet = ([20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 120])
    else:
        binArrayGen = ([0, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75,\
            80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150, 160, 170, 180, 250])
        binArrayDet = ([0, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75,\
            80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150, 160, 170, 180, 250])

    binArrayGen = array('d',binArrayGen)
    binArrayDet = array('d',binArrayDet)
    ptBinArrayDict = {'mcGen':binArrayGen, 'mcDet':binArrayDet} 

    return ptBinArrayDict
################################################################################


if __name__ == "__main__":
    unfoldingProcessForJ()



