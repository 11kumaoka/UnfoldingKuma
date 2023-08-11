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
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TMath

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

import UnfoldPlots
import genePlotSets

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

ispp = False
isClosureTest = False
isClosure  = False

lJetR = [2, 3, 4, 5]
lPtCut = [0, 3, 5, 7, 10]
# lCentDivKind = [[0, 1], [1, 2, 3], [3, 4, 5], [5, 6, 7], [7, 8, 9, 10]]
lCentDivKind = [[0,5],[5,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80]]
lCentRange = [[0,5],[5,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80],[80,90]]

TAADiv10 = [23.26, 14.4, 8.767,  5.086, 2.747, 1.352,  0.5992,  0.2385, 0.08383, 0.02527]
lPsi2Reso = [0.4513, 0.6388, 0.6801, 0.6463, 0.5609, 0.4455, 0.2860, 0.1400, 0.0477]

minPtGen = 30.
maxPtGen = 200.
minPtDet = 40. 
maxPtDet = 150.
# 0: ptMin, 1:ptMax
ptRangeDict = {'mcGen':[30,250], 'mcDet':[40,120], 'meas':[40,120], 'reported':[30,150]}
# binArrayGen = ([30, 40, 50, 60, 70, 80, 100, 120, 140, 190, 250])
# binArrayGen = ([30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130])
# binArrayDet = ([40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 120])
binArrayGen = ([30, 50, 70, 90, 110, 130])
binArrayDet = ([40, 50, 60, 70, 80, 90, 100, 110, 120])
binArrayGen = array('d',binArrayGen)
binArrayDet = array('d',binArrayDet)

jetRadius = 0.2
matchingDistance = 0.6

lEPLabel = ['OutOfPlane', 'InPlane']

errorType = RooUnfold.kCovariance

### Main START  ################################################################
def testUnfold():
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;"); #..to supress a lot of standard output

    baseLabel = 'RMresuts'

    #Things to set for
    useRMflexibleR = False  
    
    regParam = 4
    powerLawOffset = 0
    
    
    inputDir = '../'
    inputFileName = 'EmbedPtHardScaledResults.root'
    inputFile = ROOT.TFile(inputDir + inputFileName, "READ")

    mainTree = inputFile.Get('mainTree')

    outputDir = './'
    outputFileName = 'UnfoldedPtDists.root'
    outputFile = ROOT.TFile(outputDir + outputFileName, "RECREATE")
    lOMainTree = ROOT.TList()
    
    label = ''
    label = baseLabel
    for centBin in range(0, 5):
        lOEachCentJetTree = ROOT.TList()
        lOEachCentJetTree.SetName('jetTree_Cent{}'.format(centBin))
        lOMainTree.Add(lOEachCentJetTree)
    

    lHOPlaneJetMeasPerBin = list()
    lHOPlaneJetTrueUncutPerBin = list()
    lHOPlaneJetUnfoldedPerGeV = list()
    lHIPlaneJetMeasPerBin = list()
    lHIPlaneJetTrueUncutPerBin = list()
    lHIPlaneJetUnfoldedPerGeV = list()
    lJetV2 = list()

    epLabel = lEPLabel[0]
    IORMMainTree = mainTree.FindObject(epLabel)
    for centBin in range(0, 5):
        hOPlaneJetMeasPerBin, hOPlaneJetTrueUncutPerBin, hOPlaneJetUnfoldedPerGeV = getSpectrums(centBin, epLabel, mainTree, IORMMainTree)
        outputFile.cd()
        lOMainTree[centBin].Add(hOPlaneJetMeasPerBin)
        lOMainTree[centBin].Add(hOPlaneJetTrueUncutPerBin)
        lOMainTree[centBin].Add(hOPlaneJetUnfoldedPerGeV)
        lHOPlaneJetMeasPerBin.append(hOPlaneJetMeasPerBin)
        lHOPlaneJetTrueUncutPerBin.append(hOPlaneJetTrueUncutPerBin)
        lHOPlaneJetUnfoldedPerGeV.append(hOPlaneJetUnfoldedPerGeV)

    inputFile.cd()
    epLabel = lEPLabel[1]
    IORMMainTree = mainTree.FindObject(epLabel)
    for centBin in range(0, 5):
        hIPlaneJetMeasPerBin, hIPlaneJetTrueUncutPerBin, hIPlaneJetUnfoldedPerGeV = getSpectrums(centBin, epLabel, mainTree, IORMMainTree)

        outputFile.cd()
        lOMainTree[centBin].Add(hIPlaneJetMeasPerBin)
        lOMainTree[centBin].Add(hIPlaneJetTrueUncutPerBin)
        lOMainTree[centBin].Add(hIPlaneJetUnfoldedPerGeV)
        lHIPlaneJetMeasPerBin.append(hIPlaneJetMeasPerBin)
        lHIPlaneJetTrueUncutPerBin.append(hIPlaneJetTrueUncutPerBin)
        lHIPlaneJetUnfoldedPerGeV.append(hIPlaneJetUnfoldedPerGeV)

    for centBin in range(0, 5):
        hJetV2N = lHIPlaneJetUnfoldedPerGeV[centBin].Clone()
        histname = "hJetV2Cent"+str(centBin)
        hJetV2N.SetName(histname)
        hJetV2D = lHIPlaneJetUnfoldedPerGeV[centBin].Clone()
        histname = "hIJetDV2Cent"+str(centBin)
        hJetV2D.SetName(histname)

        hOJetNV2 = lHOPlaneJetUnfoldedPerGeV[centBin].Clone()
        histname = "hOJetNV2Cent"+str(centBin)
        hOJetNV2.SetName(histname)
        hOJetDV2 = lHOPlaneJetUnfoldedPerGeV[centBin].Clone()
        histname = "hOJetDV2Cent"+str(centBin)
        hOJetDV2.SetName(histname)
        
        hOJetNV2.Scale(-1)
        hJetV2N.Add(hOJetNV2)
        hJetV2D.Add(hOJetDV2)
        hJetV2N.Divide(hJetV2D)
        hJetV2N.Scale(math.pi/4)
        hJetV2N.Scale(1/lPsi2Reso[centBin])

        lJetV2.append(hJetV2N)

        outputFile.cd()
        lOMainTree[centBin].Add(hJetV2N)
    
    
    outputFile.cd()
    lOMainTree.Write('mainTree', 1)
    outputFile.Close()

    print(outputDir)
### Main END  ##################################################################

################################################################################
# Prepare to unfold jet spectrum                                              ##
# Returns RooUnfoldResponse object.                                           ##
################################################################################
def prepareScaleSpectra(mainTree, jetSpectDict, scalPraDict, centBin):

    # Normalization (after re-binning has been done)
    if jetSpectDict['measPerBin']:
        # Normalize the data spectrum by Nevents, and by TAA
        # Taa = 23.07        # mb^-1 (value taken from http://cds.cern.ch/record/2636623)
        # jetSpectDict['measPerBin'].Scale(1./Taa)
        # print("Scale the PbPb spectrum by TAA ({})".format(Taa))

        jetSpectDict['measPerBin'].Scale(1./scalPraDict['nEventsData'])
        print("Scale the Data spectrum by nEvents ({})".format(scalPraDict['nEventsData']))

    # Normalize the truth spectrum by avg Nevents per pT hard bin (in the embedded data), 
    # to form the pp cross-section
    print("Scale MC Particle Level spectrum by nEventsResponse ({})"\
        .format(scalPraDict['nEventsResponse']))
    # jetSpectDict['mcGenJetPerFineBin'].Scale(1./scalPraDict['nEventsResponse'])

    # Normalize MC det spectrum by avg Nevents per pT hard bin (in the embedded data)
    # jetSpectDict['mcDetPerBin'].Scale(1./scalPraDict['nEventsResponse'])
    print("Scale MC Detector Level spectrum by nEventsResponse ({})"\
        .format(scalPraDict['nEventsResponse']))
################################################################################


################################################################################
def eventInfoExtract(mainTree, centBin):
    # Get the number of accepted events in data in the specified centrality range
    # hCentEvents = mainTree.FindObject("Centrality_rawData")
    hCentEvents = mainTree.FindObject("Centrality_raw")
    hCentEventsClone = hCentEvents.Clone()
    hCentEventsClone.SetName("Centrality_raw"+'_Cent{}'.format(centBin))
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
    
    # Assume 0 errors on response matrix
    for bin in range(1, hRM_rebin.GetNcells() + 1):
        hRM_rebin.SetBinError(bin, 0)

    return hRM_rebin
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


def getSpectrums(centBin, epLabel, mainTree, IORMMainTree):
    label = '_cent{}'.format(centBin)
    lOEachCentJetTree = ROOT.TList()
    lOEachCentJetTree.SetName('jetTree_Cent{}'.format(centBin))
    # lOMainTree.Add(lOEachCentJetTree)

    nEventsData, nEventsResponse, vertexEfficiency = eventInfoExtract(mainTree, centBin)
    scalPraDict = {'nEventsData':nEventsData, 'nEventsResponse':nEventsResponse,\
        'vertexEfficiency':vertexEfficiency}

    treeName = 'eachCentHists_Cent'+ str(centBin)
    mainCentTree = IORMMainTree.FindObject(treeName)
    hOrigRawJetName = 'hJetCorrPtLocal_' + str(centBin)
    hOrigMeaseardJetPerBin = mainCentTree.FindObject(hOrigRawJetName)
    hOrigMeaseardJetPerBin.GetXaxis().SetRangeUser(minPtDet, maxPtDet)
    histname =  hOrigMeaseardJetPerBin.GetName() +epLabel+ '_rebin'+'Cent'+str(centBin)
    nBinsDet = len(binArrayDet) - 1
    hJetSpectrumMeasuredPerBin = hOrigMeaseardJetPerBin.Rebin(nBinsDet, histname, binArrayDet)
    
    hJetSpectrumMeasuredPerGeV = hJetSpectrumMeasuredPerBin.Clone()
    histname = 'hJetSpectrumMeasuredPerGeV' +epLabel+'Cent'+str(centBin)
    hJetSpectrumMeasuredPerGeV.SetName(histname)
    hJetSpectrumMeasuredPerGeV.Scale(1., "width")

    baseMatchedJetHists = IORMMainTree.FindObject('hResponseMatrixDiff')
    
    hMatchingDistance0 = baseMatchedJetHists.Clone()
    histName = 'RMHists0'+epLabel+'Cent'+str(centBin)
    hMatchingDistance0.SetName(histName)
    hMatchingDistance0.GetAxis(0).SetRangeUser(minPtGen, maxPtGen)
    hMatchingDistance0.GetAxis(1).SetRangeUser(minPtDet, maxPtDet)
    hMatchingDistance0.GetAxis(2).SetRangeUser(-1,jetRadius*matchingDistance)
    hMatchingDistance0.GetAxis(5).SetRangeUser(lCentRange[centBin][0], lCentRange[centBin][1])
    hRMFineBinned = hMatchingDistance0.Projection(0, 1, "")
    histname = 'hRMFineBinned'+epLabel+'Cent'+str(centBin)
    hRMFineBinned.SetName(histName)

    nBinsDet = len(binArrayDet) - 1
    nBinsTruth = len(binArrayGen) - 1
    hRM = ROOT.TH2D(histname, histname, nBinsDet, binArrayDet, nBinsTruth, binArrayGen)
    for ibin in range(1, hRMFineBinned.GetNbinsX() + 1):
        for jbin in range(1, hRMFineBinned.GetNbinsY() + 1):
            oldContent = hRMFineBinned.GetBinContent(ibin, jbin)
            x = hRMFineBinned.GetXaxis().GetBinCenter(ibin)
            y = hRMFineBinned.GetYaxis().GetBinCenter(jbin)
            hRM.Fill(x, y, oldContent)
    
    for bin in range(1, hRM.GetNcells() + 1):
        hRM.SetBinError(bin, 0)
    
    
    hMatchingDistance1 = baseMatchedJetHists.Clone()
    histName = 'RMHists1'+epLabel+'Cent'+str(centBin)
    hMatchingDistance1.SetName(histName)
    hMatchingDistance1.GetAxis(0).SetRangeUser(0, maxPtGen)
    # hMatchingDistance1.GetAxis(1).SetRangeUser(minPtDet, maxPtDet)
    hMatchingDistance1.GetAxis(2).SetRangeUser(-1,jetRadius*matchingDistance)
    hMatchingDistance1.GetAxis(5).SetRangeUser(lCentRange[centBin][0], lCentRange[centBin][1])
    hRMUncut = hMatchingDistance1.Projection(0, 1, "")
    
    hJetSpectrumTrueUncutPerBin = hRMUncut.ProjectionY()
    histName = 'hJetGenUncutHistsRebin'+epLabel+'Cent'+str(centBin)
    hJetSpectrumTrueUncutPerBin = hJetSpectrumTrueUncutPerBin.Rebin(len(binArrayGen)-1, histName, binArrayGen)

    histname = "hRMMain" +epLabel+'Cent'+str(centBin)
    response = RooUnfoldResponse(0, hJetSpectrumTrueUncutPerBin, hRM, histname, histname)
    response.UseOverflow(False)
    unfoldBayes = RooUnfoldBayes(response, hJetSpectrumMeasuredPerBin, 4)
    hJetSpectrumUnfoldedPerGeV = unfoldBayes.Hreco(errorType)

    label = '_Cent' + str(centBin)
    plotResultFolded(nEventsData, nEventsResponse, response, \
        hJetSpectrumTrueUncutPerBin, hJetSpectrumMeasuredPerGeV, 4, "nIter", label)
    :
    hRM1 = hRM.Clone()
    hisname = 'hRM1' +epLabel+'Cent'+str(centBin)
    hRM1.SetName(hisname)
    response1 = RooUnfoldResponse(0, 0, hRM1, "hResponseMatrix1", "hResponseMatrix1")
    hRM2 = hRM.Clone()
    hisname = 'hRM2' +epLabel+'Cent'+str(centBin)
    hRM2.SetName(hisname)
    response2 = RooUnfoldResponse(0, 0, hRM2, "hResponseMatrix2", "hResponseMatrix2")
    plotRefoldingTest(nEventsData, response1, response2, 
        hJetSpectrumMeasuredPerBin, hJetSpectrumMeasuredPerGeV, 4, "nIter", label)

    return  hJetSpectrumMeasuredPerBin, hJetSpectrumTrueUncutPerBin, hJetSpectrumUnfoldedPerGeV


#################################################################################
# Apply RM to unfolded result,                                               ###
# and check that I obtain measured spectrum (simple technical check)         ###
################################################################################
def plotResultFolded(nEventsData, nEventsResponse,response, 
    hJetSpectrumTrueUncutPerBin, hJetSpectrumMeasuredPerGeV, i, regularizationParamName, label):
    # Produces folded distribution PerBin 
    # (unfolded spectrum is also PerBin at the moment, despite its name)
    hJetMeasuredPerGeVCP = hJetSpectrumMeasuredPerGeV.Clone()
    histname = hJetMeasuredPerGeVCP.GetName() + '_clone'
    hJetMeasuredPerGeVCP.SetName(histname)
    hJetMeasuredPerGeVCP.Scale(1/nEventsData)
    hJetMeasuredPerGeVCP.Scale(1/1000) # ??????????????????
    hFoldedTruthPerBin = response.ApplyToTruth(hJetSpectrumTrueUncutPerBin)
    hFoldedTruthPerGeV = hFoldedTruthPerBin.Clone()
    hFoldedTruthPerGeV.SetName("hFoldedPerGeV")
    hFoldedTruthPerGeV.Scale(1., "width")
    xRangeMin = 0
    xRangeMax = 250
    
    yAxisTitle = "#frac{1}{T_{AA}}#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [mb (GeV/c)^{-1}]"
    legendTitle = ""
    h1LegendLabel = "Folded truth, {} = {}".format(regularizationParamName,i)
    h2LegendLabel = "Measured Pb-Pb"
    ratioYAxisTitle = "Folded truth / Measured"
    outputDirFoldedPbPbTruth = "./FoldedPbPbTruth/"
    if not os.path.exists(outputDirFoldedPbPbTruth):
        os.makedirs(outputDirFoldedPbPbTruth)
    histName = 'hJetSpectraFoldedPbPbTruthCent'+label+'.root'
    outputFilename = os.path.join(outputDirFoldedPbPbTruth, histName)
    genePlotSets.plotSpectra(hFoldedTruthPerGeV, hJetMeasuredPerGeVCP, "", 1., \
        xRangeMin, xRangeMax, yAxisTitle, ratioYAxisTitle, outputFilename, \
            "", legendTitle, h1LegendLabel, h2LegendLabel)
################################################################################

################################################################################
# Refolding test                                                              ##
################################################################################
def plotRefoldingTest(nEventsData, response1, response2, \
    hJetSpectrumMeasuredPerBin, hJetSpectrumMeasuredPerGeV, \
        i, regularizationParamName, label):
    hJetMeasuredPerGeVCP = hJetSpectrumMeasuredPerGeV.Clone()
    histname = hJetMeasuredPerGeVCP.GetName() + '_cloneForRefo'
    hJetMeasuredPerGeVCP.SetName(histname)
    hJetMeasuredPerGeVCP.Scale(1/nEventsData)

    # unfold measured spectrum with response1, then apply response2 to unfolded result.
    unfold1 = RooUnfoldBayes(response1, hJetSpectrumMeasuredPerBin, i)
    # Produces the truth distribution, with errors, PerBin (will scale by bin width below)
    hJetSpectrumUnfolded1PerGeV = unfold1.Hreco() 
    # Produces folded distribution PerBin 
    # (unfolded spectrum is also PerBin at the moment, despite its name)
    hFoldedPbPbTruth1PerGeV = response2.ApplyToTruth(hJetSpectrumUnfolded1PerGeV)
    hFoldedPbPbTruth1PerGeV.Scale(1/nEventsData)
    hFoldedPbPbTruth1PerGeV.Scale(1., "width") # Divide by bin width to create per GeV spectrum
    xRangeMin = 10
    xRangeMax = 250
    yAxisTitle = "#frac{1}{T_{AA}}#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [mb (GeV/c)^{-1}]"
    h1LegendLabel = "Refolded Pb-Pb, {} = {}".format(regularizationParamName,i)
    h2LegendLabel = "Measured Pb-Pb"
    legendTitle = ""

    ratioYAxisTitle = "Refolded / Measured"
    outputDirRefoldingTest = "./RefoldingTest/"
    if not os.path.exists(outputDirRefoldingTest):
        os.makedirs(outputDirRefoldingTest)
    histName = 'hJetSpectraRefoldedCent'+label+'.root'
    outputFilename = os.path.join(outputDirRefoldingTest, histName)
    genePlotSets.plotSpectra(hFoldedPbPbTruth1PerGeV, hJetMeasuredPerGeVCP, "", 1., \
        xRangeMin, xRangeMax, yAxisTitle, ratioYAxisTitle, outputFilename, "", \
            legendTitle, h1LegendLabel, h2LegendLabel,"",2.2)
#################################################################################


if __name__ == "__main__":
    testUnfold()


