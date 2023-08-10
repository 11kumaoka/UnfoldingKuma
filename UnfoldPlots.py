#! /usr/bin/env python
# unfoldPlots


import argparse
import ctypes
import os
import gc
from array import array
import numpy as np
import math

import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F

# Load the RooUnfold library
#ROOT.gSystem.Load("$ALICE_WORK_DIR/osx_x86-64/RooUnfold/latest/lib/libRooUnfold.dylib")
# ROOT.gSystem.Load("$ALIBUILD_WORK_DIR/slc7_x86-64/RooUnfold/latest/lib/libRooUnfold.so")
ROOT.gSystem.Load("/persist/sw/slc7_x86-64/RooUnfold/latest/lib/libRooUnfold.so")

from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
from ROOT import RooUnfoldSvd
from ROOT import TCanvas
from ROOT import gROOT

import genePlotSets

# Set plotting options
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

fileFormat = '.root'


################################################################################
# Prepare to unfold jet spectrum                        ##
# Returns RooUnfoldResponse object.                     ##
################################################################################
def plotHistBeforUnfold(responsDict, jetSpectDict, label, outputDir):

    # Set up the RooUnfoldResponse object
    # One can supply the measured and truth distributions is to incorporate fakes and inefficiency
    # For the truth-level, we pass in the projection 
    # before the RM was cut to a specific pT-Measured range,
    # in order to account for the kinematic efficiency


    # Plot Pb-Pb data vs. pp+Pb-Pb matched Measured-level MC ("Measured-level Raa")
    hJetSpectMeasPerGeV = jetSpectDict['measPerGeV']
    hJetSpectMCDetPerGeV = jetSpectDict['mcGenPerGeV']

    xRangeMin = 0
    xRangeMax = 200
    yAxisTitle = "#frac{1}{T_{AA}}#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [GeV^{-1}]"
    legendTitle = ""
    h1LegendLabel = "Pb-Pb data"
    h2LegendLabel = "pp+Pb-Pb, det-level"
    ratioYAxisTitle = "R_{AA} det-level"

    outputRatioDir = outputDir + '/' + label + '/RatioPlot'
    if not os.path.exists(outputRatioDir):
        os.makedirs(outputRatioDir)
    outputFilename = os.path.join(outputRatioDir, "hDetLevelRaa" + fileFormat)
    #plotSpectra(hJetSpectrumMeasuredPerGeV, hJetSpectrumMCDetPerGeV, "", 1., \
    # xRangeMin, xRangeMax, yAxisTitle, ratioYAxisTitle, outputFilename, \
    # "", legendTitle, h1LegendLabel, h2LegendLabel)

    # Apply response matrix to truth-level MC distribution, 
    # and compare it to det-level MC distribution
    # (sanity check for exact agreement -- if kinematic efficiency is not used)
    applyRMToTruth = False
    if applyRMToTruth:
        # Since truth spectrum and response matrix are per-bin, so is the folded result
        hJetSpectMCTruPerBin = jetSpectDict['mcGenPerBin'].Clone()
        hJetSpectMCTruPerBin.SetName(jetSpectDict['mcGenPerBin'].GetName()+'_clone')
        hFoldedPerBin = response['main'].ApplyToTruth(hJetSpectMCTruPerBin) # TH1

        hFoldedPerGeV = hFoldedPerBin.Clone() # Scale by bin width to transform to per-GeV result
        hFoldedPerGeV.SetName("hFoldedPerGeV")
        hFoldedPerGeV.Scale(1., "width")

        yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [GeV^{-1}]"
        legendTitle = ""
        h1LegendLabel = "Folded pp truth"
        h2LegendLabel = "pp+Pb-Pb, det-level"
        ratioYAxisTitle = "Folded/Det-level"
        outputFilename = os.path.join(outputDir, "hJetSpectraFoldedMCTruth" + fileFormat)

        genePlotSets.plotSpectra(hFoldedPerGeV, hJetSpectrumMCDetPerGeV, "", 1., xRangeMin, \
            xRangeMax, yAxisTitle, ratioYAxisTitle, outputFilename, "", legendTitle, \
                h1LegendLabel, h2LegendLabel)
################################################################################


################################################################################
# Plot SVD D-vector plots                                                     ##
################################################################################
def svdDVectorPlot(unfoldSVD, regPara, label,outputDir):
    # Plot d-vector, as a function of k
    svdUnfoldObject = unfoldSVD.Impl() # Get TSVDUnfold object
    hDVector = svdUnfoldObject.GetD()
    hDVector.GetXaxis().SetTitle("k")
    hDVector.GetYaxis().SetTitle("d_{k}")
    hDUsed = hDVector.Clone("hDUsed")
    for bin in range(1, hDVector.GetNbinsX() + 1):
        if bin!=regPara: hDUsed.SetBinContent(bin,0)
    
    outputDirDVector = outputDir  + '/' + label + "/DVector/"
    if not os.path.exists(outputDirDVector):
        os.makedirs(outputDirDVector)
    outputFilename = os.path.join(outputDirDVector, "hDVector{}".format(fileFormat))
    genePlotSets.plotHist(hDVector, hDUsed, outputFilename, "", True, False, "", False, True)
################################################################################


################################################################################
# Plot JES shift distribution for various fixed pT-gen                        ##
################################################################################
def buildAngularityWeightedResponse(hist3D, label, factorHigh, factorLow, outputDir):

    hAgularvsPtTruth = hist3D.Project3D("zy")
    hAgularvsPtTruth.SetName("histhAgularVsPtTruth")

    #Do the reweighting in 5GeV slizes of pt_Truth
    #Determine mean if angularity in this pt_Truth slize
    histRebinned = hAgularvsPtTruth.Clone("rebinnedForProfile")
    histRebinned.RebinX(5)

    histAngularityProf = histRebinned.ProfileX()
    histAngularityProf.GetYaxis().SetTitle("mean angularity")

    meanAngulArray=[]
    #Fill array with lower bin edges of data histogram
    for bin in range(1, histAngularityProf.GetNbinsX()+1):
        mean    = histAngularityProf.GetBinContent(bin)
        meanBin = hAgularvsPtTruth.GetYaxis().FindBin(mean)
        meanAngulArray.append(meanBin)
    meanAngulArray2 = array('i',meanAngulArray)

    #Project the 3D to the RM_2D
    hResponseMatrixFineBinned = hist3D.Project3D("yx")

    histRebinnedEmptyLow  = hResponseMatrixFineBinned.Clone("rebinnedForAngularityLow")
    histRebinnedEmptyLow.Reset("ICESM")
    histRebinnedEmptyHigh = hResponseMatrixFineBinned.Clone("rebinnedForAngularityHigh")
    histRebinnedEmptyHigh.Reset("ICESM")
    histRebinnedSum = hResponseMatrixFineBinned.Clone("rebinnedForAngularitySum")
    histRebinnedSum.Reset("ICESM")

    aguBins=0
    #loop over ptTruth (1GeV binning), select angularity, project
    for yBin in range (1, hist3D.GetNbinsY()+1):
        arrayBin = yBin%5
        if yBin%5==0:
            aguBins=aguBins+1
        hist3D.GetYaxis().SetRangeUser(yBin,yBin)#set a single pT truth bin

        hist3D.GetZaxis().SetRange(-1,meanAngulArray2[arrayBin])#Set Angularity Range
        #hist3D.GetZaxis().SetRange(-1,hist3D.GetNbinsZ()+1)#Set complete Angularity Range for testing
        hResponseMatrixFineBinnedLowAngu = hist3D.Project3D("yxo")
        histRebinnedEmptyLow.Add(hResponseMatrixFineBinnedLowAngu)
        histRebinnedEmptyLow.SetName("LowAnguRM")
        #
        hist3D.GetZaxis().SetRange(meanAngulArray2[arrayBin]+1,hist3D.GetNbinsZ()+1) #Set Angularity Range
        hResponseMatrixFineBinnedHighAngu = hist3D.Project3D("yxo")
        histRebinnedEmptyHigh.Add(hResponseMatrixFineBinnedHighAngu)
        histRebinnedEmptyHigh.SetName("HighAnguRM")

    histRebinnedSum.Add(histRebinnedEmptyHigh,factorHigh)
    histRebinnedSum.Add(histRebinnedEmptyLow,factorLow)

    histRebinnedEmptyHigh.SetMinimum(1e-7)
    histRebinnedEmptyHigh.SetMaximum(1e9)
    outputFilename = os.path.join(outputDir, "hResponseMatrixHighAngul_{}.png".format(label))
    genePlotSets.plotHist(histRebinnedEmptyHigh, '', outputFilename, "colz", False, True)

    histRebinnedEmptyLow.SetMinimum(1e-7)
    histRebinnedEmptyLow.SetMaximum(1e9)
    outputFilename = os.path.join(outputDir, "hResponseMatrixLowAngul_{}.png".format(label))
    genePlotSets.plotHist(histRebinnedEmptyLow, '', outputFilename, "colz", False, True)
    
    return histRebinnedSum
################################################################################


################################################################################
# Re-weight the pythia truth spectrum to reduce tension with data at low pT  ###
################################################################################
def reweightTruthSpectrum(hFoldedPythiaTruth, hJetSpectMeasPerBin,outputDir):

    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    c.SetLogy()
    ROOT.gPad.SetLeftMargin(0.15)

    leg = ROOT.TLegend(0.6,0.6,0.88,0.83,"")
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    
    hFoldedPythiaTruth.GetXaxis().SetRangeUser(0, 250)
    hFoldedPythiaTruth.GetYaxis().SetRangeUser(2e-7,2e-2)
    
    hFoldedPythiaTruth.SetMarkerColor(2)
    hFoldedPythiaTruth.SetMarkerStyle(21)
    leg.AddEntry(hFoldedPythiaTruth, "Pythia", "P")
    hFoldedPythiaTruth.Draw("EP")
    
    hJetSpectMeasPerBin.SetMarkerColor(4)
    hJetSpectMeasPerBin.SetMarkerStyle(20)
    leg.AddEntry(hJetSpectMeasPerBin, "Pb-Pb data", "P")
    hJetSpectMeasPerBin.Draw("EP same")
    
    leg.Draw("same")

    outputFilename = os.path.join(outputDir, "hPythiaTruth{}".format(fileFormat))
    c.SaveAs(outputFilename)

    hRatio = hJetSpectMeasPerBin.Clone()
    hRatio.SetName("tgsgthyH")
    hRatio.Divide(hRatio, hFoldedPythiaTruth, 1., 1., "B")
    outputFilename = os.path.join(outputDir, "hPythiaTruthRatio{}".format(fileFormat))
    genePlotSets.plotHist(hRatio, outputFilename)
################################################################################


################################################################################
# Plot Pearson correlation coefficients from the covariance matrix           ###
# The correlation coefficent C_xy is related to Cov(x,y)                     ###
# by C_xy = Cov(x,y)/(sigma_x * sigma_y)                                     ###
################################################################################
def plotCorrelationCoefficients(covarianceMatrix, i, label,outputDir):
    nBinsX = covarianceMatrix.GetNrows()
    nBinsY = covarianceMatrix.GetNcols()
    
    histName = 'correlationCoefficientMatrix'
    correlationCoeffMatrix = ROOT.TH2D(histName, histName, nBinsX, 0, nBinsX, nBinsY, 0, nBinsY)

    for xbin in range(0, nBinsX):
        varianceX = covarianceMatrix(xbin, xbin)
        sigmaX = math.sqrt(varianceX)

        for ybin in range(0, nBinsY):
            varianceY = covarianceMatrix(ybin, ybin)
            sigmaY = math.sqrt(varianceY)

            covXY = covarianceMatrix(xbin, ybin)
            if sigmaX > 0 and sigmaY > 0:
                    Cxy = covXY / (sigmaX * sigmaY)
                    correlationCoeffMatrix.SetBinContent(xbin+1, ybin+1, Cxy)

                    #print "sigma x: {}, sigmay: {}".format(sigmaX, sigmaY)
                    #print "cov (x,y) = {}".format(covXY)
                    #print "Cxy = {}".format(Cxy)

    outputDirCovCoeff = outputDir+ '/' + label + "/CorrelationCoefficients/"
    if not os.path.exists(outputDirCovCoeff): os.makedirs(outputDirCovCoeff)
    outputFilename = os.path.join(outputDirCovCoeff, \
        "hCorrelationCoefficientMatrix{}{}".format(i, fileFormat))
    genePlotSets.plotHist(correlationCoeffMatrix, '', outputFilename, "colz")
################################################################################


#########################################################################################
# Plot kinematic efficiency                                                            ##
# (i.e. (pT-truth projection of response matrix with measured pT-det range selected)   ##
# / (pT-truth projection of response matrix with full pT-det range selected)           ##
#########################################################################################
def plotKinematicEfficiency(hRMWCut, hRMUncut, jetSpectDict, ptRangeDict, ptBinArrayDict, label,outputDir):
    hRMUncutForKinEff = hRMUncut.Clone('hRMUncutForKinEff')
    hRMWithCutForKinEff = hRMWCut.Clone('hRMWithCutForKinEff')
    
    hMCGenJetUncut = hRMUncutForKinEff.ProjectionY()
    histName = 'hMCGenJetUnCut'
    hMCGenJetUncut.SetName(histName)

    # Get the truth-level jet spectrum (matched) from response matrix (already re-binned)   
    # Do exclude under and overflow bins
    hMCGenJetWithCut = hRMWithCutForKinEff.ProjectionY()
    histName = 'hMCGenJetWithCut'
    hMCGenJetWithCut.SetName(histName)

    # Plot the ratio of the spectra
    hKinematicEfficiency = hMCGenJetWithCut.Clone()
    hKinematicEfficiency.SetName("hKinematicEfficiency")
    hKinematicEfficiency.Divide(hMCGenJetWithCut, hMCGenJetUncut, 1., 1., "B")
    
    hKinematicEfficiency.GetXaxis().SetTitle("#it{p}_{T}^{gen}")
    hKinematicEfficiency.GetYaxis().SetTitle("Kinematic Efficiency")
    hKinematicEfficiency.SetMarkerStyle(21)
    hKinematicEfficiency.SetMarkerColor(2)
    
    hinputRange = hKinematicEfficiency.Clone("inputRange")
    #set bins to zero outside the needed range
    for bin in range(1, hKinematicEfficiency.GetNcells() + 1):
        #if bin<hinputRange.FindBin(minPtReported) or bin>hinputRange.FindBin(maxPtReported)+1:
        if bin<hinputRange.FindBin(ptRangeDict['reported'][1]) \
            or bin>=hinputRange.FindBin(ptRangeDict['reported'][1]):
            hinputRange.SetBinContent(bin,0)
    text = "p_{T}^{det} #in [%d, %d] GeV" % (ptBinArrayDict['mcDet'][0], ptBinArrayDict['mcDet'][-1])
    
    outputDirKinEff = outputDir + '/' + label + "/KinEff/"
    if not os.path.exists(outputDirKinEff):
        os.makedirs(outputDirKinEff)
    outputFilename = os.path.join(outputDirKinEff, "hKinematicEfficiency_{}_{}{}"\
        .format(ptBinArrayDict['mcDet'][0], ptBinArrayDict['mcDet'][-1], fileFormat))
        
    genePlotSets.plotHist(hKinematicEfficiency, hinputRange, outputFilename, \
        "P E", False, False, text, False, True)

    return hKinematicEfficiency
################################################################################


#################################################################################
# Apply RM to unfolded result,                                               ###
# and check that I obtain measured spectrum (simple technical check)         ###
################################################################################
def plotFoldedTest(response, hJetSpectrumUnfoldedPerGeV, hJetSpectrumMeasuredPerGeV, \
    i, regularizationParamName, type, ptRangeDict, ptBinArrayDict, label,outputDir):
    histName = hJetSpectrumMeasuredPerGeV.GetName() + '_ForFoldTest'
    hJetSpectrumMeasuredPerGeVCP = hJetSpectrumMeasuredPerGeV.Clone()
    # nBins = len(ptDetBinArray) - 1
    # hJetSpectrumMeasuredPerGeVCP = hJetSpectrumMeasuredPerGeVCP.Rebin(nBins, histName, ptDetBinArray)
    # hJetSpectrumMeasuredPerGeVCP.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
    hJetSpectrumMeasuredPerGeVCP.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])
    
    # Produces folded distribution PerBin 
    # (unfolded spectrum is also PerBin at the moment, despite its name)
    hFoldedTruthPerGeV = response.ApplyToTruth(hJetSpectrumUnfoldedPerGeV) 
    hFoldedTruthPerGeV.Scale(1., "width") # Divide by bin width to create per GeV spectrum
    hFoldedTruthPerGeV.Scale(9.) # tuning ??????????????????????????????????????

    # nBins = len(ptDetBinArray) - 1
    # hFoldedTruthPerGeV = hFoldedTruthPerGeV.Rebin(nBins, histName, ptDetBinArray)
    # hFoldedTruthPerGeV.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
    hFoldedTruthPerGeV.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])

    yAxisTitle = "#frac{1}{T_{AA}}#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [mb (GeV/c)^{-1}]"
    legendTitle = ""
    h1LegendLabel = "Folded truth, {} = {}".format(regularizationParamName,i)
    h2LegendLabel = "Measured Pb-Pb"
    ratioYAxisTitle = "Folded truth / Measured Pb-Pb"
    outputDirFoldedPbPbTruth = outputDir + '/' + label + "/FoldedPbPbTruth/"
    if not os.path.exists(outputDirFoldedPbPbTruth): os.makedirs(outputDirFoldedPbPbTruth)
    outputFilename = os.path.join(outputDirFoldedPbPbTruth, \
        "hJetSpectraFoldedPbPbTruth{}_{}{}".format(type, i, fileFormat))
    genePlotSets.plotSpectra(hFoldedTruthPerGeV, hJetSpectrumMeasuredPerGeV, "", 1., \
        ptRangeDict['reported'][0], ptRangeDict['reported'][1], yAxisTitle, ratioYAxisTitle, outputFilename, \
            "", legendTitle, h1LegendLabel, h2LegendLabel)
################################################################################


################################################################################
# Refolding test                                                              ##
################################################################################
def plotRefoldingTest(response1, response2, hJetSpectrumMeasuredPerBin, hJetSpectrumMeasuredPerGeV, hKinEff,\
        i, regularizationParamName, type, ptRangeDict, ptBinArrayDict, label,outputDir):
    histName = hJetSpectrumMeasuredPerGeV.GetName() + '_ForReFoldTest'
    hJetSpectrumMeasuredPerGeVCP = hJetSpectrumMeasuredPerGeV.Clone()
    # nBins = len(ptDetBinArray) - 1
    # hJetSpectrumMeasuredPerGeVCP = hJetSpectrumMeasuredPerGeVCP.Rebin(nBins, histName, ptDetBinArray)
    # hJetSpectrumMeasuredPerGeVCP.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
    hJetSpectrumMeasuredPerGeVCP.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])
    
    # unfold measured spectrum with response1, then apply response2 to unfolded result.
    unfold1 = None
    if "SVD" in type: unfold1 = RooUnfoldSvd(response1, hJetSpectrumMeasuredPerBin, i)
    elif "Bayes" in type: unfold1 = RooUnfoldBayes(response1, hJetSpectrumMeasuredPerBin, i)
    

    # Produces the truth distribution, with errors, PerBin (will scale by bin width below)
    hJetSpectrumUnfolded1PerGeV = unfold1.Hreco() 
    # applyKinEff(hJetSpectrumUnfolded1PerGeV, hKinEff)
    # Produces folded distribution PerBin 
    # (unfolded spectrum is also PerBin at the moment, despite its name)
    hFoldedPbPbTruth1PerGeV = response2.ApplyToTruth(hJetSpectrumUnfolded1PerGeV) 
    
    # Then compare the refolded result to the measured spectrum
    hFoldedPbPbTruth1PerGeV.Scale(1., "width") # Divide by bin width to create per GeV spectrum
    # nBins = len(ptDetBinArray) - 1
    # hFoldedPbPbTruth1PerGeV = hFoldedPbPbTruth1PerGeV.Rebin(nBins, histName, ptDetBinArray)
    # hFoldedPbPbTruth1PerGeV.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
    hFoldedPbPbTruth1PerGeV.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])
    
    yAxisTitle = "#frac{1}{T_{AA}}#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [mb (GeV/c)^{-1}]"
    h1LegendLabel = "Refolded Pb-Pb, {} = {}".format(regularizationParamName,i)
    h2LegendLabel = "Measured Pb-Pb"
    legendTitle = ""

    ratioYAxisTitle = "#frac{Refolded}{Measured}"
    outputDirRefoldingTest = outputDir + '/' + label + "/RefoldingTest/"
    if not os.path.exists(outputDirRefoldingTest):
        os.makedirs(outputDirRefoldingTest)
    outputFilename = os.path.join(outputDirRefoldingTest, \
        "hJetSpectraRefoldingTest_{}{}".format( i, fileFormat))
    genePlotSets.plotSpectra(hFoldedPbPbTruth1PerGeV, hJetSpectrumMeasuredPerGeVCP, "", 1., \
        ptRangeDict['reported'][0], ptRangeDict['reported'][1], \
            yAxisTitle, ratioYAxisTitle, outputFilename, "", \
            legendTitle, h1LegendLabel, h2LegendLabel,"",2.2)
#################################################################################


################################################################################
# Closure test                                                                ##
################################################################################
def plotClosureTest(responseOri, hJetSpectrumMCDetPerBin, hJetSpectrumTruePerBin, hKinEff, \
    i, regularizationParamName, type, ptRangeDict, ptBinArrayDict, label,outputDir):

    histName = 'hMCGenYeild_ForClosureTest'
    hJetSpectrumTruePerGeVCP = hJetSpectrumTruePerBin.Clone(histName)
    # hJetSpectrumTruePerGeVCP.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
    hJetSpectrumTruePerGeVCP.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])
    hJetSpectrumTruePerGeVCP.Scale(1., "width") # Divide by bin width to create per GeV spectrum

    # Unfold smeared det-level spectrum with RM
    histName = responseOri.GetName() + '_ForClosureTest'
    response = responseOri.Clone(histName)
    unfold2 = None
    if "SVD" in type: unfold2 = RooUnfoldSvd(response, hJetSpectrumMCDetPerBin, i)
    elif "Bayes" in type: unfold2 = RooUnfoldBayes(response, hJetSpectrumMCDetPerBin, i)
    hJetSpectrumUnfolded2PerGeV = unfold2.Hreco()
    # applyKinEff(hJetSpectrumUnfolded2PerGeV, hKinEff)
    # nBins = len(ptDetBinArray) - 1
    histName = 'hUnfoldedMCDetYeild_ForClosureTest'
    # hJetSpectrumUnfolded2PerGeV = hJetSpectrumUnfolded2PerGeV.Rebin(nBins, histName, ptDetBinArray)
    # hJetSpectrumUnfolded2PerGeV.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
    hJetSpectrumUnfolded2PerGeV.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])
    hJetSpectrumUnfolded2PerGeV.Scale(1., "width") # Divide by bin width to create per GeV spectrum
    
    yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [GeV^{-1}]"
    legendTitle = ""
    h1LegendLabel = "Unfolded MC det, {} = {}".format(regularizationParamName,i)
    h2LegendLabel = "MC Truth"
    ratioYAxisTitle = "#frac{Unfolded MC det}{Truth}"
    outputDirClosureTest = outputDir + '/' + label + "/ClosureTest/"
    if not os.path.exists(outputDirClosureTest): os.makedirs(outputDirClosureTest)
    outputFilename = os.path.join(outputDirClosureTest, \
        "hJetSpectraClosureTest_{}{}".format(i, fileFormat))
    genePlotSets.plotSpectra(hJetSpectrumUnfolded2PerGeV, hJetSpectrumTruePerGeVCP, "", 1., \
        ptRangeDict['reported'][0], ptRangeDict['reported'][1], yAxisTitle, ratioYAxisTitle, outputFilename, \
            "", legendTitle, h1LegendLabel, h2LegendLabel)
################################################################################


####################################################################################
# Plot unfolded results and measured spectrum, to see the effect of unfolding   ####
####################################################################################
def plotEffectOfUnfolding(hUFJetsPerGeVOrig, hMeasJetsPerGeVOrig, \
    i, regParaName,type, ptRangeDict, ptBinArrayDict,label,outputDir):
    histName = hUFJetsPerGeVOrig.GetName() + '_ForUFEff'
    hUFJetsPerGeV = hUFJetsPerGeVOrig.Clone(histName)
    # nBins = len(ptDetBinArray) - 1
    # hUFJetsPerGeV = hUFJetsPerGeV.Rebin(nBins, histName, ptDetBinArray)
    # hUFJetsPerGeV.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
    hUFJetsPerGeV.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])
    histName = hMeasJetsPerGeVOrig.GetName() + '_ForUFEff'
    hMeasJetsPerGeV = hMeasJetsPerGeVOrig.Clone(histName)
    # nBins = len(ptDetBinArray) - 1
    # hMeasJetsPerGeV = hMeasJetsPerGeV.Rebin(nBins, histName, ptDetBinArray)
    # hMeasJetsPerGeV.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
    hMeasJetsPerGeV.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])
    
    yAxisTitle = "#frac{1}{T_{AA}}#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [mb (GeV/c)^{-1}]"
    legendTitle = ""
    h1LegendLabel = "Unfolded spectrum, {} = {}".format(regParaName,i)
    h2LegendLabel = "Measured spectrum"
    ratioYAxisTitle = "Unfolded/Measured"
    outputDirUnfoldedMeasuredRatio = outputDir + '/' + label + "/UnfoldedMeasuredRatio/"
    if not os.path.exists(outputDirUnfoldedMeasuredRatio):
        os.makedirs(outputDirUnfoldedMeasuredRatio)
    outputFilename = os.path.join(outputDirUnfoldedMeasuredRatio, \
        "hJetSpectraUnfoldedMeasuredRatio{}_{}{}".format(type, i, fileFormat))
    genePlotSets.plotSpectra(hUFJetsPerGeV, hMeasJetsPerGeV, \
        "", 1., ptRangeDict['reported'][0], ptRangeDict['reported'][1],\
            yAxisTitle, ratioYAxisTitle, outputFilename, \
                "", legendTitle, h1LegendLabel, h2LegendLabel,"",3.5)
####################################################################################


################################################################################
# Plot unfolded k developed result                                        　　###
################################################################################
def plotUnfoldedKDevelopedSpectra(FoldType, lhUnfKDevelop, regPara, ptRangeDict,\
    ptBinArrayDict, label,lColor, outputDir):

    cAll = ROOT.TCanvas("cAll","cAll: pT",800,850)
    cAll.cd()
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad1.SetTopMargin(0.05)
    pad1.SetLogy()
    pad1.Draw()
    pad1.cd()
    
    hMainResult = lhUnfKDevelop[regPara-1].Clone('hNominal_kIteration'+str(regPara))
    hMainResult.SetLineColor(1)
    hMainResult.SetLineWidth(2)
    hMainResult.SetLineStyle(1)
    hMainResult.GetYaxis().SetTitle("counts")
    hMainResult.GetYaxis().SetTitleSize(0.06)
    # hMainResult.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
    hMainResult.Rebin(len(ptBinArrayDict['reported'])-1, 'hNominal_kIteration'+str(regPara), ptBinArrayDict['reported'])
    hMainResult.GetYaxis().SetRangeUser(2e-10,2e-3)
    hMainResult.GetYaxis().SetLabelFont(43)
    hMainResult.GetYaxis().SetLabelSize(20)
    #hMainResult.GetXaxis().SetTitle("")
    hMainResult.Draw("hist E")

    TitleName = ''
    if FoldType == 0: TitleName = '{} Unfolding'.format(type)
    elif FoldType == 1: TitleName = '{} Refolding'.format(type)
    leg2 = ROOT.TLegend(0.7,0.6,0.88,0.93,TitleName)
    leg2.SetFillColor(10)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextSize(0.04)

    for kDevLoop in range(0, len(lhUnfKDevelop)):
        hUnfKDevelop = lhUnfKDevelop[kDevLoop]
        if kDevLoop == regPara-1: hUnfKDevelop.SetLineColor(1)
        else :hUnfKDevelop.SetLineColor(lColor[kDevLoop])
        hUnfKDevelop.SetLineWidth(3)
        hUnfKDevelop.DrawCopy("hist same E")
        leg2.AddEntry(hUnfKDevelop, "k={}".format(kDevLoop+1), "l")

    leg2.Draw("same")

    #- - - - - - - -
    cAll.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.15)
    pad2.SetRightMargin(0.05)
    pad2.Draw()
    pad2.cd()

    histName = ''
    if FoldType == 0: histName = 'hUnfoldedSpectraRatio_k=1/k='+str(regPara)
    elif FoldType == 1: histName = 'hRefoldedSpectraRatio_k=1/k='+str(regPara)
    hUnfKDeveRatio = lhUnfKDevelop[0].Clone(histName)
    hUnfKDeveRatio.Divide(hMainResult)
    hUnfKDeveRatio.GetYaxis().SetLabelSize(0.1)
    hUnfKDeveRatio.GetXaxis().SetLabelSize(0.1)
    hUnfKDeveRatio.GetYaxis().SetTitleOffset(0.4)
    hUnfKDeveRatio.GetXaxis().SetTitleOffset(0.9)
    hUnfKDeveRatio.GetYaxis().SetTitleSize(0.1)
    hUnfKDeveRatio.GetXaxis().SetTitleSize(0.15)
    hUnfKDeveRatio.GetYaxis().SetTitle("k=n/k={}".format(regPara))
    hUnfKDeveRatio.GetXaxis().SetTitle("#it{p}_{T,jet} (GeV/#it{c})")
    hUnfKDeveRatio.GetYaxis().SetNdivisions(505)
    # hUnfKDeveRatio.GetXaxis().SetRangeUser(ptRangeDict['reported'][0], ptRangeDict['reported'][1])
    hUnfKDeveRatio.Rebin(len(ptBinArrayDict['reported'])-1, histName, ptBinArrayDict['reported'])
    hUnfKDeveRatio.GetYaxis().SetRangeUser(0.85, 1.15)
    hUnfKDeveRatio.SetLineColor(lColor[0])
    hUnfKDeveRatio.Draw()

    # plot ratio h/h2
    hUnfKDeveRatio.DrawCopy("P E")
    for kDevLoop in range(1, len(lhUnfKDevelop)):
        histName = ''
        if FoldType == 0: histName = 'hUnfoldedSpectraRatio_k'+str(kDevLoop+1)\
            +'/k='+str(regPara)
        elif FoldType == 1: histName = 'hRefoldedSpectraRatio_k'+str(kDevLoop+1)\
            +'/k='+str(regPara)
        tempHUnfKDevelop = lhUnfKDevelop[kDevLoop].Clone(histName)
        if kDevLoop == regPara-1: tempHUnfKDevelop.SetLineColor(1)
        else :tempHUnfKDevelop.SetLineColor(lColor[kDevLoop])
        tempHUnfKDevelop.Divide(hMainResult)
        tempHUnfKDevelop.DrawCopy("same P E")

    outputRatioDir = outputDir + '/'  + label + '/RatioPlot'
    if not os.path.exists(outputRatioDir):
        os.makedirs(outputRatioDir)
    fileName = ''
    if FoldType == 0: fileName = 'hJetSpectraUnfoldedRatioAll'
    elif FoldType == 1: fileName = 'hJetSpectraRefoldedRatioAll'
    outputFilename = os.path.join(outputRatioDir, fileName + fileFormat)
    cAll.SaveAs(outputFilename)
    cAll.Close()
    ##------


def applyKinEff(hUFJet, hKinEff):
    nBins = hUFJet.GetNbinsX()
    for ptBin in range(1,nBins+1):
        xPtVal = hUFJet.GetXaxis().GetBinCenter(ptBin)
        kinEffBin = hKinEff.FindBin(xPtVal)
        kinEffVal = hKinEff.GetBinContent(kinEffBin)
        ufYVal = hUFJet.GetBinContent(ptBin)
        # corrUFYVal = ufYVal*kinEffVal
        corrUFYVal = ufYVal/kinEffVal
        
        hUFJet.SetBinContent(ptBin, corrUFYVal)
        ufYErr = hUFJet.GetBinError(ptBin)
        hUFJet.SetBinError(ptBin, ufYErr)


def plotRMRatio(lOMainTree, numOfCentBin):
    for centBin in range(0, numOfCentBin):
        histName = 'hRM_nomali' + '_Cent' + str(centBin)
        hRMOEP = lOMainTree[0][centBin][1].FindObject(histName)
        hRMIEP = lOMainTree[1][centBin][1].FindObject(histName)

        histName = 'hRMOIEPRatio' + '_Cent' + str(centBin)
        hRMOIEPRatio = hRMOEP.Clone(histName)
        hRMOIEPRatio.Divide(hRMIEP)

        lOMainTree[2].Add(hRMOIEPRatio)

