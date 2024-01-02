#! /usr/bin/env python
import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F
import argparse
import ctypes
import gc
from array import array
import numpy as np

import histSetting
import PtRangeList

################################################################################
# Plot JER                                                                    ##
################################################################################
def getJESshiftProf(hJESshift2D, ptBinArray, name):
    nBin = len(ptBinArray)-1
    # histJER = ROOT.TProfile('temp'+name, 'temp'+name, nBin, ptBinArray)
    histJER = ROOT.TH1D('temp'+name, 'temp'+name, nBin, ptBinArray)
    histName = hJESshift2D.GetName() + '_reBin'
    
    for iBin in range(1, histJER.GetNbinsX()):
        histName = "hTempProjJES_pTbin" + str(iBin)
        lPtVal = ptBinArray[iBin-1]
        uPtVal = ptBinArray[iBin]
        hJESshift2D_Rebin = hJESshift2D.Clone(histName)
        hTempProjJES = hJESshift2D_Rebin.GetXaxis().SetRangeUser(lPtVal, uPtVal)
        hTempProjJES = hJESshift2D_Rebin.ProjectionY()
        # JER = hTempProjJES.GetRMS()
        JER = hTempProjJES.GetStdDev()
        JERError = hTempProjJES.GetRMSError()
        # histJER.Fill(lPtVal, JER)
        histJER.SetBinContent(iBin, JER)
        histJER.SetBinError(iBin, JERError)
    
    return histJER


################################################################################
# Plot JER                                                                    ##
################################################################################
def getJER1(hJESshift2D, ptBinArray, name):
    nBin = len(ptBinArray)-1
    # histJER = ROOT.TProfile('temp'+name, 'temp'+name, nBin, ptBinArray)
    histJER = ROOT.TH1D('temp'+name, 'temp'+name, nBin, ptBinArray)
    histName = hJESshift2D.GetName() + '_reBin'
    
    for iBin in range(1, histJER.GetNbinsX()):
        histName = "hTempProjJES_pTbin" + str(iBin)
        lPtVal = ptBinArray[iBin-1]
        uPtVal = ptBinArray[iBin]
        hJESshift2D_Rebin = hJESshift2D.Clone(histName)
        hTempProjJES = hJESshift2D_Rebin.GetXaxis().SetRangeUser(lPtVal, uPtVal)
        hTempProjJES = hJESshift2D_Rebin.ProjectionY()
        # JER = hTempProjJES.GetRMS()
        JER = hTempProjJES.GetStdDev()
        JERError = hTempProjJES.GetRMSError()
        # histJER.Fill(lPtVal, JER)
        histJER.SetBinContent(iBin, JER)
        histJER.SetBinError(iBin, JERError)
    
    return histJER


def getJER2(histResponseMatrix, name):
    # For each pT^gen, compute the standard deviation of the pT^det distribution
    # Get the pT^gen profile, with errors as standard deviation of pT^det distribution
    histPtGenProf = histResponseMatrix.ProfileY("histPtGenProf" + name, 1, -1, "s")
    histPtGenProf.SetName(name)
    
    # Create histo to be used to fill JER values
    histJER = ROOT.TProfile(name, name, 100, 0., 200.) # same binning for pT^gen as in task
    
    # Loop through the bins, and fill the JER
    for bin in range(1,histPtGenProf.GetNbinsX()+1):
        
        sigma = histPtGenProf.GetBinError(bin)
        pTgen = histPtGenProf.GetXaxis().GetBinCenter(bin)
        JER = sigma/pTgen
        histJER.Fill(pTgen, JER)
    
        return histJER


################################################################################
# Plot Jet Reconstruction Efficiency                                          ##
################################################################################
def jetRecoEffPlot(lMainTree):
    lGenJetHists = list()
    lMatchedGenJetHists = list()
    lGenJetEffBaseHists = list()
    lGenJetEffBaseHists.append(lMainTree[0][1].Clone())
    for centBin in range(0,5):
        lGenEachCentEffBaseHists = list()
        lGenEachCentEffBaseHists.append(lMainTree[centBin][1].Clone())
        lGenEachCentEffBaseHists.append(lMainTree[centBin][3].Clone())
        lGenJetEffBaseHists.append(lMainTree[centBin][3].Clone())

        # hMatchGenJetPt.Scale(1/hMatchGenJetPt.GetEntries())
        # hMatchDetJetPt.Scale(1/hMatchDetJetPt.GetEntries())

        # effCanvas = lMainTree.FindObject('jetRecoEff_Cent'+str(centBin))
        effCanvas = lMainTree[centBin][0]
        ratioHistoYTitle = 'Jet reco efficiency'
        hXRange = [0, 250]
        hYRange = [0, 50000000]
        histSetting.ratioMergeCanvas(lMainTree, lGenEachCentEffBaseHists, \
            effCanvas, ratioHistoYTitle, hXRange, hYRange, centBin)

        hEachCentJetRecoEff = lMainTree[centBin][3].Clone()
        hEachCentJetRecoEff.SetName('jetRecoEffMatch_cent'+str(centBin))
        hEachCentJetRecoEff.Divide(lMainTree[centBin][3], lMainTree[centBin][1] ,1.,1.,"B")
        lMainTree[centBin].Add(hEachCentJetRecoEff)


####################################################################################################
# Plot kinematic efficiency                                                                       ##
# (i.e. (pT-truth projection of response matrix with measured pT-det range selected)              ##
# / (pT-truth projection of response matrix with full pT-det range selected)                      ##
####################################################################################################
def plotKinematicEfficiency1(useRMflexibleR, jetRadius, matchingDistance, \
    outputListResponse, minCent, maxCent, minPtDet, maxPtDet, minPtGen, maxPtGen, \
        minPtReported, maxPtReported, outputDir, fileFormat):
    
    # Get fine-binned response matrix (Measured, True),
    #  with pT-det range cut to desired range, and project truth distribution
    hResponseMatrixCut = getResponseMatrix(ispp, useRMflexibleR, jetRadius, \
        matchingDistance, outputListResponse, minCent, maxCent, \
            minPtDet, maxPtDet, 0, maxPtGen, 0, 0, 0, 0, "Cut",outputDir)
    hJetSpectrumTrueCutPerBin = hResponseMatrixCut.ProjectionY()
    hJetSpectrumTrueCutPerBin.SetName("hJetSpectrumTrueCutPerBin")
    rebinVal = 5
    hJetSpectrumTrueCutPerBin.Rebin(5)
    hJetSpectrumTrueCutPerBin.Scale(1., "width")
    
    # Get fine-binned response matrix (Measured, True), with full pT-det range,
    #  and project truth distribution
    fMinPt = -100 # Min value of pTdet in response matrix in JetPerformance task
    hResponseMatrixUncut = getResponseMatrix(ispp, useRMflexibleR, jetRadius, \
        matchingDistance, outputListResponse, minCent, maxCent, fMinPt, \
            maxPtGen, 0, maxPtGen, 0, 0, 0, 0, "Uncut",outputDir)
    hJetSpectrumTrueUncutPerBin = hResponseMatrixUncut.ProjectionY()
    hJetSpectrumTrueUncutPerBin.SetName("hJetSpectrumTrueUncutPerBin_KinEff")
    hJetSpectrumTrueUncutPerBin.Rebin(5)
    hJetSpectrumTrueUncutPerBin.Scale(1., "width")

    # Plot the ratio of the spectra
    hKinematicEfficiency = hJetSpectrumTrueCutPerBin.Clone()
    hKinematicEfficiency.SetName("hKinematicEfficiency")
    hKinematicEfficiency.Divide(hJetSpectrumTrueCutPerBin, hJetSpectrumTrueUncutPerBin, 1., 1., "B")
    
    hKinematicEfficiency.GetXaxis().SetTitle("#it{p}_{T}^{gen}")
    hKinematicEfficiency.GetYaxis().SetTitle("Kinematic Efficiency")
    hKinematicEfficiency.SetMarkerStyle(21)
    hKinematicEfficiency.SetMarkerColor(2)
    
    hinputRange = hKinematicEfficiency.Clone("inputRange")
    #set bins to zero outside the needed range
    for bin in range(1, hKinematicEfficiency.GetNcells() + 1):
        #if bin<hinputRange.FindBin(minPtReported) or bin>hinputRange.FindBin(maxPtReported)+1:
        if bin<hinputRange.FindBin(minPtReported) or bin>=hinputRange.FindBin(maxPtReported):
            hinputRange.SetBinContent(bin,0)
    text = "p_{T}^{det} #in [%d, %d] GeV" % (minPtDet, maxPtDet)
    
    outputFilename = os.path.join(outputDir, "hKinematicEfficiency_{}_{}{}".format(minPtDet, maxPtDet, fileFormat))
    plotHistKinEff(hKinematicEfficiency, hinputRange, outputFilename, "P E", False, False, text, ispp)


####################################################################################################
# Plot kinematic efficiency                                                                       ##
# (i.e. (pT-truth projection of response matrix with measured pT-det range selected)              ##
# / (pT-truth projection of response matrix with full pT-det range selected)                      ##
####################################################################################################
def plotKinematicEfficiency2(hRMUncut, minPtDet, maxPtDet, name):
    histName = hRMUncut.GetName() + "_UncutCP"
    hRMUncutCP = hRMUncut.Clone(histName)

    histName = hRMUncut.GetName() + "_CutCP"
    hRMCutCP = hRMUncut.Clone(histName)
    hRMCutCP.GetXaxis().SetRangeUser(minPtDet, maxPtDet)

    hKinDeno = hRMUncutCP.ProjectionY()
    hKinEff = hRMCutCP.ProjectionY()
    hKinEff.SetName(name)
    hKinEff.Divide(hKinEff, hKinDeno, 1., 1., "B")
    
    return hKinEff

