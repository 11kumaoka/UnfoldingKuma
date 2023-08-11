#! /usr/bin/env python

import argparse
import ctypes
import os
import gc
from array import array
import numpy as np

import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F


def main(outList, hOriginRM, label):
    ispp = True
    isClosure = False
    variation = 1
    minPtGen = 5.
    maxPtGen = 250.
    minPtDet = 5. 
    maxPtDet = 250.

    nGenBins, genBinArray, nDetBins, detBinArray \
        = returnGenDetPtBins(label, minPtGen, maxPtGen, minPtDet, maxPtDet, \
            ispp, isClosure, variation)
    hRM_rebin = rebinRM(hOriginRM, nGenBins, genBinArray, nDetBins, detBinArray)
    hRM_norm = normalizeRM(hRM_rebin, outList, minPtGen, maxPtGen, minPtDet, maxPtDet)
    
    return hRM_norm


################################################################################
# Rebin the response matrix to have variable binning                          ##
################################################################################
def rebinRM(hOriginRM, nGenBins, genBinArray, nDetBins, detBinArray):
    histname = "{}NewRebinned".format(hOriginRM.GetName())
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

            #print "Adding value {} from bin ({},{}) = ({},{})".format(oldContent, ibin, jbin, x, y)
            #newBin = hRM_rebin.FindBin(x,y)
            #print "New bin content: {}".format(hRM_rebin.GetBinContent(newBin))
    
    # Assume 0 errors on response matrix
    for bin in range(1, hRM_rebin.GetNcells() + 1):
        hRM_rebin.SetBinError(bin, 0)

    return hRM_rebin


################################################################################
# Normalize response matrix                                                   ##
# Normalize the pT-truth projection to 1                                      ##
################################################################################
def normalizeRM(hRM_rebin, outputDir, minPtGen, maxPtGen, minPtDet, maxPtDet):

    # Make projection onto pT-true axis (y-axis), and scale appropriately
    # Do exclude under and overflow bins
    hGenProjBefore = hRM_rebin.ProjectionY("_py",1,hRM_rebin.GetNbinsX()) 
    hGenProjBefore.SetName("hGenProjectionBefore")

    # Loop through truth-level bins, and apply normalization factor to all bins.
    nBinsY = hRM_rebin.GetNbinsY() # pT-gen
    nBinsX = hRM_rebin.GetNbinsX() # pT-det
    for genBin in range(1,nBinsY+1):
        normFactor = hGenProjBefore.GetBinContent(genBin)
        if normFactor > 0:
            genBinCenter = hGenProjBefore.GetXaxis().GetBinCenter(genBin)

            for detBin in range(1,nBinsX+1):
                binContent = hRM_rebin.GetBinContent(detBin, genBin)
                hRM_rebin.SetBinContent(detBin, genBin, binContent/normFactor)
        
    # Plot response matrix
    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    c.cd().SetLeftMargin(0.15)

    hRM_rebin.Draw("colz")
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
    
    c.SaveAs('hoge.root')
    c.Close()

    return hRM_rebin

def returnGenDetPtBins(label, minPtGen, maxPtGen, minPtDet, maxPtDet, ispp, isClosure, variation):
    
    binArrayDet = None
    binArrayTruth = None
    if ispp and not isClosure:
        #For PbPb comparision
        #binArrayTruth = ([minPtGen, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 140, 190, maxPtGen]) 
        binArrayTruth = ([minPtGen, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 140, 190, maxPtGen]) 
        #- - - - - - - - - - - - - - - - - - - - - - -
        #biased - mainyly for comparision to PbPb
        if "02" in label:
            binArrayDet   = ([minPtDet, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, \
                65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, maxPtDet])
        elif "04" in label:
            if "5GeV" in label:
                #binArrayDet = ([minPtDet, 20, 25, 30, 35, 40, 45, 50, 55, 60, \
                # 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, maxPtDet])
                binArrayDet = ([minPtDet, 20, 25, 30, 35, 40, 45, 50, 55, 60, \
                    65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 120, maxPtDet])
            if "7GeV" in label:
                binArrayDet = ([minPtDet, 20, 25, 30, 35, 40, 45, 50, 55, 60, \
                    65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, maxPtDet])
        #--  --  --  --  --  --  --  --  --  --  --
        #unbiased for pure pp analysis
        if "0GeV" in label:
            binArrayTruth = ([minPtGen, 10, 20, 30, 40, 50, 60, 70, 80, \
                100, 120, 140, 190, maxPtGen])
            if "05" in label:
                #binArrayTruth = ([minPtGen, 10, 20, 30, 40, 50, 60, 70, 80, \
                # 100, 120, 130, 190, maxPtGen])
                binArrayTruth = ([minPtGen, 10, 20, 30, 40, 50, 60, 70, 80, \
                    100, 120, 140, 190, maxPtGen])
        #- - - - - - - - - - - - - - - - - - - - - - -
        #Detector level input - measured spectrum
        if minPtDet<9:
            binArrayDet = ([minPtDet, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, \
                65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 120, maxPtDet])
        elif minPtDet<15:
            binArrayDet = ([minPtDet, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, \
                65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 120, maxPtDet])
        elif minPtDet<20:
            binArrayDet = ([minPtDet, 20, 25, 30, 35, 40, 45, 50, 55, 60, \
                65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 120, maxPtDet])
        else:
            binArrayDet = ([minPtDet, 25, 30, 35, 40, 45, 50, 55, 60, \
                65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 120, maxPtDet])
        if maxPtDet<=120:
            if minPtDet<9:
                #binArrayDet = ([minPtDet, 10, 15, 20, 25, 30, 35, 40, 45, 50, \
                # 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, maxPtDet])
                binArrayDet = ([minPtDet, 10, 15, 20, 25, 30, 35, 40, 45, 50, \
                    55, 60, 65, 70, 75, 80, 85, 90, 95, 105, maxPtDet])
            elif minPtDet<15:
                binArrayDet = ([minPtDet, 15, 20, 25, 30, 35, 40, 45, 50, \
                    55, 60, 65, 70, 75, 80, 85, 90, 95, 105, maxPtDet])
                #binArrayDet = ([minPtDet, 15, 20, 25, 30, 35, 40, 45, 50, \
                # 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, maxPtDet])
            elif minPtDet<20:
                binArrayDet = ([minPtDet, 20, 25, 30, 35, 40, 45, 50, \
                    55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, maxPtDet])
            else:
                binArrayDet = ([minPtDet, 25, 30, 35, 40, 45, 50, \
                    55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, maxPtDet])
        if maxPtDet<=90:
            if minPtDet<9:
                binArrayDet = ([minPtDet, 10, 15, 20, 25, 30, 35, 40, 45, 50, \
                    55, 60, 65, 70, 75, 80, maxPtDet])
            elif minPtDet<15:
                binArrayDet = ([minPtDet, 15, 20, 25, 30, 35, 40, 45, 50, \
                    55, 60, 65, 70, 75, 80, maxPtDet])
            elif minPtDet<20:
                binArrayDet = ([minPtDet, 20, 25, 30, 35, 40, 45, 50, \
                    55, 60, 65, 70, 75, 80, maxPtDet])
            else:
                binArrayDet = ([minPtDet, 25, 30, 35, 40, 45, 50, \
                    55, 60, 65, 70, 75, 80, maxPtDet])
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #For the PbPb case:
    else:
        if "02" in label:
            if variation==2: # rangeVar = -5
                binArrayDet = ([minPtDet, 20, 25, 30, 35, 40, 45, 50, 55, 60, \
                    65, 70, 75, 80, 85, 90, 95, 100, 105, 110, maxPtDet])
            elif variation==3: # rangeVar = +5
                binArrayDet = ([minPtDet, 30, 35, 40, 45, 50, 55, 60, \
                    65, 70, 75, 80, 85, 90, 95, 100, 105, 110, maxPtDet])
            else:
                binArrayDet = ([minPtDet, 25, 30, 35, 40, 45, 50, 55, 60, \
                    65, 70, 75, 80, 85, 90, 95, 100, 105, 110, maxPtDet])
                binArrayTruth = ([minPtGen, 10, 20, 30, 40, 50, 60, \
                    70, 80, 100, 120, 140, 190, maxPtGen])
        elif "04" in label:
            if variation==2: # rangeVar = -5
                binArrayDet = ([minPtDet, 35, 40, 45, 50, 55, 60, 65, 70, 75, \
                    80, 85, 90, 95, 100, 105, 110, 115, maxPtDet])
            elif variation==3: # rangeVar = +5
                binArrayDet = ([minPtDet, 45, 50, 55, 60, 65, 70, 75, \
                    80, 85, 90, 95, 100, 105, 110, 115, maxPtDet])
        else:
            binArrayDet = ([minPtDet, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, \
                90, 95, 100, 105, 110, 115, maxPtDet])
            binArrayTruth = ([minPtGen, 10, 20, 30, 40, 50, 60, 70, 80, 100, \
                120, 140, 190, maxPtGen])

    ## test array  ##########################################################################
    binArrayDet = ([minPtDet, 10, 15, 20, 25, 30, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, \
                90, 95, 100, 105, 110, 115, 120, 130, 140, 150, 170, 180, 190, 200, maxPtDet])
    binArrayTruth = ([minPtDet, 10, 15, 20, 25, 30, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, \
                90, 95, 100, 105, 110, 115, 120, 130, 140, 150, 170, 180, 190, 200, maxPtDet])
    ## test array  ##########################################################################

    nDetBins = len(binArrayDet) - 1
    detBinArray = array('d',binArrayDet)
    nGenBins = len(binArrayTruth) - 1
    genBinArray = array('d',binArrayTruth)

    return nGenBins, genBinArray, nDetBins, detBinArray