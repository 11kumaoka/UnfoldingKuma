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


###################################################################################
def JESshiftCheck():
    # inFileDir = "./"
    inFileDir = '~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18qr/pass3/Ch/Embedding/'
    inFileName = inFileDir + 'EmbedPtHardScaledResults_TrackPtCut5__Ver1.root'
    # outFileDir = './'
    outFileDir = '~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18qr/pass3/Ch/Embedding/JESshift/'
    outFile = outFileDir + 'JESShiftCheck.root'

    outRootFile = ROOT.TFile(outFile, 'RECREATE')
    lOMainTree = ROOT.TList()
    
    f = ROOT.TFile(inFileName, "READ")
    mainTaskName = 'mainTree'
    lMainTree = f.Get(mainTaskName)
    f.ls()
    lOEPTree = lMainTree.FindObject('OutOfPlane')
    lIEPTree = lMainTree.FindObject('InPlane')
    
    hJESshiftOEP = lOEPTree.FindObject('hJESshift')
    hJESshiftIEP = lIEPTree.FindObject('hJESshift')
    # hJESshiftOEP = lOEPTree.FindObject('hJESshiftHybDet')
    # hJESshiftIEP = lIEPTree.FindObject('hJESshiftHybDet')
    
    for centBin in range(0, numOfCentBin):
        ptRangeDict, ptBinArrayDict = PtRangeList.eachJetPtBinDef(0,centBin)
        
        histName = 'hJESshiftOEPCent'+str(centBin)
        hJESshiftOEP_CP = hJESshiftOEP.Clone(histName)
        histName = 'hJESshiftIEPCent'+str(centBin)
        hJESshiftIEP_CP = hJESshiftIEP.Clone(histName)

        hJESshiftOEP_CP.GetXaxis().SetRangeUser(centRangeList[centBin][0], centRangeList[centBin][1])
        hJESshiftIEP_CP.GetXaxis().SetRangeUser(centRangeList[centBin][0], centRangeList[centBin][1])
        
        hJESshift2DOEP = hJESshiftOEP_CP.Project3D("zyeo")
        hJESshift2DIEP = hJESshiftIEP_CP.Project3D("zyeo")
        
        hJESshiftProfOEP = hJESshift2DOEP.ProfileX()
        hJESshiftProfIEP = hJESshift2DIEP.ProfileX()

        nBins = len(ptBinArrayDict['mcDet']) - 1
        binArray = ptBinArrayDict['mcDet']
        
        hJESshiftProfOEP.Sumw2()
        histName = 'hJESshiftProfOEPCent'+str(centBin)
        hJESshiftProfOEP = hJESshiftProfOEP.Rebin(nBins, histName, binArray)
        hJESshiftProfIEP.Sumw2()
        histName = 'hJESshiftProfIEPCent'+str(centBin)
        hJESshiftProfIEP = hJESshiftProfIEP.Rebin(nBins, histName, binArray)

        # = s =    Plot spectra and ratio of h (and h3, if supplied) to h2           ###
        c = ROOT.TCanvas("c","c: pT",800,850)
        c.cd()
        pad1 = ROOT.TPad("pad1", "pad1", 0,0.3,1,1)
        pad1.SetBottomMargin(0)
        pad1.SetLeftMargin(0.15)
        pad1.SetRightMargin(0.05)
        pad1.SetTopMargin(0.05)
        pad1.SetLogy()
        pad1.Draw()
        pad1.cd()
        
        yTitle = "#frac{#it{p}_{T}^{det} - #it{p}_{T}^{particle}}{#it{p}_{T}^{particle}}"
        hJESshiftProfOEP.GetYaxis().SetTitle(yTitle)
        hJESshiftProfOEP.GetYaxis().SetTitleSize(0.06)
        # hJESshiftProfOEP.GetYaxis().SetRangeUser(2e-10,2e-3)
        hJESshiftProfOEP.GetYaxis().SetLabelFont(43)
        hJESshiftProfOEP.GetYaxis().SetLabelSize(20)
        xTitle = '#it{p}_{T, particle}^{jet}'
        hJESshiftProfOEP.GetXaxis().SetTitle(xTitle)

        hJESshiftProfOEP.Draw("hist PE")
        hJESshiftProfIEP.Draw("hist same PE")

        legendTitle = 'JESshift'
        leg1 = ROOT.TLegend(0.7,0.7,0.93,0.93,legendTitle)
        leg1.SetFillColor(10)
        leg1.SetBorderSize(0)
        leg1.SetFillStyle(0)
        leg1.SetTextSize(0.05)
        leg1.AddEntry(hJESshiftProfOEP, 'JESshift Out-of-plane', "l")
        leg1.AddEntry(hJESshiftProfIEP, 'JESshift In-plane', "l")
        leg1.Draw("same")

        c.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.02, 1, 0.3)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.)
        pad2.SetLeftMargin(0.15)
        pad2.SetRightMargin(0.05)
        # pad2.SetLogy()
        pad2.Draw()
        pad2.cd()

        # plot ratio h/h2 
        hRatioName = 'hRatioJESshiftCent'+str(centBin)
        hRatio = hJESshiftProfOEP.Clone(hRatioName)
        hRatio.Divide(hJESshiftProfIEP)
        
        ratioYAxisTitle = '#frac{Out-of-plane}{In-Plane}'
        hRatio.GetYaxis().SetTitle(ratioYAxisTitle)
        hRatio.GetYaxis().SetTitleSize(15)
        hRatio.GetYaxis().SetTitleFont(43)
        hRatio.GetYaxis().SetTitleOffset(2.2)
        hRatio.GetYaxis().SetLabelFont(43)
        hRatio.GetYaxis().SetLabelSize(20)
        hRatio.GetYaxis().SetNdivisions(505)

        #automatic zoom-in for a very small scatter of the points
        hRatio.GetYaxis().SetRangeUser(0.2,2.5)
        
        hRatio.Draw("P E")
        
        pad1.cd()
        
        outputFilename = outFileDir + 'JESshiftProfCent'+str(centBin)+'.root'
        c.SaveAs(outputFilename)
        c.Close()

        
        lOMainTree.Add(hJESshiftProfOEP)
        lOMainTree.Add(hJESshiftProfIEP)
        lOMainTree.Add(hRatio)
        
        # = e =    Plot spectra and ratio of h (and h3, if supplied) to h2           ###
    
    # outRootFile = ROOT.TFile(outFile, 'RECREATE')
    # lMainTree = ROOT.TList()
    # OFileStracture(lMainTree, numOfCentBin)
    
    outRootFile.cd()
    lOMainTree.Write('mainTree', 1)
    outRootFile.Close()
    print('root ' + outFile)


if __name__ == "__main__":
    JESshiftCheck()


