import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TMath
import argparse
import ctypes
import os
import gc
from array import array
import numpy as np

import sys

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

numOfCentBin = 5

lEPLabel = ['OutOfPlane', 'InPlane', 'Inclusive']

###################################################################################
# Main function (0/1/2/3/4, "LHC18q/LHC18r", 0/5/7, "98%: ''/94%: 'TrackEff094'"
def HAddOutputScaledFiles(LHCPeriod, leadingTrackPtCut,  diffSys):
    # outEmbFileDir = './'
    outEmbFileDir = '~/ALICE/cernbox/SWAN_projects/outputFiles/'+LHCPeriod+'/pass3/Ch/Embedding/'
    outEmbFile = outEmbFileDir + 'EmbedPtHardScaledResults'\
        +'_TrackPtCut'+str(leadingTrackPtCut)+'_'+diffSys+'.root'

    outRootFile = ROOT.TFile(outEmbFile, 'RECREATE')
    lMainTree = ROOT.TList()
    OFileStracture(lMainTree, numOfCentBin)


    for centBin in range(0, numOfCentBin):
        for epBin in range(0, 2):
            listName = lEPLabel[epBin]
            inEmbFileDir = '~/ALICE/cernbox/SWAN_projects/outputFiles/'\
                +LHCPeriod+'/pass3/Ch/Embedding/'
            inEmbFileName = inEmbFileDir + 'EmbedPtHardScaledResults'\
                +'_TrackPtCut'+str(leadingTrackPtCut)+'_'+diffSys+'_CentBin'+str(centBin)+'.root'
            print('Open: '+inEmbFileName)
            inputFile = ROOT.TFile(inEmbFileName, "READ")
            iMainTree = inputFile.Get('mainTree')
            iLEBList = iMainTree.FindObject(lEPLabel[epBin])
            print(' Add: '+listName)
            listName = 'lCent'+str(centBin)
            iLCentList = iLEBList.FindObject(listName)
            lMainTree[epBin].Add(iLCentList)

            if centBin == numOfCentBin-1:
                for histBin in range(numOfCentBin, len(iLEBList)):
                    lMainTree[epBin].Add(iLEBList[histBin])
                # print(len(iMainTree))
                for histBin in range(3, len(iMainTree)):
                    lMainTree.Add(iMainTree[histBin])


    for centBin in range(0, numOfCentBin):
        lTempCent = ROOT.TList()
        lTempCent.SetName('lCent'+str(centBin))
        lMainTree[2].Add(lTempCent)
        for histBin in range(0, len(lMainTree[0][centBin])):
            tempHist = lMainTree[0][centBin][histBin].Clone()
            lMainTree[2][centBin].Add(tempHist)
            lMainTree[2][centBin][histBin].Add(lMainTree[1][centBin][histBin])
    for histBin in range(numOfCentBin, len(lMainTree[0])):
        tempHist = lMainTree[0][histBin].Clone()
        # lMainTree[2].Add(lMainTree[0][histBin])
        lMainTree[2].Add(tempHist)
        lMainTree[2][histBin].Add(lMainTree[1][histBin])

    outRootFile.cd()
    lMainTree.Write('mainTree', 1)
    outRootFile.Close()
    print('root ' + outEmbFile)

################################################################################
def OFileStracture(lMainTree, numOfCentBin):
    for epBin in range(0, 3):
        listName = lEPLabel[epBin]
        lEPList = ROOT.TList()
        lEPList.SetName(listName)
            
        lMainTree.Add(lEPList)
################################################################################

if __name__ == "__main__":
    args = sys.argv
    HAddOutputScaledFiles(args[1], int(args[2]),args[3])

