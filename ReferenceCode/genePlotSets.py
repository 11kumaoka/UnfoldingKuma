
#! /usr/bin/env python
# unfoldPlots

import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F
import argparse
import ctypes
import os
import gc
from array import array
import numpy as np

import math

# Set plotting options
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


# = s = Plot basic histogram  ##################################################
def plotHist(h, h2, outputFilename, drawOptions = "", \
    setLogy = False, setLogz = False, text = "", ispp = False, secondColor = False):
    
    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    c.cd().SetLeftMargin(0.15)

    if setLogy: c.SetLogy()
    if setLogz: c.SetLogz()
    h.Draw(drawOptions)

    if h2 and secondColor:
        #h2.SetFillColorAlpha(600-6, 0.3)
        h2.SetFillColorAlpha(632, 0.4)
        #h2.DrawCopy("same E3")
        h2.DrawCopy("same hist")

    if text:
        textFit = ROOT.TLatex()
        textFit.SetTextSize(0.04)
        textFit.SetNDC()
        if ispp:
            textFit.DrawLatex(0.6,0.8,text)
        else:
            textFit.DrawLatex(0.4,0.6,text)

    c.SaveAs(outputFilename)
    c.Close()
# = e = Plot basic histogram  ##################################################



# = s =    Plot spectra and ratio of h (and h3, if supplied) to h2           ###
def plotSpectra(h, h2, h3, nEvents, xRangeMin, xRangeMax, yAxisTitle, ratioYAxisTitle, \
    outputFilename, scalingOptions = "", legendTitle = "",\
        hLegendLabel = "", h2LegendLabel = "", h3LegendLabel = "", yRatioMax = 2.2):
    
    c = ROOT.TCanvas("c","c: pT",800,850)
    c.cd()
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad1.SetTopMargin(0.05)
    pad1.SetLogy()
    pad1.Draw()
    pad1.cd()
    
    h.SetLineColor(1)
    h.SetLineWidth(2)
    h.SetLineStyle(1)

    h.Scale(1./nEvents, scalingOptions)
    h.GetYaxis().SetTitle(yAxisTitle)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetXaxis().SetRangeUser(xRangeMin, xRangeMax)
    h.GetYaxis().SetRangeUser(2e-10,2e-3)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(20)
    xAxisTitle = h.GetXaxis().GetTitle()
    h.GetXaxis().SetTitle("")

    h2.SetLineColor(4)
    h2.SetLineWidth(2)
    h2.SetLineStyle(1)

    h.Draw("hist E")
    h2.Draw("hist same E")
    
    if h3:
        h3.SetLineColor(2)
        h3.SetLineWidth(2)
        h3.SetLineStyle(1)
        h3.Scale(1./nEvents, scalingOptions)
        h3.Draw("hist same")

    c.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.15)
    pad2.SetRightMargin(0.05)
    pad2.Draw()
    pad2.cd()

    # plot ratio h/h2
    hRatio = h.Clone()
    hRatio.Divide(h2)
    hRatio.SetMarkerStyle(21)
    hRatio.SetMarkerColor(1)

    hRatio.GetXaxis().SetRangeUser(xRangeMin,xRangeMax)
    hRatio.GetXaxis().SetTitleSize(30)
    hRatio.GetXaxis().SetTitleFont(43)
    hRatio.GetXaxis().SetTitleOffset(4.)
    hRatio.GetXaxis().SetLabelFont(43)
    hRatio.GetXaxis().SetLabelSize(20)
    hRatio.GetXaxis().SetTitle(xAxisTitle)
    
    hRatio.GetYaxis().SetTitle(ratioYAxisTitle)
    hRatio.GetYaxis().SetTitleSize(20)
    hRatio.GetYaxis().SetTitleFont(43)
    hRatio.GetYaxis().SetTitleOffset(2.2)
    hRatio.GetYaxis().SetLabelFont(43)
    hRatio.GetYaxis().SetLabelSize(20)
    hRatio.GetYaxis().SetNdivisions(505)

    min= hRatio.GetBinContent(hRatio.GetMinimumBin())
    max= hRatio.GetBinContent(hRatio.GetMaximumBin())
    #automatic zoom-in for a very small scatter of the points
    if min>0.5 and max<1.5:
        hRatio.GetYaxis().SetRangeUser(0.5,1.5)
    elif yRatioMax>2:
        hRatio.GetYaxis().SetRangeUser(0,yRatioMax)
    else:
        hRatio.GetYaxis().SetRangeUser(2-yRatioMax,yRatioMax)

    hRatio.Draw("P E")
    
    # plot ratio h3/h2
    if h3:
        hRatio3 = h3.Clone()
        hRatio3.Divide(h2)
        hRatio3.SetMarkerStyle(21)
        hRatio3.SetMarkerColor(2)
        hRatio3.Draw("P E same")

    line = ROOT.TLine(xRangeMin,1,xRangeMax,1)
    line.SetLineColor(920+2)
    line.SetLineStyle(2)
    line.Draw()
    
    pad1.cd()

    if nEvents > 2:
        textNEvents = ROOT.TLatex()
        textNEvents.SetNDC()
        textNEvents.DrawLatex(0.55,0.6,"#it{N}_{events} = %d" % nEvents)

    leg2 = ROOT.TLegend(0.3,0.7,0.88,0.93,legendTitle)
    leg2.SetFillColor(10)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextSize(0.04)
    leg2.AddEntry(h, hLegendLabel, "l")
    if h3:
        leg2.AddEntry(h3, h3LegendLabel, "l")
    if h2:
        leg2.AddEntry(h2, h2LegendLabel, "l")
    leg2.Draw("same")
    
    c.SaveAs(outputFilename)
    c.Close()
# = e =    Plot spectra and ratio of h (and h3, if supplied) to h2           ###



#  = s =  Plot three different histograms     ##################################
def plotHistCompare(plotRatioToo, h1, h2, h3, outputFilename, drawOptions = "", \
    setLogy = False, setLogz = False, text = ""):

    if plotRatioToo:
        #c = ROOT.TCanvas("c","c: hist",600,450)
        c = ROOT.TCanvas("cSpectrum{}".format(plotRatioToo),"cSpectrum: hist",800,850)
    else:
        c = ROOT.TCanvas("cSpectrum{}".format(plotRatioToo),"cSpectrum: hist",600,600)
    c.Draw()
    c.cd()
    
    # Set pad and histo arrangement
    if plotRatioToo:
        myPad = ROOT.TPad("myPad", "The pad",0,0.3,1,1)
    else:
        myPad = ROOT.TPad("myPad", "The pad",0,0,1,1)
    myPad.SetLeftMargin(0.22)
    myPad.SetTopMargin(0.04)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.15)
    if plotRatioToo:
        myPad.SetBottomMargin(0.0)
    myPad.Draw()
    myPad.SetLogy()
    myPad.cd()

    if setLogy: c.SetLogy()
    if setLogz: c.SetLogz()

    h1.SetMarkerStyle(21)
    h1.SetMarkerColor(4)
    h1.GetYaxis().SetTitleOffset(1.8)

    h1.SetNdivisions(505)
    h1.SetXTitle("#it{p}_{T,jet} (GeV/#it{c})")
    h1.GetYaxis().SetTitleOffset(2.2)
    h1.SetYTitle("Normalized counts (Integral=1)")
    #h1.SetMaximum(8e-4)
    #h1.SetMinimum(2e-8)
    h1.GetYaxis().SetTitleSize(0.06)
    h1.GetYaxis().SetTitleOffset(1.6)
    h1.GetYaxis().SetLabelSize(0.06)
    h1.GetXaxis().SetRangeUser(10,250)
    h1.DrawCopy(drawOptions)

    if h2:
        h2.SetMarkerStyle(20)
        h2.SetMarkerColor(2)
        h2.DrawCopy(drawOptions + "same")

    if h3:
        h3.SetMarkerStyle(22)
        h3.SetMarkerColor(8)
        h3.DrawCopy(drawOptions + "same")
    #This is to be used for JES comparision - currently unused
    h1LegendLabel = "R = 0.2"
    h2LegendLabel = "R = 0.3"
    h3LegendLabel = "R = 0.4"
    #This is for prior comparision
    h1LegendLabel = "Particle level proj. of the response (cut)"
    h2LegendLabel = "POWHEG or unfolded spectrum as prior"
    h3LegendLabel = "smeared spectrum (A) or uncut projection (B)"
    leg = ROOT.TLegend(0.45,0.7,0.6,0.85,"")
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(h1, h1LegendLabel, "P")
    if h2:
        leg.AddEntry(h2, h2LegendLabel, "P")
    if h3:
        leg.AddEntry(h3, h3LegendLabel, "P")
    leg.Draw("same")

    if text:
        textFit = ROOT.TLatex()
        textFit.SetTextSize(0.04)
        textFit.SetNDC()
        textFit.DrawLatex(0.6,0.8,text)

    if plotRatioToo:
        c.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.02, 1, 0.3)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.4)
        pad2.SetLeftMargin(0.22)
        pad2.SetRightMargin(0.04)
        pad2.Draw()
        pad2.cd()

        if h2:
            #h2copy = h2
            h2copy = h2.Clone("h2_CopyForDivision")
            h2copy.Divide(h1)

            h2copy.SetYTitle("#frac{New Prior}{Orig. Prior}")
            h2copy.SetXTitle("#it{p}_{T,jet} (GeV/#it{c})")
            h2copy.GetXaxis().SetTitleSize(30)
            h2copy.GetXaxis().SetTitleFont(43)
            h2copy.GetXaxis().SetTitleOffset(4.)

            h2copy.GetXaxis().SetLabelFont(43)
            #h2copy.GetXaxis().SetLabelOffset(0.2)
            h2copy.GetXaxis().SetLabelSize(30)

            h2copy.GetYaxis().SetTitleSize(20)
            h2copy.GetYaxis().SetTitleFont(43)
            h2copy.GetYaxis().SetTitleOffset(2.2)
            h2copy.GetYaxis().SetLabelFont(43)
            h2copy.GetYaxis().SetLabelSize(20)
            h2copy.GetYaxis().SetNdivisions(505)
            h2copy.GetXaxis().SetRangeUser(10,250)
            h2copy.GetYaxis().SetRangeUser(0.5,2.5)
            #h2copy.GetYaxis().SetRangeUser(0.1,1.2)
            h2copy.DrawCopy("E")
        if h3:
            #h3copy = h3
            h3copy = h3.Clone("h3_CopyForDivision")
            h3copy.Divide(h1)
            h3copy.DrawCopy("same E")

    c.SaveAs(outputFilename)
    c.Close()
    #  = e =  Plot three different histograms     ##############################



