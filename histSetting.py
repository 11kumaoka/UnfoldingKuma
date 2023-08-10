import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F

lFillMerker = [20, 21, 22, 23, 34, 29, 33, 31, 24, 24, 26, 32, 20, 21, 22, 23, 34, 29]
lNoFillMerker = [24, 25, 26, 32, 28, 30, 27, 24, 25, 26]

#[kBrack, kRed, kBlue, kGreen+3, kOrange+8, kAzure+10, kSpring+4, kViolet+1, kPink+9, kSpring-9, kOrange+8, kAzure+1, kTeal+4, kRed-7, kBlue-7, kGreen-2, kOrange-2, kViolet+7, kSpring+4, kTeal-5]
lFillColor = [1, 632, 600, 416+3, 800+8, 860+10, 900+10, 820+4, 880+1, 900+9, 820-9, 800+1, 860+1, 840+4, 632-7, 600-7, 416-2, 800-2, 880+7, 820+4, 840-5]

def histLabelSetting(outputHist, histTitle, xTitle, yTitle, histNumber, fillType, Rebin):
    outputHist.SetTitle(histTitle)
    outputHist.SetXTitle(xTitle)
    outputHist.SetYTitle(yTitle)
    
    outputHist.GetXaxis().SetTitleSize(20)
    outputHist.GetXaxis().SetTitleFont(43)
    outputHist.GetXaxis().SetTitleOffset(2.)

    outputHist.GetYaxis().SetTitleSize(20)
    outputHist.GetYaxis().SetTitleFont(43)
    outputHist.GetYaxis().SetTitleOffset(2.)
    
    outputHist.SetLineColor(lFillColor[histNumber])
    outputHist.SetMarkerColor(lFillColor[histNumber])
    if(fillType): outputHist.SetMarkerStyle(lFillMerker[histNumber])
    else: outputHist.SetMarkerStyle(lNoFillMerker[histNumber])
    outputHist.SetMarkerSize(1.8)
    outputHist.SetLineWidth(5)

    #outputHist.SetMinimum(0.000000000000000000000001)
    #outputHist.SetMaximum(50)

    outputHist.Rebin(Rebin)


def ratioMergeCanvas(outputBranch, lBaseHist, canvas, ratioHistoYTitle, hXRange, hYRange, kind):

    lRatioHist = list() 
    ratioHistXTitle = lBaseHist[0].GetXaxis().GetTitle()
    lBaseHist[0].GetYaxis().SetTitleOffset(1.8)

    for histNum in range(1, len(lBaseHist)):
        ratioHistName = "ratioHisKind" + str(kind)
        tempRatioHist = ROOT.TH1F()
        tempRatioHist = lBaseHist[histNum].Clone(ratioHistName)

        # tempRatioHist.GetXaxis().SetRangeUser(-50., 250.)
        # tempRatioHist.GetYaxis().SetRangeUser(1.*10e-9, 50.)

        tempRatioHist.GetYaxis().SetNdivisions(505)
        # tempRatioHist.Sumw2()
        tempRatioHist.SetStats(0)
        tempRatioHist.Divide(tempRatioHist,lBaseHist[0],1.,1.,"B")
        # tempRatioHist.Write("eff_cent"+str(kind))
        lRatioHist.append(tempRatioHist)

    # Draw Part Start  ===============================================
    # Setting Pad Parameters start  ##################################
    gap = 0
    margin = 0.025 # for bottom y label (like 0)
    uHeight = 0.7
    dHeight = 1.0 - uHeight
    
    # Setting Pad Parameters end   ###################################
    # canvas = ROOT.TCanvas(ratioCanvasName, ratioCanvasName, 800, 800)
    canvas.cd()

    uPad = ROOT.TPad("uPad", "uPad", 0, 1.0-uHeight-margin, 1, 1.0)
    dPad = ROOT.TPad("dPad", "dPad", 0, 0.00, 1, dHeight)
    uPad.SetGridx()
    dPad.SetGridx()

    uPad.SetLeftMargin(0.25)
    dPad.SetLeftMargin(0.25)
    #important setting for deviding plot #############################
    uPad.SetBottomMargin(margin / uHeight)
    dPad.SetTopMargin(gap * margin / dHeight)
    dPad.SetBottomMargin(0.30)
    dPad.SetFillColorAlpha(10,0.0)
    #important setting for deviding plot #############################

    #draw graphics start  ############################################
    canvas.cd()
    uPad.Draw()
    uPad.cd() # from canvas to uPad
    
    # lBaseHist[0].GetXaxis().SetRangeUser(-50., 250.)
    # lBaseHist[0].GetYaxis().SetRangeUser(1.*10e-9, 50.)
    lBaseHist[0].SetStats(0)
    lBaseHist[0].GetXaxis().SetLabelSize(0.)
    lBaseHist[0].Draw()
    for histNum in range(1, len(lBaseHist)):
        lBaseHist[histNum].Draw("same E")

    uPad.SetLogy(1)
    canvas.cd() # from uPad to canvas
    dPad.Draw()
    dPad.cd() # from canvas to dPad
    # to erase the top value on y axis of lower plot
    lRatioHist[0].SetMaximum(1.40 - 0.01*(1-gap))
    lRatioHist[0].GetXaxis().SetLabelSize(0.09)
    lRatioHist[0].GetXaxis().SetTitleOffset(3.)
    lRatioHist[0].GetYaxis().SetLabelSize(0.09)
    # lRatioHist[0].GetYaxis().SetTitleOffset(3.)
    lRatioHist[0].Draw()
    for histNum in range(0, len(lRatioHist)):
        lRatioHist[histNum].Draw("same E")
        

    # draw graphics end ##############################################
    # Draw Part End  ================================================
    outputBranch[kind].Add(canvas)
