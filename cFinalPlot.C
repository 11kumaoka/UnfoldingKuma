
#include <TRandom3.h>
#include <TMath.h>
#include <TChain.h>


#include <TClonesArray.h>
#include <TFile.h>
#include <THashList.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THnSparse.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TList.h>

#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>

#include "TGrid.h"

void overwrightSomePlots(TH1D* hBlank, TH1D* lHists[5], \
    TLegend *legALICE,TLegend* legJet,TLegend* leg, Int_t colorList[5], Int_t styleList[5], \
        TString label, TString addDrawOption, TString outputDir);
void overwrightSomePlotsWithErr(TH1D* hBlank,TH1D* lHists[5], TGraphErrors* TGroughlHErr[5], \
    Int_t listSize, TLegend* legALICE,TLegend* legJet,TLegend* leg,\
        Int_t colorList[5], Int_t styleList[5], TString label, TString addDrawOption,\
            TString outputDir);
void setupPadMargin(TPad *pad);
void addSomeHistsLegend(TLegend *leg, std::vector<TH1D* > lHists, std::vector<TString > lHistTags, \
    Double_t textSize);
void addSomeHistsLegend(TLegend *leg, TH1D* lHists[5], TString lHistTags[5], \
    Int_t listSize, Double_t textSize);
void addALICELegend(TH1D* blankObj, TLegend* leg, \
    TString centLabel, TString addLabel, Double_t textSize);
void addJetLegend(TH1D *blankObj, TLegend* leg, Double_t Reso, Int_t ptCut, Double_t textSize);
void setHistLooks(TH1D* h, Int_t style, Int_t color);

Int_t numOfCentBin = 5;
TString lCentLabel[5] = {"0-5","5-10","10-30","30-50","50-80"};
// std::vector<TH1D*> lHChJetV2;
// std::vector<TGraphErrors*> lHChJetV2SysErr;
TH1D* lHChJetV2[5];
TGraphErrors* lHChJetV2SysErr[5];
Int_t leadingTrackPtCut = 5;
Double_t JetReso = 0.2;

TString iFileDir = "~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18qr/pass3/Ch/FinalPlots/";

// TString outEmbFileDir = "./";
TString oFileDir = "~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18qr/pass3/Ch/FinalPlots/";

int main(){
    TString iFileName = TString::Format("%sFinalResultsWithSys.root", iFileDir.Data());
    TFile *iFile = new TFile(iFileName.Data(), "READ");
    
    TString mainTreeName = "mainTree";
    TList *lMainTree = (TList *) iFile->Get(mainTreeName.Data());

    for(Int_t centBin=0; centBin<numOfCentBin; centBin++){
        TString centTreeName = TString::Format("lCent%d", centBin).Data();
        TList *lCentTree = (TList *) lMainTree->FindObject(centTreeName);
        TList *lJetTree = (TList *) lCentTree->FindObject("lChJetV2");
        TString chJetV2Name = TString::Format("hChJetV2_Cent%d_Nominal", centBin).Data();
        TH1D *hChJetV2 = (TH1D *)lJetTree->FindObject(chJetV2Name);
        // lHChJetV2.back_push(hChJetV2);
        lHChJetV2[centBin] = hChJetV2;
        
        TList *lSysErrTree = (TList*) lCentTree->FindObject("lSysUncErrGrough");
        TGraphErrors* sysErrChJetV2 = (TGraphErrors*) lSysErrTree->FindObject("Graph");
        // lHChJetV2SysErr.back_push(sysErrChJetV2);
        lHChJetV2SysErr[centBin] = sysErrChJetV2;


        TH1D *tempHChJetV2[5];
        tempHChJetV2[0] = hChJetV2;
        TGraphErrors *tempHChJetV2SysErr[5];
        tempHChJetV2SysErr[0] = sysErrChJetV2;
        TH1D *hBlank = new TH1D();
        TLegend *legALICE = new TLegend(0.3,0.7,0.88,0.93,"");
        addALICELegend(hBlank, legALICE, "", "Work in progress", 0.03);
        TLegend *legJet = new TLegend(0.3,0.7,0.88,0.93,"");
        addJetLegend(hBlank, legJet, 0.2, leadingTrackPtCut, 0.03);
        TLegend *legChJetV2AllCent = new TLegend(0.3,0.7,0.88,0.93,"Centrality[%]");
        addSomeHistsLegend(legChJetV2AllCent, tempHChJetV2, lCentLabel, 5, 0.03);

        TString label = TString::Format("FinalPlotSingleCent%dSysUncError_TrackPtCut%d", \
            centBin, leadingTrackPtCut);
        Int_t lPlotColor[5] = {1, 600, 419, 909, 808};
        Int_t lPlotStyle[5] = {0,0,0,0,0};
        overwrightSomePlotsWithErr(hBlank, tempHChJetV2, tempHChJetV2SysErr, 1, \
            legALICE,legJet,legChJetV2AllCent, lPlotColor, lPlotStyle, label, "", oFileDir);
    }


    TH1D *hBlank = new TH1D();
    TLegend *legALICE = new TLegend(0.3,0.7,0.88,0.93,"");
    addALICELegend(hBlank, legALICE, "", "Work in progress", 0.03);
    TLegend *legJet = new TLegend(0.3,0.7,0.88,0.93,"");
    addJetLegend(hBlank, legJet, 0.2, leadingTrackPtCut, 0.03);
    TLegend *legChJetV2AllCent = new TLegend(0.3,0.7,0.88,0.93,"Centrality[%]");
    addSomeHistsLegend(legChJetV2AllCent, lHChJetV2, lCentLabel, 5, 0.03);

    TString label = TString::Format("FinalPlotSysUncError_TrackPtCut%d", leadingTrackPtCut);
    Int_t lPlotColor[5] = {1, 600, 419, 909, 808};
    Int_t lPlotStyle[5] = {0,0,0,0,0};
    overwrightSomePlotsWithErr(hBlank, lHChJetV2, lHChJetV2SysErr, 5, \
        legALICE,legJet,legChJetV2AllCent, lPlotColor, lPlotStyle, label, "", oFileDir);
    
    std::cout << "root " <<  oFileDir << std::endl;

}


// # = s = Overwright Some Plots on One canvas ####################################
void overwrightSomePlots(TH1D* hBlank, TH1D* lHists[5], Int_t listSize,\
    TLegend *legALICE,TLegend* legJet,TLegend* leg, Int_t colorList[5], Int_t styleList[5], \
        TString label, TString addDrawOption, TString outputDir){
    TString canvasName = "canvas_" + label;
    TCanvas* c = new TCanvas(canvasName.Data(),canvasName.Data(), 800,800);
    c->cd();
    TString padName = "pad_" + label;
    TPad* pad = new TPad(padName.Data(), padName.Data(), 0, 0.3, 1, 1.0);
    setupPadMargin(pad);
    // # pad.SetLogy();
    pad->Draw();
    pad->cd();

    TString drawOption = addDrawOption + "same";
    for(Int_t hBin=0;hBin<listSize; hBin++){
        setHistLooks(lHists[hBin], styleList[hBin], colorList[hBin]);
        lHists[hBin]->Draw(drawOption.Data());
    }
    hBlank->Draw("same");
    legALICE->Draw("same");
    legJet->Draw("same");
    leg->Draw("same");

    TString oFileName = outputDir + "canv" + label + ".root";
    c->SaveAs(oFileName.Data());
    c->Close();
}
// # = e = Overwright Some Plots on One canvas ####################################


// # = s = Overwright Some Plots on One canvas ####################################
void overwrightSomePlotsWithErr(TH1D* hBlank,TH1D* lHists[5], TGraphErrors* lHErr[5], \
    Int_t listSize, TLegend* legALICE,TLegend* legJet,TLegend* leg,\
        Int_t colorList[5], Int_t styleList[5], TString label, TString addDrawOption,\
            TString outputDir){
    TString canvasName = "canvas_" + label;
    TCanvas* c = new TCanvas(canvasName.Data(),canvasName.Data(),800,800);
    c->cd();
    TString padName = "pad_" + label;
    TPad *pad = new TPad(padName.Data(), padName.Data(), 0, 0.3, 1, 1.0);
    setupPadMargin(pad);
    // # pad.SetLogy()
    pad->Draw();
    pad->cd();
    
    TString drawOption = addDrawOption + "same";
    for(Int_t hBin=0;hBin<listSize;hBin++){
        setHistLooks(lHists[hBin], styleList[hBin], colorList[hBin]);
        lHists[hBin]->SetMarkerSize(0.8);
        lHists[hBin]->SetMarkerStyle(20);
        lHists[hBin]->SetMarkerColor(colorList[hBin]);
        lHists[hBin]->Draw(drawOption.Data());
        lHErr[hBin]->SetLineColor(colorList[hBin]);
        // lHErr[hBin]->SetFillColorAlpha(colorList[hBin], 0.8);
        lHErr[hBin]->SetFillColor(colorList[hBin]);
        lHErr[hBin]->SetMarkerColor(colorList[hBin]);
        lHErr[hBin]->SetMarkerSize(0.8);
        lHErr[hBin]->SetMarkerStyle(20);
        lHErr[hBin]->Draw("epz 5 same");
    }
        
    hBlank->Draw("same");
    legALICE->Draw("same");
    legJet->Draw("same");
    leg->Draw("same");

    TString oFileName = outputDir + "canv" + label + ".root";
    c->SaveAs(oFileName);
    c->Close();
}
// # = e = Overwright Some Plots on One canvas ####################################


// # == s == Set legend parameters  #############################################
void setupPadMargin(TPad *pad){
    pad->SetBottomMargin(0.1);
    pad->SetLeftMargin(0.15);
    pad->SetRightMargin(0.15);
    pad->SetTopMargin(0.08);
}
// # == s == Set legend parameters  #############################################

// # == s == Make Some Hists        #############################################
void addSomeHistsLegend(TLegend *leg, std::vector<TH1D*> lHists, std::vector<TString> lHistTags, Double_t textSize){
    for(Int_t histKinds=0; histKinds<lHists.size(); histKinds++){
        leg->AddEntry(lHists.at(histKinds), lHistTags.at(histKinds), "plf");
    }

    leg->SetFillColorAlpha(0, 0.0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);
    leg->SetMargin(0.25);
    leg->SetTextSize(textSize);
    leg->SetEntrySeparation(0.5);
}
// # == e == Make Some Hists        #############################################


// # == s == Make Some Hists        #############################################
void addSomeHistsLegend(TLegend *leg, TH1D* lHists[5], TString lHistTags[5], \
    Int_t listSize, Double_t textSize){
    for(Int_t histKinds=0; histKinds<listSize; histKinds++){
        leg->AddEntry(lHists[histKinds], lHistTags[histKinds].Data(), "plf");
    }

    leg->SetFillColorAlpha(0, 0.0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);
    leg->SetMargin(0.25);
    leg->SetTextSize(textSize);
    leg->SetEntrySeparation(0.5);
}
// # == e == Make Some Hists        #############################################

// # == s == Make ALICE Legend       #############################################
void addALICELegend(TH1D* blankObj, TLegend* leg, \
    TString centLabel, TString addLabel, Double_t textSize){
    if(addLabel.Data()!="") leg->AddEntry(blankObj, addLabel, "");
    leg->AddEntry(blankObj, "ALICE Pb-Pb #sqrt{#it{s}_{NN}}=5.02 TeV/#it{c}", "");
    TString label = TString::Format("Centrality: %s%", centLabel.Data());
    if(centLabel.Data()!="") leg->AddEntry(blankObj, label.Data(), "");

    leg->SetFillColorAlpha(0, 0.0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);
    leg->SetMargin(0.25);
    leg->SetTextSize(textSize);
    leg->SetEntrySeparation(0.05);
}
// # == e == Make ALICE Legend       #############################################


// # == s == Make Jet Legend       #############################################
void addJetLegend(TH1D *blankObj, TLegend* leg, Double_t Reso, Int_t ptCut, Double_t textSize){
    TString label = "";
    leg->AddEntry(blankObj, "Anti-#it{k}_{T}, charged jet", "");
    label = TString::Format("#it{R} = %f, |#it{#eta}_{jet} < %f",Reso,0.9-Reso);
    leg->AddEntry(blankObj,label.Data(), "");
    label = TString::Format("#it{p}_{T}^{leading track} > %d GeV/#it{c} ",ptCut);
    leg->AddEntry(blankObj, label.Data(), "");

    leg->SetFillColorAlpha(0, 0.0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);
    leg->SetMargin(0.25);
    leg->SetTextSize(textSize);
    leg->SetEntrySeparation(0.05);
}
// # == e == Make ALICE Legend       #############################################

// # == s == Set legend parameters  #############################################
void setupLegend(TLegend *leg, Double_t textSize){
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);
    leg->SetMargin(0.25);
    leg->SetTextSize(textSize);
    leg->SetEntrySeparation(0.5);
}
// # == s == Set legend parameters  #############################################

// # == e == Set histgram looks  #################################################
void setHistLooks(TH1D* h, Int_t style, Int_t color){
    h->SetLineColor(color);
    h->SetLineWidth(4);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(style);
    h->SetMarkerSize(1.5);
}
// # == e == Set histgram looks  #################################################




