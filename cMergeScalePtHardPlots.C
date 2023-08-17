
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
#include <TProfile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TString.h>

#include "TGrid.h"


// import histSetting
// import plotPerformanceHists
// import PtRangeList

Int_t lFillMerker[18] = {20, 21, 22, 23, 34, 29, 33, 31, 24, 24, 26, 32, 20, 21, 22, 23, 34, 29};
Int_t lNoFillMerker[10] = {24, 25, 26, 32, 28, 30, 27, 24, 25, 26};

// #[kBrack, kRed, kBlue, kGreen+3, kOrange+8, kAzure+10, kSpring+4, kViolet+1, kPink+9, kSpring-9, kOrange+8, kAzure+1, kTeal+4, kRed-7, kBlue-7, kGreen-2, kOrange-2, kViolet+7, kSpring+4, kTeal-5]
Int_t lFillColor[21] = {1, 632, 600, 416+3, 800+8, 860+10, 900+10, 820+4, 880+1, 900+9, 820-9, 800+1, 860+1, 840+4, 632-7, 600-7, 416-2, 800-2, 880+7, 820+4, 840-5};

Bool_t plotQAHists = false; // true, false

// Int_t numOfCentBin = 10;
// Double_t centRangeList[10][2] = {{0,5},{5,10},{10,30},{30,50},{50,90}};
Int_t numOfCentBin = 4;
Int_t iniCentBin = 3;
Double_t centRangeList[5][2] = {{0,5},{5,10},{10,30},{30,50},{50,90}};

Double_t ptRangeDictGen[2] = {15.,200.};
Double_t ptRangeDictDet[2] = {20.,150.};
Double_t ptRangeDictReport[2] = {40.,120.};

TString lEPLabel[3] = {"OutOfPlane", "InPlane", "Inclusive"};

TString LHCPeriod = "LHC18q";
Int_t mCentBin = 0;
Int_t leadingTrackPtCut = 5; // 0, 5, 7 GeV/c
// // 'BKGNoFit', 'BKGV2', 'V0C', 'V0A'
TString diffSys = ""; // 98%: "", 94%: "TrackEff094"

Int_t PtHardBins = 21;
// Int_t PtHardBins = 2;


Int_t variation = 1;
Double_t minPtGen = 10.;
Double_t maxPtGen = 200.;
Double_t minPtDet = 20.; 
Double_t maxPtDet = 150.;

//LHC18q
// Double_t lrefPtHardScalFactor[20] = {2.0890671432532143e-07, 6.483648558567263e-08, \
// 2.942434838768743e-08, 1.0116785397551089e-08, 3.4497027354228536e-09, 3.5551552007287986e-08, \
// 4.308919172443723e-10, 1.3882748477994637e-10, 1.6081466071541884e-09, 1.8538095794454472e-11, \
// 6.760726612338514e-12, 2.3879130596397466e-12, 1.140758647139845e-12, 5.295523483080566e-13, \
// 2.593055941557842e-13, 2.569815638073137e-13, 7.278621817851156e-14, 1.4427633953104294e-14, \
// 4.122334085737104e-14, 1.1094286846225388e-13};


//LHC18r
// Double_t lScaleFactor[20] = \
//     {2.456280125649383e-07, 7.61698249503861e-08, 3.22581122679987e-08, 1.1812252486320773e-08,\
//         4.049254578128142e-09, 1.6319721765403964e-09, 4.65763556935813e-10, 1.6123117645106004e-10, \
//         6.392822378921412e-11, 2.2751243702969846e-11, 8.132893034348747e-12, 2.9160375006758812e-12,\
//         1.483583932062433e-12, 7.274705162337206e-13, 3.486126907680189e-13, 1.1722360279760496e-12, \
//         9.228858393189566e-14, 4.5892010095811983e-14, 2.908187970843322e-14, 2.538214345770791e-13};
// Double_t lScaleFactor[20] = \
//     {2.446091855262593e-07,7.314975239763274e-08,3.362231880594244e-08,1.1895882634924634e-08,\
//     4.040376437588744e-09, 8.080455194155069e-09, 4.750910913607525e-10, 1.566097040702298e-10,\
//     3.4940423913705284e-10, 2.118136932704947e-11, 8.003818748583827e-12, 2.8357605981006887e-12,\
//     1.3703850028545194e-12, 6.368826505468165e-13, 3.1137657459266875e-13, 2.221801809420425e-13,\
//     8.647496708263527e-14, 2.7811870500499066e-14, 3.492944863521209e-14, 6.539761692063409e-14};
// Double_t lPtHardEventNum[20] = \
//     {298012544.0, 286666656.0, 300784576.0, 299847456.0, 298053152.0, 285056032.0,\
//         292419200.0, 288783872.0, 293465568.0, 284348704.0, 297642656.0, 296179264.0, \
//             290634720.0, 283121344.0, 285828544.0, 114284304.0, 290955040.0, 298294720.0,\
//                 275080032.0, 104309320.0};


Double_t lScaleFactor[20] = {};
Double_t lPtHardEventNum[20] = {};

Double_t ptHardL[20] = {5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235};
Double_t ptHardH[20] = {7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, -1};

void SetScaleFactor(Double_t* lScaleFactor, Double_t* lPtHardEventNum, TString LHCPeriod, TString diffSys);
Int_t GetNEvents(TString inEmbFileBaseName, TString taskName, Int_t bin, TH1D *hNEvents, Bool_t bAccEvents);
// Int_t GetNEvents(TString inEmbFileBaseName, TString taskName, Int_t bin, Bool_t bAccEvents);
void eachPtHardBinScaleFact(TString inEmbFileBaseName, TString taskName, Int_t PtHardBins, \
    TH1D *hNEventsAcc, TH1D *hNEventsTot, Int_t nEventsAccAvg, Int_t nNEventsTotAvg,\
        Double_t *lScaleFactor, Double_t *lPtHardEventNum, TList *lMainTree);
void ExtractDataRawJet(TString inRawJetFile, Int_t leadingTrackPtCut, TString diffSys, Int_t mCentBin, TList *lMainTree);
void ScaleJetHists(TList * lMainTree, Double_t lScaleFactor[20], Int_t PtHardBins, TString inEmbFileBaseName, TString jetTaskName);
void PlotPerformanceHistEachCent(TList *lMainTree);

Double_t movePhaseShift(Double_t angle);

TProfile* getJER1(TProfile* hJESshift, TString name);
TProfile* getJER2(TH2D* histResponseMatrix, TString name);

void addHistIntoList(TList *lTree, TH1D* inputHist, TString histName, Double_t scaleFactor, Int_t ptHardBin);
void addHistIntoList(TList *lTree, TH2D* inputHist, TString histName, Double_t scaleFactor, Int_t ptHardBin);
void addHistIntoList(TList *lTree, TH3D* inputHist, TString histName, Double_t scaleFactor, Int_t ptHardBin);
void addHistIntoList(TList *lTree, TProfile* inputHist, TString histName, Double_t scaleFactor, Int_t ptHardBin);
void addHistIntoList(TList *lTree, THnSparse* inputHist, TString histName, Double_t scaleFactor, Int_t ptHardBin);

void OFileStracture(TList *lMainTree, Int_t numOfCentBin);

void histLabelSetting(TH1D *outputHist, TString histTitle, TString xTitle, TString yTitle, \
    Int_t histNumber, Int_t fillType, Double_t Rebin);
void histLabelSetting(TH2D *outputHist, TString histTitle, TString xTitle, TString yTitle, \
    Int_t histNumber, Int_t fillType, Double_t Rebin);
void histLabelSetting(TH3D *outputHist, TString histTitle, TString xTitle, TString yTitle, \
    Int_t histNumber, Int_t fillType, Double_t Rebin);
void histLabelSetting(TProfile *outputHist, TString histTitle,TString xTitle,TString yTitle, \
    Int_t histNumber, Int_t fillType, Double_t Rebin);

// TString inRawJetFile = "~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/RawJet/TrainAnaResTreeBKGWay.root";
// TString inEmbFileDir = "./";
// TString inEmbFileDir = "/Volumes/KumaSSD/TrainOutput/LHC18q/Embedding/8870/";

TString inEmbFileDir = "~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Embedding/Train/";

TString inEmbFileBaseName = "";
// TString outEmbFileDir = "./";
TString outEmbFileDir = "~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Embedding/";


// ###################################################################################
// # Main function
int main(int argc, char *argv[]){
    std::cout << "input : " << argv[1] << ", "<< argv[2] << ", "<< argv[3] << ", "<< argv[4] << std::endl;
    
    mCentBin = atoi(argv[1]);
    LHCPeriod = argv[2];
    leadingTrackPtCut = atoi(argv[3]);
    diffSys = argv[4];

    TString inRawJetFile = TString::Format("~/ALICE/cernbox/SWAN_projects/outputFiles/%s/pass3/Ch/RawJet/", LHCPeriod.Data());
    if((diffSys != "V0M") && (diffSys != "V0C") && (diffSys != "V0A")) {
        inRawJetFile += "TrainAnaResTreeBKGWay.root";
    }else  inRawJetFile += "TrainAnaResTreeBKGWayDiffV0.root";
    
    if((diffSys != "V0M") && (diffSys != "V0C") && (diffSys != "V0A")) {
        inEmbFileDir += "8870/";
    }else inEmbFileDir += "8922/";
    inEmbFileBaseName = inEmbFileDir + "AnalysisResultsPtHard";

    TString outEmbFile = TString::Format("%sEmbedPtHardScaledResults_TrackPtCut%d_%s_CentBin%d.root"\
        , outEmbFileDir.Data(), leadingTrackPtCut, diffSys.Data(), mCentBin);

    TString embHelperTaskName = "AliAnalysisTaskEmcalEmbeddingHelper_histos";
    TString tempDiffSys = diffSys;
    if((diffSys == "Norm")||(diffSys == "V0M")) tempDiffSys = "";
    TString jetTaskName = TString::Format("AliAnalysisTaskEmbeddingJetWithEP_R02PtCut%d%s_histos", leadingTrackPtCut, tempDiffSys.Data());
    
    // # Create histogram of NEvents accepted and NEvents acc+rej, as a function of pT-hard bin
    TH1D *hNEventsAcc = new TH1D("hNEventsAcc", "hNEventsAccepted", PtHardBins+1, 0, PtHardBins+1);
    TH1D *hNEventsTot = new TH1D("hNEventsTot", "hNEventsTotal", PtHardBins+1, 0, PtHardBins+1);
    Int_t nEventsAccSum = 0;
    Int_t nEventsTotSum = 0;
    for (Int_t bin =1;bin< PtHardBins; bin++){
        hNEventsAcc->GetXaxis()->SetBinLabel(bin, TString::Format("%f-%f", ptHardL[bin-1], ptHardH[bin-1]));
        hNEventsTot->GetXaxis()->SetBinLabel(bin, TString::Format("%f-%f", ptHardL[bin-1], ptHardH[bin-1]));
    }

    for (Int_t bin =1;bin< PtHardBins; bin++){
        Int_t nNEventsTot = GetNEvents(inEmbFileBaseName, embHelperTaskName,\
            bin, hNEventsTot, false);
        Int_t nEventsAcc  = GetNEvents(inEmbFileBaseName, embHelperTaskName,\
            bin, hNEventsAcc, true);
        std::cout << "nEventsAcc ******* = " << nEventsAcc << std::endl;
        
        nEventsAccSum += nEventsAcc;
        nEventsTotSum += nNEventsTot;
    }
    std::cout << "nEventsAccSum ******* = " << nEventsAccSum << std::endl;
    Double_t nEventsAccAvg = nEventsAccSum/PtHardBins;
    Double_t nNEventsTotAvg= nEventsTotSum/PtHardBins;
    
    Double_t lScaleFactor[20];
    Double_t lPtHardEventNum[20];

    TFile *outRootFile = new TFile(outEmbFile, "RECREATE");
    TList *lMainTree = new TList();
    OFileStracture(lMainTree, 5);
    
    eachPtHardBinScaleFact(inEmbFileBaseName, embHelperTaskName, PtHardBins, \
        hNEventsAcc, hNEventsTot, nEventsAccAvg, nNEventsTotAvg,\
            lScaleFactor, lPtHardEventNum, lMainTree);
    
    ExtractDataRawJet(inRawJetFile, leadingTrackPtCut, diffSys, mCentBin, lMainTree);
    
    SetScaleFactor(lScaleFactor, lPtHardEventNum, LHCPeriod, diffSys);
    std::cout << "lScaleFactor = [" << "\n";
    for (int i = 0; i < 20; i++) std::cout << lScaleFactor[i] << ", " << "\n"; 
    ScaleJetHists(lMainTree, lScaleFactor, PtHardBins, inEmbFileBaseName, jetTaskName);
    
    outRootFile->cd();
    lMainTree->Write("mainTree", 1);
    outRootFile->Close();
    std::cout << "root " << outEmbFile << std::endl;
}

// ###################################################################################
void ExtractDataRawJet(TString inRawJetFile, Int_t leadingTrackPtCut, TString diffSys, Int_t mCentBin, TList *lMainTree){
    TString histName = ""; 
    TString histTitle = "";
    TString listName = "";

    TFile *f = new TFile(inRawJetFile.Data(), "READ");
    
    TString diffSysName = "";    
    if((diffSys == "Norm")||(diffSys == "V0M")) diffSysName = "";
    else if((diffSys == "V0C") || (diffSys == "V0A")) diffSysName = diffSys;
    else if(diffSys == "BKGNoFit") diffSysName = "kNoFit";
    else if(diffSys == "BKGV2") diffSysName = "kV2";
    TString mainJetTaskName \
        = TString::Format("AliAnalysisTaskRawJetWithEP_R02PtCut%d%s_histos",\
            leadingTrackPtCut,diffSysName.Data());
    TList *lMainHJetDists = (TList *) f->Get(mainJetTaskName.Data());
    
    TString jetTaskName = "Jet_AKTChargedR020_tracks_pT0150_pt_scheme";
    TList *lHJetDists = (TList *) lMainHJetDists->FindObject(jetTaskName.Data());
    f->Close();

    if(!lHJetDists) std::cout << "ERROR no " << jetTaskName << " found" << std::endl;

    TH1D *fNormalisationHist = (TH1D *) lMainHJetDists->FindObject("fNormalisationHist");
    lMainTree->Add(fNormalisationHist);
    TH1D *hCentrality_raw = (TH1D *) lMainHJetDists->FindObject("Centrality_raw");
    lMainTree->Add(hCentrality_raw);

    TList *lHJetDistsInc = (TList *)lHJetDists->FindObject("Inclusive");
    if(!lHJetDistsInc) std::cout << "ERROR no Inclusive found" << std::endl;
    
    TList* tempLMainInc = (TList *) lMainTree->FindObject("Inclusive");
    TString jetHistName = "hJetCorrPtVsEP2AngleVsCent";
    TH3D *hRawJet3D = (TH3D *)lHJetDistsInc->FindObject(jetHistName.Data());
    if(!hRawJet3D) std::cout << "ERROR no hJetCorrPtVsEP2AngleVsCent found" << std::endl;
    
    // ## S ## Centrality Loop #################
    // for(Int_t mCentBin=iniCentBin; mCentBin<numOfCentBin; mCentBin++){
    listName = TString::Format("lCent%d",mCentBin);
    TList* tempLMainIncCent = (TList *) tempLMainInc->FindObject(listName.Data());
    histName = TString::Format("hRawJet3D_Cent%d",mCentBin);
    TH3D *hRawJet3DCP = (TH3D *) hRawJet3D->Clone(histName.Data());
    hRawJet3DCP->GetZaxis()->SetRangeUser(centRangeList[mCentBin][0], centRangeList[mCentBin][1]);
    TH2D *hRawJet2D = (TH2D *) hRawJet3DCP->Project3D("xy");
    
    hRawJet2D->SaveAs("hRawJet2D.root");
    histName = TString::Format("hRawJet2D_Cent%d",mCentBin);
    histTitle = TString::Format("hRawJet2D_Cent%d; #phi^{jet}-#Psi_{EP, 2}; #it{p}_{T,jet}^{corr} [GeV/#it{c}];",mCentBin);
    TH2D* hRawJet2DShiftPhase = new TH2D(histName, histTitle,  6, 0., TMath::Pi()/2, 300, -50., 250.);
    
    Int_t origAllNBin = hRawJet2D->GetNcells();
    for(Int_t allBin=1; allBin<origAllNBin+1; allBin++){
        Int_t phiBin = 0;
        Int_t ptBin = 0;
        Int_t zBin = 0;
        hRawJet2D->GetBinXYZ(allBin, phiBin, ptBin, zBin);
        // Double_t fillContent = hRawJet2D->GetBinContent(phiBin, ptBin);
        Double_t fillContent = hRawJet2D->GetBinContent(allBin);
        Double_t jetAngle = hRawJet2D->GetXaxis()->GetBinCenter(phiBin);
        Double_t jetPt = hRawJet2D->GetYaxis()->GetBinCenter(ptBin);
        // std::cout << "content:: phi, pt = " << fillContent << ", " << jetAngle  << ", " << jetPt << std::endl;
        jetAngle = movePhaseShift(jetAngle);
        hRawJet2DShiftPhase->Fill(jetAngle, jetPt, fillContent);
    }
    
    tempLMainIncCent->Add(hRawJet2DShiftPhase);
    // } // ## E ## Centrality Loop #################
    
    // for(Int_t epBin=0; epBin < 2; epBin++){
    //     TList* tempLMainEP = (TList *) lMainTree->FindObject(lEPLabel[epBin].Data());
    //     for(Int_t centBin=iniCentBin; centBin<numOfCentBin; centBin++){
    //         listName = TString::Format("lCent%d",centBin);
    //         TList* tempLMainEPCent = (TList *) tempLMainEP->FindObject(listName.Data());
    // 
    //         listName = lEPLabel[epBin];
    //         // lHJetDistsWEP = lHJetDists->FindObject(listName);
    //         // jetHistName = "hJetCorrPtLocal_" + str(centBin);
    //         // hRawJet = lHJetDistsWEP->FindObject(jetHistName);
    //         // tempLMainEPCent->Add(hRawJet);
    // 
    //         TString histName = TString::Format("hRawJet3D_Cent%d_EP%d",centBin,epBin);
    //         TH3D *hRawJet3DCP = (TH3D *) hRawJet3D->Clone(histName.Data());
    //         hRawJet3DCP->GetZaxis()->SetRangeUser(centRangeList[centBin][0], centRangeList[centBin][1]);
    //         TH2D *hRawJet2D = (TH2D *) hRawJet3DCP->Project3D("xy");
    // 
    //         TString histName = TString::Format("hRawJet3D_Cent%d_EP%d",centBin,epBin);
    //         TH2D* hRawJet2DShiftPhase = new TH2D();
    //         tempLMainEPCent->Add(hRawJet2D);
    // 
            // if(epBin==0){
            //     histName = TString::Format("hRawJet2DCP1_Cent%d_EP%d", centBin, epBin);
            //     TH2D *hRawJet2DCP1 = (TH2D *) hRawJet2D->Clone(histName.Data());
            //     hRawJet2DCP1->GetXaxis()->SetRangeUser(0, TMath::Pi()/4);
            //     TH1D *hRawJet1 = hRawJet2DCP1->ProjectionY();
            //     histName = TString::Format("hRawJet2DCP2_Cent%d_EP%d", centBin, epBin);
            //     TH2D* hRawJet2DCP2 = (TH2D* )hRawJet2D->Clone(histName.Data());
            //     hRawJet2DCP2->GetXaxis()->SetRangeUser(TMath::Pi()*3/4, TMath::Pi()*5/4);
            //     TH1D* hRawJet2 = hRawJet2DCP2->ProjectionY();
            //     histName = TString::Format("hRawJet2DCP3_Cent%d_EP%d", centBin, epBin);
            //     TH2D* hRawJet2DCP3 = (TH2D* )hRawJet2D->Clone(histName.Data());
            //     hRawJet2DCP3->GetXaxis()->SetRangeUser(TMath::Pi()*7/4, TMath::Pi()*2);
            //     TH1D* hRawJet3 = hRawJet2DCP3->ProjectionY();
            //
            //     hRawJet1->Add(hRawJet2);
            //     hRawJet1->Add(hRawJet3);
            //
            //     histName = TString::Format("hJetCorrPtLocal_%d",centBin);
            //     hRawJet1->SetName(histName.Data());
            //
            //     tempLMainEPCent->Add(hRawJet1);
            // }else{
            //     histName = TString::Format("hRawJet2DCP1_Cent%d_EP%d", centBin, epBin);
            //     TH2D* hRawJet2DCP1 = (TH2D* ) hRawJet2D->Clone(histName.Data());
            //     TH1D* hRawJet1 = hRawJet2DCP1->ProjectionY();
            //     hRawJet2DCP1->GetXaxis()->SetRangeUser(TMath::Pi()/4, TMath::Pi()*3/4);
            //     histName = TString::Format("hRawJet2DCP2_Cent%d_EP%d", centBin, epBin);
            //     TH2D* hRawJet2DCP2 = (TH2D* )hRawJet2D->Clone(histName.Data());
            //     hRawJet2DCP2->GetXaxis()->SetRangeUser(TMath::Pi()*3/4, TMath::Pi()*7/4);
            //     TH1D* hRawJet2 = hRawJet2DCP2->ProjectionY();
            // 
            //     histName = TString::Format("hJetCorrPtLocal_%d",centBin);
            //     hRawJet1->SetName(histName.Data());
            //     hRawJet1->Add(hRawJet2);
            // 
            //     tempLMainEPCent->Add(hRawJet1);
            // }
        // }
    // }
}

// ###################################################################################
// # Given event list name eventListName, pT-hard bin number, and histogram hNEvents of appropriate form, 
// # this function fills the number of events 
// # (accepted events only if bAcceptedEventsOnly=True, otherwise all events)
Int_t GetNEvents(TString inEmbFileBaseName, TString taskName, Int_t bin, TH1D *hNEvents, Bool_t bAccEvents){
// Int_t GetNEvents(TString inEmbFileBaseName, TString taskName, Int_t bin, Bool_t bAccEvents){
    TString inEmbFile = TString::Format("%s%d.root", inEmbFileBaseName.Data(), bin);
    TFile *f = new TFile(inEmbFile.Data(), "READ");
    
    TList *lHEventInfo = (TList *) f->Get(taskName.Data());
    f->Close();
    
    Int_t nEvents = 0;
    if(!lHEventInfo) std::cout << "ERROR no lHEventInfo found" << std::endl;
    
    TList *eventCutList = (TList *) lHEventInfo->FindObject("EventCuts");
    
    Bool_t bCutHist = false;
    if(eventCutList) bCutHist = true;
    
    TH1D *hNEventsPtHard = new TH1D();
    Int_t histBin = 0;
    if(bCutHist){
        hNEventsPtHard = (TH1D *)eventCutList->FindObject("fCutStats");
        if(bAccEvents) histBin = hNEventsPtHard->GetXaxis()->FindBin("All cuts");
        else histBin = hNEventsPtHard->GetXaxis()->FindBin("No cuts");
    }
    else{
        hNEventsPtHard = (TH1D *)lHEventInfo->FindObject("fHistEventCount");
        histBin = 1;
    }
    
    nEvents = hNEventsPtHard->GetBinContent(histBin);
    if((!bAccEvents)&&(!bCutHist)) nEvents += hNEventsPtHard->GetBinContent(2);
    
    hNEvents->Fill(bin-0.5, nEvents);

    return nEvents;
}

void eachPtHardBinScaleFact(TString inEmbFileBaseName, TString taskName, Int_t PtHardBins, \
        TH1D *hNEventsAcc, TH1D *hNEventsTot, Int_t nEventsAccAvg, Int_t nNEventsTotAvg,\
        Double_t *lScaleFactor, Double_t *lPtHardEventNum, TList *lMainTree){
    std::cout << "ooo Determine Scale Factors from single pT Hard Bins" << std::endl;

    TH1D* hXSecPerEvent = new TH1D("hXSecPerEvent", "hXSecPerEvent", PtHardBins+1, 0, PtHardBins+1);
    TH1D* hNTrialsPerEvent = new TH1D("hNTrialsPerEvent", "hNTrialsPerEvent", PtHardBins+1, 0, PtHardBins+1);
    TH1D* hScaleFactor = new TH1D("hScaleFactor", "hScaleFactor", PtHardBins+1, 0, PtHardBins+1);
    TH1D* hPtHardSpectOrig = new TH1D();
    TH1D* hPtHardSpectScaled = new TH1D();

    for(Int_t bin=1; bin < PtHardBins; bin++){
        // # Label histograms
        TString histLabel = TString::Format("%f-%f",ptHardL[bin-1], ptHardH[bin-1]);
        hXSecPerEvent->GetXaxis()->SetBinLabel(bin, histLabel);
        hNTrialsPerEvent->GetXaxis()->SetBinLabel(bin, histLabel);
        hScaleFactor->GetXaxis()->SetBinLabel(bin, histLabel);
    }

    // # Extract cross sections from pT Hard bins
    for(Int_t bin=1;bin < PtHardBins; bin++){
        // # Open input file and get relevant lists
        TString inEmbFile = TString::Format("%s%d.root", inEmbFileBaseName.Data(), bin);
        TFile *f = new TFile(inEmbFile.Data(), "READ");
        std::cout << TString::Format("ooo Pt-hard bin %d", bin) << std::endl;
        std::cout << TString::Format("  o File: %s", inEmbFile.Data()) << std::endl;
        
        TList *lHEventInfo = (TList *)f->Get(taskName.Data());
        f->Close();
        
        std::cout << TString::Format("ooo Computing scaling factors with list: %s", lHEventInfo->GetName()) << std::endl;

        TH1D *hNEvent = (TH1D* ) lHEventInfo->FindObject("fHistEventCount");
        // # hXsecPtHard      = lHEventInfo->FindObject("hXsec")
        // # hTrialsPtHard    = lHEventInfo->FindObject("hNtrials")
        TH1D *hXsecPtHard   = (TH1D* ) lHEventInfo->FindObject("fHistXsection");
        TH1D *hTrialsPtHard = (TH1D* ) lHEventInfo->FindObject("fHistTrials");
        
        // # Compute: scale factor = xsec per event / trials per event
        Double_t nEvent = hNEvent->GetBinContent(1);
        Int_t nEventsTot = hNEventsTot->GetBinContent(bin);
        Int_t nEventsAcc = hNEventsAcc->GetBinContent(bin);

        // # xsec = hXsecPtHard->GetBinContent(1) / hXsecPtHard->GetEntries() 
        // # trials = 1.*hTrialsPtHard->GetBinContent(1) / nEventsTot
        
        std::cout << "entry = " << hXsecPtHard->GetEntries() << std::endl;
        Double_t xsec = hXsecPtHard->GetBinContent(bin+1);
        Double_t xsecEntry = hXsecPtHard->GetEntries();
        Double_t xsecScale = xsec * xsecEntry;
        // xsec/=xsecEntry;
        // Double_t trials = hTrialsPtHard->GetBinContent(bin+1) / nEventsTot;
        Double_t trials = hTrialsPtHard->GetBinContent(bin+1);
        Double_t scaleFactor = xsec/trials; // # Check ????
        // Double_t scaleFactor = lrefPtHardScalFactor[bin-1] / lrefPtHardNumOfEvent[bin-1];
        // # scaleFactor = 0.5
        // # scaleFactor = 1.

        // # also scale to account that there are different number of events in each Pt-hard bin
        Double_t eventScaleFactorAcc = 0;
        Double_t eventScaleFactorTot = 0;
        std::cout << "nEventsAccAvg, nEventsAcc = " << nEventsAccAvg<< ", "<<nEventsAcc << std::endl;
        if(nEventsAcc>0) eventScaleFactorAcc = nEventsAccAvg/nEventsAcc;
        else eventScaleFactorAcc = 0;
        if(nEventsTot>0) eventScaleFactorTot = nNEventsTotAvg/nEventsTot;
        else eventScaleFactorTot = 0;

        hXSecPerEvent->Fill(bin-0.5, xsec);
        hNTrialsPerEvent->Fill(bin-0.5, trials);
        
        // hScaleFactor->Fill(bin-0.5, eventScaleFactorAcc*scaleFactor);
        hScaleFactor->Fill(bin-0.5, lScaleFactor[bin-1]);//?????????????????????

        std::cout << TString::Format("xSec: %f, trials: %f, histVal:  %f, nEventsTot: %d, nEventsTot: %d, eventScalFactor: %f, scaleFactor: %f"\
            ,xsec, trials, hTrialsPtHard->GetBinContent(1), nEventsTot, nEventsAcc, eventScaleFactorAcc, scaleFactor) << std::endl;

        lScaleFactor[bin-1] = eventScaleFactorAcc*scaleFactor;
        // lScaleFactor[bin-1] = 1;
        lPtHardEventNum[bin-1] = nEventsTot;

        TH1D *origHPtHardSpectScaled = (TH1D *) lHEventInfo->FindObject("fHistPtHard");

        TH1D *cpHPtHardSpectOrig = (TH1D *) origHPtHardSpectScaled->Clone();
        cpHPtHardSpectOrig->SetName(TString::Format("hPtHardSpectOrig%d",bin));
        if(bin == 1) hPtHardSpectOrig = cpHPtHardSpectOrig;
        else hPtHardSpectOrig->Add(cpHPtHardSpectOrig);

        TH1D *cpHPtHardSpectScaled = (TH1D *) origHPtHardSpectScaled->Clone();
        cpHPtHardSpectScaled->SetName(TString::Format("hPtHardSpectScaled%d",bin));
        cpHPtHardSpectScaled->Scale(scaleFactor);
        if(bin == 1) hPtHardSpectScaled = cpHPtHardSpectScaled;
        else hPtHardSpectScaled->Add(cpHPtHardSpectScaled);
        
        f->Close();
    }
    
    lMainTree->Add(hXSecPerEvent);
    lMainTree->Add(hNTrialsPerEvent);
    lMainTree->Add(hScaleFactor);
    lMainTree->Add(hNEventsAcc);
    lMainTree->Add(hNEventsTot);
    lMainTree->Add(hPtHardSpectOrig);
    lMainTree->Add(hPtHardSpectScaled);
}

void ScaleJetHists(TList * lMainTree, Double_t lScaleFactor[20], Int_t PtHardBins, TString inEmbFileBaseName, TString jetTaskName){
    TString histName = "";
    TString listName = "";
    
    std::cout << "000000 s 00000 START Scaling 00000000000000000 " << std::endl;
    for(Int_t ptHardBin=1;ptHardBin< PtHardBins;ptHardBin++){
        std::cout << "  &&&& s &&&& START pT hard bin: " << ptHardBin << "   &&&& " << std::endl;
        TString inEmbFile = TString::Format("%s%d.root",inEmbFileBaseName.Data(),ptHardBin);
        TFile *f = new TFile(inEmbFile.Data(), "READ");
        TList *lHEventInfo = (TList *) f->Get(jetTaskName.Data());
        f->Close();
        
        if(!lHEventInfo) std::cout << "ERROR no " << jetTaskName << " found"<<std::endl;

        Double_t ptHardScaleFactor = lScaleFactor[ptHardBin-1];

        // == s ==  Excute 4D RM =============
        TString titleFourDRM = "hRMWithEP2Angle;#it{p}_{T}^{truth} (GeV/#it{c}); #it{p}_{T,corr}^{det} (GeV/#it{c}); #phi^{truth jet}-#Psi_{EP, 3}; #phi^{det jet}-#Psi_{EP, 3}; Centrality (%)";
        Int_t nbinsFourDRM[4]  = {50, 50, 6, 6};
        Double_t minValFourDRM[4] = {0., 0., 0., 0.};
        Double_t maxValFourDRM[4] = {250, 250, TMath::Pi()/2, TMath::Pi()/2};
        THnSparse* hRMWithEP2 = new THnSparseD("hRMWithEP2Angle", titleFourDRM, 4, nbinsFourDRM, minValFourDRM, maxValFourDRM);
        
        TList* tempLMainInc = (TList *) lMainTree->FindObject("Inclusive");
        // ## S ## Centrality Loop #################
        // for(Int_t mCentBin=iniCentBin; mCentBin < numOfCentBin; mCentBin++){
        std::cout<< "          ## s 11 START centrality : "<< mCentBin <<"  ##"<<std::endl;
        listName = TString::Format("lCent%d", mCentBin);
        TList* tempLMainIncCent = (TList *) tempLMainInc->FindObject(listName.Data());
        
        histName = "hRMWithEP2Angle";
        THnSparse *hTempMD = (THnSparse *) lHEventInfo->FindObject(histName.Data());
        
        Int_t nBins0 = hTempMD->GetAxis(0)->GetNbins();
        Int_t nBins1 = hTempMD->GetAxis(1)->GetNbins();
        Int_t nBins2 = hTempMD->GetAxis(2)->GetNbins();
        Int_t nBins3 = hTempMD->GetAxis(3)->GetNbins();
        Int_t nBins4 = hTempMD->GetAxis(4)->GetNbins();
        
        Int_t printCount = 0;
        Int_t binList[4] = {0,0,0,0};
        Double_t binValList[4] = {0.,0.,0.,0.};
        Double_t tempBinValList[5] = {0.,0.,0.,0.,0.};
        Int_t embCentBinNunRange[5][2] = {{1,1},{2,2},{3,6},{7,10},{11,17}};
        
        hTempMD->GetAxis(4)->SetRange(embCentBinNunRange[mCentBin][0], embCentBinNunRange[mCentBin][1]);
        Int_t dimX[4] = {0,1,2,3};
        THnSparse *hTemp4D = hTempMD->Projection (4, dimX);
        Long64_t allNBinId = hTemp4D->GetNbins();
        for(Long64_t allBin = 1; allBin < allNBinId+1; allBin++){
            Long64_t totalBinNum = allBin;
            Int_t binList2[4] = {0,0,0,0};
            Double_t fillBinContent = hTemp4D->GetBinContent(totalBinNum, binList2);
            Double_t fillBinErr = hTemp4D->GetBinError2(totalBinNum);
            
            Int_t binContent = hTemp4D->GetBinContent(binList);
            Double_t genJetPt = hTemp4D->GetAxis(0)->GetBinCenter(binList2[0]);
            Double_t detJetPt = hTemp4D->GetAxis(1)->GetBinCenter(binList2[1]);
            Double_t genJetAngle = hTemp4D->GetAxis(2)->GetBinCenter(binList2[2]);
            Double_t detJetAngle = hTemp4D->GetAxis(3)->GetBinCenter(binList2[3]);
            genJetAngle = movePhaseShift(genJetAngle);
            detJetAngle = movePhaseShift(detJetAngle);
            
            binValList[0] = genJetPt;
            binValList[1] = detJetPt;
            binValList[2] = genJetAngle;
            binValList[3] = detJetAngle;
            
            Int_t fillBinList[4] = {0,0,0,0};
            for(Int_t i = 0; i<4; i++) fillBinList[i] = binList2[i];
            Int_t fillTotalBinNum = hRMWithEP2->GetBin(binValList);
            
            if(genJetAngle!=TMath::Pi()/2 || detJetAngle!=TMath::Pi()/2){
                hRMWithEP2->Fill(binValList, fillBinContent);
            }
            
            
            // if(printCount%100010==0){
            // if(printCount>0){
            if(0){
                std::cout << "fill count2 : binerror :: binId : genPt : detPt : genPhi : detPhi = " << \
                fillBinContent << " : "<< fillBinErr << " :: " \
                << totalBinNum << " : " << genJetPt << " : "<< detJetPt << " : " \
                << genJetAngle << " : " << detJetAngle << std::endl;
            }
            printCount+= 1;
        }
        
        hRMWithEP2->Sumw2();
        hRMWithEP2->Scale(ptHardScaleFactor);
        if(ptHardBin == 1) tempLMainIncCent->Add(hRMWithEP2);
        else{
            THnSparse* tempHist = (THnSparse*) tempLMainIncCent->FindObject(histName.Data());
            tempHist->Add(hRMWithEP2);
        }
        
        // } // ## E ## Centrality Loop #################
        // == e ==  Excute 4D RM =============
        
        // # Look for the EventCutOutput from AliEventCuts, and if it doesn't exist, look for histo fHistEventCount
        // TH1D* lHybRawJet = (TH1D *)lHEventInfo->FindObject("hybridRawJet");
        // TH1D* lParRawJet = (TH1D *)lHEventInfo->FindObject("particleRawJet");
        
        // for(Int_t epBin=0;epBin<2;epBin++){
        //     std::cout << "      %%% s %%% START ep : " << lEPLabel[epBin] << "  %%% " << std::endl;
        //     TString epLabel = lEPLabel[epBin];
        //     TList* tempLMainEP = (TList *) lMainTree->FindObject(epLabel.Data());
        //     TH1D *lHybRawJetEP = (TH1D *)lHybRawJet->FindObject(epLabel.Data());
        //     TH1D *lParRawJetEP = (TH1D *)lParRawJet->FindObject(epLabel.Data());
            // 
        //     Double_t reBin = 1.;
        //     // for(Int_t centBin=iniCentBin;centBin<numOfCentBin;centBin++){
        //     //     std::cout<< "          ## s 11 START centrality : "<< centBin <<"  ##"<<std::endl;
        //     //     listName = TString::Format("lCent%d",centBin);
        //     //     TList* tempLMainEPCent = (TList *) tempLMainEP->FindObject(listName.Data());
        //     // 
        //     //     histName = TString::Format("hJetCorrPtLocal_%d",centBin);
        //     //     TH1D* hTempHybRawJetEP = (TH1D *)lHybRawJetEP->FindObject(histName.Data());
        //     //     TH1D* hHybRawJetEP = (TH1D *)hTempHybRawJetEP->Clone();
        //     //     histName = TString::Format("hHybJetCorrPtLocal_%d",centBin);
        //     //     hHybRawJetEP->SetName(histName.Data());
        //     //     hHybRawJetEP->GetXaxis()->SetRangeUser(0, 250);
        //     // 
        //     //     reBin = 1.;
        //     //     histLabelSetting(hHybRawJetEP, histName, \
        //     //         "#it{p}_{T, det}^{jet} [GeV/#it{c}]", "count", centBin, 1, reBin);
        //     //     addHistIntoList(tempLMainEPCent, hHybRawJetEP, histName, \
        //     //         ptHardScaleFactor, ptHardBin);
        //     // 
        //     //     histName = TString::Format("hJetCorrPtLocal_%d",centBin);
        //     //     TH1D *hTempParRawJetEP = (TH1D *)lParRawJetEP->FindObject(histName.Data());
        //     //     TH1D *hParRawJetEP = (TH1D *) hTempParRawJetEP->Clone();
        //     //     histName = TString::Format("hParJetCorrPtLocal_%d",centBin);
        //     //     hParRawJetEP->SetName(histName.Data());
        //     //     hParRawJetEP->GetXaxis()->SetRangeUser(0, 250);
        //     //     reBin = 1.;
        //     //     histLabelSetting(hParRawJetEP, histName, \
        //     //         "#it{p}_{T, gen}^{jet} [GeV/#it{c}]", "count", centBin, 1, reBin);
        //     //     addHistIntoList(tempLMainEPCent, hParRawJetEP, histName, \
        //     //         ptHardScaleFactor, ptHardBin);
        //     // }
            // 
            // 
        //     // TH1D* hTemp1D = new TH1D();
        //     // TH2D* hTemp2D = new TH2D();
        //     // TH3D* hTemp3D = new TH3D();
        //     // 
        //     // listName = "MatchedJetHisto_" + epLabel;
        //     // TList *lMatchingJet = (TList *)lHEventInfo->FindObject(listName.Data());
        //     // histName = "hResponseMatrixDiff";
        //     // THnSparse *hTempMD = (THnSparse *) lMatchingJet->FindObject(histName.Data());
        //     // hTempMD->Sumw2();
        //     // hTempMD->Scale(ptHardScaleFactor);
        //     // if(ptHardBin == 1) tempLMainEP->Add(hTempMD);
        //     // else{
        //     //     THnSparse* tempHist = (THnSparse*) tempLMainEP->FindObject(histName.Data());
        //     //     tempHist->Add(hTempMD);
        //     // }
        //     // 
        //     // if(plotQAHists){
        //     //     histName = "hJESshift";
        //     //     hTemp3D = (TH3D *) lMatchingJet->FindObject(histName.Data());
        //     //     addHistIntoList(tempLMainEP, hTemp3D, histName, 1, ptHardBin);
        //     //     histName = "hEmbDeltaPt";
        //     //     hTemp3D = (TH3D *) lMatchingJet->FindObject(histName.Data());
        //     //     addHistIntoList(tempLMainEP, hTemp3D, histName, 1, ptHardBin);
        //     //     histName = "hJESshiftHybDet";
        //     //     hTemp3D = (TH3D *) lMatchingJet->FindObject(histName.Data());
        //     //     addHistIntoList(tempLMainEP, hTemp3D, histName, 1, ptHardBin);
        //     //     histName = "hEmbDeltaPtHybDet";
        //     //     hTemp3D = (TH3D *) lMatchingJet->FindObject(histName.Data());
        //     //     addHistIntoList(tempLMainEP, hTemp3D, histName, 1, ptHardBin);
        //     //     histName = "hJESshiftDetPar";
        //     //     hTemp3D = (TH3D *) lMatchingJet->FindObject(histName.Data());
        //     //     addHistIntoList(tempLMainEP, hTemp3D, histName, 1, ptHardBin);
        //     //     histName = "hEmbDeltaPtDetPar";
        //     //     hTemp3D = (TH3D *) lMatchingJet->FindObject(histName.Data());
        //     //     addHistIntoList(tempLMainEP, hTemp3D, histName, 1, ptHardBin);
        //     // 
        //     //     histName = "hNEFVsPt";
        //     //     hTemp3D = (TH3D *) lMatchingJet->FindObject(histName.Data());
        //     //     addHistIntoList(tempLMainEP, hTemp3D, histName, ptHardScaleFactor, ptHardBin);
        //     //     histName = "hZLeadingVsPt";
        //     //     hTemp3D = (TH3D *) lMatchingJet->FindObject(histName.Data());
        //     //     addHistIntoList(tempLMainEP, hTemp3D, histName, ptHardScaleFactor, ptHardBin);
        //     //     histName = "hMatchingDistance";
        //     //     hTemp3D = (TH3D *) lMatchingJet->FindObject(histName.Data());
        //     //     addHistIntoList(tempLMainEP, hTemp3D, histName, ptHardScaleFactor, ptHardBin);
        //     //     // histName = "fHistJetMatchingQA"
        //     //     // # hTemp1D = (TH1D *) lMatchingJet->FindObject(histName.Data());
        //     //     // # addHistIntoList(tempLMainEP, hTemp, histName, 1, ptHardBin);
        //     // }
        // 
        // }
    }
}

void PlotPerformanceHistEachCent(TList *lMainTree){
    TString histName = "";
    TString listName = "";
    Double_t reBin = 1.;
    std::cout << "=== s === Add Performance Hists  #############################" << std::endl;
    THnSparse *baseMatchedJetHists = (THnSparse *)lMainTree->FindObject("hResponseMatrixDiff");
    // TH3D *baseHJESshift = (TH3D *) lMainTree->FindObject("hJESshift");
    // TH3D *baseHEmbDeltaPt = (TH3D *) lMainTree->FindObject("hEmbDeltaPt");
    // TH3D *baseHJESshiftHybDet = (TH3D *) lMainTree->FindObject("hJESshiftHybDet");
    // TH3D *baseHEmbDeltaPtHybDet = (TH3D *) lMainTree->FindObject("hEmbDeltaPtHybDet");
    // TH3D *baseHJESshiftDetPar = (TH3D *) lMainTree->FindObject("hJESshiftDetPar");
    // TH3D *baseHEmbDeltaPtDetPar = (TH3D *) lMainTree->FindObject("hEmbDeltaPtDetPar");
    // TH3D *baseHNEFVsPt = (TH3D *) lMainTree->FindObject("hNEFVsPt");
    // TH3D *baseHZLeadingVsPt = (TH3D *) lMainTree->FindObject("hZLeadingVsPt");
    // TH3D *baseHMatchingDistance = (TH3D *) lMainTree->FindObject("hMatchingDistance");

    for(Int_t centBin=iniCentBin;centBin< numOfCentBin; centBin++){
        listName = TString::Format("lCent%d",centBin);
        TList* tempLMainEPCent = (TList *) lMainTree->FindObject(listName.Data());

        // ########################################################################
        std::cout << "1. Set Up RM to Each Centrality for all histograms  #####" << std::endl;
        THnSparse *tempOriginMatchedJetHists = (THnSparse *)baseMatchedJetHists->Clone();
        histName = TString::Format("RMHistsCent%d",centBin);
        tempOriginMatchedJetHists->SetName(histName.Data());
        tempOriginMatchedJetHists->GetAxis(5)->SetRangeUser(\
            centRangeList[centBin][0], centRangeList[centBin][1]);
        // #######################################################################

        if(plotQAHists){
            // ###  Add Matched Gen Jet hist  #########################################
            std::cout << "2. Add Matched Particle Level Jet hist Centrality   #####" << std::endl;
            histName = TString::Format("hMatchGenJetPt_Cent%d",centBin);
            THnSparse* hTempMatchedGenJet = (THnSparse* )tempOriginMatchedJetHists->Clone();
            TH1D* hMatchGenJetPt = (TH1D* )hTempMatchedGenJet->Projection(0, histName.Data());
            hMatchGenJetPt->SetName(histName.Data());
            reBin = 1;
            histLabelSetting(hMatchGenJetPt, histName, \
                "#it{p}_{T, gen}^{Matched jet} [GeV/#it{c}]", "count", centBin, 1, reBin);
            hMatchGenJetPt->Draw("AP");
            tempLMainEPCent->Add(hMatchGenJetPt);
            // ########################################################################

            // ###  Add Matched Det Jet hist  #########################################
            std::cout << "3. Add Matched Detector Level Jet hist Centrality   #####" << std::endl;
            histName = TString::Format("hMatchDetJetPt_Cent%d",centBin);
            THnSparse* hTempMatchedDetJet = (THnSparse* )tempOriginMatchedJetHists->Clone();
            TH1D *hMatchDetJetPt = (TH1D *) hTempMatchedDetJet->Projection(1, histName.Data());
            hMatchDetJetPt->SetName(histName.Data());
            reBin = 1;
            histLabelSetting(hMatchDetJetPt, histName, \
                "#it{p}_{T, det}^{Matched jet} [GeV/#it{c}]", "count", centBin, 1, reBin);
            hMatchDetJetPt->Draw("AP");
            tempLMainEPCent->Add(hMatchDetJetPt);
            // ##############################################################################

            // ########################################################################
            std::cout << "4. Add Matched Jets Delta R hist Centrality       #####" << std::endl;
            histName = TString::Format("hMatchJetDeltaR_Cent%d",centBin);
            THnSparse* hTempMatchJetDeltaR = (THnSparse* )tempOriginMatchedJetHists->Clone();
            TH1D *hMatchJetDeltaR = (TH1D *)hTempMatchJetDeltaR->Projection(2, histName);
            hMatchJetDeltaR->SetName(histName.Data());
            histLabelSetting(hMatchJetDeltaR, histName, "#Delta R", "count", centBin, 1, 1);
            hMatchJetDeltaR->Draw("AP");
            tempLMainEPCent->Add(hMatchJetDeltaR);
            // ########################################################################

            // ########################################################################
            std::cout << "5. Add Matched Jets Delta Angle hist Centrality   #####" << std::endl;
            histName = TString::Format("hMatchJetDAngle_Cent%d",centBin);
            THnSparse* hTemphMatchJetDAngle = (THnSparse* ) tempOriginMatchedJetHists->Clone();
            TH1D *hMatchJetDAngle = (TH1D *) hTemphMatchJetDAngle->Projection(4, histName);
            hMatchJetDAngle->SetName(histName.Data());
            histLabelSetting(hMatchJetDAngle, histName, "angle", "count", centBin, 1, 1);
            hMatchJetDAngle->Draw("AP");
            tempLMainEPCent->Add(hMatchJetDAngle);
            // ########################################################################
        }

        // ###  Add RM hist  ######################################################
        std::cout << "6. Add Original RM                   #####" << std::endl;
        histName = TString::Format("hRM_Cent%d",centBin);
        THnSparse* hTempEachCentRM = (THnSparse* ) tempOriginMatchedJetHists->Clone();
        TH2D *hEachCentRM = (TH2D *) hTempEachCentRM->Projection(0, 1, "");
        hEachCentRM->SetName(histName.Data());
        tempLMainEPCent->Add(hEachCentRM);
        // ########################################################################

        // ###  Add RM hist  ###################################################
        std::cout << "7. Add rebin and normalize RM        #####" << std::endl;
        histName = TString::Format("hRM_forUnfold_Cent%d",centBin);
        THnSparse* hTempEachCentRM_forUnfold = (THnSparse* )tempOriginMatchedJetHists->Clone();
        TH2D *hEachCentRM_forUnfold = hTempEachCentRM_forUnfold->Projection(0, 1, "");
        hEachCentRM_forUnfold->SetName(histName.Data());
        tempLMainEPCent->Add(hEachCentRM_forUnfold);
        // ########################################################################

        // if(plotQAHists){
        //     // ###  Add hJESshiftEMCal hist   #########################################
        //     std::cout << "8. Add Jet Energy Scale Shift hists  #####" << std::endl;
        //     Double_t lJESshiftPtRange[3][2] = {{20,30}, {50,70}, {100,120}};
        //     TList *lJESshiftDistHists = new TList();
        //     lJESshiftDistHists->SetName("hJESshiftList");
        // 
        //     histName = TString::Format("hJESshift_Cent%d",centBin);
        //     TH3D* hJESshift = (TH3D* )baseHJESshift->Clone();
        //     hJESshift->SetName(histName.Data());
        //     hJESshift->GetXaxis()->SetRangeUser(centRangeList[centBin][0], centRangeList[centBin][1]);
        // 
        //     histName = TString::Format("hJESshiftForProj_Cent%d",centBin);
        //     TH3D* hJESshift_clone = (TH3D* )hJESshift->Clone();
        //     hJESshift_clone->SetName(histName.Data());
        //     TH2D *hJESshift_proj = (TH2D* ) hJESshift_clone->Project3D("zyeo");
        //     TProfile *hJESshift_prof = (TProfile *)hJESshift_proj->ProfileX();
        // 
        //     TString xTitle = "#it{p}_{T, particle}^{jet}";
        //     TString yTitle = "#frac{#it{p}_{T}^{det} - #it{p}_{T}^{particle}}{#it{p}_{T}^{particle}}";
        //     histLabelSetting(hJESshift_prof, histName, xTitle, yTitle, centBin, 1, reBin);
        //     tempLMainEPCent->Add(hJESshift_prof);
        // 
        //     std::cout << "8.1 Add Jet Energy Scale Shift projected hists for three jet pt range   #####" << std::endl;
        //     Int_t lPtRangeColor[3] = {1, 632, 600};
        //     for(Int_t ptRangeKind=0;ptRangeKind<2;ptRangeKind++){
        //         histName = TString::Format("hJESshiftDist_Cent%d_ptRange%d",centBin,ptRangeKind);
        //         TH3D* hJESshif_clone = (TH3D* ) hJESshift->Clone();
        //         hJESshif_clone->SetName(histName.Data());
        //         hJESshif_clone->GetYaxis()->SetRangeUser(\
        //             lJESshiftPtRange[ptRangeKind][0], lJESshiftPtRange[ptRangeKind][1]);
        //         TH1D* hJESshif_clone_proj = (TH1D *) hJESshif_clone->Project3D("ze");
        // 
        //         xTitle = "#frac{#it{p}_{T}^{det} - #it{p}_{T}^{gen}}{#it{p}_{T}^{gen}}";
        //         yTitle = "probability";
        //         histLabelSetting(hJESshif_clone_proj, histName, xTitle, yTitle, \
        //             ptRangeKind, 1, reBin);      
        //         hJESshif_clone_proj->SetLineColor(lPtRangeColor[ptRangeKind]);
        //         hJESshif_clone_proj->SetMarkerColor(lPtRangeColor[ptRangeKind]);
        //         tempLMainEPCent->Add(hJESshif_clone_proj);
        //     }
        //     // ########################################################################
        // 
        //     // ###  Add JER1 hist  ####################################################
        //     std::cout << "9. Add Jet Energy Resolution hists 1 way   #####" << std::endl;
        //     histName = TString::Format("hJESForJER1_Cent%d",centBin);
        //     TProfile* hTempJESforJER = (TProfile* ) hJESshift_prof->Clone();
        //     hTempJESforJER->SetName(histName.Data());
        //     TProfile* hEachCentJER1 = getJER1(hTempJESforJER, TString::Format("JER1_cent%d",centBin));
        // 
        //     xTitle = "#it{p}_{T, gen}^{jet} [GeV/#it{c}]";
        //     yTitle = "#sigma(#frac{#it{p}_{T}^{det} - #it{p}_{T}^{particle}}{#it{p}_{T}^{particle}})";
        //     histLabelSetting(hEachCentJER1, histName, xTitle, yTitle, centBin, 1, reBin);
        //     tempLMainEPCent->Add(hEachCentJER1);
        //     // ########################################################################
        // 
        //     // ###  Add JER2 hist  ####################################################
        //     std::cout << "10. Add Jet Energy Resolution hists 2 way  #####" << std::endl;
        //     histName = TString::Format("hRM_forJER2_Cent%d",centBin);
        //     THnSparse* hTempEachCentRM_forJER = (THnSparse* )tempOriginMatchedJetHists->Clone();
        //     TH2D* hEachCentRM_forJER = (TH2D* ) hTempEachCentRM_forJER->Projection(0, 1, "");
        //     hEachCentRM_forJER->SetName(histName.Data());
        //     TProfile* hEachCentJER2 = getJER2(hEachCentRM_forJER,TString::Format("JER2_cent%d",centBin));
        // 
        //     xTitle = "#it{p}_{T, gen}^{jet} [GeV/#it{c}]";
        //     yTitle = "#frac{#sigma(#it{p}_{T}^{gen})}{#it{p}_{T}^{gen}}";
        //     histLabelSetting(hEachCentJER2, histName, xTitle, yTitle, centBin, 1, reBin);
        //     tempLMainEPCent->Add(hEachCentJER2);
        //     // ########################################################################
        // }
        
        // // ###  Add hNEFVsPt hist   ###############################################
        // std::cout << "10. Add NEFVsPt                       #####" << std::endl;
        // histName = TString::Format("hNEFVsPt_Cent%d",centBin);
        // hNEFVsPt = baseHNEFVsPt->Clone();
        // hNEFVsPt->SetName(histName.Data());
        // hNEFVsPt->GetXaxis()->SetRangeUser(centRangeList[centBin][0], centRangeList[centBin][1]);
        // hNEFVsPt_proj = hNEFVsPt->Project3D("zyeo");
        // hNEFVsPt_profX = hNEFVsPt_proj->ProfileX();
        // 
        // xTitle = "#it{p}_{T, gen}^{jet}";
        // yTitle = "NEF";
        // histLabelSetting(hNEFVsPt_profX, histName, xTitle, yTitle, centBin, 1, reBin);
        // tempLMainEPCent->Add(hNEFVsPt_profX);
        // // ##############################################################################
        // 
        // ###  Add hZLeadingVsPt hist   ################################################
        // histName = TString::Format("hZLeadingVsPt_Cent%d",centBin);
        // hZLeadingVsPt = baseHZLeadingVsPt->Clone();
        // hZLeadingVsPt->SetName(histName.Data());
        // hZLeadingVsPt->GetXaxis()->SetRangeUser(centRangeList[centBin][0], centRangeList[centBin][1]);
        // hZLeadingVsPt_proj = hZLeadingVsPt->Project3D("zyeo");
        // hZLeadingVsPt_profX = hZLeadingVsPt_proj->ProfileX();
        // 
        // xTitle = "#it{p}_{T, gen}^{jet} [GeV/#it{c}]";
        // yTitle = "zLeading #it{p}_{T}";
        // histLabelSetting(hZLeadingVsPt_profX, histName, xTitle, yTitle, centBin, 1, reBin);
        // tempLMainEPCent->Add(hZLeadingVsPt_profX);
        // ##############################################################################
    }
}

// def EditRMForUF(outList, hOriginRM, label):
//     ptBinArrayDict = PtRangeList->eachJetPtBinDef(1, ptRangeDict)
//     hRM_rebin = rebinRM(outList, hOriginRM, ptBinArrayDict)
//     hRM_norm = normalizeRM(outList,hRM_rebin, outList, ptBinArrayDict)
    // 
//     return hRM_norm

// ################################################################################
// # Rebin the response matrix to have variable binning                          ##
// ################################################################################
// def rebinRM(outList, hOriginRM, ptBinArrayDict):
//     genBinArray = ptBinArrayDict["mcGen"]
//     nGenBins = len(genBinArray )
//     detBinArray = ptBinArrayDict["mcDet"]
//     nDetBins = len(detBinArray )
//     histname = "{}_Rebinned"->format(hOriginRM->GetName())
//     title = histname + ";#it{p}_{T,corr}^{det} (GeV/#it{c});#it{p}_{T}^{truth} (GeV/#it{c})"
//     hRM_rebin = ROOT->TH2D(histname, title, nDetBins, detBinArray, nGenBins, genBinArray)

//     # Loop over all bins in fine-binned response matrix, 
//     # and fill appropriate bin in new response matrix
//     # Assume that the bin edges overlap appropriately
//     for detBin in range(1, hOriginRM->GetNbinsX() + 1):
//         for genBin in range(1, hOriginRM->GetNbinsY() + 1):
//             oldContent = hOriginRM->GetBinContent(detBin, genBin)

//             # Find the bin that should be filled in the new histogram, and fill it
//             # Need to get (x,y) location from bins (ibin, jbin)
//             x = hOriginRM->GetXaxis()->GetBinCenter(detBin)
//             y = hOriginRM->GetYaxis()->GetBinCenter(genBin)
//             hRM_rebin->Fill(x, y, oldContent)
//             
//             #print "Adding value {} from bin ({},{}) = ({},{})"->format(oldContent, ibin, jbin, x, y)
//             #newBin = hRM_rebin->FindBin(x,y)
//             #print "New bin content: {}"->format(hRM_rebin->GetBinContent(newBin))
//  
//     # Assume 0 errors on response matrix
//     for bin in range(1, hRM_rebin->GetNcells() + 1):
//         hRM_rebin->SetBinError(bin, 0)
// 
//     outList->Add(hRM_rebin)
// 
//     return hRM_rebin

// ################################################################################
// # Normalize response matrix                                                   ##
// # Normalize the pT-truth projection to 1                                      ##
// ################################################################################
// def normalizeRM(outList,hRM_rebin, outputDir, ptBinArrayDict):
//     genBinArray = ptBinArrayDict["mcGen"]
//     nGenBins = len(genBinArray )
//     detBinArray = ptBinArrayDict["mcDet"]
//     nDetBins = len(detBinArray )
//     # Make projection onto pT-true axis (y-axis), and scale appropriately
//     # Do exclude under and overflow bins
//     hGenProjBefore = hRM_rebin->ProjectionY("_py",1,hRM_rebin->GetNbinsX()) 
//     hGenProjBefore->SetName("hGenProjectionBefore");
// 
//     histname = "{}Normed"->format(hRM_rebin->GetName())
//     title = histname + ";#it{p}_{T,corr}^{det} (GeV/#it{c});#it{p}_{T}^{truth} (GeV/#it{c})"
//     hRM_norm = ROOT->TH2D(histname, title, nDetBins, detBinArray, nGenBins, genBinArray)
// 
//     # Loop through truth-level bins, and apply normalization factor to all bins->
//     nBinsY = hRM_rebin->GetNbinsY() # pT-gen
//     nBinsX = hRM_rebin->GetNbinsX() # pT-det
//     for genBin in range(1,nBinsY+1):
//         normFactor = hGenProjBefore->GetBinContent(genBin)
//         if normFactor > 0:
//             genBinCenter = hGenProjBefore->GetXaxis()->GetBinCenter(genBin)
// 
//             for detBin in range(1,nBinsX+1):
//                 binContent = hRM_rebin->GetBinContent(detBin, genBin)
//                 # hRM_rebin->SetBinContent(detBin, genBin, binContent/normFactor)
//                 hRM_norm->SetBinContent(detBin, genBin, binContent/normFactor)
// 
//     # Plot response matrix
//     c = ROOT->TCanvas("c","c: hist",600,450)
//     c->cd()
//     c->cd()->SetLeftMargin(0.15)
// 
//     # hRM_rebin->Draw("colz");
//     hRM_norm->Draw("colz");
//     line = ROOT->TLine(minPtDet,0,minPtDet,250)
//     line->SetLineColor(0)
//     line->SetLineStyle(2)
//     line->Draw("same");
//     line2 = ROOT->TLine(maxPtDet,0,maxPtDet,250)
//     line2->SetLineColor(0)
//     line2->SetLineStyle(2)
//     line2->Draw("same");
//     line3 = ROOT->TLine(0,minPtGen,100,minPtGen)
//     line3->SetLineColor(0)
//     line3->SetLineStyle(2)
//     line3->Draw("same");
//     line4 = ROOT->TLine(0,maxPtGen,100,maxPtGen)
//     line4->SetLineColor(0)
//     line4->SetLineStyle(2)
//     line4->Draw("same");
// 
//     # c->SaveAs("hoge->root");
//     c->Close()
// 
//     outList->Add(hRM_norm)
// 
//     return hRM_norm

Double_t movePhaseShift(Double_t jetAngle){
    if(jetAngle > TMath::Pi())   jetAngle -= TMath::Pi();
    if(jetAngle > TMath::Pi()/2) jetAngle = TMath::Pi() - jetAngle;
    
    return jetAngle;
}

// ################################################################################
// # Plot JER                                                                    ##
// ################################################################################
TProfile* getJER1(TProfile* hJESshift, TString name){
    TProfile* histJER = new TProfile(name.Data(), name.Data(), 100, 0., 200.);

    for(Int_t iBin=1; iBin < hJESshift->GetNbinsX() + 1;iBin++){
        Double_t JER = hJESshift->GetBinError(iBin);
        Double_t partJetPt = hJESshift->GetXaxis()->GetBinCenter(iBin);
        histJER->Fill(partJetPt, JER);
    }

    return histJER;
}

TProfile* getJER2(TH2D* histResponseMatrix, TString name){
    // # For each pT^gen, compute the standard deviation of the pT^det distribution
    // # Get the pT^gen profile, with errors as standard deviation of pT^det distribution
    TString histName = "histPtGenProf" +name;
    TProfile* histPtGenProf = (TProfile* ) histResponseMatrix->ProfileY(histName.Data(), 1, -1, "s");
    histPtGenProf->SetName(name.Data());
    
    // # Create histo to be used to fill JER values
    TProfile* histJER = new TProfile(name.Data(), name.Data(), 100, 0., 200.);
    
    // # Loop through the bins, and fill the JER
    for(Int_t bin=1;bin<histPtGenProf->GetNbinsX()+1;bin++){
        Double_t sigma = histPtGenProf->GetBinError(bin);
        Double_t pTgen = histPtGenProf->GetXaxis()->GetBinCenter(bin);
        Double_t JER = sigma/pTgen;
        histJER->Fill(pTgen, JER);
    }

    return histJER;
}

void addHistIntoList(TList *lTree, TH1D* inputHist, TString histName, Double_t scaleFactor, Int_t ptHardBin){
    inputHist->Sumw2();
    inputHist->Scale(scaleFactor);

    if(ptHardBin == 1) lTree->Add(inputHist);
    else{
        TH1D* tempHist = (TH1D*) lTree->FindObject(histName.Data());
        tempHist->Add(inputHist);
    }
}
void addHistIntoList(TList *lTree, TH2D* inputHist, TString histName, Double_t scaleFactor, Int_t ptHardBin){
    inputHist->Sumw2();
    inputHist->Scale(scaleFactor);

    if(ptHardBin == 1) lTree->Add(inputHist);
    else{
        TH2D* tempHist = (TH2D*) lTree->FindObject(histName.Data());
        tempHist->Add(inputHist);
    }
}
void addHistIntoList(TList *lTree, TH3D* inputHist, TString histName, Double_t scaleFactor, Int_t ptHardBin){
    inputHist->Sumw2();
    inputHist->Scale(scaleFactor);

    if(ptHardBin == 1) lTree->Add(inputHist);
    else{
        TH3D* tempHist = (TH3D*) lTree->FindObject(histName.Data());
        tempHist->Add(inputHist);
    }
}
void addHistIntoList(TList *lTree, TProfile* inputHist, TString histName, Double_t scaleFactor, Int_t ptHardBin){
    inputHist->Sumw2();
    inputHist->Scale(scaleFactor);

    if(ptHardBin == 1) lTree->Add(inputHist);
    else{
        TProfile * tempHist = (TProfile *) lTree->FindObject(histName.Data());
        tempHist->Add(inputHist);
    }
}
void addHistIntoList(TList *lTree, THnSparse* inputHist, TString histName, Double_t scaleFactor, Int_t ptHardBin){
    inputHist->Sumw2();
    inputHist->Scale(scaleFactor);

    if(ptHardBin == 1) lTree->Add(inputHist);
    else{
        THnSparse* tempHist = (THnSparse*) lTree->FindObject(histName.Data());
        tempHist->Add(inputHist);
    }
}

void histLabelSetting(TH1D *outputHist, TString histTitle, TString xTitle, TString yTitle, \
    Int_t histNumber, Int_t fillType, Double_t Rebin){
    outputHist->SetTitle(histTitle.Data());
    outputHist->SetXTitle(xTitle.Data());
    outputHist->SetYTitle(yTitle.Data());
    
    outputHist->GetXaxis()->SetTitleSize(20);
    outputHist->GetXaxis()->SetTitleFont(43);
    outputHist->GetXaxis()->SetTitleOffset(2.);

    outputHist->GetYaxis()->SetTitleSize(20);
    outputHist->GetYaxis()->SetTitleFont(43);
    outputHist->GetYaxis()->SetTitleOffset(2.);
    
    outputHist->SetLineColor(lFillColor[histNumber]);
    outputHist->SetMarkerColor(lFillColor[histNumber]);
    if(fillType) outputHist->SetMarkerStyle(lFillMerker[histNumber]);
    else outputHist->SetMarkerStyle(lNoFillMerker[histNumber]);
    outputHist->SetMarkerSize(1.8);
    outputHist->SetLineWidth(5);

    outputHist->Rebin(Rebin);
}
void histLabelSetting(TH2D *outputHist, TString histTitle, TString xTitle, TString yTitle, \
    Int_t histNumber, Int_t fillType, Double_t Rebin){
    outputHist->SetTitle(histTitle.Data());
    outputHist->SetXTitle(xTitle.Data());
    outputHist->SetYTitle(yTitle.Data());
    
    outputHist->GetXaxis()->SetTitleSize(20);
    outputHist->GetXaxis()->SetTitleFont(43);
    outputHist->GetXaxis()->SetTitleOffset(2.);

    outputHist->GetYaxis()->SetTitleSize(20);
    outputHist->GetYaxis()->SetTitleFont(43);
    outputHist->GetYaxis()->SetTitleOffset(2.);
    
    outputHist->SetLineColor(lFillColor[histNumber]);
    outputHist->SetMarkerColor(lFillColor[histNumber]);
    if(fillType) outputHist->SetMarkerStyle(lFillMerker[histNumber]);
    else outputHist->SetMarkerStyle(lNoFillMerker[histNumber]);
    outputHist->SetMarkerSize(1.8);
    outputHist->SetLineWidth(5);

    outputHist->Rebin(Rebin);
}
void histLabelSetting(TH3D *outputHist, TString histTitle, TString xTitle, TString yTitle, \
    Int_t histNumber, Int_t fillType, Double_t Rebin){
    outputHist->SetTitle(histTitle.Data());
    outputHist->SetXTitle(xTitle.Data());
    outputHist->SetYTitle(yTitle.Data());
    
    outputHist->GetXaxis()->SetTitleSize(20);
    outputHist->GetXaxis()->SetTitleFont(43);
    outputHist->GetXaxis()->SetTitleOffset(2.);

    outputHist->GetYaxis()->SetTitleSize(20);
    outputHist->GetYaxis()->SetTitleFont(43);
    outputHist->GetYaxis()->SetTitleOffset(2.);
    
    outputHist->SetLineColor(lFillColor[histNumber]);
    outputHist->SetMarkerColor(lFillColor[histNumber]);
    if(fillType) outputHist->SetMarkerStyle(lFillMerker[histNumber]);
    else outputHist->SetMarkerStyle(lNoFillMerker[histNumber]);
    outputHist->SetMarkerSize(1.8);
    outputHist->SetLineWidth(5);

    outputHist->Rebin(Rebin);
}
void histLabelSetting(TProfile *outputHist, TString histTitle, TString xTitle, TString yTitle, \
    Int_t histNumber, Int_t fillType, Double_t Rebin){
    outputHist->SetTitle(histTitle.Data());
    outputHist->SetXTitle(xTitle.Data());
    outputHist->SetYTitle(yTitle.Data());
    
    outputHist->GetXaxis()->SetTitleSize(20);
    outputHist->GetXaxis()->SetTitleFont(43);
    outputHist->GetXaxis()->SetTitleOffset(2.);

    outputHist->GetYaxis()->SetTitleSize(20);
    outputHist->GetYaxis()->SetTitleFont(43);
    outputHist->GetYaxis()->SetTitleOffset(2.);
    
    outputHist->SetLineColor(lFillColor[histNumber]);
    outputHist->SetMarkerColor(lFillColor[histNumber]);
    if(fillType) outputHist->SetMarkerStyle(lFillMerker[histNumber]);
    else outputHist->SetMarkerStyle(lNoFillMerker[histNumber]);
    outputHist->SetMarkerSize(1.8);
    outputHist->SetLineWidth(5);

    outputHist->Rebin(Rebin);
}


void SetScaleFactor(Double_t* lScaleFactor, Double_t* lPtHardEventNum, TString LHCPeriod, TString diffSys){

    Double_t lrefPtHardScalFactorLHC18q[20] = {9.177734972060997e-08, \
        2.8349028721358537e-08, 1.2906816933214331e-08, 4.438265121771655e-09, \
        1.5143096874822496e-09, 1.481292114033884e-08, 1.8896079634207416e-10, \
        6.094080087747223e-11, 6.719599044722054e-10, 8.133222071600738e-12, \
        2.9663998343658287e-12, 1.0475069333113453e-12, 5.007576907062431e-13, \
        2.325499458678617e-13, 1.1388060504181516e-13, 1.1133849925716706e-13, \
        3.1958597696234174e-14, 6.339668437202558e-15, 1.7674154782042322e-14, \
        4.734806373943197e-14};
    Double_t lrefPtHardNumOfEventLHC18q[20] = {130929720.0, 126156976.0, \
        127772880.0, 131483728.0, 130965840.0, 25415294.0, 123289592.0, 126142400.0, \
        24295196.0, 127763496.0, 132379808.0, 132791088.0, 134328048.0, 134483312.0, \
        134274000.0,96676568.0, 132849256.0, 215552448.0, 94747208.0, 65914416.0};

    Double_t lrefPtHardScalFactorLHC18r[20] = {1.0311858837350005e-07, \
        3.1965864804665285e-08, 1.3544813757767156e-08, 4.960101754275058e-09, \
        1.7000300954266733e-09, 6.850195617723509e-10, 1.9551046046071345e-10, \
        6.767757025526565e-11, 2.6830176915083785e-11, 9.551208018080331e-12, \
        3.415149269564962e-12, 1.22434827333091e-12, 6.227175402055881e-13, \
        3.053814558072003e-13, 1.4631121940665949e-13, 4.824270054036721e-13, \
        3.8743223134692397e-14, 1.9262639859449602e-14, 1.2205335164987764e-14, \
        1.076972637929458e-13};
    Double_t lrefPtHardNumOfEventLHC18r[20] = {119152816.0, 114575400.0, 120281920.0, \
        119913736.0, 119175216.0, 113954416.0, 116901768.0, 115446128.0, 117300200.0, \
        113688128.0, 119033856.0, 118434224.0, 116181552.0, 113190656.0, 114248560.0, \
        44793380.0, 116328024.0, 119243576.0, 109950448.0, 42151224.0};

    // == s == V0 root file ==
    Double_t lrefPtHardScalFactorLHC18qV0[20] = {1.0321989824871018e-07, 2.93822098900898e-08,\
        1.4154134786746665e-08, 5.193423001593628e-09, 1.7056757467304182e-09, \
        6.225948475912749e-10, 1.845875082499011e-10, 6.201911334052037e-11, \
        2.5421819691359822e-11, 8.421122396813474e-12, 3.2963503770197606e-12, \
        1.165303758790266e-12, 5.759455256796066e-13, 2.7035640908113497e-13, \
        1.3222542894576327e-13, 4.830849232450998e-14, 3.690801210396656e-14, \
        1.940868489423863e-14, 1.0371460954031887e-14, 1.7789622079360143e-14};
    Double_t lrefPtHardNumOfEventLHC18qV0[20] = {130869888.0, 131426952.0, 129399008.0, \
        128898584.0, 130849480.0, 131401312.0, 132276952.0, 132600872.0, 132539584.0, \
        133147424.0, 133166312.0, 133498968.0, 132816752.0, 132287688.0, 132129824.0, \
        155653296.0, 131128264.0, 130593248.0, 131126784.0, 114060416.0};

    Double_t lrefPtHardScalFactorLHC18rV0[20] = {1.0854657739954753e-07, \
        3.125689558557477e-08, 1.4795268318798147e-08, 5.361323264319431e-09, \
        1.820748781478633e-09, 6.573623469384724e-10, 1.9574099151899413e-10, \
        6.716684207816774e-11, 2.7645351065349576e-11, 9.16027296183672e-12, \
        3.5620883973149024e-12, 1.2764713296291094e-12, 6.216401942811304e-13, \
        4.763351619556121e-13, 1.417433574153358e-13, 1.199755862199768e-13, \
        3.906567932266237e-14, 1.9893576310816953e-14, 1.0690590074901366e-14, \
        1.3872420490884595e-14};
    Double_t lrefPtHardNumOfEventLHC18rV0[20] = {120807904.0, 120575528.0, 119736904.0,\
        120069312.0, 119875776.0, 121040976.0, 121558056.0, 120573872.0, 120247048.0, \
        120768424.0, 121215752.0, 120690928.0, 120959304.0, 94321376.0, 120782072.0, \
        93542600.0, 120635328.0, 122086672.0, 122201992.0, 122211872.0};
    // == e == V0 root file ==


    for(Int_t ptHardBin=0; ptHardBin<20; ptHardBin++){
        if(LHCPeriod=="LHC18q"){
            if((diffSys != "V0M") && (diffSys != "V0C") && (diffSys != "V0A")) {
                lScaleFactor[ptHardBin] = lrefPtHardScalFactorLHC18q[ptHardBin];
                lPtHardEventNum[ptHardBin] = lrefPtHardNumOfEventLHC18q[ptHardBin];
            }else {
                lScaleFactor[ptHardBin] = lrefPtHardScalFactorLHC18qV0[ptHardBin];
                lPtHardEventNum[ptHardBin] = lrefPtHardNumOfEventLHC18qV0[ptHardBin];
            }
        }else if(LHCPeriod=="LHC18r"){
            if((diffSys != "V0M") && (diffSys != "V0C") && (diffSys != "V0A")) {
                lScaleFactor[ptHardBin] = lrefPtHardScalFactorLHC18r[ptHardBin];
                lPtHardEventNum[ptHardBin] = lrefPtHardNumOfEventLHC18r[ptHardBin];
            }else {
                lScaleFactor[ptHardBin] = lrefPtHardScalFactorLHC18rV0[ptHardBin];
                lPtHardEventNum[ptHardBin] = lrefPtHardNumOfEventLHC18rV0[ptHardBin];
            }
        }

    }
}

// ################################################################################
void OFileStracture(TList *lMainTree, Int_t numOfCentBin){
    for(Int_t epBin = 0; epBin< 3; epBin++){
        TString listName = lEPLabel[epBin];
        TList *lEPList = new TList();
        lEPList->SetName(listName.Data());
        for (Int_t centBin = 0; centBin < numOfCentBin; centBin++){
            TList *lEachCentHists = new TList();
            TString listName = TString::Format("lCent%d",centBin);
            lEachCentHists->SetName(listName);
            lEPList->Add(lEachCentHists);
        }
        
        lMainTree->Add(lEPList);
    }
    
}
// ################################################################################
