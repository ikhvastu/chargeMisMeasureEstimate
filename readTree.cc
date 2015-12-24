#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>

#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TSystem.h"
#include "TVector.h"
#include "TText.h"
#include "TCOlor.h"

#include "TFrame.h"

#include "TLatex.h"

#include "fstream"
#include "sstream"  

#include "TApplication.h"   

//#include "mt2_bisect.h"

TString lumi_13TeV = "2.11 fb^{-1}";

using namespace std;

namespace{
  const int bands = 999;
  int rainbow[bands];
  int patriotic[bands];
}

int _multiWP[3] = {1,0,0};//Medium - ele, Loose - mu

void showHist(TVirtualPad*, TH1D *, THStack *, string, string, string, double, TLegend *);
double calculateZbi(double signal, double bkg, double unc);
double errCalc(int, int);

TH2D TranslateHisto(const TH2 &);
void Print2D(TH2 const * const h_data_in, TH2 const * const h_mc_in, const TString &ext);
void PrintTable(TH2 const * const histo, const TString &ext);
void PrintLine(ofstream &file, TH2 const * const histo, int bin, const TString &label);
void GetPatrioticPalette();
void GetRainbowPalette();

void readTree()
{
    
    const int nLeptonsMax = 10;

    int fontToUse = 42;
    gStyle->SetPaintTextFormat("1.5f");
    gStyle->SetOptFit(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetNdivisions(505,"XY");
     /*
    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);
    
    gStyle->SetLabelFont(fontToUse,"XYZ");
    gStyle->SetLabelSize(0.05,"XYZ");
    gStyle->SetLabelOffset(0.001,"X");
    
    gStyle->SetTitleFont(fontToUse,"");
    gStyle->SetTitleFont(fontToUse,"XYZ");
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetTitleSize(0.06,"XY");
    //gStyle->SetTitleXOffset(1.0);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYOffset(1.2);
    
    gStyle->SetErrorX(0.);
    
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadRightMargin(0.035);
    gStyle->SetPadLeftMargin(0.15);
    
    gStyle->SetStatFont(fontToUse);
    gStyle->SetStatColor(10);
    gStyle->SetStatFontSize(0.08);
    gStyle->SetTitleFont(fontToUse);
    gStyle->SetTitleFontSize(0.08);
    
    gStyle->SetMarkerSize(1.);
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerColor(gStyle->GetHistLineColor());
    
    gStyle->SetPalette(1,0);
    
    gStyle->SetFuncColor(kRed);
    */

    // List of branches
    TBranch        *b__eventNb;   //!
    TBranch        *b__runNb;   //!
    TBranch        *b__lumiBlock;   //!
    TBranch        *b__weight;   //!
    TBranch        *b__genqpt;   //!
    TBranch        *b__nLeptons;   //!
    TBranch        *b__lPt;   //!
    TBranch        *b__lEta;   //!
    TBranch        *b__lPhi;   //!
    TBranch        *b__lE;   //!

    TBranch        *b__lPtmc;
    TBranch        *b__lchargemc;

    TBranch        *b__nEle;   //!
    TBranch        *b__nMu;   //!
    TBranch        *b__nTau;   //!
    TBranch        *b__flavors;   //!
    TBranch        *b__charges;   //!
    TBranch        *b__indeces;   //!

    TBranch        *b__isolation;   //!
    TBranch        *b__miniisolation;   //!
    TBranch        *b__multiisolation;   //!

    TBranch        *b__ptrel;   //!
    TBranch        *b__ptratio;   //!

    TBranch        *b__origin;   //!
    TBranch        *b__originReduced;   //!

    TBranch        *b__PVchi2;   //!
    TBranch        *b__PVerr;   //!
    TBranch        *b__ipPV;   //!
    TBranch        *b__ipPVerr;   //!
    TBranch        *b__ipZPV;   //!
    TBranch        *b__ipZPVerr;   //!
    TBranch        *b__ipPVmc;   //!
    TBranch        *b__3dIP;   //!
    TBranch        *b__3dIPerr;   //!
    TBranch        *b__3dIPsig;   //!

    TBranch        *b__mt;   //!

    TBranch        *b__isloose;   //!
    TBranch        *b__istight;   //!
    TBranch        *b__chargeConst; 
    TBranch        *b__istightID;   //!

    TBranch        *b__met;   //!
    TBranch        *b__met_phi;   //!
    TBranch        *b_HT;   //!
    TBranch        *b__genmet; 
    TBranch        *b__genmet_phi;
    TBranch        *b__n_bJets;   //!
    TBranch        *b__n_Jets;   //!
    TBranch        *b__bTagged;   //!
    TBranch        *b__jetEta;   //!
    TBranch        *b__jetPhi;   //!
    TBranch        *b__jetPt;   //!
    TBranch        *b__jetE; 
    TBranch        *b__jetFlavour;
    TBranch        *b__csv;   //!
    TBranch        *b__jetDeltaR;   //!

    
    // Declaration of leaf types
    ULong64_t       _eventNb;
    ULong64_t       _runNb;
    ULong64_t       _lumiBlock;
    Double_t        _weight;
    Double_t        _genqpt;
    Int_t           _nLeptons;
    Double_t        _lPt[10];
    Double_t        _lEta[10];
    Double_t        _lPhi[10];
    Double_t        _lE[10];

    Double_t        _lPtmc[10];
    Double_t        _lchargemc[10];

    Int_t           _nEle;
    Int_t           _nMu;
    Int_t           _nTau;
    Int_t           _flavors[10];
    Double_t        _charges[10];
    Int_t           _indeces[10];

    Double_t        _isolation[10];
    Double_t        _miniisolation[10][2];
    Bool_t          _multiisolation[10][3];
    
    Double_t        _ptrel[10];
    Double_t        _ptratio[10];

    Int_t           _origin[10];
    Int_t           _originReduced[10];

    Double_t        _PVchi2;
    Double_t        _PVerr[3];
    Double_t        _ipPV[10];
    Double_t        _ipPVerr[10];
    Double_t        _ipZPV[10];
    Double_t        _ipZPVerr[10];
    Double_t        _ipPVmc[10];
    Double_t        _3dIP[10];
    Double_t        _3dIPerr[10];
    Double_t        _3dIPsig[10];

    Double_t        _mt[10];

    Bool_t          _isloose[10];
    Bool_t          _istight[10];
    Bool_t          _istightID[10];
    Bool_t          _chargeConst[10];

    Double_t        _met;
    Double_t        _met_phi;
    Double_t        HT;
    Double_t        _genmet;
    Double_t        _genmet_phi;

    Int_t           _n_bJets;
    Int_t           _n_Jets;
    Bool_t          _bTagged[20];
    Double_t        _jetEta[20];
    Double_t        _jetPhi[20];
    Double_t        _jetPt[20];
    Double_t        _jetE[20];
    Int_t           _jetFlavour[20];
    Double_t        _csv[20];
    Double_t        _jetDeltaR[20][10];

    Int_t           _n_PV;
    Double_t        _n_MCTruth_PV;

    double _mvaValue[nLeptonsMax];

    int _hitsNumber[nLeptonsMax];

    bool _vtxFitConversion[nLeptonsMax];

    
    TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
    
    const int nSamples = 17;
    
    
    //TString fileList[nSamples] = {"ttbar_2l.root"};
    TString fileList[nSamples] = {"dy50_OS.root",
        "ttbar_2l_OS.root",
        //"dy50_13.root", 
        "dy50_OS.root", 
        "wjets_13.root", "wjets_1_13.root", "wjets_2_13.root", "wjets_3_13.root", "wjets_4_13.root", "wjets_5_13.root", "wjets_6_13.root", "wjets_7_13.root",
        "WZ13.root", "ttZ13.root", "ttH13.root", "tZq13.root",
        "ttW13.root",
        //"data_ee.root"};
        //"data_mumu.root"};
        "data_ee_OS.root"};
    

    double xSections[nSamples] = {
        6025.2, 
        87.35, 6025.2,
        61526, 1347, 360, 48.9, 12.8, 5.26, 1.33, 0.03089,
        4.42965, 0.2529, 0.2, 0.0758,
        0.2043,
        1
    };

    Color_t colsStack[nSamples] = {kBlue, kRed, kGreen, kYellow, kYellow, kYellow, kYellow, kYellow, kYellow, kYellow, kYellow, kBlue, kBlue, kBlue, kBlue, kBlue, kBlack};
    //Color_t colsStack[nSamples] = {kBlue, kRed, kGreen, kBlack};
    int stylesStack[nSamples] = {20, 21, 22, 23};
    

    TFile *hfile[nSamples];
    TTree* inputTreeFake[nSamples];
    
    //TString names[nSamples] = {"t#bar{t}","dy_M-5_50","dy_M-10_50_1","dy_M-10_50_2","dy_M-10_50_3","dy_M-10_50_4","dy_M-50","dy_M-50_1","dy_M-50_2","dy_M-50_3","dy_M-50_4","WZ","t#bar{t}Z","t#bar{t}W"};
    TString names[nSamples] = {"DY", "t#bar{t}", "DY Z-window","WJets", "WJets_1_13", "WJets_2_13", "WJets_3_13", "WJets_4_13", "WJets_5_13", "WJets_6_13", "WJets_7_13","WZ","t#bar{t}Z", "t#bar{t}H", "tZq","t#bar{t}W","data"};

    const int numberBinspt = 7;
    double ptb23[numberBinspt] = {15., 40., 60., 80., 100., 200., 300.}; 

    const int numberBinseta = 4;
    double etab23[numberBinspt] = {0., 0.8, 1.479, 2.5}; 

    TH1D* h_leadpt[nSamples][2];
    TH1D* h_2ndleadpt[nSamples];
    TH1D* h_trailpt[nSamples];
    TH1D* h_sumpt[nSamples];
    TH1D* h_mll[nSamples];
    
    TH1D* h_njets[nSamples];
    
    TH1D* h_SR[nSamples];
    
    TH1D* h_dilep[nSamples];

    TH1D* h_ele_pt[nSamples];
    TH1D* h_ele_ptc[nSamples];

    TH1D* h_mu_pt[nSamples];
    TH1D* h_mu_ptc[nSamples];

    TH2D* h_ele_pt_eta[nSamples][2];

    TH1D* h_deltaPt[nSamples][2];

    TH1D* h_Mee[nSamples];
    TH1D* h_eta[nSamples];

    TGraphAsymmErrors *gh_ele_ptc[nSamples];
    
    for (int sample=0; sample!=nSamples; ++sample)  {
        TString name;
        //for(int j = 0; j < 4; j++){
            for(int i = 0; i < 2; i++){
                name = Form("lead_leptons_1D_%d_%d",sample, i);
                h_leadpt[sample][i] = new TH1D(name, "Leading lepton p_{T} "+names[sample]+";Leading lepton p_{T} [GeV]; events / 10 GeV", numberBinspt-1, ptb23);
                h_leadpt[sample][i]->SetLineColor(colsStack[sample]);
                h_leadpt[sample][i]->SetFillColor(colsStack[sample]);
                h_leadpt[sample][i]->SetFillStyle(1001);
                h_leadpt[sample][i]->SetMarkerColor(colsStack[sample]);
                h_leadpt[sample][i]->SetMarkerStyle(stylesStack[sample]);
                h_leadpt[sample][i]->Sumw2();
            }
       // }

        for(int i = 0; i < 2; i++){
            name = Form("lead_leptons_%d_%d",sample, i);
            h_ele_pt_eta[sample][i] = new TH2D(name, "Q mis-ID "+names[sample]+";Leading lepton p_{T} [GeV]; #eta", numberBinspt - 1, ptb23, numberBinseta -1, etab23);
            h_ele_pt_eta[sample][i]->SetLineColor(colsStack[sample]);
            h_ele_pt_eta[sample][i]->SetFillColor(colsStack[sample]);
            h_ele_pt_eta[sample][i]->SetFillStyle(1001);
            h_ele_pt_eta[sample][i]->SetMarkerColor(colsStack[sample]);
            h_ele_pt_eta[sample][i]->Sumw2();
        

            name = Form("deltaPt_%d_%d",sample,i);
            h_deltaPt[sample][i] = new TH1D(name, "Delta p_{T} "+names[sample]+";Delta lepton p_{T} [GeV]; events", 100, 0, 2);
            h_deltaPt[sample][i]->SetLineColor(colsStack[sample]);
            h_deltaPt[sample][i]->SetFillColor(colsStack[sample]);
            h_deltaPt[sample][i]->SetFillStyle(1001);
            h_deltaPt[sample][i]->SetMarkerColor(colsStack[sample]);
            h_deltaPt[sample][i]->Sumw2();
        }
        
        name = Form("M_ee_%d",sample);
        h_Mee[sample] = new TH1D(name, "M_{ee} "+names[sample]+";M_{ee} [GeV]; events / 5 GeV", 50, 70, 120);
        h_Mee[sample]->SetLineColor(colsStack[sample]);
        h_Mee[sample]->SetFillColor(colsStack[sample]);
        h_Mee[sample]->SetFillStyle(1001);
        h_Mee[sample]->SetMarkerColor(colsStack[sample]);
        h_Mee[sample]->Sumw2();

        name = Form("Eta_%d",sample);
        h_eta[sample] = new TH1D(name, "#eta "+names[sample]+";#eta; events", 10, -2.5, 2.5);
        h_eta[sample]->SetLineColor(colsStack[sample]);
        h_eta[sample]->SetFillColor(colsStack[sample]);
        h_eta[sample]->SetFillStyle(1001);
        h_eta[sample]->SetMarkerColor(colsStack[sample]);
        h_eta[sample]->Sumw2();
        

        name = Form("2ndlead_leptons_%d",sample);
        h_2ndleadpt[sample] = new TH1D(name, "2nd Leading lepton p_{T} "+names[sample]+";2nd Leading lepton p_{T} [GeV]; events / 10 GeV", 10, 0, 100);
        h_2ndleadpt[sample]->SetLineColor(colsStack[sample]);
        h_2ndleadpt[sample]->SetFillColor(colsStack[sample]);
        h_2ndleadpt[sample]->SetFillStyle(1001);
        h_2ndleadpt[sample]->SetMarkerColor(colsStack[sample]);
        //h_leadpt[sample]->SetMarkerStyle(20+sample*5);
        h_2ndleadpt[sample]->Sumw2();
        
        name = Form("trail_leptons_%d",sample);
        h_trailpt[sample] = new TH1D(name, "Trailing lepton p_{T} "+names[sample]+";Trailing lepton p_{T} [GeV]; events / 10 GeV", 10, 0, 100);
        h_trailpt[sample]->SetLineColor(colsStack[sample]);
        h_trailpt[sample]->SetFillColor(colsStack[sample]);
        h_trailpt[sample]->SetFillStyle(1001);
        h_trailpt[sample]->SetMarkerColor(colsStack[sample]);
        //h_trailpt[sample]->SetMarkerStyle(20+sample*5);
        h_trailpt[sample]->Sumw2();
        
        name = Form("sumpt_leptons_%d",sample);
        h_sumpt[sample] = new TH1D(name, "Sum lepton p_{T} "+names[sample]+";Sum lepton p_{T} [GeV]; events / 10 GeV", 40, 0, 400);
        h_sumpt[sample]->SetLineColor(colsStack[sample]);
        h_sumpt[sample]->SetFillColor(colsStack[sample]);
        h_sumpt[sample]->SetFillStyle(1001);
        h_sumpt[sample]->SetMarkerColor(colsStack[sample]);
        //h_trailpt[sample]->SetMarkerStyle(20+sample*5);
        h_sumpt[sample]->Sumw2();
        
        name = Form("mll_%d",sample);
        h_mll[sample] = new TH1D(name, "Invariant mass of 2 lepton "+names[sample]+";Invariant mass of 2 lepton [GeV]; events / 10 GeV", 12, 0, 120);
        h_mll[sample]->SetLineColor(colsStack[sample]);
        h_mll[sample]->SetFillColor(colsStack[sample]);
        h_mll[sample]->SetFillStyle(1001);
        h_mll[sample]->SetMarkerColor(colsStack[sample]);
        h_mll[sample]->Sumw2();
        

        name = Form("h_njets_%d",sample);
        h_njets[sample] = new TH1D(name, "N_{jets} "+names[sample]+";N_{jets}; events / 1", 10, 0, 10);
        h_njets[sample]->SetLineColor(colsStack[sample]);
        h_njets[sample]->SetFillColor(colsStack[sample]);
        h_njets[sample]->SetFillStyle(1001);
        h_njets[sample]->SetMarkerColor(colsStack[sample]);
        
        h_njets[sample]->Sumw2();


        name = Form("h_SR_%d",sample);
        h_SR[sample] = new TH1D(name, "SR "+names[sample]+";SR; events / 1", 400, 0, 400);
        h_SR[sample]->SetLineColor(colsStack[sample]);
        h_SR[sample]->SetFillColor(colsStack[sample]);
        h_SR[sample]->SetFillStyle(1001);
        h_SR[sample]->SetMarkerColor(colsStack[sample]);

        
        name = Form("h_each_dilepton_%d",sample);
        h_dilep[sample] = new TH1D(name, "dilep "+names[sample]+";dilep; events / 1", 7, -0.5, 6.5);
        h_dilep[sample]->SetLineColor(colsStack[sample]);
        h_dilep[sample]->SetFillColor(colsStack[sample]);
        h_dilep[sample]->SetFillStyle(1001);
        h_dilep[sample]->SetMarkerColor(colsStack[sample]);
        h_dilep[sample]->Sumw2();
        

        name = Form("electron_pt_%d", sample);
        h_ele_pt[sample] = new TH1D(name, "Electron p_{T} "+names[sample]+";Electron p_{T} (GeV); events", numberBinspt - 1, ptb23);
        h_ele_pt[sample]->SetLineColor(colsStack[sample]);
        h_ele_pt[sample]->SetFillColor(colsStack[sample]);
        h_ele_pt[sample]->SetMarkerColor(colsStack[sample]);
        h_ele_pt[sample]->Sumw2();

        name = Form("electron_pt_cut_%d", sample);
        h_ele_ptc[sample] = new TH1D(name, "Electron p_{T} cut "+names[sample]+";Electron p_{T} (GeV); events", numberBinspt - 1, ptb23);
        h_ele_ptc[sample]->SetLineColor(colsStack[sample]);
        h_ele_ptc[sample]->SetFillColor(colsStack[sample]);
        h_ele_ptc[sample]->SetMarkerColor(colsStack[sample]);
        h_ele_ptc[sample]->Sumw2();

        name = Form("muon_pt_%d", sample);
        h_mu_pt[sample] = new TH1D(name, "Muon p_{T} "+names[sample]+";Muon p_{T} (GeV); events", numberBinspt - 1, ptb23);
        h_mu_pt[sample]->SetLineColor(colsStack[sample]);
        h_mu_pt[sample]->SetFillColor(colsStack[sample]);
        h_mu_pt[sample]->SetMarkerColor(colsStack[sample]);
        h_mu_pt[sample]->Sumw2();

        name = Form("muon_pt_cut_%d", sample);
        h_mu_ptc[sample] = new TH1D(name, "Muon p_{T} cut"+names[sample]+";Muon p_{T} (GeV); events", numberBinspt - 1, ptb23);
        h_mu_ptc[sample]->SetLineColor(colsStack[sample]);
        h_mu_ptc[sample]->SetFillColor(colsStack[sample]);
        h_mu_ptc[sample]->SetMarkerColor(colsStack[sample]);
        h_mu_ptc[sample]->Sumw2();

    }
    
    THStack* st_leadpt = new THStack("st_leadpt","Leading lepton p_{T}");
    THStack* st_2ndleadpt = new THStack("st_2ndleadpt","Leading lepton p_{T}; Leading lepton p_{T}, [GeV]; events");
    THStack* st_trailpt = new THStack("st_trailpt","Trailing lepton p_{T}; Trailing lepton p_{T}, [GeV]; events");
    THStack* st_sumpt = new THStack("st_sumpt","Sum lepton p_{T}");
    THStack* st_mll = new THStack("st_mll","Invariant mass of 2 leptons");
    THStack* st_njets = new THStack("st_njets","N_{jets}");
    THStack* st_SR = new THStack("st_SR","SR");
    THStack* st_dilep = new THStack("st_dilep","dilep");
    THStack* st_Mee = new THStack("st_mee","Invariant mass of 2 OS electrons; M_{e^{-}e^{+}}, [GeV];events");
    THStack* st_eta = new THStack("st_eta","#eta;#eta;events");
    
    for (int i=0; i!=nSamples-2; ++i) {
        //st_leadpt->Add(h_leadpt[i]);
        //
        
        //st_sumpt->Add(h_sumpt[i]);
        //st_mll->Add(h_mll[i]);
        
        st_SR->Add(h_SR[i]);
        st_dilep->Add(h_dilep[i]);
        if(i == 2) continue;
        st_njets->Add(h_njets[i]);
        st_Mee->Add(h_Mee[i]);
        st_eta->Add(h_eta[i]);
        st_2ndleadpt->Add(h_2ndleadpt[i]);
        st_trailpt->Add(h_trailpt[i]);
    }
    
    const int nVars  = 18;
    TH1D* distribs[nVars][nSamples];
    THStack* distribsST[nVars];
    TString varN[nVars] = {"p_{T}^{leading} [GeV]","p_{T}^{trailing} [GeV]",
        "|#eta|^{max}","M_{ll}(OSSF)","M_{T} [GeV]","M_{ee} [GeV]","R_{mini}",
        "E_{T}^{miss}","H_{T}","N_{jets}","N_{jets40}","N_{b jets}","M_{T2}","M_{T2}^{true}","#DeltaM(OSSF,Z)", "M_{T2ll}", "M_{T2blbl}",
        "NPV"};
    double varMin[nVars] = {0,0,
        0,0,0,0,0,
        0,0,0,0,0,0,0,0, 0, 0,
        0.};
    
    double varMax[nVars] = {200,100,
        2.4,200,200,200,10,
        500,800,10,10,5,200,200,100, 200, 400, 
        50.
    };
    
    int nBins[nVars] = {10,20,
        12,50,40,40,50,
        10,16,10,10,5,40,40,20, 20, 40,
        50
    };
    
    for (int i=0; i!=nVars;++i) {
        TString name = Form("varST_%d",i);
        distribsST[i] = new THStack(name,varN[i]);
        for (int j=nSamples-1; j!=-1; --j) {
            name = Form("var_%d_%d",i,j);
            distribs[i][j] = new TH1D(name,name+";"+varN[i],nBins[i],varMin[i],varMax[i]);
            distribs[i][j]->SetLineColor(colsStack[j]);
            if (j < nSamples - 1)
                distribs[i][j]->SetFillColor(colsStack[j]);
            distribs[i][j]->SetMarkerColor(colsStack[j]);
            distribs[i][j]->SetMarkerStyle(1);
            distribs[i][j]->Sumw2();
            if (j < nSamples - 1)
                distribsST[i]->Add(distribs[i][j]);
        }
        
    }

    //0 - vloose, 1 - vloose FOIDEmu, 2 - vlooseFOIDISOEMU, 3 - tight
    double valuesMVA[3][4];
    valuesMVA[0][0] = -0.16;
    valuesMVA[1][0] = -0.65;
    valuesMVA[2][0] = -0.74;

    valuesMVA[0][1] = -0.70;
    valuesMVA[1][1] = -0.83; 
    valuesMVA[2][1] = -0.92; 

    valuesMVA[0][2] = -0.155;
    valuesMVA[1][2] = -0.56;
    valuesMVA[2][2] = -0.76; 

    valuesMVA[0][3] = 0.87;
    valuesMVA[1][3] = 0.60;
    valuesMVA[2][3] = 0.17;

    double _stiso_el[2] = {0.0766, 0.0678};
    //double _stiso_el[2] = {0.4, 0.4};

    int numberOS[nSamples][4] = {0};
    int numberSS[nSamples][4] = {0};

    int numberOSTotal = 0;
    int numberSSTotal = 0;

    ofstream outFile;
    outFile.open("chargemisID.txt");

    double stiso[2] = {0.0354, 0.0646};

    TFile *file_dataMC = TFile::Open("/Users/ikhvastu/Desktop/CERN/pileupCorr/dataMC.root","READ");
    TH1D *h_dataMC = (TH1D*)file_dataMC->Get("h3");    
    
    for (int sam = 0; sam != nSamples; ++sam) {
        //if(sam == 2 || sam == 3) continue;
        //if(sam == 1) continue;
        if(sam != 0 && sam != 1 && sam != 16) continue;

        //if(sam == 16)
        //    hfile[sam] = new TFile("data/"+fileList[sam],"read");
        //else
        //hfile[sam] = new TFile("data/"+fileList[sam],"read");
        hfile[sam] = new TFile("/Users/ikhvastu/Desktop/CERN/MCsamples/"+fileList[sam],"read");
        
        hfile[sam]->cd("FakeElectrons");
        inputTreeFake[sam] = static_cast<TTree*>(hfile[sam]->Get("FakeElectrons/fakeTree"));
        
        inputTreeFake[sam]->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
        inputTreeFake[sam]->SetBranchAddress("_runNb", &_runNb, &b__runNb);
        inputTreeFake[sam]->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
        inputTreeFake[sam]->SetBranchAddress("_weight", &_weight, &b__weight);
        inputTreeFake[sam]->SetBranchAddress("_genqpt", &_genqpt, &b__genqpt);
        inputTreeFake[sam]->SetBranchAddress("_nLeptons", &_nLeptons, &b__nLeptons);
        inputTreeFake[sam]->SetBranchAddress("_lPt", _lPt, &b__lPt);
        inputTreeFake[sam]->SetBranchAddress("_lEta", _lEta, &b__lEta);
        inputTreeFake[sam]->SetBranchAddress("_lPhi", _lPhi, &b__lPhi);
        inputTreeFake[sam]->SetBranchAddress("_lE", _lE, &b__lE);

        inputTreeFake[sam]->SetBranchAddress("_lPtmc", _lPtmc, &b__lPtmc);
        inputTreeFake[sam]->SetBranchAddress("_lchargemc", _lchargemc, &b__lchargemc);

        inputTreeFake[sam]->SetBranchAddress("_nEle", &_nEle, &b__nEle);
        inputTreeFake[sam]->SetBranchAddress("_nMu", &_nMu, &b__nMu);
        inputTreeFake[sam]->SetBranchAddress("_nTau", &_nTau, &b__nTau);
        inputTreeFake[sam]->SetBranchAddress("_flavors", _flavors, &b__flavors);
        inputTreeFake[sam]->SetBranchAddress("_charges", _charges, &b__charges);
        inputTreeFake[sam]->SetBranchAddress("_indeces", _indeces, &b__indeces);

        inputTreeFake[sam]->SetBranchAddress("_isolation", _isolation, &b__isolation);
        inputTreeFake[sam]->SetBranchAddress("_miniisolation", _miniisolation, &b__miniisolation);
        inputTreeFake[sam]->SetBranchAddress("_multiisolation", _multiisolation, &b__multiisolation);

        inputTreeFake[sam]->SetBranchAddress("_ptrel", _ptrel, &b__ptrel);
        inputTreeFake[sam]->SetBranchAddress("_ptratio", _ptratio, &b__ptratio);

        inputTreeFake[sam]->SetBranchAddress("_origin", _origin, &b__origin);
        inputTreeFake[sam]->SetBranchAddress("_originReduced", _originReduced, &b__originReduced);

        inputTreeFake[sam]->SetBranchAddress("_PVchi2", &_PVchi2, &b__PVchi2);
        inputTreeFake[sam]->SetBranchAddress("_PVerr", _PVerr, &b__PVerr);
        inputTreeFake[sam]->SetBranchAddress("_ipPV", _ipPV, &b__ipPV);
        inputTreeFake[sam]->SetBranchAddress("_ipPVerr", _ipPVerr, &b__ipPVerr);
        inputTreeFake[sam]->SetBranchAddress("_ipZPV", _ipZPV, &b__ipZPV);
        inputTreeFake[sam]->SetBranchAddress("_ipZPVerr", _ipZPVerr, &b__ipZPVerr);
        inputTreeFake[sam]->SetBranchAddress("_ipPVmc", _ipPVmc, &b__ipPVmc);
        inputTreeFake[sam]->SetBranchAddress("_3dIP", _3dIP, &b__3dIP);
        inputTreeFake[sam]->SetBranchAddress("_3dIPerr", _3dIPerr, &b__3dIPerr);
        inputTreeFake[sam]->SetBranchAddress("_3dIPsig", _3dIPsig, &b__3dIPsig);

        inputTreeFake[sam]->SetBranchAddress("_mt", _mt, &b__mt);

        inputTreeFake[sam]->SetBranchAddress("_isloose", _isloose, &b__isloose);
        inputTreeFake[sam]->SetBranchAddress("_istight", _istight, &b__istight);
        inputTreeFake[sam]->SetBranchAddress("_chargeConst", _chargeConst, &b__chargeConst);
        
        inputTreeFake[sam]->SetBranchAddress("_istightID", _istightID, &b__istightID);

        inputTreeFake[sam]->SetBranchAddress("_met", &_met, &b__met);
        inputTreeFake[sam]->SetBranchAddress("_met_phi", &_met_phi, &b__met_phi);

        inputTreeFake[sam]->SetBranchAddress("HT", &HT, &b_HT);
        inputTreeFake[sam]->SetBranchAddress("_genmet", &_genmet, &b__genmet);
        inputTreeFake[sam]->SetBranchAddress("_genmet_phi", &_genmet_phi, &b__genmet_phi);

        inputTreeFake[sam]->SetBranchAddress("_n_bJets", &_n_bJets, &b__n_bJets);
        inputTreeFake[sam]->SetBranchAddress("_n_Jets", &_n_Jets, &b__n_Jets);
        inputTreeFake[sam]->SetBranchAddress("_bTagged", _bTagged, &b__bTagged);
        inputTreeFake[sam]->SetBranchAddress("_jetEta", _jetEta, &b__jetEta);
        inputTreeFake[sam]->SetBranchAddress("_jetPhi", _jetPhi, &b__jetPhi);
        inputTreeFake[sam]->SetBranchAddress("_jetPt", _jetPt, &b__jetPt);
        inputTreeFake[sam]->SetBranchAddress("_jetE", _jetE, &b__jetE);
        //inputTreeFake[sam]->SetBranchAddress("_jetFlavour", _jetFlavour, &b__jetFlavour);
        inputTreeFake[sam]->SetBranchAddress("_csv", _csv, &b__csv);
        inputTreeFake[sam]->SetBranchAddress("_jetDeltaR", _jetDeltaR, &b__jetDeltaR);

        inputTreeFake[sam]->SetBranchAddress("_mvaValue", &_mvaValue);

        inputTreeFake[sam]->SetBranchAddress("_hitsNumber", &_hitsNumber);

        inputTreeFake[sam]->SetBranchAddress("_vtxFitConversion", &_vtxFitConversion);

        inputTreeFake[sam]->SetBranchAddress("_n_PV", &_n_PV);

        if(sam < nSamples-1){
            inputTreeFake[sam]->SetBranchAddress("_n_MCTruth_PV", &_n_MCTruth_PV);
        }
        _hCounter->Read("hCounter");
        Double_t scale = xSections[sam]*2110/(_hCounter->GetBinContent(1));
        
        Long64_t nEntries = inputTreeFake[sam]->GetEntries();
        std::cout<<"Entries in "<<fileList[sam]<<" "<<nEntries<<std::endl;
        std::cout<<xSections[sam]<<" "<<_hCounter->GetBinContent(1)<<" "<<scale<<std::endl;

        int counter = 0;

        int nLoc = 0;

        int posCharge = 0;
        int negCharge = 0;

        int leptInd[6];
    
        TLorentzVector l0p4, l1p4;

        int numberEvents = 0;
        int numberEventsPassed = 0;

        //numberOS = 0;
        //numberSS = 0;

        numberOSTotal = 0;
        numberSSTotal = 0;

        double weight = 1.;


        for (Long64_t it=0; it!=nEntries; ++it) {

            inputTreeFake[sam]->GetEntry(it);

            if (it%100000 == 0)
                cout<<'.'<<flush;

            //if(it == 10000) break;

            if(sam == 16){
                scale = 1.;
                _weight = 1.;
            }

            //if(sam == 1 && _genqpt > 100) continue;
           
            if (_nLeptons < 2) continue;
            
            nLoc = 0;

            for (int i=0; i!=_nLeptons; ++i) {
                //if (_flavors[i] > 0) continue;
                /*
                bool passedMVA = false;
                if (TMath::Abs(_lEta[i]) < 0.8 ) {
                    passedMVA = _mvaValue[i]> valuesMVA[0][3];
                } else if (TMath::Abs(_lEta[i]) < 1.479 ) {
                    passedMVA = _mvaValue[i]> valuesMVA[1][3];
                } else {
                    passedMVA = _mvaValue[i]> valuesMVA[2][3];
                }
                */
                    
                if(_miniisolation[i][0] > 0.4) continue;
                if(!_chargeConst[i]) continue;
                if(_lPt[i] < 10) continue;
                if(_3dIPsig[i] > 4) continue;

                if(_flavors[i] == 0 && _istight[i] && ((_miniisolation[i][0] < 0.12) && ((_ptratio[i] > 0.80) || (_ptrel[i] > 7.2)))){
                //if(_flavors[i] == 0 && _istight[i] && _isolation[i] < stiso[TMath::Abs(_lEta[i]) > 1.479]){
                //if(_flavors[i] == 1 && ((_miniisolation[i][0] < 0.16) && ((_ptratio[i] > 0.76) || (_ptrel[i] > 7.2)))){    
                    leptInd[nLoc] = i;
                    nLoc++;
                }
                     
            }
            
            //std::cout << "Electrons number: " << nLoc << std::endl;
            if (nLoc != 2) continue;

            double maxPt = 0.;
            double maxEta = 0.;
            double minPt = 0.;
            double munEta = 0.;
            int etaInd = -1;

            if(_lPt[leptInd[0]] > _lPt[leptInd[1]]){
                maxPt = _lPt[leptInd[0]];
                maxEta = _lEta[leptInd[0]];
                minPt = _lPt[leptInd[1]];
                //if(_charges[leptInd[0]] != _charges[leptInd[1]]){
                    if(TMath::Abs(_lEta[leptInd[0]]) < 1.479 && TMath::Abs(_lEta[leptInd[1]]) < 1.479)
                        etaInd = 0;
                    if(TMath::Abs(_lEta[leptInd[0]]) > 1.479 && TMath::Abs(_lEta[leptInd[1]]) < 1.479)
                        etaInd = 1;
                    if(TMath::Abs(_lEta[leptInd[0]]) < 1.479 && TMath::Abs(_lEta[leptInd[1]]) > 1.479)
                        etaInd = 2;
                    if(TMath::Abs(_lEta[leptInd[0]]) > 1.479 && TMath::Abs(_lEta[leptInd[1]]) > 1.479)
                        etaInd = 3;
                //}
                
            }
            else{
                maxPt = _lPt[leptInd[1]];
                maxEta = _lEta[leptInd[1]];
                minPt = _lPt[leptInd[0]];
                //if(_charges[leptInd[0]] != _charges[leptInd[1]]){
                    if(TMath::Abs(_lEta[leptInd[1]]) < 1.479 && TMath::Abs(_lEta[leptInd[0]]) < 1.479)
                        etaInd = 0;
                    if(TMath::Abs(_lEta[leptInd[1]]) > 1.479 && TMath::Abs(_lEta[leptInd[0]]) < 1.479)
                        etaInd = 1;
                    if(TMath::Abs(_lEta[leptInd[1]]) < 1.479 && TMath::Abs(_lEta[leptInd[0]]) > 1.479)
                        etaInd = 2;
                    if(TMath::Abs(_lEta[leptInd[1]]) > 1.479 && TMath::Abs(_lEta[leptInd[0]]) > 1.479)
                        etaInd = 3;
                //}
            }

            double deltaMZ = 999999.;
            l0p4.SetPtEtaPhiE(_lPt[leptInd[0]],_lEta[leptInd[0]],_lPhi[leptInd[0]],_lE[leptInd[0]]);
            l1p4.SetPtEtaPhiE(_lPt[leptInd[1]],_lEta[leptInd[1]],_lPhi[leptInd[1]],_lE[leptInd[1]]);
            l1p4+=l0p4;
            double mdiL = l1p4.M();
            deltaMZ = fabs(mdiL - 91);
            if(deltaMZ > 15)
                continue;

            if(_charges[leptInd[0]] == - _charges[leptInd[1]]){
                    numberOSTotal++;
                    h_Mee[sam]->Fill(mdiL, scale*_weight);
                    h_eta[sam]->Fill(_lEta[leptInd[0]], scale*_weight);
                    h_eta[sam]->Fill(_lEta[leptInd[1]], scale*_weight);
                    h_2ndleadpt[sam]->Fill(maxPt, scale*_weight);
                    h_trailpt[sam]->Fill(minPt, scale*_weight);

                    numberOS[sam][etaInd] += 1;    
                    h_leadpt[sam][1]->Fill(maxPt, scale*_weight);
                    h_ele_pt_eta[sam][1]->Fill(maxPt, maxEta, scale*_weight);
                    double dataMC_SF = 1.;
                    /*
                    if(sam < nSamples-1)
                        dataMC_SF = h_dataMC->GetBinContent(h_dataMC->GetXaxis()->FindBin(_n_MCTruth_PV));
                    else
                        dataMC_SF = 1;
                    */
                    distribs[17][sam]->Fill(TMath::Min(double(_n_PV),varMax[17]-0.1),scale*_weight*dataMC_SF);

            }
            

            if(_charges[leptInd[0]] == _charges[leptInd[1]] && _charges[leptInd[0]] != 0){
                if(sam != 2 && sam != 16){
                    if(_charges[leptInd[0]] == - _lchargemc[leptInd[0]] || _charges[leptInd[1]] == - _lchargemc[leptInd[1]]){
                    //if(_originReduced[leptInd[0]] == 0 && _originReduced[leptInd[1]] == 0){
                        
                        numberSSTotal++;
                        numberSS[sam][etaInd] += 1;
                        h_leadpt[sam][0]->Fill(maxPt, scale*_weight);
                        h_ele_pt_eta[sam][0]->Fill(maxPt, maxEta, scale*_weight);
                        //h_Mee[sam]->Fill(mdiL, scale*_weight);

                    }

                    //}
                } 
                else{
                    numberSSTotal++;
                    numberSS[sam][etaInd] += 1;
                    h_leadpt[sam][0]->Fill(maxPt, scale*_weight);
                    h_ele_pt_eta[sam][0]->Fill(maxPt, maxEta, scale*_weight);
                    //h_Mee[sam]->Fill(mdiL, scale*_weight);
                }
            }

            weight = _weight;

               
        }
        cout<<endl;
        std::cout << "Scale: " << scale << "; Weight: " << TMath::Abs(weight) << std::endl;
        std::cout << "SS total events: " << numberSSTotal * scale * TMath::Abs(weight) << std::endl;
        std::cout << "OS total events: " << numberOSTotal * scale * TMath::Abs(weight) << std::endl;
        
        //for(int i = 0; i < 4; i++)
        gh_ele_ptc[sam] = new TGraphAsymmErrors(h_leadpt[sam][0], h_leadpt[sam][1]);
            

    }
    cout<<endl;

    std::cout<<"Done"<<std::endl;
    
    string regions[4] = {"Cent-Cent", "Forw-Cent", "Cent-Forw", "Forw-Forw"}; 

    outFile << std::fixed << setprecision(5) << "\n";
    for(int i = 0; i < 4; i++){
        outFile << regions[i];
        for(int j=0; j < 4; j++)
            outFile << " & " << 1. * numberSS[j][i] / numberOS[j][i] << "$\\pm$ " << errCalc(numberSS[j][i], numberOS[j][i]) ;
        outFile << "\\\\ \\hline \n";
    }

    const char * titles[4] = {"Central-Central", "Forward-Central", "Central-Forward", "Forward-Forward"}; 
    
    TLegend* mtleg1 = new TLegend(0.2,0.88,0.55,0.60);
    mtleg1->SetFillColor(0);
    mtleg1->SetFillStyle(0);
    mtleg1->SetBorderSize(0);
    for (int i=0; i!=nSamples-1; ++i) {
        if(i != 1) continue;
        mtleg1->AddEntry(h_Mee[i],names[i],"f");
    }

    mtleg1->AddEntry(h_Mee[nSamples-1],names[nSamples-1],"l");

    /*
    TCanvas *c1 = new TCanvas("c1", "c1");
    //c1->Divide(2,2);
    //for(int i = 0; i < 4; i++){
        //c1->cd(i+1);
        gh_ele_ptc[1]->SetTitle("Leading lepton p_{T}");
        gh_ele_ptc[1]->GetXaxis()->SetTitle("p_{T}^{lead} [GeV]");
        gh_ele_ptc[1]->GetYaxis()->SetTitle("Q Mis-ID Rate");
        //gh_ele_ptc[0]->Draw("");
        gh_ele_ptc[1]->Draw();
        //gh_ele_ptc[2]->Draw("same");
        gh_ele_ptc[16]->Draw("same");
        mtleg1->Draw("same");
    //}
        */
    

    /*
    TCanvas *c2 = new TCanvas("c2", "c2");
    c2->Divide(2,2);
    for(int i = 0; i < 4; i++){
        c2->cd(i+1);
        h_deltaPt[i][1]->SetLineColor(kBlack);
        h_deltaPt[i][1]->SetFillColor(kBlack);
        h_deltaPt[i][1]->SetMarkerColor(kBlack);
        h_deltaPt[i][1]->DrawNormalized();
        h_deltaPt[i][0]->SetFillColor(kRed);
        h_deltaPt[i][0]->SetMarkerColor(kRed);
        h_deltaPt[i][0]->SetLineColor(kRed);
        h_deltaPt[i][0]->DrawNormalized("same");
    }
    */


    TLegend* mtleg = new TLegend(0.2,0.88,0.55,0.60);
    mtleg->SetFillColor(0);
    mtleg->SetFillStyle(0);
    mtleg->SetBorderSize(0);
    for (int i=0; i!=nSamples-1; ++i) {
        if(i != 0 && i != 1) continue;
        mtleg->AddEntry(h_Mee[i],names[i],"f");
    }

    mtleg->AddEntry(h_Mee[nSamples-1],names[nSamples-1],"l");

    
    TCanvas *c3 = new TCanvas("c3", "c3");
    //c3->Divide(2,2);
    //c3->cd(1);
    //h_Mee[nSamples-1]->SetMaximum(500);
    //h_Mee[nSamples-1]->Draw("P");
    //st_Mee->SetMaximum(1000);
    //h_Mee[nSamples-1]->GetXaxis()->SetTitle("M_{e^{-}e^{+}}");
    //h_Mee[nSamples-1]->GetYaxis()->SetTitle("events");

    showHist(c3, h_Mee[nSamples-1], st_Mee, "", "M_{ll}", "events", 1.6, mtleg);

    TCanvas* plotplot = new TCanvas("ptLepplot","ptLepplot",600,400);

    showHist(plotplot, distribs[17][nSamples-1], distribsST[17], "", "NPV", "events", 1.6, mtleg);

    /*
    c3->cd(2);
    
    //TCanvas *c4 = new TCanvas("c4", "c4");
    //h_eta[nSamples-1]->SetMaximum(500);
    //h_eta[nSamples-1]->Draw("P");
    //st_eta->SetMaximum(1000);
    //st_eta->GetXaxis()->SetTitle("#eta");
    //st_eta->GetYaxis()->SetTitle("events");
    showHist(st_eta,"","","Events / 1 ",1.);
    h_eta[nSamples-1]->SetLineWidth(2.);
    h_eta[nSamples-1]->SetFillColor(0);
    h_eta[nSamples-1]->SetLineColor(1);
    h_eta[nSamples-1]->SetLineStyle(1);
    h_eta[nSamples-1]->Draw("Psame");

    c3->cd(3);
    //st_2ndleadpt->GetXaxis()->SetTitle("leading lepton p_{T}, [GeV]");
    //st_2ndleadpt->GetYaxis()->SetTitle("events");
    showHist(st_2ndleadpt,"","","Events / 1 ",1.);
    h_2ndleadpt[nSamples-1]->SetLineWidth(2.);
    h_2ndleadpt[nSamples-1]->SetFillColor(0);
    h_2ndleadpt[nSamples-1]->SetLineColor(1);
    h_2ndleadpt[nSamples-1]->SetLineStyle(1);
    h_2ndleadpt[nSamples-1]->Draw("Psame");

    c3->cd(4);
    //st_trailpt->GetXaxis()->SetTitle("trailing lepton p_{T}, [GeV]");
    //st_trailpt->GetYaxis()->SetTitle("events");
    showHist(st_trailpt,"","","Events / 1 ",1.);
    h_trailpt[nSamples-1]->SetLineWidth(2.);
    h_trailpt[nSamples-1]->SetFillColor(0);
    h_trailpt[nSamples-1]->SetLineColor(1);
    h_trailpt[nSamples-1]->SetLineStyle(1);
    h_trailpt[nSamples-1]->Draw("Psame");
    */

    //mtleg->Draw("same");
    

    /*
    Print2D(h_ele_pt_eta[0][0], h_ele_pt_eta[0][1], "ttbar");
    Print2D(h_ele_pt_eta[1][0], h_ele_pt_eta[1][1], "DY");
    Print2D(h_ele_pt_eta[2][0], h_ele_pt_eta[2][1], "DY_onZ");
    Print2D(h_ele_pt_eta[3][0], h_ele_pt_eta[3][1], "data");
    */
    
   
    std::cout<<"Done 2"<<std::endl;
    
}


void showHist(TVirtualPad* c1, TH1D *hist, THStack *stack, string title, string titleX, string titleY, double num, TLegend *leg){   

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetTopMargin(0.1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    
    double xmin = hist->GetXaxis()->GetXmin();
    double xmax = hist->GetXaxis()->GetXmax();
    //pad1->DrawFrame(xmin, -0.1, xmax, 1.1);
    
    hist->SetTitle(title.c_str());
    //hist->GetXaxis()->SetTitle(titleX.c_str());
    hist->GetYaxis()->SetTitle(titleY.c_str());
    hist->SetMaximum(hist->GetMaximum() * num);
    hist->SetLineWidth(1);
    hist->SetFillColor(0);
    hist->SetLineColor(1);
    hist->SetLineStyle(1);
    hist->SetMarkerSize(0.5);
    hist->SetMarkerStyle(20);
    //hist->GetXaxis()->SetLabelSize(0.15);
    //hist->GetXaxis()->SetTitleSize(0.08);
    hist->GetYaxis()->SetLabelSize(0.08);
    hist->GetYaxis()->SetTitleSize(0.08);
    hist->GetYaxis()->SetTitleOffset(0.6);
    hist->Draw("e");
    stack->Draw("histsame");
    hist->Draw("esame");
    /*
    pad1->cd();
    pad1->Update();
    pad1->RedrawAxis();
    pad1->GetFrame()->Draw();
    */
    leg->Draw("same");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);  

    latex.SetTextFont(42);
    latex.SetTextAlign(31); 
    latex.SetTextSize(0.07);    
    latex.DrawLatex(0.98, 0.93,lumi_13TeV + "(13 TeV)");

    TLatex cmsText;
    cmsText.SetNDC();
    cmsText.SetTextAngle(0);
    cmsText.SetTextColor(kBlack);  

    cmsText.SetTextFont(61);
    cmsText.SetTextAlign(31); 
    cmsText.SetTextSize(0.08);  
    cmsText.DrawLatex(0.2, 0.93,"CMS");


    TLatex extraText;
    extraText.SetNDC();
    extraText.SetTextAngle(0);
    extraText.SetTextColor(kBlack);  

    extraText.SetTextFont(52);
    extraText.SetTextAlign(31); 
    extraText.SetTextSize(0.07);  
    extraText.DrawLatex(0.35, 0.93,"Preliminary");

    

    c1->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.4);
    pad2->Draw();
    pad2->cd();
    TH1D * histcopy = (TH1D*)hist->Clone("histcopy");
    //TH1D * histo_stack = (TH1D*)stack->GetHistogram();
    //histcopy->Sumw2();
    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);
    histcopy->SetStats(0);

    histcopy->SetTitle("");
    histcopy->GetXaxis()->SetTitle(titleX.c_str());
    histcopy->GetYaxis()->SetTitle("data/pred");

    histcopy->GetXaxis()->SetLabelSize(0.1);
    histcopy->GetYaxis()->SetLabelSize(0.12);
    histcopy->GetXaxis()->SetTitleSize(0.1);
    histcopy->GetYaxis()->SetTitleSize(0.15);

    histcopy->GetYaxis()->SetTitleOffset(0.2);

    histcopy->SetMaximum(1.5);
    histcopy->SetMinimum(0.5);
    //
    histcopy->SetMarkerStyle(1);

    TList *histos = stack->GetHists();
    TIter next(histos);
    TH1D *histo_stack = (TH1D*)stack->GetHistogram();;
    while ((hist =(TH1D*)next())) {
      //cout << "Adding " << hist->GetName() << endl;
      histo_stack->Add(hist);
    }

    histcopy->Divide(histo_stack);

    histcopy->Draw("e");
    line->Draw("same");
    /*
    hist2->SetLineWidth(2.);
    hist2->SetFillColor(0);
    hist2->SetLineColor(1);
    hist2->SetLineStyle(3);
    hist2->Draw("hist same");
    hist2->Draw("axis same");
    */

}

int main(int argc, char *argv[]){

    //TApplication *rootapp = new TApplication("example", &argc, argv);

    readTree();

    //rootapp->Run();

    return 0;
}

double calculateZbi(double signal, double bkg, double unc){

 double n_on = signal+bkg;
 double mu_b_hat=bkg;
 double sigma_b=unc*bkg;
 double tau = mu_b_hat/(sigma_b*sigma_b);
 double n_off = tau*mu_b_hat;
 double P_Bi = TMath::BetaIncomplete(1./(1.+tau),n_on,n_off+1);
 double Z_Bi = sqrt(2)*TMath::ErfInverse(1 - 2*P_Bi);
 return Z_Bi;
 //  std::cout<<"The calculated Zbi for a signal of "<<signal<<" events and background of "<<bkg<<" events with a systematic uncertainty of "<<unc*100<<"% is "<<Z_Bi<<std::endl;

}

double errCalc(int a, int b){
    return TMath::Sqrt(a / TMath::Power(b, 2) + TMath::Power(a, 2)/ TMath::Power(b, 3));    
    //TMath::Sqrt(1. * a / (b * b) + a * a / (b * b * b)) ;
}

TH2D TranslateHisto(const TH2 &input){
  int nx = input.GetNbinsX();
  int ny = input.GetNbinsY();
  TH2D output(input.GetName(), input.GetTitle(), nx, 0.5, nx+0.5, ny, 0.5, ny+0.5);
  output.Sumw2();
  output.SetStats(false);
  output.SetMarkerSize(2);
  output.SetLabelSize(0.05,"XYZ");
  output.SetTitleSize(0.05,"XYZ");
  for(int ix = 0; ix <= nx+1; ++ix){
    for(int iy = 0; iy <= ny+1; ++iy){
      output.SetBinContent(ix, iy, input.GetBinContent(ix, iy));
      output.SetBinError(ix, iy, input.GetBinError(ix, iy));
    }
  }
  
  for(int ix = 1; ix <= nx; ++ix){
    const TAxis *iaxis = input.GetXaxis();
    TAxis *oaxis = output.GetXaxis();
    if(iaxis == NULL || oaxis == NULL) continue;
    oaxis->SetTitle(iaxis->GetTitle());
    ostringstream oss;
    oss << iaxis->GetBinLowEdge(ix) << "-" << iaxis->GetBinUpEdge(ix) << flush;
    oaxis->SetBinLabel(ix, oss.str().c_str());
  }

  for(int iy = 1; iy <= ny; ++iy){
    const TAxis *iaxis = input.GetYaxis();
    TAxis *oaxis = output.GetYaxis();
    if(iaxis == NULL || oaxis == NULL) continue;
    oaxis->SetTitle(iaxis->GetTitle());
    ostringstream oss;
    oss << iaxis->GetBinLowEdge(iy) << "-" << iaxis->GetBinUpEdge(iy) << flush;
    oaxis->SetBinLabel(iy, oss.str().c_str());
  }

  return output;
}

void Print2D(TH2 const * const h_data_in, TH2 const * const h_mc_in, const TString &ext){
  if(h_data_in == NULL || h_mc_in == NULL) return;

  TH2D h_data = TranslateHisto(*h_data_in);
  TH2D h_mc = TranslateHisto(*h_mc_in);

  TCanvas canvas;
  gStyle->SetPalette(bands, rainbow);
  h_data.SetMarkerSize(2);
  h_mc.SetMarkerSize(2);

  h_data.SetTitle("Fullsim Eff");
  h_data.Draw("colz");
  h_data.Draw("textesame");
  canvas.Print("plots/2d_full_"+ext+".pdf");
  //canvas.Print("plots/2d_data_"+ext+".png");
  //PrintTable(&h_data, "data_"+ext);
  h_mc.SetTitle("Fastsim Eff");
  h_mc.Draw("colz");
  h_mc.Draw("textesame");
  canvas.Print("plots/2d_fast_"+ext+".pdf");
  //canvas.Print("plots/2d_mc_"+ext+".png");
  //PrintTable(&h_mc, "mc_"+ext);

  gStyle->SetPalette(bands, patriotic);
  h_data.Divide(&h_mc);
  h_data.SetMinimum(0.8);
  h_data.SetMaximum(1.25);
  canvas.SetLogz();
  h_data.SetTitle("Flip rate");
  h_data.Draw("colz");
  h_data.Draw("textesame");
  canvas.Print("plots/flip_"+ext+".pdf");
  //canvas.Print("plots/sf_"+ext+".png");
  //PrintTable(&h_data, "sf_"+ext);
}

void PrintTable(TH2 const * const histo, const TString &ext){
  int eta = ext.Index("abseta");
  int pt = ext.Index("_et_");
  if(eta<0 || pt<0 || pt<eta) return;

  ofstream file("tables/"+ext+".tex");
  file << "\\documentclass{article}\n\n";

  file << "\\begin{document}\n";
  file << "\\begin{table}\n";
  file << "  \\begin{tabular}{r|rrrrr}\n";
  file << "    \\hline\\hline\n";
  file << "    & \\multicolumn{5}{c}{$p_T$ [GeV]}\\\\\n";
  file << "    $|\\eta|$ & 10-20 & 20-30 & 30-40 & 40-50 & 50-200\\\\\n";
  file << "    \\hline\n";
  PrintLine(file, histo, 1, "0-1.442");
  PrintLine(file, histo, 2, "1.442-1.556");
  PrintLine(file, histo, 3, "1.556-2.5");
  file << "    \\hline\\hline\n";
  file << "  \\end{tabular}\n";
  file << "\\end{table}\n";
  file << "\\end{document}\n";

  file.flush();
  file.close();
}

void PrintLine(ofstream &file, TH2 const * const histo, int bin, const TString &label){
  if(!file.is_open()) return;

  if(histo == NULL || histo->GetNbinsX() < bin || histo->GetNbinsY()<5) return;
  file << "    " << label;
  for(int y = 1; y <= 5; ++y){
    file << " & $"
         << fixed << setprecision(3) << histo->GetBinContent(bin, y)
         << "\\pm"
         << fixed << setprecision(3) << histo->GetBinError(bin, y)
         << "$";
  }
  file << "\\\\\n";
}

void GetPatrioticPalette(){
  const unsigned num = 3;
  double red[num] = {0., 1., 1.};
  double green[num] = {0., 1., 0.};
  double blue[num] = {1., 1., 0.};
  double stops[num] = {0., 0.5, 1.};
  int fi = TColor::CreateGradientColorTable(num,stops,red,green,blue,bands);
  for(int i = 0; i < bands; ++i){
    patriotic[i] = fi + i;
  }
}

void GetRainbowPalette(){
  const unsigned num = 6;
  double red[num] =   {1.,0.,0.,0.,1.,1.};
  double green[num] = {0.,0.,1.,1.,1.,0.};
  double blue[num] =  {1.,1.,1.,0.,0.,0.};
  double stops[num] = {0.,0.2,0.4,0.6,0.8,1.};
  int fi = TColor::CreateGradientColorTable(num,stops,red,green,blue,bands);
  for(int i = 0; i < bands; ++i){
    rainbow[i] = fi+i;
  }
}
