//
// Created by ming on 6/4/22.
//

//#include "TCanvas.h"
//#include "TH2D.h"
//#include "TTree.h"
//#include "TFile.h"
//#include "TNtuple.h"
//
//void plotNtuple(){
//
//    TCanvas* c1 = new TCanvas("c1", "G4_hw");
//    TFile f("output.root");
//    TNtuple* ntuple = (TNtuple*)f.Get("data");
//    c1->cd();
//    ntuple->Draw("x:y:edep","2D");
//}

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TF1.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "iostream"

Double_t data[1000000]={0};
void bragg()
{
    //%jsroot on
    //  TH1D *tdiff=new TH1D("tdiff","td-tu",800,-1,2);
    TGraph *tdiff = new TGraph();
    TCanvas *c1=new TCanvas("c1","c1");
    //1.打开文件，得到TTree指针
    //TFile *ipf=new TFile("tree.root");//root -l tree.root
//    TFile *ipf = new TFile("12C_10kparticles.root");
//    TFile *ipf = new TFile("gamma_10kparticles.root");
//    TFile *ipf = new TFile("neutron_10kparticles.root");
    TFile *ipf = new TFile("proton_10kparticles.root");
//    TFile *ipf = new TFile();
    ipf->ls();//ROOT 环境下 > .ls
    //观察tree的结构
    ipf->cd();
    TTree *tree=(TTree*)ipf->Get("data");//得到名字为“tree”的TTree指针
//    TNtuple* Ntuple = (Ntuple*)ipf->Get("data");
    //TNtuple *ntuple=(TNtuple*)ipf->Get("B4");

    //2. 声明tree的Branch变量
    //Double_t x;
    //Double_t y;
    Double_t edepStep;
    Double_t Egap;
    Double_t E;
    Double_t stepx;
    Double_t stepy;
    Double_t stepz;
    Double_t energy7;
    Double_t energy4;
    double e = 0;

//    tree->SetBranchAddress("edep",&edepStep);
//    tree->SetBranchAddress("x",&stepx);
//    tree->SetBranchAddress("y",&stepy);
//    tree->SetBranchAddress("z",&stepz);

    tree->SetBranchAddress("edep2",&edepStep);
    tree->SetBranchAddress("x2",&stepx);
    tree->SetBranchAddress("y2",&stepy);
    tree->SetBranchAddress("z2",&stepz);

    /*tree->SetBranchAddress("Eabs",&Eabs);
    tree->SetBranchAddress("Egap",&Egap);
    tree->SetBranchAddress("E",&E);
    tree->SetBranchAddress("xStep",&xStep);
    tree->SetBranchAddress("yStep",&yStep);
    tree->SetBranchAddress("zStep",&zStep);
    tree->SetBranchAddress("energy7",&energy7);*/
    //  tree->SetBranchAddress("energy4",&energy4);

    Long64_t nentries=tree->GetEntries();//得到事件总数
    for(Long64_t jentry = 0; jentry<nentries; jentry++) {//对每个事件进行遍历
        tree->GetEntry(jentry);
//        data[int (stepx + 10)]+=edepStep;
        data[int (stepy + 10)]+=edepStep;
//        data[int (stepz + 10)]+=edepStep;
        //tdiff->AddPoint(0,data[int(xStep)]);
    }
    for(int i=0;i<901;i++){
        tdiff->AddPoint(i - 10,data[i]);

    }
    tdiff->Draw();
    c1->Draw();


}