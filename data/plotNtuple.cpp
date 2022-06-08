//
// Created by ming on 6/7/22.
//

#include "TCanvas.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TNtuple.h"

void plotNtuple(){

    TCanvas* c1 = new TCanvas("c1", "G4_hw");
//    TFile f("12C_10kparticles.root");
//    TFile f("gamma_10particles.root");
//    TFile f("neutron_10kParticles.root");
//    TFile f("proton_10kParticles.root");
    TFile f("output.root");
    TNtuple* ntuple = (TNtuple*)f.Get("data");
    c1->cd();
//    ntuple->Draw("z:y:x:edep");
    ntuple->Draw("z2:y2:x2:edep2");
}
