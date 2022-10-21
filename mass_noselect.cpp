//get a mass histogram of three-pion triplets with only charge selection
//this motivates the selection we implemented
//does not include electrons

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <ROOT/RLogger.hxx>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "Math/Vector4D.h"
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <stdlib.h>
#include <math.h>

template <typename T> 
using Vec = ROOT::RVec<T>;
using ROOT::Math::XYZTVector;
using ROOT::Math::PtEtaPhiMVector;

float minmass = 60;
float maxmass = 100;

auto maketriplets(Vec<int> pdgids, Vec<int> charges, Vec<float> pts, Vec<float> etas, Vec<float> phis, Vec<float> masses) {
    Vec<Vec<std::size_t>> triplets; //stores all passing triplets (best one selected at the end)
    Vec<std::size_t> ix; //pion indeces

    for (std::size_t i = 0; i < pdgids.size(); ++i) { //make list of all hadrons
        if (abs(pdgids[i]) == 211) ix.push_back(i);
    }
   if (ix.size() > 2) { //if there are 3+ pions
   //make the triplets
       for (std::size_t i1 = 0; i1 < ix.size(); ++i1) {
           for (std::size_t i2 = 0; i2 < ix.size(); ++i2) {
                if (i2 <= i1) continue;
               for (std::size_t i3 = 0; i3 < ix.size(); ++i3) {
                   if (i3 <= i2) continue;
                    Vec<std::size_t> tr{ix[i1], ix[i2], ix[i3]}; //triplet of indeces

                    if (abs(charges[tr[0]]*charges[tr[1]]*charges[tr[2]]) == 1) {
                        //make Lorentz vectors for each triplet
                        ROOT::Math::PtEtaPhiMVector p1(pts[tr[0]],etas[tr[0]],phis[tr[0]],masses[tr[0]]);
                        ROOT::Math::PtEtaPhiMVector p2(pts[tr[1]],etas[tr[1]],phis[tr[1]],masses[tr[1]]);
                        ROOT::Math::PtEtaPhiMVector p3(pts[tr[2]],etas[tr[2]],phis[tr[2]],masses[tr[2]]);
                        auto mass = (p1 + p2 + p3).M();
                        if (mass >= minmass and mass <= maxmass) { //MASS test
                            triplets.push_back(tr);
                        } // mass cut  
                    } //charge cut 
               } //i3
           } //i2
       } //i1
   } //if 3 or more pions

return triplets;
}

bool notempty(Vec<Vec<std::size_t>> index) { //used to check if any triplets passed in an event
    return !index.empty();
}

auto tripletmass(Vec<Vec<std::size_t>> triplets, Vec<float> pts, Vec<float> etas, Vec<float> phis, Vec<float> masses) {
    Vec<float> massvec;
    for (Vec<std::size_t> t : triplets) {
        ROOT::Math::PtEtaPhiMVector p1(pts[t[0]],etas[t[0]],phis[t[0]],masses[t[0]]);
        ROOT::Math::PtEtaPhiMVector p2(pts[t[1]],etas[t[1]],phis[t[1]],masses[t[1]]);
        ROOT::Math::PtEtaPhiMVector p3(pts[t[2]],etas[t[2]],phis[t[2]],masses[t[2]]);
        double mass = (p1 + p2 + p3).M();
        massvec.push_back(mass);
    }
    return massvec;
}

int main() {
    //ROOT::EnableImplicitMT();
    //ROOT::RDataFrame d("Events", "/eos/cms/store/cmst3/group/l1tr/gpetrucc/w3pi/ntuples/v1/l1Nano_WTo3Pion_PU200.v1_more.root");
    ROOT::RDataFrame d("Events", "/eos/cms/store/cmst3/group/l1tr/gpetrucc/w3pi/ntuples/v1/l1Nano_WTo3Pion_PU200.v1.root");
    auto d2 = d
        .Define("Triplet_Index", maketriplets, {"L1Puppi_pdgId", "L1Puppi_charge", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"})
        .Filter(notempty,{"Triplet_Index"})
        .Define("Triplet_Mass",tripletmass,{"Triplet_Index", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"});

    auto masshist = d2.Histo1D({"masshist", "W Boson mass from selected pions; mass (GeV/c^2)", 100, 0, 100},"Triplet_Mass");
    auto c1 = new TCanvas("c1","c1",800,600);
    masshist->Draw();
    c1->Draw();
    c1->SaveAs("figures/masshist_noselect.pdf");
}