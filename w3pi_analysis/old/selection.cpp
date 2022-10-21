
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

//all our cut values and other parameters
    double minpt1 = 9;
    double minpt2 = 15;
    double minpt3 = 20;
    double mindeltar = 0.5;
    double minmass = 60;
    double maxmass = 100;
    double mindr = 0.01; //size of cone where we calculate isolation, not angular separation of pions
    double maxdr = 0.25;
    double maxiso = 0.4; 

//TEST functions
bool isolation(int pidex, Vec<double> eta, Vec<double> phi, Vec<double> pt, double mindr, double maxdr, double maxiso) {
    bool passed=false;
    double psum = 0;
    for (int j = 0u; j < pt.size(); ++j) { //loop over other particles
        if (pidex == j) continue;
        //auto dr = ROOT::VecOps::DeltaR(eta[pidex], eta[j], phi[pidex], phi[j]);
        double dr = sqrt(pow(eta[pidex]-eta[j],2)+pow(phi[pidex]-phi[j],2));
        if (dr >= mindr && dr <= maxdr) psum += pt[j];
    }
    if (psum/pt[pidex] <= maxiso) passed=true;
    return passed;
}

bool deltar(ROOT::Math::PtEtaPhiMVector p1, ROOT::Math::PtEtaPhiMVector p2, ROOT::Math::PtEtaPhiMVector p3,
              double mindeltar) {
    bool passed=true;
    Vec<ROOT::Math::PtEtaPhiMVector> pis{p1, p2, p3};

    for (int i = 0; i<3; ++i) {
        int j;
        if (i < 2) j = i + 1;
        else j = 0;
        double deltaeta = pis[i].Eta() - pis[j].Eta();
        double deltaphi = pis[i].Phi() - pis[j].Phi();
        double dr = sqrt(pow(deltaeta, 2.0) + pow(deltaphi, 2.0));
        if (dr < mindeltar) passed=false;
    }
    return passed;
} 

bool initptcut(Vec<float> pts, Vec<int> pdgids) {
    bool pass;
    int intermediatecut = 0;
    int highcut = 0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (abs(pdgids[i]) == 11 or abs(pdgids[i]) == 211) {
            if (pts[i] >= minpt2) intermediatecut += 1;
            if (pts[i] >= minpt3) highcut += 1;
        }
    }
    if (intermediatecut > 1 and highcut > 0) pass = true;
    else pass = false;
    return pass;
}

//generates combinations of pion indeces after pt lists are made
//high: list of pion indeces with pt > min largest
//inter: list of pion indeces with pt > min intermediate
//low: list of pion indeces with pt > min smallest
auto combiner(Vec<std::size_t> high, Vec<std::size_t> inter, Vec<std::size_t> low) {
    using namespace ROOT::VecOps;
    Vec<Vec<std::size_t>> triplets;
    Vec<Vec<std::size_t>> c1 = Combinations(high, inter);
    Vec<Vec<std::size_t>> c2 = Combinations(c1[0].size(),low.size());
    for (std::size_t i = 0; i < c2[0].size(); ++i) {
        Vec<std::size_t> triplet{high[c1[0][c2[0][i]]], inter[c1[1][c2[0][i]]], low[c2[1][i]]};
        if (triplet[0] != triplet[1] and triplet[1] != triplet[2] and triplet[2] != triplet[0]) {
            triplets.push_back(triplet);
        }
    }
    return triplets;
}

//SELECTION
//generates list of pion triplets passing all tests
auto maketriplets(Vec<int> pdgids, Vec<int> charges, Vec<float> pts, Vec<float> etas, Vec<float> phis, Vec<float> masses) {
    //cutoff values and other parameters
    
    Vec<Vec<std::size_t>> triplets;
    Vec<std::size_t> indeces;
    Vec<double> masspass;
    Vec<double> ptsums;
    Vec<bool> isopasses;
    Vec<std::size_t> isoindex;

    for (std::size_t i = 0; i < pdgids.size(); ++i) { //make list of all hadrons
        if ((abs(pdgids[i]) == 211 or abs(pdgids[i]) == 11) and pts[i] >= minpt1) indeces.push_back(i);
    }
   if (indeces.size() > 2) { //if there are 3+ pions
       float ptsum;
       
       //make lists of pts passing pt1, pt2, and pt3 cuts
       Vec<std::size_t> highpt;
       Vec<std::size_t> intermediatept;
       Vec<std::size_t> lowpt;
       for (std::size_t k : indeces) {
           if (pts[k] >= minpt1) lowpt.push_back(k);
           if (pts[k] >= minpt2) intermediatept.push_back(k);
           if (pts[k] >= minpt3) highpt.push_back(k);
       }
       
       if (highpt.size() > 0 and intermediatept.size() > 1 and lowpt.size() > 2) {//PT test
       
       Vec<Vec<std::size_t>> tr = combiner(highpt, intermediatept, lowpt);
       for (std::size_t i = 0; i < tr.size(); ++i) { //go through list of triplets
           Vec<std::size_t> triplet = tr[i];
           if (abs(charges[triplet[0]] + charges[triplet[1]] + charges[triplet[2]]) == 1) { //CHARGE test
                   //make Lorentz vectors for each triplet
                   ROOT::Math::PtEtaPhiMVector p1(pts[triplet[0]],etas[triplet[0]],phis[triplet[0]],masses[triplet[0]]);
                   ROOT::Math::PtEtaPhiMVector p2(pts[triplet[1]],etas[triplet[1]],phis[triplet[1]],masses[triplet[1]]);
                   ROOT::Math::PtEtaPhiMVector p3(pts[triplet[2]],etas[triplet[2]],phis[triplet[2]],masses[triplet[2]]);
                   auto mass = (p1 + p2 + p3).M();
                   if (mass >= minmass and mass <= maxmass) { //MASS test
                       if (deltar(p1, p2, p3, mindeltar) == true) { //ANGULAR SEPARATION test
                           //ISOLATION test
                           bool isop = true;
                           for (int j = 0; j < triplet.size(); ++j) {
                               auto isopos_ = find(isoindex.begin(), isoindex.end(), triplet[j]);
                               if (isopos_ != isoindex.end()) { //check if isolation has already been calculated
                                   int isopos = std::distance(isoindex.begin(), isopos_);
                                   if (isopasses[isopos] == false) {
                                       isop = false;
                                       break;
                                   }
                               }
                               else { //calculate isolation and add true or false to isopasses
                                   isoindex.push_back(triplet[j]); //add particle index to isolation index list
                                   if (isolation(triplet[j], etas, phis, pts, mindr, maxdr, maxiso) == true) {
                                       isopasses.push_back(true);
                                   }
                                   else {
                                       isopasses.push_back(false);
                                       isop = false;
                                       break;
                                   }
                               }
                           } 
                           if (isop == true) {

                            triplets.push_back(triplet);
                            masspass.push_back(mass);
                            ptsum = p1.Pt() + p2.Pt() + p3.Pt();
                            ptsums.push_back(ptsum);
                         } // iso
                        } // delta R
                   } // mass
           } // charge
       } // i (for each triplet)
       } //pt test
   } // if there are 3+ pions
   // if (choosebest == "on") {
        Vec<Vec<std::size_t>> ret;
        if (triplets.size() > 1) { //choose the best triplet for each event
            Vec<double> scores;
            int best;
            for (int i = 0; i < triplets.size(); ++i) {
                double score = masspass[i]*ptsums[i];
                scores.push_back(score);
                if (i == 0)
                if (score > scores[i-1]) best = i; //i-1 works fine bc the [-1] element is 0 and this is all positive
            }
            ret = {triplets[best]};
            return ret;
        }
    else return triplets;
   // }
}

//processing after selection: calculate mass and check the triplet isn't empty

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
    //increase verbosity to see how long this is taking
    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);

//ask for event type: single neutrino, w-> 3 pi, ttbar
std::cout << "Please select an option by entering 1, 2, or 3:" << std::endl;
std::cout << "1. single neutrino" << std::endl;
std::cout << "2. W->3 pion " << std::endl;
std::cout << "3. TTbar" << std::endl;
int n;
std::cin >> n;
if (n != 1 and n!=2 and n!=3) {
	std::cout << "Invalid option :(" << std::endl;
	return 0;
	}
    //ROOT::EnableImplicitMT();
    //W->3pi    
if (n == 2){
 ROOT::RDataFrame d("Events", "/eos/cms/store/cmst3/group/l1tr/gpetrucc/w3pi/ntuples/v1/l1Nano_WTo3Pion_PU200.v1_more.root");
    auto d2 = d.Filter(initptcut, {"L1Puppi_pt","L1Puppi_pdgId"})
        .Define("Triplet_Index", maketriplets, {"L1Puppi_pdgId", "L1Puppi_charge", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"})
        .Filter(notempty,{"Triplet_Index"})
        .Define("Triplet_Mass",tripletmass,{"Triplet_Index", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"}); //simulates analysis on 1 orbit

    auto masshist = d2.Histo1D({"masshist", "W Boson mass from selected pions; mass (GeV/c^2)", 100, 0, 100},"Triplet_Mass");
    auto c1 = new TCanvas("c1","c1",800,600);
    masshist->Draw();
    c1->Draw();
    c1->SaveAs("figures/masshist.pdf");
}
    //ttbar  
    if (n == 3) {
	ROOT::RDataFrame d("Events", "/eos/cms/store/cmst3/group/l1tr/gpetrucc/w3pi/ntuples/v1/l1Nano_TTbar_PU200.v1.root");
    auto d2 = d.Filter(initptcut, {"L1Puppi_pt","L1Puppi_pdgId"})
        .Define("Triplet_Index", maketriplets, {"L1Puppi_pdgId", "L1Puppi_charge", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"})
        .Filter(notempty,{"Triplet_Index"})
        .Define("Triplet_Mass",tripletmass,{"Triplet_Index", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"}); //simulates analysis on 1 orbit

    auto masshist = d2.Histo1D({"masshist", "W Boson mass from selected pions; mass (GeV/c^2)", 100, 0, 100},"Triplet_Mass");
    auto c1 = new TCanvas("c1","c1",800,600);
    masshist->Draw();
    c1->Draw();
    c1->SaveAs("figures/masshist.pdf");
}    
//single neutrino
    if (n == 1) {
	 ROOT::RDataFrame d("Events", "/eos/cms/store/cmst3/group/l1tr/gpetrucc/w3pi/ntuples/v1/l1Nano_SingleNeutrino_PU200.v1_more.root");
    auto d2 = d.Filter(initptcut, {"L1Puppi_pt","L1Puppi_pdgId"})
        .Define("Triplet_Index", maketriplets, {"L1Puppi_pdgId", "L1Puppi_charge", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"})
        .Filter(notempty,{"Triplet_Index"})
        .Define("Triplet_Mass",tripletmass,{"Triplet_Index", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"}); //simulates analysis on 1 orbit
    
    auto masshist = d2.Histo1D({"masshist", "W Boson mass from selected pions; mass (GeV/c^2)", 100, 0, 100},"Triplet_Mass");
    auto c1 = new TCanvas("c1","c1",800,600);
    masshist->Draw();
    c1->Draw();
    c1->SaveAs("figures/masshist.pdf");
}
}
