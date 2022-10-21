//reads files saved w/Snapshot
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
    float minpt1 = 9; //9
    float minpt2 = 15;
    float minpt3 = 20;
    float mindeltar = 0.5;
    float minmass = 60;
    float maxmass = 100;
    float mindr = 0.01; //size of cone where we calculate isolation, not angular separation of pions
    float maxdr = 0.25;
    float maxiso = 0.4; //0.4

//TEST functions
bool isolation(std::size_t pidex, Vec<float> eta, Vec<float> phi, Vec<float> pt, float mindr, float maxdr, float maxiso) {
    bool passed=false;
    float psum = 0;
    for (int j = 0u; j < pt.size(); ++j) { //loop over other particles
        if (pidex == j) continue;
        //auto dr = ROOT::VecOps::DeltaR(eta[pidex], eta[j], phi[pidex], phi[j]);
        float dr = sqrt(pow(eta[pidex]-eta[j],2)+pow(phi[pidex]-phi[j],2));
        if (dr >= mindr && dr <= maxdr) psum += pt[j];
    }
    if (psum/pt[pidex] <= maxiso) passed=true;
    return passed;
}

bool deltar(std::size_t pi1, std::size_t pi2, float eta1, float eta2, float phi1, float phi2, float mindeltar) {
    bool passed=true;
        float deta = eta1 - eta2;
        float dphi = phi1 - phi2;
        float dr2 = deta*deta + dphi*dphi;
        if (dr2 < mindeltar*mindeltar) {
            passed = false;
            return passed;
        }
    return passed;
} 

//speeds up code to have an initial filter that does not read many branches
//this eliminates most of single neutrino events
bool initptcut(Vec<float> pts, Vec<int> pdgids) {
    bool pass;
    int lowcut = 0;
    int intermediatecut = 0;
    int highcut = 0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (abs(pdgids[i]) == 11 or abs(pdgids[i]) == 211) {
            if (pts[i] >= minpt1) {
                lowcut +=1;
                if (pts[i] >= minpt2) {
                    intermediatecut += 1;
                    if (pts[i] >= minpt3) highcut += 1;
                }
            }
        }
    }
    if (lowcut > 2 and intermediatecut > 1 and highcut > 0) pass = true;
    else pass = false;
    return pass;
}

//no longer used
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
    Vec<Vec<std::size_t>> triplets; //stores all passing triplets (best one selected at the end)
    Vec<std::size_t> ix; //pion indeces
    Vec<float> masspass;
    Vec<float> ptsums;
    Vec<std::size_t> iso(pdgids.size(),0); //stores whether a particle passes isolation test so we don't calculate reliso twice

    for (std::size_t i = 0; i < pdgids.size(); ++i) { //make list of all hadrons
        if ((abs(pdgids[i]) == 211 or abs(pdgids[i]) == 11) and pts[i] >= minpt1) ix.push_back(i);
    }
   if (ix.size() > 2) { //if there are 3+ pions
       float ptsum;
       
       for (std::size_t i1 = 0; i1 < ix.size(); ++i1) {
           if (pts[ix[i1]] < minpt3) continue; //high pt cut
           if (isolation(ix[i1], etas, phis, pts, mindr, maxdr, maxiso) == 0) continue; //check iso of high pt pion
           for (std::size_t i2 = 0; i2 < ix.size(); ++i2) {
               if (i2 == i1 || pts[ix[i2]] < minpt2) continue;
               if (pts[ix[i2]] > pts[ix[i1]] || (pts[ix[i2]] == pts[ix[i1]] and i2 < i1)) continue; //intermediate pt cut
               if (deltar(ix[i1], ix[i2], etas[ix[i1]], etas[ix[i2]], phis[ix[i1]], phis[ix[i2]], mindeltar) == false) continue; //angular sep of top 2 pions
               for (std::size_t i3 = 0; i3 < ix.size(); ++i3) {
                   if (i3 == i1 or i3 == i2) continue;
                   if (pts[ix[i2]] < minpt1) continue; //low pt cut
                   if (pts[ix[i3]] > pts[ix[i1]] || (pts[ix[i3]] == pts[ix[i1]] and i3 < i1)) continue;
                   if (pts[ix[i3]] > pts[ix[i2]] || (pts[ix[i3]] == pts[ix[i2]] and i3 < i2)) continue;
                    Vec<std::size_t> tr{ix[i1], ix[i2], ix[i3]}; //triplet of indeces

                    if (abs(charges[tr[0]] + charges[tr[1]] + charges[tr[2]]) == 1) {
                        //make Lorentz vectors for each triplet
                        ROOT::Math::PtEtaPhiMVector p1(pts[tr[0]],etas[tr[0]],phis[tr[0]],masses[tr[0]]);
                        ROOT::Math::PtEtaPhiMVector p2(pts[tr[1]],etas[tr[1]],phis[tr[1]],masses[tr[1]]);
                        ROOT::Math::PtEtaPhiMVector p3(pts[tr[2]],etas[tr[2]],phis[tr[2]],masses[tr[2]]);
                        auto mass = (p1 + p2 + p3).M();
                        if (mass >= minmass and mass <= maxmass) { //MASS test
                            if (deltar(ix[i1], ix[i3], etas[ix[i1]], etas[ix[i3]], phis[ix[i1]], phis[ix[i3]], mindeltar) == true
                                and deltar(ix[i2], ix[i3], etas[ix[i2]], etas[ix[i3]], phis[ix[i2]], phis[ix[i3]], mindeltar) == true) { //ANGULAR SEPARATION test
                                //ISOLATION test for lower 2 pions
                                bool isop = true;
                                for (int j = 1; j < 3; ++j) {
                                    if (iso[tr[j]] == 0) {
                                        if (isolation(tr[j], etas, phis, pts, mindr, maxdr, maxiso) == false) {
                                            iso[tr[j]] = 2;
                                            isop = false;
                                            break;
                                        }
                                        else {
                                            iso[tr[j]] = 1;
                                            }
                                    }
                                    if (iso[tr[j]] == 2) {
					                    isop = false;
					                    break; //fail triplet if one bad isolation 
                                    } 
                                }
                                    if (isop == true) {
                                        triplets.push_back(tr);
                                        masspass.push_back(mass);
                                        ptsum = pts[tr[0]]+pts[tr[1]]+pts[tr[2]];
                                        ptsums.push_back(ptsum);
                                    } // iso
                                } // delta R
                        } // mass   
                    } //charge
               } //low pt cut
           } //intermediate pt cut
       } //high pt cut
   } //if 3 or more pions
    if (triplets.empty()) {
        Vec<std::size_t> ret;
        return ret;
    }
    if (triplets.size() == 1) return triplets[0];
    
    //if there are multiple triplets passing, choose the best
    float bestscore = 0; 
    std::size_t best = 0; //index of best triplet in triplet array
    for (std::size_t i = 0; i < triplets.size(); ++i) {
        float score = masspass[i]*ptsums[i];
        if (score > bestscore) {
            bestscore = score;
            best = i;
        }
    }
    return triplets[best];
}

//processing after selection: calculate mass and check the triplet isn't empty

bool notempty(Vec<std::size_t> index) { //used to check if any triplets passed in an event
    return !index.empty();
}

auto tripletmass(Vec<std::size_t> t, Vec<float> pts, Vec<float> etas, Vec<float> phis, Vec<float> masses) {
    ROOT::Math::PtEtaPhiMVector p1(pts[t[0]],etas[t[0]],phis[t[0]],masses[t[0]]);
    ROOT::Math::PtEtaPhiMVector p2(pts[t[1]],etas[t[1]],phis[t[1]],masses[t[1]]);
    ROOT::Math::PtEtaPhiMVector p3(pts[t[2]],etas[t[2]],phis[t[2]],masses[t[2]]);
    float mass = (p1 + p2 + p3).M();
return mass;
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
	return 1;
	}
    //ROOT::EnableImplicitMT();
std::cout << "Compression? (0/1)" <<std::endl;
int c;
std::cin >> c;
if (c != 1 and c != 0) {
    std::cout << "Invalid option :(" <<std::endl;
    return 1;
}

auto filename = "";
if (n == 1) { //single neutrino
    if (c == 1) filename = "/dev/shm/millerca/l1Nano_SingleNeutrino_PU200_iso.root";
    if (c == 0) filename = "/dev/shm/millerca/l1Nano_SingleNeutrino_PU200_iso_comp0.root";
}
if (n == 2) { //W to 3 pion
    if (c == 1) filename = "/dev/shm/millerca/l1Nano_WTo3Pion_PU200_iso.root";
    if (c == 0) filename = "/dev/shm/millerca/l1Nano_WTo3Pion_PU200_iso_comp0.root";
}
if (n == 3) { //TTbar
    if (c == 1) filename = "/dev/shm/millerca/l1Nano_TTbar_PU200_iso.root";
    if (c == 0) filename = "/dev/shm/millerca/l1Nano_TTbar_PU200_iso_comp0.root";
}

 ROOT::RDataFrame d("Events", filename);
    auto d2 = d.Filter(initptcut, {"L1Puppi_pt","L1Puppi_pdgId"}, "Initial Cut")
        .Define("Triplet_Index", maketriplets, {"L1Puppi_pdgId", "L1Puppi_charge", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"})
        .Filter(notempty,{"Triplet_Index"})
        .Define("Triplet_Mass",tripletmass,{"Triplet_Index", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"});

    auto masshist = d2.Histo1D({"masshist", "W Boson mass from selected pions; mass (GeV/c^2)", 100, 0, 100},"Triplet_Mass");
    auto c1 = new TCanvas("c1","c1",800,600);
    masshist->Draw();
    c1->Draw();
    c1->SaveAs("figures/masshist.pdf");
    //auto cutReport = d2.Report();
    //cutReport->Print();
return 0;
}
