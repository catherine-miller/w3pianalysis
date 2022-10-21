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

auto maketriplets(Vec<int> pdgids, Vec<int> charges, Vec<float> pts, Vec<float> etas, Vec<float> phis, Vec<float> masses) {
    Vec<std::size_t> ret;
    if (pdgids.size() >= 3) {
        for (std::size_t i = 0; i < 3; ++i) {
            ret.push_back(i);
        }
    }
    return ret;
}

bool notempty(Vec<std::size_t> index) { //used to check if any triplets passed in an event
    if (index.size() == 3) return false; //normally we would check if there is a triplet
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
    //W->3pi    
if (n == 2){
 ROOT::RDataFrame d("Events", "/eos/cms/store/cmst3/group/l1tr/gpetrucc/w3pi/ntuples/v1/l1Nano_WTo3Pion_PU200.v1_more.root");
    auto d2 = d
        .Define("Triplet_Index", maketriplets, {"L1Puppi_pdgId", "L1Puppi_charge", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"})
        .Filter(notempty,{"Triplet_Index"})
        .Define("Triplet_Mass",tripletmass,{"Triplet_Index", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"});

    auto masshist = d2.Histo1D({"masshist", "W Boson mass from selected pions; mass (GeV/c^2)", 100, 0, 100},"Triplet_Mass");
    auto c1 = new TCanvas("c1","c1",800,600);
    masshist->Draw();
    c1->Draw();
    c1->SaveAs("figures/masshist.pdf");
}
    //ttbar  
    if (n == 3) {
	ROOT::RDataFrame d("Events", "/eos/cms/store/cmst3/group/l1tr/gpetrucc/w3pi/ntuples/v1/l1Nano_TTbar_PU200.v1.root");
    auto d2 = d
        .Define("Triplet_Index", maketriplets, {"L1Puppi_pdgId", "L1Puppi_charge", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"})
        .Filter(notempty,{"Triplet_Index"})
        .Define("Triplet_Mass",tripletmass,{"Triplet_Index", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"});

    auto masshist = d2.Histo1D({"masshist", "W Boson mass from selected pions; mass (GeV/c^2)", 100, 0, 100},"Triplet_Mass");
    auto c1 = new TCanvas("c1","c1",800,600);
    masshist->Draw();
    c1->Draw();
    c1->SaveAs("figures/masshist.pdf");
}    
//single neutrino
    if (n == 1) {
	 ROOT::RDataFrame d("Events", "/eos/cms/store/cmst3/group/l1tr/gpetrucc/w3pi/ntuples/v1/l1Nano_SingleNeutrino_PU200.v1_more.root");
    auto d2 = d
        .Define("Triplet_Index", maketriplets, {"L1Puppi_pdgId", "L1Puppi_charge", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"})
        .Filter(notempty,{"Triplet_Index"})
        .Define("Triplet_Mass",tripletmass,{"Triplet_Index", "L1Puppi_pt", "L1Puppi_eta", "L1Puppi_phi", "L1Puppi_mass"});
    
    auto masshist = d2.Histo1D({"masshist", "W Boson mass from selected pions; mass (GeV/c^2)", 100, 0, 100},"Triplet_Mass");
    auto c1 = new TCanvas("c1","c1",800,600);
    masshist->Draw();
    c1->Draw();
    c1->SaveAs("figures/masshist.pdf");
}
return 0;
}
