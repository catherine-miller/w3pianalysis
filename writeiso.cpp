//writes isolation into trees and saves

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
//#include <RSnapshotOptions.hxx>

template <typename T> using Vec = ROOT::RVec<T>; 
using ROOT::Math::XYZTVector;
using ROOT::Math::PtEtaPhiMVector;

//parameters for dr cone
double mindr = 0.01;
double maxdr = 0.25;

auto isolation(Vec<float> eta, Vec<float> phi, Vec<float> pt) {
    Vec<float> ret;
    for (std::size_t i = 0; i < eta.size(); ++i) { //loop over particles whose isolation to calculate
        double psum = 0;
        for (std::size_t j = 0; j < eta.size(); ++j) { //loop over other particles
            if (i == j) continue;
            double dr = sqrt(pow(eta[i]-eta[j],2)+pow(phi[i]-phi[j],2));
            if (dr >= mindr && dr <= maxdr) psum += pt[j];
        }
        float iso = psum/pt[i];
        ret.push_back(iso);
    }
    return ret;
}

int main() {
//this time set compression to 0
std::unique_ptr<TFile> f1_( TFile::Open("/eos/home-m/millerca/SWAN_projects/W_3pi/data/l1Nano_WTo3Pion_PU200_iso_comp0.root", "RECREATE") );
std::unique_ptr<TFile> f2_( TFile::Open("/eos/home-m/millerca/SWAN_projects/W_3pi/data/l1Nano_SingleNeutrino_PU200_iso_comp0.root", "RECREATE") );
std::unique_ptr<TFile> f3_( TFile::Open("/eos/home-m/millerca/SWAN_projects/W_3pi/data/l1Nano_TTbar_PU200_iso_comp0.root", "RECREATE") );

ROOT::EnableImplicitMT();
//W->3pi 
ROOT::RDF::RSnapshotOptions opts;
opts.fCompressionLevel = 0;
ROOT::RDataFrame dw("Events", "/eos/cms/store/cmst3/group/l1tr/gpetrucc/w3pi/ntuples/v1/l1Nano_WTo3Pion_PU200.v1_more.root");
auto dw20 = dw.Define("RelIso", isolation, {"L1Puppi_eta", "L1Puppi_phi", "L1Puppi_pt"});
dw20.Snapshot("Events","/eos/home-m/millerca/SWAN_projects/W_3pi/data/l1Nano_WTo3Pion_PU200_iso_comp0.root", 
              {"run","luminosityBlock","event","nGenPi","GenPi_eta","GenPi_mass","GenPi_phi","GenPi_pt",
               "GenPi_vz","GenPi_charge","GenPi_pdgId","GenPi_prompt","nGenW","GenW_eta","GenW_mass","GenW_phi",
               "GenW_pt","GenW_vz","GenW_charge","GenW_pdgId","nL1PF","L1PF_eta","L1PF_mass","L1PF_phi","L1PF_pt",
               "L1PF_vz","L1PF_charge","L1PF_pdgId",
                  "RelIso","nL1Puppi","L1Puppi_eta","L1Puppi_mass","L1Puppi_phi","L1Puppi_pt","L1Puppi_vz","L1Puppi_charge",
              "L1Puppi_pdgId","L1PF_GenPiIdx","L1PF_GenPiFlav","L1Puppi_GenPiIdx","L1Puppi_GenPiFlav"},opts);
//single neutrino
ROOT::RDataFrame dsn_("Events", "/eos/cms/store/cmst3/group/l1tr/gpetrucc/w3pi/ntuples/v1/l1Nano_SingleNeutrino_PU200.v1_more.root");
auto dsn2_ = dsn_.Define("RelIso", isolation, {"L1Puppi_eta", "L1Puppi_phi", "L1Puppi_pt"});
auto colnames2 = dsn2_.GetColumnNames();
dsn2_.Snapshot("Events","/eos/home-m/millerca/SWAN_projects/W_3pi/data/l1Nano_SingleNeutrino_PU200_iso_comp0.root", colnames2, opts);

//ttbar
ROOT::RDataFrame dt_("Events", "/eos/cms/store/cmst3/group/l1tr/gpetrucc/w3pi/ntuples/v1/l1Nano_TTbar_PU200.v1.root");
auto dt2_ = dt_.Define("RelIso", isolation, {"L1Puppi_eta", "L1Puppi_phi", "L1Puppi_pt"});
auto colnames3 = dt2_.GetColumnNames();
dt2_.Snapshot("Events","/eos/home-m/millerca/SWAN_projects/W_3pi/data/l1Nano_TTbar_PU200_iso_comp0.root",colnames3, opts);

f1_->Write();
f1_->Close();
f2_->Write();
f2_->Close();
f3_->Write();
f3_->Close();

return 0;
}