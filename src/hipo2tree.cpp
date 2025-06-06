#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <string>
#include <TLorentzVector.h>

#include "clas12reader.h"
#include "HipoChain.h"
#include "Kinematics.h"     

const float PI=3.14159265;
const double Mp = 0.938272;
const double Me = 0.000511;

// A minimal struct for reconstructed and MC particles
struct SimplePart {
  int pid;
  double px, py, pz;
  double E, p, theta, phi;
};

int hipo2tree(
    const char* hipoFile,
    const char* outputFile = "dihadron.root",
    const double beamE   = 10.6041
) {
  // Prepare output
  TFile * fOut = new TFile(outputFile, "RECREATE");
  TTree * tree = new TTree("dih", "Dihadron Tree");

  // --- Branch variables ---
  double x, Q2, y, W, nu;
  int    hel;

  double px_e,  py_e,  pz_e,  E_e,  p_e,  theta_e,  phi_e;
  double px_pip,py_pip,pz_pip,E_pip,p_pip,theta_pip,phi_pip;
  double px_pim,py_pim,pz_pim,E_pim,p_pim,theta_pim,phi_pim;
  double Mh, z1, z2, z, dihadron_th, phi_h, phi_R, pT_lab;

  tree->Branch("x",    &x,    "x/D");
  tree->Branch("Q2",   &Q2,   "Q2/D");
  tree->Branch("y",    &y,    "y/D");
  tree->Branch("W",    &W,    "W/D");
  tree->Branch("nu",   &nu,   "nu/D");
  tree->Branch("hel",  &hel,  "hel/I");

  tree->Branch("px_e",  &px_e,  "px_e/D");
  tree->Branch("py_e",  &py_e,  "py_e/D");
  tree->Branch("pz_e",  &pz_e,  "pz_e/D");
  tree->Branch("E_e",   &E_e,   "E_e/D");
  tree->Branch("p_e",   &p_e,   "p_e/D");
  tree->Branch("theta_e", &theta_e, "theta_e/D");
  tree->Branch("phi_e",   &phi_e,   "phi_e/D");

  tree->Branch("px_pip",  &px_pip,  "px_pip/D");
  tree->Branch("py_pip",  &py_pip,  "py_pip/D");
  tree->Branch("pz_pip",  &pz_pip,  "pz_pip/D");
  tree->Branch("E_pip",   &E_pip,   "E_pip/D");
  tree->Branch("p_pip",   &p_pip,   "p_pip/D");
  tree->Branch("theta_pip", &theta_pip, "theta_pip/D");
  tree->Branch("phi_pip",   &phi_pip,   "phi_pip/D");

  tree->Branch("px_pim",  &px_pim,  "px_pim/D");
  tree->Branch("py_pim",  &py_pim,  "py_pim/D");
  tree->Branch("pz_pim",  &pz_pim,  "pz_pim/D");
  tree->Branch("E_pim",   &E_pim,   "E_pim/D");
  tree->Branch("p_pim",   &p_pim,   "p_pim/D");
  tree->Branch("theta_pim", &theta_pim, "theta_pim/D");
  tree->Branch("phi_pim",   &phi_pim,   "phi_pim/D");

  tree->Branch("Mh", &Mh, "Mh/D");
  tree->Branch("z1", &z1, "z1/D");
  tree->Branch("z2", &z2, "z2/D");
  tree->Branch("z", &z, "z/D");
  tree->Branch("phi_h", &phi_h, "phi_h/D");
  tree->Branch("phi_R1", &phi_R, "phi_R1/D");
  tree->Branch("th", &dihadron_th, "th/D");
  tree->Branch("pT_lab", &pT_lab, "pT_lab/D");
  // --- Open HIPO chain ---
  clas12root::HipoChain chain;
  chain.Add(hipoFile);
  auto c12 = chain.GetC12Reader();
  auto &_c12= chain.C12ref();
  // Require at least 1 electron, 1 pi+, 1 pi-
  c12->addAtLeastPid(11,  1);
  c12->addAtLeastPid(211, 1);
  c12->addAtLeastPid(-211,1);

  // Determine helicity from filename pattern
  std::string fname(hipoFile);
  if      (fname.find("_LU_p_")!=std::string::npos) hel = -1; // electron helicity should be opposite quark helicity
  else if (fname.find("_LU_n_")!=std::string::npos) hel = +1; // electron helicity should be opposite quark helicity
  else                                                 hel =  0;

  Kinematics kin;

  // Event loop
  while (chain.Next()) {
    auto event = _c12->event();
    // --- Read REC::Particle bank ---
    std::vector<SimplePart> recParts;
    auto particles=_c12->getDetParticles();
    for(unsigned int idx = 0 ; idx < particles.size() ; idx++){
      auto particle = particles.at(idx);
      int    pid = particle->getPid(); 
      double theta = particle->getTheta();
      double phi   = particle->getPhi();
      double p     = particle->getP();
      double mass = 0.0;
      if(pid!=22)
        mass = particle->getPdgMass();
      else
        mass = 0;
      double px = p * sin(theta) * cos(phi);
      double py = p * sin(theta) * sin(phi);
      double pz = p * cos(theta);
        
      double E    = std::sqrt(px*px + py*py + pz*pz + mass*mass);
      recParts.push_back({pid, px, py, pz, E, p, theta, phi});
    }

    // --- Find highest-energy electron ---
    int ie = -1;
    double maxE = -1;
    for (size_t i = 0; i < recParts.size(); ++i){
      std::cout << recParts[i].pid << std::endl;
      if (recParts[i].pid==11 && recParts[i].E>maxE) {
        maxE = recParts[i].E;
        ie   = i;
      }
    }
    if (ie<0) continue;
    // --- Collect pi+ and pi- indices ---
    std::vector<int> idx_pip, idx_pim;
    for (size_t i = 0; i < recParts.size(); ++i){
      if      (recParts[i].pid==211)  idx_pip.push_back(i);
      else if (recParts[i].pid==-211) idx_pim.push_back(i);
    }
    if (idx_pip.empty() || idx_pim.empty()) continue;

    // --- Loop over all pi+ / pi- pairs ---
    for (int ip : idx_pip) {
      for (int im : idx_pim) {
        auto &e  = recParts[ie];
        auto &pip= recParts[ip];
        auto &pim= recParts[im];

        TLorentzVector lv_target(0, 0, 0, Mp);                    // Target nucleon at rest
        TLorentzVector lv_beam(0, 0, beamE, beamE);               // Incoming electron (massless)
        TLorentzVector lv_scattered(px_e, py_e, pz_e, E_e);       // Scattered electron
        TLorentzVector lv_q = lv_beam - lv_scattered;             // Virtual photon q = k - k'
        // Construct Lorentz vectors for pi+ and pi-
        TLorentzVector lv_pip(px_pip, py_pip, pz_pip, E_pip);
        TLorentzVector lv_pim(px_pim, py_pim, pz_pim, E_pim);
        // Dihadron Lorentz vector and invariant mass
        TLorentzVector lv_dihadron = lv_pip + lv_pim;
          
        
        // Electron
        px_e   = e.px;  py_e   = e.py;  pz_e   = e.pz;
        E_e    = e.E;   p_e    = e.p;
        theta_e= e.theta; phi_e = e.phi;

        // pi+
        px_pip   = pip.px;  py_pip   = pip.py;  pz_pip   = pip.pz;
        E_pip    = pip.E;   p_pip    = pip.p;
        theta_pip= pip.theta; phi_pip = pip.phi;

        // pi-
        px_pim   = pim.px;  py_pim   = pim.py;  pz_pim   = pim.pz;
        E_pim    = pim.E;   p_pim    = pim.p;
        theta_pim= pim.theta; phi_pim = pim.phi;

        // Kinematics
        Q2 = kin.Q2(beamE, E_e, std::cos(theta_e));
        y  = kin.y(beamE, E_e);
        nu = kin.nu(beamE, E_e);
        x  = kin.x(Q2,  lv_q,  lv_target);
        W  = kin.W(Q2,  Mp,  nu);
        z1 = kin.z(lv_target, lv_pip, lv_q);
        z2 = kin.z(lv_target, lv_pim, lv_q);
        z  = kin.z(lv_target, lv_dihadron, lv_q);
        Mh = lv_dihadron.M();  // Invariant mass of the pion pair
        pT_lab = lv_dihadron.Pt();
        phi_h = kin.phi_h(lv_q, lv_scattered, lv_pip, lv_pim);
        phi_R = kin.phi_R(lv_q, lv_scattered, lv_pip, lv_pim);
        dihadron_th = kin.com_th(lv_pip,lv_pim);
        tree->Fill();
      }
    }
  }

  // Write and close
  fOut->Write();
  fOut->Close();
  return 0;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input.hipo> [output.root]" << std::endl;
    return 1;
  }
  const char* input = argv[1];
  const char* output = (argc >= 3) ? argv[2] : "dihadron.root";
  return hipo2tree(input, output);
}
