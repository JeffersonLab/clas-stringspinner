#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <string>
#include <TLorentzVector.h>

#include "clas12reader.h"
#include "HipoChain.h"
#include "Kinematics.h"

const float PI = 3.14159265358979323846;
const double Mp = 0.938272;
const double Me = 0.000511;

// A minimal struct for reconstructed particles
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
  // -----------------------------------
  // 1) Prepare output TTrees
  // -----------------------------------
  TFile* fOut = new TFile(outputFile, "RECREATE");

  // We will create 10 TTrees: original + acceptance for each of 5 channels
  TTree* tree_ppim        = new TTree("dih_piplus_piminus",            "π⁺–π⁻ dihadron");
  TTree* tree_ppim_acc    = new TTree("dih_piplus_piminus_acceptance", "π⁺–π⁻ dihadron (acceptance)");
  TTree* tree_ppiz        = new TTree("dih_piplus_pi0",                 "π⁺–π⁰ dihadron");
  TTree* tree_ppiz_acc    = new TTree("dih_piplus_pi0_acceptance",      "π⁺–π⁰ dihadron (acceptance)");
  TTree* tree_pmiz        = new TTree("dih_piminus_pi0",                "π⁻–π⁰ dihadron");
  TTree* tree_pmiz_acc    = new TTree("dih_piminus_pi0_acceptance",     "π⁻–π⁰ dihadron (acceptance)");
  TTree* tree_pppp        = new TTree("dih_piplus_piplus",              "π⁺–π⁺ dihadron");
  TTree* tree_pppp_acc    = new TTree("dih_piplus_piplus_acceptance",   "π⁺–π⁺ dihadron (acceptance)");
  TTree* tree_pimpimp     = new TTree("dih_piminus_piminus",            "π⁻–π⁻ dihadron");
  TTree* tree_pimpimp_acc = new TTree("dih_piminus_piminus_acceptance", "π⁻–π⁻ dihadron (acceptance)");

  // --- Common branch variables for all trees ---
  double x, Q2, y, W, nu;
  int    hel;

  double px_e,   py_e,   pz_e,   E_e,   p_e,   theta_e,   phi_e;
  double px_pip, py_pip, pz_pip, E_pip, p_pip, theta_pip, phi_pip;
  double px_pim, py_pim, pz_pim, E_pim, p_pim, theta_pim, phi_pim;
  double Mh, z1, z2, z, phi_h, phi_R, dihadron_th, pT_lab, Mx, xF1, xF2;

  // Helper lambda to set up identical branches on any TTree pointer:
  auto setupBranches = [&](TTree* t) {
    t->Branch("x",    &x,    "x/D");
    t->Branch("Q2",   &Q2,   "Q2/D");
    t->Branch("y",    &y,    "y/D");
    t->Branch("W",    &W,    "W/D");
    t->Branch("nu",   &nu,   "nu/D");
    t->Branch("hel",  &hel,  "hel/I");

    t->Branch("px_e",   &px_e,   "px_e/D");
    t->Branch("py_e",   &py_e,   "py_e/D");
    t->Branch("pz_e",   &pz_e,   "pz_e/D");
    t->Branch("E_e",    &E_e,    "E_e/D");
    t->Branch("p_e",    &p_e,    "p_e/D");
    t->Branch("theta_e",&theta_e,"theta_e/D");
    t->Branch("phi_e",  &phi_e,  "phi_e/D");

    t->Branch("px_pip",   &px_pip,   "px_pip/D");
    t->Branch("py_pip",   &py_pip,   "py_pip/D");
    t->Branch("pz_pip",   &pz_pip,   "pz_pip/D");
    t->Branch("E_pip",    &E_pip,    "E_pip/D");
    t->Branch("p_pip",    &p_pip,    "p_pip/D");
    t->Branch("theta_pip",&theta_pip,"theta_pip/D");
    t->Branch("phi_pip",  &phi_pip,  "phi_pip/D");

    t->Branch("px_pim",   &px_pim,   "px_pim/D");
    t->Branch("py_pim",   &py_pim,   "py_pim/D");
    t->Branch("pz_pim",   &pz_pim,   "pz_pim/D");
    t->Branch("E_pim",    &E_pim,    "E_pim/D");
    t->Branch("p_pim",    &p_pim,    "p_pim/D");
    t->Branch("theta_pim",&theta_pim,"theta_pim/D");
    t->Branch("phi_pim",  &phi_pim,  "phi_pim/D");

    t->Branch("Mh",      &Mh,      "Mh/D");
    t->Branch("Mx",      &Mx,      "Mx/D");
    t->Branch("xF1",     &xF1,     "xF1/D");
    t->Branch("xF2",     &xF2,     "xF2/D");
    t->Branch("z1",      &z1,      "z1/D");
    t->Branch("z2",      &z2,      "z2/D");
    t->Branch("z",       &z,       "z/D");
    t->Branch("phi_h",   &phi_h,   "phi_h/D");
    t->Branch("phi_R1",  &phi_R,   "phi_R1/D");
    t->Branch("th",      &dihadron_th, "th/D");
    t->Branch("pT_lab",  &pT_lab,  "pT_lab/D");
  };

  // Apply the same branches to all 10 TTrees
  for (auto tptr : { tree_ppim, tree_ppim_acc,
                     tree_ppiz, tree_ppiz_acc,
                     tree_pmiz, tree_pmiz_acc,
                     tree_pppp, tree_pppp_acc,
                     tree_pimpimp, tree_pimpimp_acc })
  {
    setupBranches(tptr);
  }

  // -----------------------------------
  // 2) Open HIPO chain and set up clas12reader
  // -----------------------------------
  clas12root::HipoChain chain;
  chain.Add(hipoFile);
  auto c12 = chain.GetC12Reader();
  auto& _c12 = chain.C12ref();

  // Require at least 1 electron, 1 pi+, 1 pi- in the event
  c12->addAtLeastPid(11,  1);

  // Determine helicity from filename pattern
  std::string fname(hipoFile);
  if      (fname.find("_LU_p_") != std::string::npos) hel = -1;
  else if (fname.find("_LU_n_") != std::string::npos) hel = +1;
  else                                                  hel =  0;

  // Instantiate kinematics helper
  Kinematics kin;

  // Acceptance angular cuts (in radians)
  const double theta_min = 5. * PI / 180.0;   // 5°
  const double theta_max = 35. * PI / 180.0;  // 35°

  // -----------------------------------
  // 3) Event loop
  // -----------------------------------
  while (chain.Next()) {
    auto event = _c12->event();

    // --- 3.1) Read all reconstructed particles into recParts[] ---
    std::vector<SimplePart> recParts;
    auto particles = _c12->getDetParticles();
    recParts.reserve(particles.size());
    for (unsigned int idx = 0; idx < particles.size(); ++idx) {
      auto particle = particles.at(idx);
      int    pid   = particle->getPid();
      double theta = particle->getTheta();
      double phi   = particle->getPhi();
      double p     = particle->getP();
      double mass  = 0.0;
      if (pid != 22)             // photon has no mass, otherwise use PDG mass
        mass = particle->getPdgMass();
      double px = p * std::sin(theta) * std::cos(phi);
      double py = p * std::sin(theta) * std::sin(phi);
      double pz = p * std::cos(theta);
      double E  = std::sqrt(px*px + py*py + pz*pz + mass*mass);
      recParts.push_back({ pid, px, py, pz, E, p, theta, phi });
    }

    // --- 3.2) Find the highest‐energy final‐state electron (pid == 11) ---
    int ie = -1;
    double maxE = -1.;
    for (size_t i = 0; i < recParts.size(); ++i) {
      if (recParts[i].pid == 11 && recParts[i].E > maxE) {
        maxE = recParts[i].E;
        ie   = i;
      }
    }
    if (ie < 0) continue;  // no final‐state electron → skip event

    // Build TLorentzVectors for beam, target, and scattered electron
    TLorentzVector lv_target(0, 0, 0, Mp);
    TLorentzVector lv_beam(  0, 0, beamE, beamE );
    auto& e = recParts[ie];
    TLorentzVector lv_scattered(e.px, e.py, e.pz, e.E);
    TLorentzVector lv_q = lv_beam - lv_scattered;  // virtual photon

    // --- 3.3) Collect charged‐pion indices ---
    std::vector<int> idx_pip, idx_pim, idx_pho;
    for (size_t i = 0; i < recParts.size(); ++i) {
      if      (recParts[i].pid == 211)  idx_pip.push_back(i);
      else if (recParts[i].pid == -211) idx_pim.push_back(i);
      else if (recParts[i].pid == 22)   idx_pho.push_back(i);  // photon
    }
    if (idx_pip.empty() || idx_pim.empty()) {
      // We still may want π⁺–π⁺ or π⁻–π⁻, but if there is no π⁺ or no π⁻, 
      // skip π⁺–π⁻ and π⁺–π⁰ / π⁻–π⁰ if photon combos are needed.
    }

    // Precompute electron branch variables once:
    px_e    = e.px;   py_e    = e.py;   pz_e    = e.pz;   E_e    = e.E;
    p_e     = e.p;    theta_e = e.theta; phi_e  = e.phi;

    // Precompute lv_target, lv_beam, lv_scattered, lv_q above.

    // ----------------------------------
    // 3.4) π⁺–π⁻ channel
    // ----------------------------------
    if (!idx_pip.empty() && !idx_pim.empty()) {
      for (int ip : idx_pip) {
        for (int im : idx_pim) {
          auto& pip = recParts[ip];
          auto& pim = recParts[im];
          TLorentzVector lv_pip(pip.px, pip.py, pip.pz, pip.E);
          TLorentzVector lv_pim(pim.px, pim.py, pim.pz, pim.E);

          // 3.4.1) Fill original π⁺–π⁻ tree
          {
            // fill electron branches (precomputed) and pion branches:
            px_pip    = pip.px;   py_pip    = pip.py;   pz_pip    = pip.pz;
            E_pip     = pip.E;    p_pip     = pip.p;    theta_pip = pip.theta;
            phi_pip   = pip.phi;

            px_pim    = pim.px;   py_pim    = pim.py;   pz_pim    = pim.pz;
            E_pim     = pim.E;    p_pim     = pim.p;    theta_pim = pim.theta;
            phi_pim   = pim.phi;

            TLorentzVector lv_dihadron = lv_pip + lv_pim;

            // Compute all kinematics (un‐cut):
            Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
            y     = kin.y(  beamE, E_e );
            nu    = kin.nu( beamE, E_e );
            x     = kin.x(    Q2, lv_q,   lv_target );
            W     = kin.W(    Q2, Mp,     nu );
            z1    = kin.z(    lv_target, lv_pip, lv_q );
            z2    = kin.z(    lv_target, lv_pim, lv_q );
            z     = kin.z(    lv_target, lv_dihadron, lv_q );
            Mh    = lv_dihadron.M();
            pT_lab= lv_dihadron.Pt();
            phi_h = kin.phi_h(  lv_q, lv_scattered, lv_pip, lv_pim );
            phi_R = kin.phi_R(  lv_q, lv_scattered, lv_pip, lv_pim );
            dihadron_th = kin.com_th(lv_pip, lv_pim);
            xF1   = kin.xF( lv_q, lv_pip,  lv_target, W );
            xF2   = kin.xF( lv_q, lv_pim,  lv_target, W );
            Mx    = (lv_beam + lv_target - lv_scattered - lv_pip - lv_pim).M();

            tree_ppim->Fill();
          }

          // 3.4.2) Now check acceptance angles & kinematic cuts for π⁺–π⁻
          bool goodAngles =
            (theta_e        >= theta_min && theta_e        <= theta_max) &&
            (pip.theta      >= theta_min && pip.theta      <= theta_max) &&
            (pim.theta      >= theta_min && pim.theta      <= theta_max);
          if (goodAngles) {
            // Recompute kinematics (we already stored px_e, px_pip/pim, etc.)
            TLorentzVector lv_dihadron = lv_pip + lv_pim;

            Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
            y     = kin.y(  beamE, E_e );
            nu    = kin.nu( beamE, E_e );
            x     = kin.x(    Q2, lv_q,   lv_target );
            W     = kin.W(    Q2, Mp,     nu );
            z1    = kin.z(    lv_target, lv_pip, lv_q );
            z2    = kin.z(    lv_target, lv_pim, lv_q );
            z     = kin.z(    lv_target, lv_dihadron, lv_q );
            Mh    = lv_dihadron.M();
            pT_lab= lv_dihadron.Pt();
            phi_h = kin.phi_h(  lv_q, lv_scattered, lv_pip, lv_pim );
            phi_R = kin.phi_R(  lv_q, lv_scattered, lv_pip, lv_pim );
            dihadron_th = kin.com_th(lv_pip, lv_pim);
            xF1   = kin.xF( lv_q, lv_pip,  lv_target, W );
            xF2   = kin.xF( lv_q, lv_pim,  lv_target, W );
            Mx    = (lv_beam + lv_target - lv_scattered - lv_pip - lv_pim).M();

            // π⁺–π⁻ acceptance cuts:
            bool goodKinCuts =
              (Mx     > 1.5 ) &&
              (xF1    > 0   ) &&
              (xF2    > 0   ) &&
              (z      < 0.95) &&
              (W      > 2   ) &&
              (Q2     > 1   ) &&
              (p_pip  > 1.25) &&
              (p_pim  > 1.25);

            if (goodKinCuts) {
              tree_ppim_acc->Fill();
            }
          }
        }
      }
    }

    // ----------------------------------
    // 3.5) π⁺–π⁰ channel
    // ----------------------------------
    // We form every pair of photons (pid == 22), require E_γ ≥ 0.2 GeV, then pair with each π⁺
    if (!idx_pip.empty() && idx_pho.size() >= 2) {
      for (int ip : idx_pip) {
        auto& pip = recParts[ip];
        TLorentzVector lv_pip(pip.px, pip.py, pip.pz, pip.E);

        for (size_t i1 = 0; i1 < idx_pho.size(); ++i1) {
          for (size_t i2 = i1 + 1; i2 < idx_pho.size(); ++i2) {
            auto& ph1 = recParts[idx_pho[i1]];
            auto& ph2 = recParts[idx_pho[i2]];

            // Photon‐energy cut:
            if (ph1.E < 0.2 || ph2.E < 0.2) continue;

            TLorentzVector lv_ph1(ph1.px, ph1.py, ph1.pz, ph1.E);
            TLorentzVector lv_ph2(ph2.px, ph2.py, ph2.pz, ph2.E);
            TLorentzVector lv_pi0 = lv_ph1 + lv_ph2;
            double p_pi0 = lv_pi0.P();
            double theta_pi0 = lv_pi0.Theta();

            // 3.5.1) Fill original π⁺–π⁰ tree:
            {
              px_pip    = pip.px;   py_pip    = pip.py;   pz_pip    = pip.pz;
              E_pip     = pip.E;    p_pip     = pip.p;    theta_pip = pip.theta;
              phi_pip   = pip.phi;

              // For π⁰ we store its 4‐vector as “particle b”:
              px_pim    = lv_pi0.Px();   // re‐use px_pim etc. to store π0
              py_pim    = lv_pi0.Py();
              pz_pim    = lv_pi0.Pz();
              E_pim     = lv_pi0.E();
              p_pim     = p_pi0;
              theta_pim = theta_pi0;
              phi_pim   = lv_pi0.Phi();

              TLorentzVector lv_dihadron = lv_pip + lv_pi0;

              Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
              y     = kin.y(  beamE, E_e );
              nu    = kin.nu( beamE, E_e );
              x     = kin.x(    Q2, lv_q,   lv_target );
              W     = kin.W(    Q2, Mp,     nu );
              z1    = kin.z(    lv_target, lv_pip, lv_q );
              z2    = kin.z(    lv_target, lv_pi0, lv_q );
              z     = kin.z(    lv_target, lv_dihadron, lv_q );
              Mh    = lv_dihadron.M();
              pT_lab= lv_dihadron.Pt();
              phi_h = kin.phi_h(  lv_q, lv_scattered, lv_pip, lv_pi0 );
              phi_R = kin.phi_R(  lv_q, lv_scattered, lv_pip, lv_pi0 );
              dihadron_th = kin.com_th(lv_pip, lv_pi0);
              xF1   = kin.xF( lv_q, lv_pip,  lv_target, W );
              xF2   = kin.xF( lv_q, lv_pi0,  lv_target, W );
              Mx    = (lv_beam + lv_target - lv_scattered - lv_pip - lv_pi0).M();

              tree_ppiz->Fill();
            }

            // 3.5.2) Check acceptance angles & π⁺–π⁰ cuts
            bool goodAngles =
              (theta_e        >= theta_min && theta_e        <= theta_max) &&
              (pip.theta      >= theta_min && pip.theta      <= theta_max) &&
              (ph1.theta      >= theta_min && ph1.theta      <= theta_max) &&
              (ph2.theta      >= theta_min && ph2.theta      <= theta_max);
            if (goodAngles) {
              TLorentzVector lv_dihadron = lv_pip + lv_pi0;

              Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
              y     = kin.y(  beamE, E_e );
              nu    = kin.nu( beamE, E_e );
              x     = kin.x(    Q2, lv_q,   lv_target );
              W     = kin.W(    Q2, Mp,     nu );
              z1    = kin.z(    lv_target, lv_pip, lv_q );
              z2    = kin.z(    lv_target, lv_pi0, lv_q );
              z     = kin.z(    lv_target, lv_dihadron, lv_q );
              Mh    = lv_dihadron.M();
              pT_lab= lv_dihadron.Pt();
              phi_h = kin.phi_h(  lv_q, lv_scattered, lv_pip, lv_pi0 );
              phi_R = kin.phi_R(  lv_q, lv_scattered, lv_pip, lv_pi0 );
              dihadron_th = kin.com_th(lv_pip, lv_pi0);
              xF1   = kin.xF( lv_q, lv_pip,  lv_target, W );
              xF2   = kin.xF( lv_q, lv_pi0,  lv_target, W );
              Mx    = (lv_beam + lv_target - lv_scattered - lv_pip - lv_pi0).M();

              // π⁺–π⁰ acceptance cuts (no p_b required):
              bool goodKinCuts =
                (Mx     > 1.5 ) &&
                (xF1    > 0   ) &&
                (xF2    > 0   ) &&
                (z      < 0.95) &&
                (W      > 2   ) &&
                (Q2     > 1   ) &&
                (p_pip  > 1.25);
              if (goodKinCuts) {
                tree_ppiz_acc->Fill();
              }
            }
          }
        }
      }
    }

    // ----------------------------------
    // 3.6) π⁻–π⁰ channel
    // ----------------------------------
    if (!idx_pim.empty() && idx_pho.size() >= 2) {
      for (int im : idx_pim) {
        auto& pim = recParts[im];
        TLorentzVector lv_pim(pim.px, pim.py, pim.pz, pim.E);

        for (size_t i1 = 0; i1 < idx_pho.size(); ++i1) {
          for (size_t i2 = i1 + 1; i2 < idx_pho.size(); ++i2) {
            auto& ph1 = recParts[idx_pho[i1]];
            auto& ph2 = recParts[idx_pho[i2]];

            if (ph1.E < 0.2 || ph2.E < 0.2) continue;

            TLorentzVector lv_ph1(ph1.px, ph1.py, ph1.pz, ph1.E);
            TLorentzVector lv_ph2(ph2.px, ph2.py, ph2.pz, ph2.E);
            TLorentzVector lv_pi0 = lv_ph1 + lv_ph2;
            double p_pi0 = lv_pi0.P();
            double theta_pi0 = lv_pi0.Theta();

            // 3.6.1) Fill original π⁻–π⁰:
            {
              px_pim    = pim.px;   py_pim    = pim.py;   pz_pim    = pim.pz;
              E_pim     = pim.E;    p_pim     = pim.p;    theta_pim = pim.theta;
              phi_pim   = pim.phi;

              // π⁰ as “particle a”:
              px_pip    = lv_pi0.Px();  // re‐use px_pip etc. to store π0
              py_pip    = lv_pi0.Py();
              pz_pip    = lv_pi0.Pz();
              E_pip     = lv_pi0.E();
              p_pip     = p_pi0;
              theta_pip = theta_pi0;
              phi_pip   = lv_pi0.Phi();

              TLorentzVector lv_dihadron = lv_pim + lv_pi0;

              Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
              y     = kin.y(  beamE, E_e );
              nu    = kin.nu( beamE, E_e );
              x     = kin.x(    Q2, lv_q,   lv_target );
              W     = kin.W(    Q2, Mp,     nu );
              z1    = kin.z(    lv_target, lv_pim, lv_q );
              z2    = kin.z(    lv_target, lv_pi0, lv_q );
              z     = kin.z(    lv_target, lv_dihadron, lv_q );
              Mh    = lv_dihadron.M();
              pT_lab= lv_dihadron.Pt();
              phi_h = kin.phi_h(  lv_q, lv_scattered, lv_pim, lv_pi0 );
              phi_R = kin.phi_R(  lv_q, lv_scattered, lv_pim, lv_pi0 );
              dihadron_th = kin.com_th(lv_pim, lv_pi0);
              xF1   = kin.xF( lv_q, lv_pim,  lv_target, W );
              xF2   = kin.xF( lv_q, lv_pi0,  lv_target, W );
              Mx    = (lv_beam + lv_target - lv_scattered - lv_pim - lv_pi0).M();

              tree_pmiz->Fill();
            }

            // 3.6.2) Check acceptance for π⁻–π⁰
            bool goodAngles =
              (theta_e        >= theta_min && theta_e        <= theta_max) &&
              (pim.theta      >= theta_min && pim.theta      <= theta_max) &&
              (ph1.theta      >= theta_min && ph1.theta      <= theta_max) &&
              (ph2.theta      >= theta_min && ph2.theta      <= theta_max);
            if (goodAngles) {
              TLorentzVector lv_dihadron = lv_pim + lv_pi0;

              Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
              y     = kin.y(  beamE, E_e );
              nu    = kin.nu( beamE, E_e );
              x     = kin.x(    Q2, lv_q,   lv_target );
              W     = kin.W(    Q2, Mp,     nu );
              z1    = kin.z(    lv_target, lv_pim, lv_q );
              z2    = kin.z(    lv_target, lv_pi0, lv_q );
              z     = kin.z(    lv_target, lv_dihadron, lv_q );
              Mh    = lv_dihadron.M();
              pT_lab= lv_dihadron.Pt();
              phi_h = kin.phi_h(  lv_q, lv_scattered, lv_pim, lv_pi0 );
              phi_R = kin.phi_R(  lv_q, lv_scattered, lv_pim, lv_pi0 );
              dihadron_th = kin.com_th(lv_pim, lv_pi0);
              xF1   = kin.xF( lv_q, lv_pim,  lv_target, W );
              xF2   = kin.xF( lv_q, lv_pi0,  lv_target, W );
              Mx    = (lv_beam + lv_target - lv_scattered - lv_pim - lv_pi0).M();

              // π⁻–π⁰ acceptance cuts (no p_b required):
              bool goodKinCuts =
                (Mx     > 1.5 ) &&
                (xF1    > 0   ) &&
                (xF2    > 0   ) &&
                (z      < 0.95) &&
                (W      > 2   ) &&
                (Q2     > 1   ) &&
                (p_pim  > 1.25);
              if (goodKinCuts) {
                tree_pmiz_acc->Fill();
              }
            }
          }
        }
      }
    }

    // ----------------------------------
    // 3.7) π⁺–π⁺ channel
    // ----------------------------------
    if (idx_pip.size() >= 2) {
      for (size_t i1 = 0; i1 < idx_pip.size(); ++i1) {
        for (size_t i2 = i1 + 1; i2 < idx_pip.size(); ++i2) {
          auto& p1 = recParts[idx_pip[i1]];
          auto& p2 = recParts[idx_pip[i2]];
          TLorentzVector lv1(p1.px, p1.py, p1.pz, p1.E);
          TLorentzVector lv2(p2.px, p2.py, p2.pz, p2.E);

          // Determine leading/subleading by z:
          double z1_ = kin.z(lv_target, lv1, lv_q);
          double z2_ = kin.z(lv_target, lv2, lv_q);
          TLorentzVector lo = (z1_ >= z2_ ? lv1 : lv2);
          TLorentzVector hi = (z1_ >= z2_ ? lv2 : lv1);
          SimplePart const& pl = (z1_ >= z2_ ? p1 : p2);
          SimplePart const& ph = (z1_ >= z2_ ? p2 : p1);

          // 3.7.1) Fill original π⁺–π⁺
          {
            // electron is already in px_e,py_e,... 
            px_pip    = pl.px;    py_pip    = pl.py;    pz_pip    = pl.pz;
            E_pip     = pl.E;     p_pip     = pl.p;     theta_pip = pl.theta;
            phi_pip   = pl.phi;

            px_pim    = ph.px;                   // re‐use px_pim etc. to store second π⁺
            py_pim    = ph.py;
            pz_pim    = ph.pz;
            E_pim     = ph.E;
            p_pim     = ph.p;
            theta_pim = ph.theta;
            phi_pim   = ph.phi;

            TLorentzVector lv_dihadron = lo + hi;

            Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
            y     = kin.y(  beamE, E_e );
            nu    = kin.nu( beamE, E_e );
            x     = kin.x(    Q2, lv_q,   lv_target );
            W     = kin.W(    Q2, Mp,     nu );
            z1    = kin.z(    lv_target, lo,    lv_q );
            z2    = kin.z(    lv_target, hi,    lv_q );
            z     = kin.z(    lv_target, lv_dihadron, lv_q );
            Mh    = lv_dihadron.M();
            pT_lab= lv_dihadron.Pt();
            phi_h = kin.phi_h(  lv_q, lv_scattered, lo, hi );
            phi_R = kin.phi_R(  lv_q, lv_scattered, lo, hi );
            dihadron_th = kin.com_th(lo, hi);
            xF1   = kin.xF( lv_q, lo,  lv_target, W );
            xF2   = kin.xF( lv_q, hi,  lv_target, W );
            Mx    = (lv_beam + lv_target - lv_scattered - lo - hi).M();

            tree_pppp->Fill();
          }

          // 3.7.2) Check acceptance for π⁺–π⁺
          bool goodAngles =
            (theta_e      >= theta_min && theta_e      <= theta_max) &&
            (lo.Theta()   >= theta_min && lo.Theta()   <= theta_max) &&
            (hi.Theta()   >= theta_min && hi.Theta()   <= theta_max);
          if (goodAngles) {
            TLorentzVector lv_dihadron = lo + hi;

            Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
            y     = kin.y(  beamE, E_e );
            nu    = kin.nu( beamE, E_e );
            x     = kin.x(    Q2, lv_q,   lv_target );
            W     = kin.W(    Q2, Mp,     nu );
            z1    = kin.z(    lv_target, lo,    lv_q );
            z2    = kin.z(    lv_target, hi,    lv_q );
            z     = kin.z(    lv_target, lv_dihadron, lv_q );
            Mh    = lv_dihadron.M();
            pT_lab= lv_dihadron.Pt();
            phi_h = kin.phi_h(  lv_q, lv_scattered, lo, hi );
            phi_R = kin.phi_R(  lv_q, lv_scattered, lo, hi );
            dihadron_th = kin.com_th(lo, hi);
            xF1   = kin.xF( lv_q, lo,  lv_target, W );
            xF2   = kin.xF( lv_q, hi,  lv_target, W );
            Mx    = (lv_beam + lv_target - lv_scattered - lo - hi).M();

            // π⁺–π⁺ acceptance cuts (both momenta > 1.25)
            bool goodKinCuts =
              (Mx     > 1.5 ) &&
              (xF1    > 0   ) &&
              (xF2    > 0   ) &&
              (z      < 0.95) &&
              (W      > 2   ) &&
              (Q2     > 1   ) &&
              (p_pip  > 1.25) &&
              (p_pim  > 1.25);
            if (goodKinCuts) {
              tree_pppp_acc->Fill();
            }
          }
        }
      }
    }

    // ----------------------------------
    // 3.8) π⁻–π⁻ channel
    // ----------------------------------
    if (idx_pim.size() >= 2) {
      for (size_t i1 = 0; i1 < idx_pim.size(); ++i1) {
        for (size_t i2 = i1 + 1; i2 < idx_pim.size(); ++i2) {
          auto& p1 = recParts[idx_pim[i1]];
          auto& p2 = recParts[idx_pim[i2]];
          TLorentzVector lv1(p1.px, p1.py, p1.pz, p1.E);
          TLorentzVector lv2(p2.px, p2.py, p2.pz, p2.E);

          double z1_ = kin.z(lv_target, lv1, lv_q);
          double z2_ = kin.z(lv_target, lv2, lv_q);
          TLorentzVector lo = (z1_ >= z2_ ? lv1 : lv2);
          TLorentzVector hi = (z1_ >= z2_ ? lv2 : lv1);
          SimplePart const& pl = (z1_ >= z2_ ? p1 : p2);
          SimplePart const& ph = (z1_ >= z2_ ? p2 : p1);

          // 3.8.1) Fill original π⁻–π⁻
          {
            px_pim    = pl.px;    py_pim    = pl.py;    pz_pim    = pl.pz;
            E_pim     = pl.E;     p_pim     = pl.p;     theta_pim= pl.theta;
            phi_pim   = pl.phi;

            px_pip    = ph.px;                   // re‐use px_pip etc. to store second π⁻
            py_pip    = ph.py;
            pz_pip    = ph.pz;
            E_pip     = ph.E;
            p_pip     = ph.p;
            theta_pip = ph.theta;
            phi_pip   = ph.phi;

            TLorentzVector lv_dihadron = lo + hi;

            Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
            y     = kin.y(  beamE, E_e );
            nu    = kin.nu( beamE, E_e );
            x     = kin.x(    Q2, lv_q,   lv_target );
            W     = kin.W(    Q2, Mp,     nu );
            z1    = kin.z(    lv_target, lo,    lv_q );
            z2    = kin.z(    lv_target, hi,    lv_q );
            z     = kin.z(    lv_target, lv_dihadron, lv_q );
            Mh    = lv_dihadron.M();
            pT_lab= lv_dihadron.Pt();
            phi_h = kin.phi_h(  lv_q, lv_scattered, lo, hi );
            phi_R = kin.phi_R(  lv_q, lv_scattered, lo, hi );
            dihadron_th = kin.com_th(lo, hi);
            xF1   = kin.xF( lv_q, lo,  lv_target, W );
            xF2   = kin.xF( lv_q, hi,  lv_target, W );
            Mx    = (lv_beam + lv_target - lv_scattered - lo - hi).M();

            tree_pimpimp->Fill();
          }

          // 3.8.2) Check acceptance for π⁻–π⁻
          bool goodAngles =
            (theta_e      >= theta_min && theta_e      <= theta_max) &&
            (lo.Theta()   >= theta_min && lo.Theta()   <= theta_max) &&
            (hi.Theta()   >= theta_min && hi.Theta()   <= theta_max);
          if (goodAngles) {
            TLorentzVector lv_dihadron = lo + hi;

            Q2    = kin.Q2( beamE, E_e, std::cos(theta_e) );
            y     = kin.y(  beamE, E_e );
            nu    = kin.nu( beamE, E_e );
            x     = kin.x(    Q2, lv_q,   lv_target );
            W     = kin.W(    Q2, Mp,     nu );
            z1    = kin.z(    lv_target, lo,    lv_q );
            z2    = kin.z(    lv_target, hi,    lv_q );
            z     = kin.z(    lv_target, lv_dihadron, lv_q );
            Mh    = lv_dihadron.M();
            pT_lab= lv_dihadron.Pt();
            phi_h = kin.phi_h(  lv_q, lv_scattered, lo, hi );
            phi_R = kin.phi_R(  lv_q, lv_scattered, lo, hi );
            dihadron_th = kin.com_th(lo, hi);
            xF1   = kin.xF( lv_q, lo,  lv_target, W );
            xF2   = kin.xF( lv_q, hi,  lv_target, W );
            Mx    = (lv_beam + lv_target - lv_scattered - lo - hi).M();

            // π⁻–π⁻ acceptance cuts
            bool goodKinCuts =
              (Mx     > 1.5 ) &&
              (xF1    > 0   ) &&
              (xF2    > 0   ) &&
              (z      < 0.95) &&
              (W      > 2   ) &&
              (Q2     > 1   ) &&
              (p_pip  > 1.25) &&
              (p_pim  > 1.25);
            if (goodKinCuts) {
              tree_pimpimp_acc->Fill();
            }
          }
        }
      }
    }

  } // end of event loop

  // -----------------------------------
  // 4) Write out all TTrees and close
  // -----------------------------------
  fOut->Write();
  fOut->Close();
  return 0;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input.hipo> [output.root]" << std::endl;
    return 1;
  }
  const char* input  = argv[1];
  const char* output = (argc >= 3) ? argv[2] : "dihadron.root";
  return hipo2tree(input, output);
}
