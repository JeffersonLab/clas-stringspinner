#include "src/common.h"
static void config_clas12(Pythia8::Pythia& pyth) {

  // Set up incoming beams, for frame with unequal beam energies.
  // the beams are back-to-back, but with different energies
  set_config(pyth, "Beams:frameType = 2");

  // Interaction mechanism.
  // $\gamma*/Z^0$ $t$-channel exchange, with full interference between $\gamma*$ and $Z^0$
  set_config(pyth, "WeakBosonExchange:ff2ff(t:gmZ) = on");

  // Phase-space cut: minimal Q2 of process.
  set_config(pyth, "PhaseSpace:Q2Min = 1.0");

  // Go down to low x-Bjorken.
  set_config(pyth, "PhaseSpace:pTHatMinDiverge = 0.5"); // extra $p_T$ cut to avoid divergences of some processes in the $p_T \to 0$ limit
  set_config(pyth, "PhaseSpace:mHatMin = 0.");          // minimum invariant mass

  // Set dipole recoil on. Necessary for DIS + shower.
  set_config(pyth, "SpaceShower:dipoleRecoil = off");

  // QED radiation off lepton not handled yet by the new procedure.
  // these are recommended when `SpaceShower:dipoleRecoil = on`
  set_config(pyth, "PDF:lepton = off");              // do not use parton densities for lepton beams
  set_config(pyth, "TimeShower:QEDshowerByL = off"); // disallow leptons to radiate photons

  // Choice of PDF = CTEQ5L LO (pSet=2).
  set_config(pyth, "PDF:pSet = 2");
  set_config(pyth, "PDF:pHardSet = 2");

  // Switch off resonance decays, ISR, FSR, MPI and Bose-Einstein.
  set_config(pyth, "ProcessLevel:resonanceDecays = off");
  set_config(pyth, "PartonLevel:FSRinResonances = off");
  set_config(pyth, "PartonLevel:FSR = off");
  set_config(pyth, "PartonLevel:ISR = off");
  set_config(pyth, "PartonLevel:MPI = off");
  set_config(pyth, "HadronLevel:BoseEinstein = off");

  // Switches for hadron production and decay.
  // set_config(pyth, "HadronLevel:Decay = off");
  set_config(pyth, "111:onMode = off"); // pi0
  set_config(pyth, "311:onMode = off"); // K0
  set_config(pyth, "211:onMode = off"); // pi+
  // set_config(pyth, "221:onMode = off"); // eta
  // set_config(pyth, "331:onMode = off"); // eta'
  // set_config(pyth, "311:onMode = off"); // K0 decay

  // Invariant mass distribution of resonances as in the string+3P0 model.
  // particles registered as having a mass width are given a mass in the range m_min < m < m_max,
  // according to a truncated relativistic Breit-Wigner, i.e. quadratic in m.
  set_config(pyth, "ParticleData:modeBreitWigner = 3");

  // Switch off automatic event listing in favour of manual.
  set_config(pyth, "Next:numberShowInfo = 0");
  set_config(pyth, "Next:numberShowProcess = 0");
  set_config(pyth, "Next:numberShowEvent = 1");

  // Settings of string fragmentation parameters.
  set_config(pyth, "StringPT:enhancedFraction = 0.0"); // the fraction of string breaks with enhanced width.
  set_config(pyth, "StringPT:enhancedWidth = 0.0");    // the enhancement of the width in this fraction.

  // Settings from `clasdis`
  //// pythia 6 -> 8 translation from: https://skands.web.cern.ch/slides/11/11-02-skands-uemb.pdf
  //// ratios of vector mesons to pseudoscalar mesons
  set_config(pyth, "StringFlav:mesonUDvector = 0.7");          // for light (u, d) mesons (analogous to PARJ(11): fraction of $\rho / \pi$)
  set_config(pyth, "StringFlav:mesonSvector = 0.75");          // for strange mesons      (analogous to PARJ(12): fraction of $K^* / K$)
  //// momentum widths
  set_config(pyth, "StringPT:sigma = 0.5");                    // pT width of the fragmentation process (analogous to PARJ(21))
  set_config(pyth, "BeamRemnants:primordialKT = off");         // Allow or not selection of primordial kT according to the parameter values.
                                                               // TODO: why do we get NO events if this is `on`?
  set_config(pyth, "BeamRemnants:primordialKTremnant = 0.64"); // The width sigma_remn, assigned as a primordial kT to beam-remnant partons. (analogous to PARJ(99))
                                                               // NOTE: this is ignored when `primordialKT == off`
}
