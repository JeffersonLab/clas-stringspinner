#include "src/common.h"
inline void config_clas12(Pythia8::Pythia& pyth) {

  // the beams are back-to-back, but with different energies; this is the recommended setting for fixed target
  set_config(pyth, "Beams:frameType = 2");

  // interaction mechanism: $\gamma*/Z^0$ $t$-channel exchange, with full interference between $\gamma*$ and $Z^0$
  set_config(pyth, "WeakBosonExchange:ff2ff(t:gmZ) = on");

  // phase-space cuts
  set_config(pyth, "PhaseSpace:Q2Min = 1.0");           // minimum Q^2
  set_config(pyth, "PhaseSpace:pTHatMinDiverge = 0.5"); // extra $p_T$ cut to avoid divergences of some processes in the $p_T \to 0$ limit (for low-x)
  set_config(pyth, "PhaseSpace:mHatMin = 0.0");         // minimum invariant mass (for low-x)

  // set dipole recoil; turning this on doesn't appear to change the kinematic distributions noticeably
  set_config(pyth, "SpaceShower:dipoleRecoil = off"); // NOTE: From footnote in the StringSpinner paper: "We also recall that in
                                                      // StringSpinner the parton showers are switched off because presently the
                                                      // string+3 P0 model does not handle the more general string configurations
                                                      // involving multiple partons that would be produced in the showering process."

  // QED radiation off lepton not handled yet by the new procedure;
  // these are recommended when `SpaceShower:dipoleRecoil = on`
  set_config(pyth, "PDF:lepton = off");              // do not use parton densities for lepton beams
  set_config(pyth, "TimeShower:QEDshowerByL = off"); // disallow leptons to radiate photons

  // PDF model
  set_config(pyth, "PDF:pSet = 13"); // NNPDF2.3 QCD+QED LO alpha_s(M_Z) = 0.130 (the current Pythia8 default)

  // switch off resonance decays, ISR, FSR, MPI and Bose-Einstein
  set_config(pyth, "ProcessLevel:resonanceDecays = off");
  set_config(pyth, "PartonLevel:FSRinResonances = off");
  set_config(pyth, "PartonLevel:FSR = off");
  set_config(pyth, "PartonLevel:ISR = off");
  set_config(pyth, "PartonLevel:MPI = off");
  set_config(pyth, "HadronLevel:BoseEinstein = off");

  // invariant mass distribution of resonances as in the string+3P0 model;
  set_config(pyth, "ParticleData:modeBreitWigner = 3");  // particles registered as having a mass width are given a mass in
                                                         // the range m_min < m < m_max, according to a truncated relativistic
                                                         // Breit-Wigner, i.e. quadratic in m.

  // fragmentation parameters
  set_config(pyth, "StringPT:enhancedFraction = 0.0"); // the fraction of string breaks with enhanced width.
  set_config(pyth, "StringPT:enhancedWidth = 1.0");    // the enhancement of the width in this fraction.
  set_config(pyth, "StringZ:aLund = 1.2");             // parameters a and b of (1/z) * (1-z)^a * exp(-b m_T^2 / z)
  set_config(pyth, "StringZ:bLund = 0.58");
  set_config(pyth, "StringFragmentation:stopMass = 0.0"); // used to define a W_min = m_q1 + m_q2 + stopMass, where m_q1 and m_q2 are
                                                          // the masses of the two current endpoint quarks or diquarks; analogous to PARJ(33)

  // ratios of vector mesons to pseudoscalar mesons
  set_config(pyth, "StringFlav:mesonUDvector = 0.7"); // for light (u, d) mesons (analogous to PARJ(11): fraction of $\rho / \pi$)
  set_config(pyth, "StringFlav:mesonSvector = 0.75"); // for strange mesons      (analogous to PARJ(12): fraction of $K^* / K$)

  // transverse momentum
  set_config(pyth, "StringPT:sigma = 0.5");                   // pT width of the fragmentation process (analogous to PARJ(21))
  set_config(pyth, "BeamRemnants:primordialKT = on");         // allow selection of primordial kT according to the parameter values.
  set_config(pyth, "BeamRemnants:primordialKThard = 0.64");   // initial kT width, analgous to PARL(3) (???)
  set_config(pyth, "BeamRemnants:halfScaleForKT = 0.0");      // set these params to zero, to try to make kT width relatively constant
  set_config(pyth, "BeamRemnants:halfMassForKT = 0.0");
  set_config(pyth, "BeamRemnants:primordialKTremnant = 0.0");

  // StringSpinner setings
  set_config(pyth, "StringSpinner:GLGT = 1.4"); // StringSpinner free parameter |GL/GT|
                                                // NOTE: fraction of long. pol. VMs: fL = |GL/GT|^2 / ( 2 + |GL/GT|^2 )
                                                // 0 <= fL <= 1

  set_config(pyth, "StringSpinner:thetaLT = 0"); // StringSpinner free parameter arg(GL/GT)
                                                 // -PI <= thetaLT <= +PI

  // handle event printouts
  set_config(pyth, "Next:numberShowInfo = 0");
  set_config(pyth, "Next:numberShowProcess = 0");
  set_config(pyth, "Next:numberShowEvent = 1");

}
