#include "src/common.h"
static void config_beam_test(Pythia8::Pythia& pyth) {

  // the beams are back-to-back, but with different energies; this is the recommended setting for fixed target
  set_config(pyth, "Beams:frameType = 2");

  // interaction mechanism: $\gamma*/Z^0$ $t$-channel exchange, with full interference between $\gamma*$ and $Z^0$
  set_config(pyth, "WeakBosonExchange:ff2ff(t:gmZ) = on");

  // phase-space cuts
  set_config(pyth, "PhaseSpace:Q2Min = 1.0");   // minimum Q^2
  set_config(pyth, "PhaseSpace:mHatMin = 0.0"); // minimum invariant mass (for low-x)

  // set dipole recoil; turning this on doesn't appear to change the kinematic distributions noticeably
  set_config(pyth, "SpaceShower:dipoleRecoil = off");

  // handle event printouts
  set_config(pyth, "Next:numberShowInfo = 0");
  set_config(pyth, "Next:numberShowProcess = 0");
  set_config(pyth, "Next:numberShowEvent = 0");

}
