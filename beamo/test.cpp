#include <Pythia8/Pythia.h>
#include <cassert>

int const I_BEAM         = 1; // beam particle index
int const I_TARGET       = 2; // target particle index
int const MIN_EVENT_SIZE = 3;
int const PRECISION      = 10; // precision for printout

//////////////////////////////////////////////////////////////////////////////////

void print_particle(Pythia8::Particle const& par, std::string const& name) {
  std::cout << std::setprecision(PRECISION) << name << "pz = " << par.pz() << " E = " << par.e() << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  Pythia8::Pythia p;
  auto& evt  = p.event;
  auto& proc = p.process;

  // pythia parameters ================================================================

  p.readString("Beams:frameType = 2");
  p.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  p.readString("PhaseSpace:Q2Min = 0.0");
  p.readString("PhaseSpace:mHatMin = 0.0");
  p.readString("SpaceShower:dipoleRecoil = off");

  p.readString("Beams:idA = 11");
  p.readString("Beams:idB = 2212");
  p.readString("Beams:eA = 10.6");
  p.readString("Beams:eB = 0.0");
  p.readString("Random:setSeed = on");
  p.readString("Random:seed = 82");

  // event loop =======================================================================
  p.init();
  for(int i=0; i<100; i++) {
    if(!p.next()) continue;
    if(evt.size() < MIN_EVENT_SIZE) continue;

    //// boost
    Pythia8::RotBstMatrix toLab;
    Pythia8::Vec4 pLepLab = p.process[1].p();
    Pythia8::Vec4 pLepEvent = p.event[1].p();
    toLab.bst(pLepEvent, pLepLab);
    p.event.rotbst(toLab);

    //// full event printouts
    // proc.list(false, false, PRECISION);
    // evt.list(false, false, PRECISION);
    // info.list();

    //// get the beam and target particles
    auto const& evt_beam    = evt[I_BEAM];
    auto const& evt_target  = evt[I_TARGET];
    auto const& proc_beam   = proc[I_BEAM];
    auto const& proc_target = proc[I_TARGET];
    for(auto const& [beam, target] : std::vector<std::pair<decltype(evt_beam),decltype(evt_beam)>>{{evt_beam, evt_target}, {proc_beam, proc_target}}) {
      assert((beam.id() == 11 && target.id() == 2212));
      assert((beam.mother1() == 0 && beam.mother2() == 0 && target.mother1() == 0 && target.mother2() == 0));
    }

    //// print the beam and target particles
    print_particle(proc_beam,   "hard process beam:   ");
    print_particle(evt_beam,    "event record beam:   ");
    print_particle(proc_target, "hard process target: ");
    print_particle(evt_target,  "event record target: ");
    std::cout << std::endl;
  }

  return 0;
}
