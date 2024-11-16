#include <Pythia8/Pythia.h>

int main(int argc, char** argv)
{
  Pythia8::Pythia p;
  Pythia8::Event& e = p.event;

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

  p.readString("Next:numberShowInfo = 0");
  p.readString("Next:numberShowProcess = 0");
  p.readString("Next:numberShowEvent = 0");

  p.init();
  for(int i=0; i<100; i++) {
    if(!p.next()) continue;
  }

  return 0;
}
