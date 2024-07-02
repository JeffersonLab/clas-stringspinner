#include "common.h"

void set_config(Pythia8::Pythia& pyth, std::string config) {
  auto result = pyth.readString(config);
  if(!result)
    throw std::runtime_error("bad configuration: '" + config + "'");
}
