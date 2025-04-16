#include <fmt/format.h>
#include "Tools.h"

enum GenCutType {
  c_theta,
  c_z
};

struct GenCut {
  double const min;
  double const max;
  GenCutType const type;
  bool check(double val) {
    return min < val && val < max;
  }
};

using GenCutMap = std::unordered_map<int, std::vector<GenCut>>; // PID -> list of GenCuts that apply

inline void ParseGenCut(char const* optarg, std::string name, GenCutMap& cut_map, GenCutType cut_type) {
  double min, max;
  std::vector<int> pdgs;
  auto token_ftn = [&min, &max, &pdgs](auto token, auto i) {
    switch(i) {
      case 0:
        min = std::stod(token);
        break;
      case 1:
        max = std::stod(token);
        break;
      default:
        pdgs.push_back(std::stoi(token));
    }
  };
  Tokenize(optarg, token_ftn);
  if(pdgs.empty())
    throw std::runtime_error(fmt::format("value of option '--{}' does not have at least 3 arguments", name));
  if(min >= max)
    throw std::runtime_error(fmt::format("option '--{}' has MIN >= MAX", name));
  for(auto const& pdg : pdgs) {
    if(cut_map.find(pdg) == cut_map.end())
      cut_map.insert(pdg, {});
    cut_map.at(pdg).emplace_back(min, max, cut_type);
  }
}
