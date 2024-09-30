#include <getopt.h>
#include <array>
#include <functional>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fmt/os.h>
#include <stringspinner/StringSpinner.h>

// configurations
#include "config/clas12.h"
#include "config/beam_test.h"
static std::map<std::string, std::function<void(Pythia8::Pythia&)>> CONFIG_MAP = {
  {"clas12", config_clas12},
  {"beam_test", config_beam_test},
};

// constants
const int EXIT_ERROR  = 1;
const int EXIT_SYNTAX = 2;
const int SEED_MAX    = 900000000;
const int BEAM_PDG    = 11;

enum obj_enum { objBeam, objTarget, nObj };
const std::string obj_name[nObj] = { "beam", "target" };

// default option values
static unsigned long            num_events       = 10000;
static std::string              out_file         = "clas-stringspinner.dat";
static double                   beam_energy      = 10.60410;
static std::string              target_type      = "proton";
static std::string              pol_type         = "UU";
static std::string              spin_type[nObj]  = {"", ""};
static std::vector<int>         cut_inclusive    = {};
static std::vector<double>      cut_theta        = {};
static std::string              config_name      = "clas12";
static std::vector<std::string> config_overrides = {};
static int                      seed             = -1;
// default flag values
static int  flag_count_before_cuts   = 0;
static int  flag_verbose_mode        = 0;
static bool enable_count_before_cuts = false;
static bool enable_verbose_mode      = false;

//////////////////////////////////////////////////////////////////////////////////

struct LundHeader {
  int    num_particles;
  double target_mass;
  int    target_atomic_num;
  double target_spin;
  double beam_spin;
  int    beam_type;
  double beam_energy;
  int    nucleon_pdg;
  int    process_id;
  double event_weight;
};

struct LundParticle {
  int    index;
  double lifetime;
  int    status;
  int    pdg;
  int    mother1;
  int    daughter1;
  double px;
  double py;
  double pz;
  double energy;
  double mass;
  double vx;
  double vy;
  double vz;
};

//////////////////////////////////////////////////////////////////////////////////

void Usage()
{
  std::vector<std::string> config_name_list;
  for(auto const& config : CONFIG_MAP)
    config_name_list.push_back(config.first);
  fmt::print(R"(USAGE: clas-stringspinner [OPTIONS]...

  --verbose                        verbose printout
  --help                           print this usage guide
  --version                        print the version number


OUTPUT FILE CONTROL:

  --num-events NUM_EVENTS          number of events
                                   warning: do not use this option on the OSG portal;
                                            instead, use the portal's options
                                   default: {num_events}

  --count-before-cuts              if used, --num-events will be the number of events
                                   before any cuts; the Lund file will thus have less
                                   than --num-events events
                                   default: number of output events == --num-events

  --out-file OUTPUT_FILE           output file name
                                   default: {out_file:?}


BEAM AND TARGET PROPERTIES:

  --beam-energy ENERGY             electron beam energy [GeV]
                                   default: {beam_energy}

  --target-type TARGET_TYPE        target type, one of:
                                     proton
                                     neutron
                                   default: {target_type:?}

  --pol-type POLARIZATION_TYPE     beam and target polarization types
                                   - two characters: beam and target
                                   - types: 'U' = unpolarized
                                            'L' = longitudinally polarized
                                            'T' = transversely polarized
                                   default: {pol_type:?}

  --beam-spin BEAM_SPIN            the spin of the beam leptons
                                   - if longitudinally polarized ('L'):
                                     'p' = spin along +z axis
                                     'n' = spin along -z axis
                                   - if transversely polarized ('T'):
                                     'p' = spin along +y axis
                                     'n' = spin along -y axis
                                   - if unpolarized ('U'): no effect

  --target-spin TARGET_SPIN        the spin of the target nucleons
                                   - same usage as --beam-spin, applied to target


GENERATOR PARAMETERS:

  Configuration parameters may be loaded from a configuration file (.h) from:
    https://github.com/JeffersonLab/clas-stringspinner/tree/main/src/config
  Use the --config option to choose one of them, and use the options below
  to set additional specific parameters

  --config CONFIG_NAME             Pythia configuration
                                   - choose a configuration file from one of the following:
                                            {config_name_list}
                                   default: {config_name:?}

  --set PARAM=VAL                  set any Pythia parameter PARAM to the value VALUE
                                   - this option is repeatable
                                   - this will OVERRIDE anything else that sets PARAM
                                   - surround your argument in single quotes (')

  --seed SEED                      random number generator seed, where:
                                   - Pythia's default seed: -1
                                   - seed based on time:  0
                                   - fixed seed:  1 to {seed_max}
                                   warning: do not use this option on the OSG portal,
                                            since it will be set for you automatically
                                   default: {seed}


CUTS FOR EVENT SELECTION:

  --cut-inclusive PDG_CODES...     only allow events which have a least these particles
                                   - delimit by commas
                                   - repeat PDG codes to require more than one
                                   - example: 1 pi- and 2 pi+s: --cut-inclusive -211,211,211
                                   default: {cut_inclusive}

  --cut-theta MIN,MAX              if set, along with --cut-inclusive, this requires the theta
                                   of all particles used in --cut-inclusive to have
                                   MIN <= theta <= MAX, with units in degrees


OPTIONS FOR OSG COMPATIBILITY:

  --trig NUM_EVENTS                same as --num-events
  --docker                         unused
    )" + std::string("\n"),
      fmt::arg("num_events", num_events),
      fmt::arg("out_file", out_file),
      fmt::arg("beam_energy", beam_energy),
      fmt::arg("target_type", target_type),
      fmt::arg("pol_type", pol_type),
      fmt::arg("cut_inclusive", cut_inclusive.empty() ? std::string("no cut") : fmt::format("{}", fmt::join(cut_inclusive, ","))),
      fmt::arg("config_name_list", fmt::join(config_name_list, "\n                                            ")),
      fmt::arg("config_name", config_name),
      fmt::arg("seed", seed),
      fmt::arg("seed_max", SEED_MAX)
      );
}

void Verbose(std::string msg)
{
  if(enable_verbose_mode)
    fmt::print(msg + "\n");
}

int Error(std::string msg)
{
  fmt::print(stderr, "[ERROR] " + msg + "\n");
  return EXIT_ERROR;
}

void Tokenize(char const* str, std::function<void(std::string,int)> func)
{
  std::istringstream token_stream(str);
  std::string token;
  char const delim = ',';
  int i = 0;
  while(getline(token_stream, token, delim))
    func(token, i);
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{

  // parse arguments
  struct option const opts[] = {
    {"num-events",      required_argument, nullptr, 'n'},
    {"trig",            required_argument, nullptr, 'n'},
    {"docker",          no_argument,       nullptr, 'D'},
    {"out-file",        required_argument, nullptr, 'o'},
    {"beam-energy",     required_argument, nullptr, 'e'},
    {"target-type",     required_argument, nullptr, 'T'},
    {"pol-type",        required_argument, nullptr, 'p'},
    {"beam-spin",       required_argument, nullptr, 'b'},
    {"target-spin",     required_argument, nullptr, 't'},
    {"cut-inclusive",   required_argument, nullptr, 'I'},
    {"cut-theta",       required_argument, nullptr, 'A'},
    {"config",          required_argument, nullptr, 'c'},
    {"seed",            required_argument, nullptr, 's'},
    {"set",             required_argument, nullptr, 'S'},
    {"help",            no_argument,       nullptr, 'h'},
    {"version",         no_argument,       nullptr, 'V'},
    {"count-before-cuts", no_argument, &flag_count_before_cuts, 1},
    {"verbose",           no_argument, &flag_verbose_mode,      1},
    {nullptr, 0, nullptr, 0}
  };

  if(argc <= 1) {
    Usage();
    return EXIT_SYNTAX;
  };

  char opt;
  while((opt = getopt_long(argc, argv, "", opts, nullptr)) != -1) {
    switch(opt) {
      case 'n': num_events = std::stol(optarg); break;
      case 'D': break;
      case 'o': out_file = std::string(optarg); break;
      case 'e': beam_energy = std::stod(optarg); break;
      case 'T': target_type = std::string(optarg); break;
      case 'p': pol_type = std::string(optarg); break;
      case 'b': spin_type[objBeam] = std::string(optarg); break;
      case 't': spin_type[objTarget] = std::string(optarg); break;
      case 'I':
        cut_inclusive.clear();
        Tokenize(optarg, [&](auto token, auto i) { cut_inclusive.push_back(std::stoi(token)); });
        break;
      case 'A': {
        cut_theta.clear();
        Tokenize(optarg, [&](auto token, auto i) { cut_theta.push_back(std::stod(token)); });
        if(cut_theta.size() != 2)
          return Error("value of option '--cut-theta' does not have 2 arguments");
        if(cut_theta[1] <= cut_theta[0])
          return Error("value of option '--cut-theta' has MAX <= MIN");
        break;
      }
      case 'c': config_name = std::string(optarg); break;
      case 's': seed = std::stoi(optarg); break;
      case 'S': config_overrides.push_back(std::string(optarg)); break;
      case 'h':
        Usage();
        return 0;
      case 'V':
        fmt::print("{}\n", CLAS_STRINGSPINNER_VERSION);
        return 0;
      case '?':
        return EXIT_ERROR;
    }
  }

  // set boolean options
  bool enable_cut_inclusive = ! cut_inclusive.empty();
  bool enable_cut_theta     = ! cut_theta.empty();
  enable_count_before_cuts  = flag_count_before_cuts == 1;
  enable_verbose_mode       = flag_verbose_mode      == 1;

  // initialize "checklist" `cut_inclusive_found` for checking if `cut_inclusive` satisfied for an event
  std::vector<std::pair<int, bool>> cut_inclusive_found;
  for(auto pdg : cut_inclusive)
    cut_inclusive_found.push_back({pdg, false});

  // check if seed is too large; if so, % SEED_MAX
  if(seed > SEED_MAX) {
    auto new_seed = seed % SEED_MAX;
    Error(fmt::format("value of option '--seed' is too large for Pythia8: {} > {}; setting it to `seed % {}` = {}", seed, SEED_MAX, SEED_MAX, new_seed));
    seed = new_seed;
  }

  // print options
  Verbose(fmt::format("{:=^82}", " Arguments "));
  Verbose(fmt::format("{:>30} = {}", "num-events", num_events));
  Verbose(fmt::format("{:>30} = {}", "count-before-cuts", enable_count_before_cuts ? "true" : "false"));
  Verbose(fmt::format("{:>30} = {:?}", "out-file", out_file));
  Verbose(fmt::format("{:>30} = {} GeV", "beam-energy", beam_energy));
  Verbose(fmt::format("{:>30} = {:?}", "target-type", target_type));
  Verbose(fmt::format("{:>30} = {:?}", "pol-type", pol_type));
  Verbose(fmt::format("{:>30} = {:?}", "beam-spin", spin_type[objBeam]));
  Verbose(fmt::format("{:>30} = {:?}", "target-spin", spin_type[objTarget]));
  Verbose(fmt::format("{:>30} = ({}) [{}]", "cut-inclusive", fmt::join(cut_inclusive, ", "), enable_cut_inclusive ? "enabled" : "disabled"));
  Verbose(fmt::format("{:>30} = ({}) [{}]", "cut-theta", fmt::join(cut_theta, ", "), enable_cut_theta ? "enabled" : "disabled"));
  Verbose(fmt::format("{:>30} = {}", "seed", seed));
  Verbose(fmt::format("{:>30} = {}", "config", config_name));
  Verbose(fmt::format("{:-^82}", ""));
  Verbose(fmt::format("{}{}", "parameter overrides (from option '--set'):", config_overrides.empty() ? " none" : ""));
  for(auto const& config_str : config_overrides)
    Verbose(fmt::format("- {:?}", config_str));
  Verbose(fmt::format("{:=^82}", ""));

  // initialize pythia
  Pythia8::Pythia pyth;
  Pythia8::Event& evt = pyth.event;
  Pythia8::ParticleData& pdt = pyth.particleData;

  // get the configuration function
  std::function<void(Pythia8::Pythia&)> apply_config_func;
  try {
    apply_config_func = CONFIG_MAP.at(config_name);
  }
  catch(std::out_of_range const& ex) {
    Error(fmt::format("value of option '--config' is {:?}, which is not found", config_name));
    return EXIT_ERROR;
  }

  // set target PDG and mass
  int target_pdg;
  int target_atomic_num;
  if(target_type == "proton") {
    target_pdg        = 2212;
    target_atomic_num = 1;
  }
  else if(target_type == "neutron") {
    target_pdg        = 2112;
    target_atomic_num = 0;
  }
  else return Error(fmt::format("unknown '--target-type' value {:?}", target_type));
  auto target_mass = pdt.constituentMass(target_pdg);

  // parse polarization type and spins -> set `spin_vec`, the spin vector for beam and target
  double spin_num[nObj]               = {0, 0};
  std::array<double,3> spin_vec[nObj] = { {0, 0, 0}, {0, 0, 0} };
  bool obj_is_polarized[nObj] = { false, false };
  enum spin_vec_enum { eX, eY, eZ };
  if(pol_type.length() != 2)
    return Error(fmt::format("option '--pol-type' value {:?} is not 2 characters", pol_type));
  for(int obj = 0; obj < nObj; obj++) {

    // parse polarization type
    std::string pol_type_name, spin_name;
    auto pol_type_char = std::toupper(pol_type.c_str()[obj]);

    // unpolarized
    if(pol_type_char == 'U') {
      pol_type_name = "unpolarized";
      spin_name = "0";
    }

    // polarized
    else {
      obj_is_polarized[obj] = true;

      // longitudinal or transverse
      switch(pol_type_char) {
        case 'L': pol_type_name = "longitudinal"; break;
        case 'T': pol_type_name = "transverse"; break;
        default:
          return Error(fmt::format("option '--pol-type' has unknown {} polarization type {:?}", obj_name[obj], pol_type.c_str()[obj]));
      }

      // use opposite sign for beam spin, since quark momentum reversed after hard scattering
      auto spin_sign = obj == objBeam ? -1.0 : 1.0;

      // parse spin type
      if(spin_type[obj].empty())
        return Error(fmt::format("option '--{}Spin' must be set when {} polarization is {}", obj_name[obj], obj_name[obj], pol_type_name));
      if(spin_type[obj].length() > 1)
        return Error(fmt::format("option '--{}Spin' value {:?} is not 1 character", obj_name[obj], spin_type[obj]));
      switch(std::tolower(spin_type[obj].c_str()[0])) {
        case 'p':
          {
            spin_num[obj] = 1.0;
            if(pol_type_char == 'L') { // longitudinal
              spin_name = "+";
              spin_vec[obj][eZ] = spin_sign;
            }
            else { // transverse
              spin_name = "up";
              spin_vec[obj][eY] = spin_sign;
            }
            break;
          }
        case 'm':
        case 'n':
          {
            spin_num[obj] = -1.0;
            if(pol_type_char == 'L') { // longitudinal
              spin_name = "-";
              spin_vec[obj][eZ] = -spin_sign;
            }
            else { // transverse
              spin_name = "down";
              spin_vec[obj][eY] = -spin_sign;
            }
            break;
          }
        default:
          return Error(fmt::format("option '--{}Spin' has unknown value {:?}", obj_name[obj], spin_type[obj]));
      }
    }

    Verbose(fmt::format("{:>30} = {}", fmt::format("{} polarization type", obj_name[obj]), pol_type_name));
    Verbose(fmt::format("{:>30} = {}", fmt::format("{} spin", obj_name[obj]), spin_name));
    Verbose(fmt::format("{:>30} = ({})", fmt::format("{} spin vector", obj == objBeam ? "quark" : obj_name[obj]), fmt::join(spin_vec[obj], ", ")));
  }

  // configure pythia
  /// plugin stringspinner hooks
  auto fhooks = std::make_shared<Pythia8::SimpleStringSpinner>();
  fhooks->plugInto(pyth);
  /// read config file
  apply_config_func(pyth);
  //// beam and target types
  set_config(pyth, fmt::format("Beams:idA = {}", BEAM_PDG));
  set_config(pyth, fmt::format("Beams:idB = {}", target_pdg));
  set_config(pyth, fmt::format("Beams:eA = {}", beam_energy));
  set_config(pyth, fmt::format("Beams:eB = {}", 0.0));
  //// seed
  set_config(pyth, "Random:setSeed = on");
  set_config(pyth, fmt::format("Random:seed = {}", seed));
  //// beam polarization
  if(obj_is_polarized[objBeam]) {
    for(auto quark : std::vector<std::string>{"u", "d", "s", "ubar", "dbar", "sbar"})
      set_config(pyth, fmt::format("StringSpinner:{}Polarisation = {}", quark, fmt::join(spin_vec[objBeam],",")));
  }
  //// target polarization
  if(obj_is_polarized[objTarget])
    set_config(pyth, fmt::format("StringSpinner:targetPolarisation = {}", fmt::join(spin_vec[objTarget],",")));
  //// finally, set the overridden parameters
  for(auto const& config_str : config_overrides)
    set_config(pyth, config_str);

  // initialize pythia
  pyth.init();

  // start LUND file: recreate it if it already exists
  auto lund_file = fmt::output_file(out_file, fmt::file::WRONLY | fmt::file::CREATE | fmt::file::TRUNC);

  // set `LundHeader` constant variables
  LundHeader lund_header{
    .target_mass       = target_mass,
    .target_atomic_num = target_atomic_num,
    .target_spin       = spin_num[objTarget],
    .beam_spin         = spin_num[objBeam],
    .beam_type         = BEAM_PDG,
    .beam_energy       = beam_energy,
    .nucleon_pdg       = target_pdg
  };

  ////////////////////////////////////////////////////////////////////
  // EVENT LOOP
  ////////////////////////////////////////////////////////////////////
  decltype(num_events) evnum = 0;
  while(true && num_events>0) {

    // next event
    if(enable_count_before_cuts && evnum >= num_events)
      break;
    if(!pyth.next())
      continue;
    Verbose(fmt::format(">>> EVENT {} <<<", evnum));
    if(enable_count_before_cuts)
      evnum++;

    // setup inclusive cut
    bool cut_inclusive_passed = false;
    decltype(cut_inclusive)::size_type n_found = 0;
    if(enable_cut_inclusive) {
      for(auto& [pdg, found] : cut_inclusive_found)
        found = false;
    }
    else cut_inclusive_passed = true;

    ///////////////////////////////////////////////////////////////////////////////

    // find the beam and target
    int beam_idx = -1;
    int target_idx = -1;
    for(int idx = 0; idx < evt.size(); ++idx) {
      auto const& par = evt[idx];
      if(par.status() == -12) {
        switch(par.id()) {
          case 11: beam_idx = idx; break;
          case 2212: target_idx = idx; break;
        }
      }
    }
    if(beam_idx < 0 || target_idx < 0)
      continue;

    // pythia's frame beam and target momenta
    auto beam_E   = evt[beam_idx].e(); // pythia's frame
    auto beam_p   = evt[beam_idx].pz();
    auto target_E = evt[target_idx].e();
    auto target_p = evt[target_idx].pz();

    // lab frame beam and target momenta (H is energy, q is pz)
    double beam_H   = beam_energy;
    double beam_q   = std::sqrt(std::pow(beam_energy,2) - std::pow(pdt.constituentMass(11),2));
    double target_H = target_mass;
    double target_q = 0.0;

    // solve for beta and gamma of the lorentz transformation from pythia frame to lab frame
    auto get_boost = [](double E, double p, double H, double q, double& beta, double& gamma, double &a, double &b) {
      double denom = (std::pow(E,2) - std::pow(p,2));
      a = (E*H - p*q) / denom;
      b = (H*p - E*q) / denom;
      // fmt::print("a={} b={}\n", a, b);
      gamma = a;
      beta = b / a;
    };
    double beam_gamma,   beam_beta, beam_a, beam_b;
    double target_gamma, target_beta, target_a, target_b;
    get_boost(beam_E,   beam_p,   beam_H,   beam_q,   beam_gamma,   beam_beta, beam_a, beam_b);
    get_boost(target_E, target_p, target_H, target_q, target_gamma, target_beta, target_a, target_b);

    // transform the event from pythia frame to lab frame
    // evt.bst(0., 0., -target_beta, target_gamma);
    // evt.bst(0., 0., -beam_beta, beam_gamma);

    // fmt::print("target_beta={} target_gamma={}\n", target_beta, target_gamma);
    // fmt::print("beam_beta={} beam_gamma={}\n", beam_beta, beam_gamma);

    // auto apply_boost = [beta=target_beta, gamma=target_gamma](Pythia8::Vec4 v) -> Pythia8::Vec4 {
    //   return Pythia8::Vec4(
    //       v.px(),
    //       v.py(),
    //       beta * gamma * v.e() + gamma * v.pz(),
    //       gamma * v.e() + beta * gamma * v.pz()
    //       );
    // };
    auto apply_boost = [a=beam_a, b=target_b](Pythia8::Vec4 v) -> Pythia8::Vec4 {
      return Pythia8::Vec4(
          v.px(),
          v.py(),
          -b*v.e() + a*v.pz(),
          a*v.e() - b*v.pz()
          );
    };

    // auto print_vec = [](std::string const& s, Pythia8::Vec4 v) {
    //   fmt::print("{}: {:5.12}  {:5.12}  {:5.12}  {:5.12}\n", s, v.px(), v.py(), v.pz(), v.e());
    // };
    // print_vec("  beam", evt[beam_idx].p());
    // print_vec("  beam", apply_boost(evt[beam_idx].p()));
    // print_vec("target", evt[target_idx].p());
    // print_vec("target", apply_boost(evt[target_idx].p()));

    evt.list(false, false, 12);

    ///////////////////////////////////////////////////////////////////////////////

    // loop over particles
    std::vector<LundParticle> lund_particles;
    Verbose("Particles:");
    Verbose(fmt::format("  {:-^10} {:-^10} {:-^12} {:-^12} {:-^12} {:-^12}", "pdg", "status", "px", "py", "pz", "theta"));
    for(auto const& par : evt) {

      // skip the "system" particle
      if(par.id() == 90)
        continue;

      if(enable_verbose_mode)
        Verbose(fmt::format("  {:10} {:10} {:12.5g} {:12.5g} {:12.5g} {:12.5g}",
              par.id(), par.status(), par.px(), par.py(), par.pz(), par.theta() * 180.0 / M_PI));

      // check if this particle is requested by `cut_inclusive`
      if(!cut_inclusive_passed && par.isFinal()) {
        for(auto& [pdg, found] : cut_inclusive_found) { // loop over `cut_inclusive` particles
          if(!found && pdg == par.id()) { // if we haven't found this one yet:

            // check if it's in `cut_theta`
            bool cut_theta_passed = ! enable_cut_theta;
            if(enable_cut_theta) {
              auto theta = par.theta() * 180.0 / M_PI;
              cut_theta_passed = theta >= cut_theta[0] && theta <= cut_theta[1];
            }
            if(cut_theta_passed) {
              Verbose("^^ good ^^");
              found = true;
              n_found++;
              break;
            }

          }
        }
        if(n_found == cut_inclusive.size()) // if all of them have been found
          cut_inclusive_passed = true;
      }

      // set lund particle variables
      int par_index = lund_particles.size() + 1;
      lund_particles.push_back({
          .index     = par_index,
          .lifetime  = par.isFinal() ? 1.0 : 0.0, // not used in GEMC; FIXME: should actually be something like `par.tau() * 1e6 / SPEED_OF_LIGHT`
          .status    = par.isFinal() ? 1 : 0,
          .pdg       = par.id(),
          .mother1   = par.mother1(),
          .daughter1 = par.daughter1(),
          .px        = par.px(),
          .py        = par.py(),
          .pz        = par.pz(),
          .energy    = par.e(),
          .mass      = par.m(),
          .vx        = par.xProd() / 10.0, // [mm] -> [cm]
          .vy        = par.yProd() / 10.0, // [mm] -> [cm]
          .vz        = par.zProd() / 10.0, // [mm] -> [cm]
          });
    }

    // apply inclusive cut
    if(!cut_inclusive_passed) {
      Verbose("cut '--cut-inclusive' did not pass");
      continue;
    }

    Verbose("All cuts PASSED");

    // set non-constant lund header variables
    lund_header.num_particles = evt.size() - 1; // one less than `evt.size()`, since PDG == 90 (entry 0) represents the system
    lund_header.process_id    = pyth.info.code();
    lund_header.event_weight  = pyth.info.weight();

    // stream to lund file
    lund_file.print("{} {:.5} {:} {:} {:} {:} {:.5} {:} {:} {:.5}\n",
      lund_header.num_particles,
      lund_header.target_mass,
      lund_header.target_atomic_num,
      lund_header.target_spin,
      lund_header.beam_spin,
      lund_header.beam_type,
      lund_header.beam_energy,
      lund_header.nucleon_pdg,
      lund_header.process_id,
      lund_header.event_weight
      );
    for(auto const& lund_particle : lund_particles)
      lund_file.print("{:} {:.5} {:} {:} {:} {:} {:.5} {:.5} {:.5} {:.5} {:.5} {:.5} {:.5} {:.5}\n",
          lund_particle.index,
          lund_particle.lifetime,
          lund_particle.status,
          lund_particle.pdg,
          lund_particle.mother1,
          lund_particle.daughter1,
          lund_particle.px,
          lund_particle.py,
          lund_particle.pz,
          lund_particle.energy,
          lund_particle.mass,
          lund_particle.vx,
          lund_particle.vy,
          lund_particle.vz
          );

    // finalize
    if(!enable_count_before_cuts) {
      if(++evnum >= num_events)
        break;
    }


    //////////////////////////////////////////////////////////////////////////////////

    /*
    for(int idx = 0; idx < evt.size(); ++idx) {
      auto const& par = evt[idx];
      if(par.isFinal()) {
        auto p = std::hypot(par.px(),par.py(),par.pz());
        auto mass = std::sqrt( std::pow(par.e(),2) - std::pow(p,2) );
        fmt::print("TEST: pdg={} p={} e()={} m()={} mass={} diff={}\n", par.id(), p, par.e(), par.m(), mass, par.m() - mass);
      }
    }
    */

    // find the scattered electron
    int ele_idx = -1;
    double ele_p = -10000;
    for(int idx = 0; idx < evt.size(); ++idx) {
      auto const& par = evt[idx];
      if(par.isFinal() && par.id() == 11) {
        auto this_p = std::hypot(par.px(), par.py(), par.pz());
        if(this_p > ele_p) {
          ele_p = par.e();
          ele_idx = idx;
        }
      }
    }
    // fmt::print("ele_idx = {}\n", ele_idx);
    if(ele_idx < 0) continue;

    // loop over pi+ pi- dihadrons
    for(int idxA = 0; idxA < evt.size(); ++idxA) {
      auto const& parA = evt[idxA];
      if(parA.isFinal() && parA.id() == 211) {
        for(int idxB = 0; idxB < evt.size(); ++idxB) {
          auto const& parB = evt[idxB];
          if(parB.isFinal() && parB.id() == -211) {

            // auto p_beam = evt[beam_idx].p();
            // auto p_target = evt[target_idx].p();
            // auto p_ele = evt[ele_idx].p();
            // auto p_pip = parA.p();
            // auto p_pim = parB.p();

            auto p_beam = apply_boost(evt[beam_idx].p());
            auto p_target = apply_boost(evt[target_idx].p());
            auto p_ele = apply_boost(evt[ele_idx].p());
            auto p_pip = apply_boost(parA.p());
            auto p_pim = apply_boost(parB.p());

            // calculate missing mass
            auto vecW = p_beam + p_target - p_ele;
            auto vecPh = p_pip + p_pim;
            auto M_X = (vecW - vecPh).mCalc();
            fmt::print("M_X {}\n", M_X);
          }
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////////


  } // end EVENT LOOP

  fmt::print("GENERATED LUND FILE: {}\n", out_file);
  fmt::print("   NUMBER OF EVENTS: {}\n", num_events);
  return 0;
}
