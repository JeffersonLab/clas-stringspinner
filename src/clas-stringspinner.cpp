#include <getopt.h>
#include <array>
#include <functional>
#include <optional>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <stringspinner/StringSpinner.h>

#include "CheckList.h"
#include "Lund.h"

//////////////////////////////////////////////////////////////////////////////////

// configurations
#include "config/clas12.h"
static std::map<std::string, std::function<void(Pythia8::Pythia&)>> CONFIG_MAP = {
  {"clas12", config_clas12}
};

//////////////////////////////////////////////////////////////////////////////////

// constants
int const SEED_MAX    = 900000000;
int const BEAM_PDG    = 11;
int const BEAM_ROW    = 1;
int const TARGET_ROW  = 2;

enum obj_enum { objBeam, objTarget, nObj };
std::string const obj_name[nObj] = { "beam", "target" };

// default option values
static unsigned long            num_events               = 10000;
static std::string              out_file                 = "clas-stringspinner.dat";
static double                   beam_energy              = 10.60410;
static std::string              target_type              = "proton";
static std::string              pol_type                 = "UU";
static std::string              spin_type[nObj]          = {"", ""};
static std::string              patch_boost              = "none";
static std::string              config_name              = "clas12";
static std::vector<std::string> config_overrides         = {};
static int                      seed                     = -1;
static bool                     enable_count_before_cuts = false;
static bool                     enable_patch_boost       = false;

// cut checklists
clas::CheckList cut_inclusive{"cut-inclusive"};
clas::CheckList cut_theta{"cut-theta"};
clas::CheckList cut_z{"cut-z"};

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

  --cut-inclusive PDG...           if set, event must include all particles with these
                                   PDG codes
                                   - PDG... is delimited by commas; no spaces
                                   - repeat PDG codes to require more than one
                                   - example: 1 pi- and 2 pi+s:
                                       --cut-inclusive -211,211,211

  --cut-theta MIN,MAX,PDG...       MIN <= theta <= MAX, for all particles in PDG...
                                   - example: charged pions in 10-30 degrees:
                                       --cut-theta 10,30,211,-211

  --cut-z MIN,MAX,PDG...           MIN <= hadron z <= MAX, for all hadrons in PDG...

MISCELLANEOUS OPTIONS:

  --patch-boost                    temporary patch for boost issue:
                                   https://gitlab.com/Pythia8/releases/-/issues/529
                                   this option ensures the event record is boosted
                                   back to the fixed-target rest frame
                                   - needed for Pythia v8.312, and possibly earlier
                                   - available choices:
                                     'beam'   = derive boost from beam momentum
                                     'target' = derive boost from target momentum
                                     'none'   = do not apply any boost
                                   default: {patch_boost:?}

OPTIONS FOR OSG COMPATIBILITY:

  --trig NUM_EVENTS                same as --num-events
  --ebeam ENERGY                   same as --beam-energy
  --docker                         unused
    )" + std::string("\n"),
      fmt::arg("patch_boost", patch_boost),
      fmt::arg("num_events", num_events),
      fmt::arg("out_file", out_file),
      fmt::arg("beam_energy", beam_energy),
      fmt::arg("target_type", target_type),
      fmt::arg("pol_type", pol_type),
      fmt::arg("config_name_list", fmt::join(config_name_list, "\n                                            ")),
      fmt::arg("config_name", config_name),
      fmt::arg("seed", seed),
      fmt::arg("seed_max", SEED_MAX)
      );
}

//////////////////////////////////////////////////////////////////////////////////

/// @returns the scattered lepton index, if found
/// @param evt the pythia event
std::optional<int> FindScatteredLepton(Pythia8::Event const& evt)
{
  for(auto const& par : evt) {
    if(par.id() == BEAM_PDG && par.status() == 23) { // hardest subprocess particle, outgoing
      for(auto const& mom_idx : std::vector<int>{par.mother1(), par.mother2()}) {
        auto const& mom = evt.at(mom_idx);
        // if status is subprocess incoming, and mother is beam
        if(mom.id() == BEAM_PDG && mom.status() == 21 && (mom.mother1() == BEAM_ROW || mom.mother2() == BEAM_ROW))
          return par.index();
      }
    }
  }
  return std::nullopt;
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{

  // parse arguments
  enum options_enum {
    opt_num_events,
    opt_docker,
    opt_out_file,
    opt_beam_energy,
    opt_target_type,
    opt_pol_type,
    opt_beam_spin,
    opt_target_spin,
    opt_cut_inclusive,
    opt_cut_theta,
    opt_cut_z,
    opt_config,
    opt_seed,
    opt_set,
    opt_patch_boost,
    opt_help,
    opt_version,
    opt_count_before_cuts,
    opt_verbose
  };
  struct option const opts[] = {
    {"num-events",        required_argument, nullptr, opt_num_events},
    {"trig",              required_argument, nullptr, opt_num_events},
    {"docker",            no_argument,       nullptr, opt_docker},
    {"out-file",          required_argument, nullptr, opt_out_file},
    {"beam-energy",       required_argument, nullptr, opt_beam_energy},
    {"ebeam",             required_argument, nullptr, opt_beam_energy},
    {"target-type",       required_argument, nullptr, opt_target_type},
    {"pol-type",          required_argument, nullptr, opt_pol_type},
    {"beam-spin",         required_argument, nullptr, opt_beam_spin},
    {"target-spin",       required_argument, nullptr, opt_target_spin},
    {"cut-inclusive",     required_argument, nullptr, opt_cut_inclusive},
    {"cut-theta",         required_argument, nullptr, opt_cut_theta},
    {"cut-z",             required_argument, nullptr, opt_cut_z},
    {"config",            required_argument, nullptr, opt_config},
    {"seed",              required_argument, nullptr, opt_seed},
    {"set",               required_argument, nullptr, opt_set},
    {"patch-boost",       required_argument, nullptr, opt_patch_boost},
    {"help",              no_argument,       nullptr, opt_help},
    {"version",           no_argument,       nullptr, opt_version},
    {"count-before-cuts", no_argument,       nullptr, opt_count_before_cuts},
    {"verbose",           no_argument,       nullptr, opt_verbose},
    {nullptr, 0, nullptr, 0}
  };

  if(argc <= 1) {
    Usage();
    return clas::EXIT_SYNTAX;
  };

  char opt;
  while((opt = getopt_long(argc, argv, "", opts, nullptr)) != -1) {
    switch(opt) {
      case opt_num_events: num_events = std::stol(optarg); break;
      case opt_out_file: out_file = std::string(optarg); break;
      case opt_beam_energy: beam_energy = std::stod(optarg); break;
      case opt_target_type: target_type = std::string(optarg); break;
      case opt_pol_type: pol_type = std::string(optarg); break;
      case opt_beam_spin: spin_type[objBeam] = std::string(optarg); break;
      case opt_target_spin: spin_type[objTarget] = std::string(optarg); break;
      case opt_cut_inclusive: cut_inclusive.Setup(optarg, false); break;
      case opt_cut_theta: cut_theta.Setup(optarg); break;
      case opt_cut_z: cut_z.Setup(optarg); break;
      case opt_config: config_name = std::string(optarg); break;
      case opt_seed: seed = std::stoi(optarg); break;
      case opt_set: config_overrides.push_back(std::string(optarg)); break;
      case opt_patch_boost: patch_boost = std::string(optarg); break;
      case opt_count_before_cuts: enable_count_before_cuts = true; break;
      case opt_verbose: clas::enable_verbose_mode = true; break;
      case opt_help:
        Usage();
        return 0;
      case opt_version:
        fmt::print("{}\n", CLAS_STRINGSPINNER_VERSION);
        return 0;
      case '?':
        return clas::EXIT_ERROR;
    }
  }

  // check if seed is too large; if so, % SEED_MAX
  if(seed > SEED_MAX) {
    auto new_seed = seed % SEED_MAX;
    clas::Error(fmt::format("value of option '--seed' is too large for Pythia8: {} > {}; setting it to `seed % {}` = {}", seed, SEED_MAX, SEED_MAX, new_seed));
    seed = new_seed;
  }

  // print options
  clas::Verbose(fmt::format("{:=^82}", " Arguments "));
  clas::Verbose(fmt::format("{:>30} = {}", "num-events", num_events));
  clas::Verbose(fmt::format("{:>30} = {}", "count-before-cuts", enable_count_before_cuts ? "true" : "false"));
  clas::Verbose(fmt::format("{:>30} = {:?}", "out-file", out_file));
  clas::Verbose(fmt::format("{:>30} = {} GeV", "beam-energy", beam_energy));
  clas::Verbose(fmt::format("{:>30} = {:?}", "target-type", target_type));
  clas::Verbose(fmt::format("{:>30} = {:?}", "pol-type", pol_type));
  clas::Verbose(fmt::format("{:>30} = {:?}", "beam-spin", spin_type[objBeam]));
  clas::Verbose(fmt::format("{:>30} = {:?}", "target-spin", spin_type[objTarget]));
  clas::Verbose(fmt::format("{:>30} = {}", "cut-inclusive", cut_inclusive.GetInfoString()));
  clas::Verbose(fmt::format("{:>30} = {}", "cut-theta", cut_theta.GetInfoString()));
  clas::Verbose(fmt::format("{:>30} = {}", "cut-z", cut_z.GetInfoString()));
  clas::Verbose(fmt::format("{:>30} = {:?}", "patch-boost", patch_boost));
  clas::Verbose(fmt::format("{:>30} = {}", "seed", seed));
  clas::Verbose(fmt::format("{:>30} = {}", "config", config_name));
  clas::Verbose(fmt::format("{:-^82}", ""));
  clas::Verbose(fmt::format("{}{}", "parameter overrides (from option '--set'):", config_overrides.empty() ? " none" : ""));
  for(auto const& config_str : config_overrides)
    clas::Verbose(fmt::format("- {:?}", config_str));
  clas::Verbose(fmt::format("{:=^82}", ""));

  // initialize pythia
  Pythia8::Pythia pyth;
  auto& evt  = pyth.event;
  auto& proc = pyth.process;
  auto& pdt  = pyth.particleData;

  // get the configuration function
  std::function<void(Pythia8::Pythia&)> apply_config_func;
  try {
    apply_config_func = CONFIG_MAP.at(config_name);
  }
  catch(std::out_of_range const& ex) {
    clas::Error(fmt::format("value of option '--config' is {:?}, which is not found", config_name));
    return clas::EXIT_ERROR;
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
  else return clas::Error(fmt::format("unknown '--target-type' value {:?}", target_type));
  auto target_mass = pdt.constituentMass(target_pdg);

  // parse polarization type and spins -> set `spin_vec`, the spin vector for beam and target
  double spin_num[nObj]               = {0, 0};
  std::array<double,3> spin_vec[nObj] = { {0, 0, 0}, {0, 0, 0} };
  bool obj_is_polarized[nObj] = { false, false };
  enum spin_vec_enum { eX, eY, eZ };
  if(pol_type.length() != 2)
    return clas::Error(fmt::format("option '--pol-type' value {:?} is not 2 characters", pol_type));
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
          return clas::Error(fmt::format("option '--pol-type' has unknown {} polarization type {:?}", obj_name[obj], pol_type.c_str()[obj]));
      }

      // use opposite sign for beam spin, since quark momentum reversed after hard scattering
      auto spin_sign = obj == objBeam ? -1.0 : 1.0;

      // parse spin type
      if(spin_type[obj].empty())
        return clas::Error(fmt::format("option '--{}Spin' must be set when {} polarization is {}", obj_name[obj], obj_name[obj], pol_type_name));
      if(spin_type[obj].length() > 1)
        return clas::Error(fmt::format("option '--{}Spin' value {:?} is not 1 character", obj_name[obj], spin_type[obj]));
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
          return clas::Error(fmt::format("option '--{}Spin' has unknown value {:?}", obj_name[obj], spin_type[obj]));
      }
    }

    clas::Verbose(fmt::format("{:>30} = {}", fmt::format("{} polarization type", obj_name[obj]), pol_type_name));
    clas::Verbose(fmt::format("{:>30} = {}", fmt::format("{} spin", obj_name[obj]), spin_name));
    clas::Verbose(fmt::format("{:>30} = ({})", fmt::format("{} spin vector", obj == objBeam ? "quark" : obj_name[obj]), fmt::join(spin_vec[obj], ", ")));
  }

  // settings for boost patch, for boosting the Pythia Event record frame back to the lab frame (fixed-target rest frame)
  int patch_boost_particle_row = -1;
  int patch_boost_particle_pdg = 0;
  if(patch_boost == "beam") {
    enable_patch_boost = true;
    patch_boost_particle_row = BEAM_ROW;
    patch_boost_particle_pdg = BEAM_PDG;
  } else if(patch_boost == "target") {
    enable_patch_boost = true;
    patch_boost_particle_row = TARGET_ROW;
    patch_boost_particle_pdg = target_pdg;
  } else if(patch_boost == "none") {
    enable_patch_boost = false;
  } else {
    return clas::Error(fmt::format("option '--patch-boost' has unknown value {:?}", patch_boost));
  }

  // configure pythia
  //// plugin stringspinner hooks
  auto fhooks = std::make_shared<Pythia8::SimpleStringSpinner>();
  fhooks->plugInto(pyth);
  //// read config file
  apply_config_func(pyth);
  //// set verbosity
  set_config(pyth, fmt::format("Next:numberShowEvent = {}", clas::enable_verbose_mode ? num_events : 0));
  // set_config(pyth, fmt::format("Next:numberShowProcess = {}", clas::enable_verbose_mode ? num_events : 0));
  // set_config(pyth, fmt::format("Next:numberShowInfo = {}", clas::enable_verbose_mode ? num_events : 0));
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
  clas::LundHeader lund_header{
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
    clas::Verbose(fmt::format("\n>>> EVENT {} ======================================================================", evnum));
    if(!pyth.next())
      continue;
    if(enable_count_before_cuts)
      evnum++;

    // boost event record back to lab frame
    // see <https://gitlab.com/Pythia8/releases/-/issues/529> for details
    if(enable_patch_boost) {
      auto const& par__evt  = evt[patch_boost_particle_row];  // beam (or target) momentum in event frame
      auto const& par__proc = proc[patch_boost_particle_row]; // beam (or target) momentum in hard-process frame, which is assumed to be the lab frame
      // check that we are using the correct beam (or target) particle
      if(par__evt.status() != -12) {
        clas::EventError("patch-boost particle is not an incoming beam (or target) particle");
        continue;
      }
      if(par__evt.id() != patch_boost_particle_pdg) {
        clas::EventError(fmt::format("patch-boost particle does not have the expected PDG: {} != {}", par__evt.id(), patch_boost_particle_pdg));
        continue;
      }
      if(par__evt.id() != par__proc.id()) {
        clas::EventError("patch-boost particle PDG mismatch between event record and hard-process record");
        continue;
      }
      // perform the boost
      Pythia8::RotBstMatrix boost_to_lab;
      boost_to_lab.bst(par__evt.p(), par__proc.p());
      evt.rotbst(boost_to_lab);
    }
    // check that the event-record frame matches the hard-process frame, which is assumed to be the lab frame
    // for(auto const& [name, row] : std::vector<std::pair<std::string,int>>{{"beam", BEAM_ROW}, {"target", TARGET_ROW}}) {
    //   auto diff = std::max(
    //       std::abs(evt[row].pz() - proc[row].pz()),
    //       std::abs(evt[row].e()  - proc[row].e())
    //       );
    //   if(diff > 0.0001)
    //     EventError(fmt::format("mismatch of event-frame and hard-process-frame {} momentum; use '--verbose' for details', and consider changing the value of the '--patch-boost' option", name));
    //   Verbose(fmt::format("hard process {:<8} pz = {:<20.10}  E = {:<20.10}", name, proc[row].pz(), proc[row].e()));
    //   Verbose(fmt::format("event record {:<8} pz = {:<20.10}  E = {:<20.10}", name, evt[row].pz(),  evt[row].e()));
    // }

    // check required inclusive particles
    if(!cut_inclusive.Check(evt))
      continue;

    // check theta cuts
    auto get_theta = [](Pythia8::Particle const& par) {
      return par.theta() * 180.0 / M_PI;
    };
    if(!cut_theta.Check(evt, get_theta))
      continue;

    // check z cuts
    std::map<int,double> z_vals;
    if(cut_z.Enabled() || clas::enable_verbose_mode) {
      // find scattered lepton
      auto const lepton_idx = FindScatteredLepton(evt);
      if(!lepton_idx.has_value()) { // no scattered lepton -> skip event
        clas::Verbose("no scattered lepton found");
        continue;
      }
      // virtaul photon momentum
      auto const vec_q = evt.at(BEAM_PDG).p() - evt.at(lepton_idx.value()).p();
      // calculate z
      auto const vec_target = evt.at(TARGET_ROW).p();
      auto get_z = [&vec_target, &vec_q] (Pythia8::Particle const& par) {
        return (vec_target * par.p()) / (vec_target * vec_q);
      };
      // check z cuts
      if(!cut_z.Check(evt, get_z))
        continue;
      // fill `z_vals` for verbose printout
      if(clas::enable_verbose_mode) {
        for(auto const& par : evt)
          z_vals.insert({par.index(), par.isFinal() ? get_z(par) : -1});
      }
    }

    // verbose printout
    if(clas::enable_verbose_mode) {
      clas::Verbose("\nAll cuts passed!\nParticle Kinematics:");
      clas::Verbose(fmt::format("  {:-^6} {:-^10} {:-^8} {:-^12} {:-^12}", "id", "pdg", "status", "theta", "z"));
      for(auto const& par : evt) {
        clas::Verbose(fmt::format("  {:6} {:10} {:8} {:12.5g} {:12.5g}",
              par.index(),
              par.id(),
              par.status(),
              get_theta(par),
              z_vals.at(par.index())
              ));
      }
    }


    // event passed all cuts -> write to output file
    std::vector<clas::LundParticle> lund_particles;
    for(auto const& par : evt) {

      // skip the "system" particle
      if(par.id() == 90)
        continue;

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

    // set non-constant lund header variables
    lund_header.num_particles = evt.size() - 1; // one less than `evt.size()`, since PDG == 90 (entry 0) represents the system
    lund_header.process_id    = pyth.info.code();
    lund_header.event_weight  = pyth.info.weight();

    // stream to lund file
    lund_header.Stream(lund_file);
    for(auto const& lund_particle : lund_particles)
      lund_particle.Stream(lund_file);

    // finalize
    if(!enable_count_before_cuts) {
      if(++evnum >= num_events)
        break;
    }

  } // end EVENT LOOP

  fmt::print("GENERATED LUND FILE: {}\n", out_file);
  fmt::print("   NUMBER OF EVENTS: {}\n", num_events);
  return 0;
}
