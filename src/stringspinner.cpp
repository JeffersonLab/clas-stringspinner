#include <getopt.h>
#include <array>
#include <filesystem>
#include <functional>
#include <fmt/format.h>
#include <fmt/os.h>
#include <Pythia8/Pythia.h>
#include <stringspinner/StringSpinner.h>

const int EXIT_ERROR = 1;
const int EXIT_SYNTAX = 2;

const int BEAM_PDG = 11;

enum obj_enum { objBeam, objTarget, nObj };
const std::string obj_name[nObj] = { "beam", "target" };

static unsigned long num_events = 10000;

static std::string      out_file             = "out.lund";
static int              verbose_mode         = 0;
static double           beam_energy          = 10.60410;
static std::string      target_type          = "proton";
static std::string      pol_type             = "UU";
static std::string      spin_type[nObj]      = {"", ""};
static double           glgt_mag             = 0.2;
static double           glgt_arg             = 0.0;
static std::vector<int> cut_string           = {2,  2101};
static std::vector<int> cut_inclusive        = {};
static bool             enable_cut_string    = false;
static bool             enable_cut_inclusive = false;
static std::string      config_file          = "clas12.cmnd";
static int              seed                 = -1;
static int              float_precision      = 5;

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
  fmt::print("USAGE: stringspinner [OPTIONS]...\n\n");
  fmt::print("  --numEvents NUM_EVENTS           number of events\n");
  fmt::print("                                   - if you apply cuts, such as `--cutString`, the real number\n");
  fmt::print("                                     of events that survive the cuts will be smaller than `NUM_EVENTS`\n");
  fmt::print("                                   default: {}\n\n", num_events);
  fmt::print("  --outFile OUTPUT_FILE            output file name\n");
  fmt::print("                                   default: {:?}\n\n", out_file);
  fmt::print("  --beamEnergy ENERGY              electron beam energy [GeV]\n");
  fmt::print("                                   default: {}\n\n", beam_energy);
  fmt::print("  --targetType TARGET_TYPE         target type, one of:\n");
  fmt::print("                                     proton\n");
  fmt::print("                                     neutron\n");
  fmt::print("                                   default: {:?}\n\n", target_type);
  fmt::print("  --polType POLARIZATION_TYPE      beam and target polarization types\n");
  fmt::print("                                   - two characters: beam and target\n");
  fmt::print("                                   - types: 'U' = unpolarized\n");
  fmt::print("                                            'L' = longitudinally polarized\n");
  fmt::print("                                            'T' = transversely polarized\n");
  fmt::print("                                   default: {:?}\n\n", pol_type);
  fmt::print("  --beamSpin BEAM_SPIN             the spin of the beam leptons\n");
  fmt::print("                                   - if longitudinally polarized ('L'):\n");
  fmt::print("                                     'p' = spin along +z axis\n");
  fmt::print("                                     'm' = spin along -z axis\n");
  fmt::print("                                   - if transversely polarized ('T'):\n");
  fmt::print("                                     'p' = spin along +y axis\n");
  fmt::print("                                     'm' = spin along -y axis\n");
  fmt::print("                                   - if unpolarized ('U'): no effect\n\n");
  fmt::print("  --targetSpin TARGET_SPIN         the spin of the target nucleons\n");
  fmt::print("                                   - same usage as --beamSpin, applied to target\n\n");
  fmt::print("  --glgtMag GLGT_MAGNITUDE         StringSpinner parameter |G_L/G_T|\n");
  fmt::print("                                   - fraction of longitudinally polarized vector mesons:\n");
  fmt::print("                                       f_L = |G_L/G_T|^2 / ( 2 + |G_L/G_T|^2 )\n");
  fmt::print("                                       with 0 <= f_L <= 1\n");
  fmt::print("                                   default: {}\n\n", glgt_mag);
  fmt::print("  --glgtArg GLGT_ARGUMENT          StringSpinner parameter theta_{{LT}} = arg(G_L/G_T)\n");
  fmt::print("                                   - related to vector meson oblique polarization\n");
  fmt::print("                                   - range: -PI <= theta_{{LT}} <= +PI\n");
  fmt::print("                                   default: {}\n\n", glgt_arg);
  fmt::print("  --cutString OBJ1,OBJ2            filter by strings, where OBJ1 and OBJ2 are PDG codes of quarks or diquarks;\n");
  fmt::print("                                   - PDG codes must be separated by a comma, with no spaces\n");
  fmt::print("                                   - examples:\n");
  fmt::print("                                       --cutString 2,2101  # selects 'u === (ud)_0' strings\n");
  fmt::print("                                       --cutString 0,0     # disable string selection\n");
  fmt::print("                                   default: {},{}\n\n", cut_string[0], cut_string[1]);
  fmt::print("  --cutInclusive PDG_CODES...      only allow events which have a least these particles\n");
  fmt::print("                                   - delimit by commas\n");
  fmt::print("                                   - repeat PDG codes to require more than one\n");
  fmt::print("                                   - example: 1 pi- and 2 pi+s: --cutInclusive -211,211,211\n\n");
  fmt::print("  --config CONFIG_FILE             choose a configuration file from one of the following:\n");
  for(auto const& entry : std::filesystem::directory_iterator(STRINGSPINNER_ETC))
    fmt::print("                                       {}\n", entry.path().filename().string());
  fmt::print("                                   default: {:?}\n\n", config_file);
  fmt::print("  --seed SEED                      random number generator seed:\n");
  fmt::print("                                       default seed: -1\n");
  fmt::print("                                      based on time:  0\n");
  fmt::print("                                         fixed seed:  1 to 900_000_000\n\n");
  fmt::print("  --floatPrecision PRECISION       floating point numerical precision for output files\n");
  fmt::print("                                   default: {}\n\n", float_precision);
  fmt::print("  --verbose                        verbose printout\n\n");
  fmt::print("  --help                           print this usage guide\n\n");
  fmt::print("NOTES:\n\n");
  fmt::print("  - view configuration files in {}\n", STRINGSPINNER_ETC);
  fmt::print("\n");
}

void Verbose(std::string msg)
{
  if(verbose_mode==1)
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
    {"numEvents",      required_argument, nullptr,       'n'},
    {"outFile",        required_argument, nullptr,       'o'},
    {"beamEnergy",     required_argument, nullptr,       'e'},
    {"targetType",     required_argument, nullptr,       'T'},
    {"polType",        required_argument, nullptr,       'p'},
    {"beamSpin",       required_argument, nullptr,       'b'},
    {"targetSpin",     required_argument, nullptr,       't'},
    {"glgtMag",        required_argument, nullptr,       'm'},
    {"glgtArg",        required_argument, nullptr,       'a'},
    {"cutString",      required_argument, nullptr,       'q'},
    {"cutInclusive",   required_argument, nullptr,       'I'},
    {"config",         required_argument, nullptr,       'c'},
    {"seed",           required_argument, nullptr,       's'},
    {"floatPrecision", required_argument, nullptr,       'f'},
    {"verbose",        no_argument,       &verbose_mode, 1},
    {"help",           no_argument,       nullptr,       'h'},
    {nullptr,          0,                 nullptr,       0}
  };

  if(argc <= 1) {
    Usage();
    return EXIT_SYNTAX;
  };

  char opt;
  while((opt = getopt_long(argc, argv, "", opts, nullptr)) != -1) {
    switch(opt) {
      case 'n': num_events = std::stol(optarg); break;
      case 'o': out_file = std::string(optarg); break;
      case 'e': beam_energy = std::stod(optarg); break;
      case 'T': target_type = std::string(optarg); break;
      case 'p': pol_type = std::string(optarg); break;
      case 'b': spin_type[objBeam] = std::string(optarg); break;
      case 't': spin_type[objTarget] = std::string(optarg); break;
      case 'm': glgt_mag = std::stod(optarg); break;
      case 'a': glgt_arg = std::stod(optarg); break;
      case 'q': {
        cut_string.clear();
        Tokenize(optarg, [&](auto token, auto i) { cut_string.push_back(std::stoi(token)); });
        if(cut_string.size() != 2)
          return Error("value of option '--cutString' does not have 2 arguments");
        break;
      }
      case 'I':
        Tokenize(optarg, [&](auto token, auto i) { cut_inclusive.push_back(std::stoi(token)); });
        break;
      case 'c': config_file = std::string(optarg); break;
      case 's': seed = std::stoi(optarg); break;
      case 'f': float_precision = std::stoi(optarg); break;
      case 'h':
        Usage();
        return EXIT_SYNTAX;
      case '?':
        return EXIT_ERROR;
    }
  }

  enable_cut_string    = ! (cut_string[0] == 0 && cut_string[1] == 0);
  enable_cut_inclusive = ! cut_inclusive.empty();
  std::vector<std::pair<int, bool>> cut_inclusive_found;
  for(auto pdg : cut_inclusive)
    cut_inclusive_found.push_back({pdg, false});

  Verbose(fmt::format("{:=^82}", " Arguments "));
  Verbose(fmt::format("{:>30} = {}", "numEvents", num_events));
  Verbose(fmt::format("{:>30} = {:?}", "outFile", out_file));
  Verbose(fmt::format("{:>30} = {} GeV", "beamEnergy", beam_energy));
  Verbose(fmt::format("{:>30} = {:?}", "targetType", target_type));
  Verbose(fmt::format("{:>30} = {:?}", "polType", pol_type));
  Verbose(fmt::format("{:>30} = {:?}", "beamSpin", spin_type[objBeam]));
  Verbose(fmt::format("{:>30} = {:?}", "targetSpin", spin_type[objTarget]));
  Verbose(fmt::format("{:>30} = {}", "|G_L/G_T|", glgt_mag));
  Verbose(fmt::format("{:>30} = {}", "arg(G_L/G_T)", glgt_arg));
  Verbose(fmt::format("{:>30} = ({})===({})  [{}]", "cutString", cut_string[0], cut_string[1], enable_cut_string ? "enabled" : "disabled"));
  Verbose(fmt::format("{:>30} = ({}) [{}]", "cutInclusive", fmt::join(cut_inclusive, ", "), enable_cut_inclusive ? "enabled" : "disabled"));
  Verbose(fmt::format("{:>30} = {}", "seed", seed));
  Verbose(fmt::format("{:>30} = {}", "config", config_file));
  Verbose(fmt::format("{:=^82}", ""));

  // initialize pythia
  Pythia8::Pythia pyth;
  Pythia8::Event& evt = pyth.event;
  Pythia8::ParticleData& pdt = pyth.particleData;

  // get path to configuration file
  // - must be installed in `STRINGSPINNER_ETC`
  // - take only `filename()` from user specified argument, to prevent them from using `../` to
  //   leave the `STRINGSPINNER_ETC` directory
  auto config_file_path = std::string(STRINGSPINNER_ETC) + "/" + std::filesystem::path{config_file}.filename().string();
  Verbose(fmt::format("config file path: {}", config_file_path));
  Verbose(fmt::format("{:=^82}", ""));

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
  else return Error(fmt::format("unknown '--targetType' value {:?}", target_type));
  auto target_mass = pdt.constituentMass(target_pdg);

  // parse polarization type and spins -> set `spin_vec`, the spin vector for beam and target
  double spin_num[nObj]               = {0, 0};
  std::array<double,3> spin_vec[nObj] = { {0, 0, 0}, {0, 0, 0} };
  bool obj_is_polarized[nObj] = { false, false };
  enum spin_vec_enum { eX, eY, eZ };
  if(pol_type.length() != 2)
    return Error(fmt::format("option '--polType' value {:?} is not 2 characters", pol_type));
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
          return Error(fmt::format("option '--polType' has unknown {} polarization type {:?}", obj_name[obj], pol_type.c_str()[obj]));
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
  pyth.readFile(config_file_path);
  //// beam and target types
  pyth.readString(fmt::format("Beams:idA = {}", BEAM_PDG));
  pyth.readString(fmt::format("Beams:idB = {}", target_pdg));
  pyth.readString(fmt::format("Beams:eA = {}", beam_energy));
  pyth.readString(fmt::format("Beams:eB = {}", target_mass));
  //// seed
  pyth.readString("Random:setSeed = on");
  pyth.readString(fmt::format("Random:seed = {}", seed));
  //// beam polarization
  if(obj_is_polarized[objBeam]) {
    for(auto quark : std::vector<std::string>{"u", "d", "s", "ubar", "dbar", "sbar"})
      pyth.readString(fmt::format("StringSpinner:{}Polarisation = {}", quark, fmt::join(spin_vec[objBeam],",")));
  }
  //// target polarization
  if(obj_is_polarized[objTarget])
    pyth.readString(fmt::format("StringSpinner:targetPolarisation = {}", fmt::join(spin_vec[objTarget],",")));
  //// stringspinner free parameters
  pyth.readString(fmt::format("StringSpinner:GLGT = {}", glgt_arg));
  pyth.readString(fmt::format("StringSpinner:thetaLT = {}", glgt_mag));

  // initialize pythia
  pyth.init();

  // start LUND file: recreate it if it already exists
  auto lund_file = fmt::output_file(out_file, fmt::file::WRONLY | fmt::file::CREATE | fmt::file::TRUNC);

  ////////////////////////////////////////////////////////////////////
  // EVENT LOOP
  ////////////////////////////////////////////////////////////////////
  for (decltype(num_events) e = 0; e < num_events; e++) {
    Verbose(fmt::format(">>> EVENT {} <<<", e));
    if(!pyth.next())
      continue;

    // string cut
    if(enable_cut_string && (evt[7].id() != cut_string[0] || evt[8].id() != cut_string[1])) {
      Verbose("cutString did not pass");
      continue;
    }

    // setup inclusive cut
    bool cut_inclusive_passed = false;
    decltype(cut_inclusive)::size_type n_found = 0;
    if(enable_cut_inclusive) {
      for(auto& [pdg, found] : cut_inclusive_found)
        found = false;
    }
    else cut_inclusive_passed = true;

    // loop over particles
    std::vector<LundParticle> lund_particles;
    Verbose("Particles:");
    Verbose(fmt::format("  {:>10} {:>10} {:>20} {:>20} {:>20}", "pdg", "status", "px", "py", "pz"));
    for(auto const& par : evt) {

      // skip the "system" particle
      if(par.id() == 90)
        continue;

      Verbose(fmt::format("  {:10} {:10} {:20.5g} {:20.5g} {:20.5g}", par.id(), par.status(), par.px(), par.py(), par.pz()));

      // check if this particle is requested by `cut_inclusive`
      if(!cut_inclusive_passed && par.isFinal()) {
        for(auto& [pdg, found] : cut_inclusive_found) {
          if(!found && pdg == par.id()) {
            found = true;
            n_found++;
            break;
          }
        }
        if(n_found == cut_inclusive.size())
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
      Verbose("cutInclusive did not pass");
      continue;
    }

    Verbose("All cuts PASSED");

    // set lund header variables
    LundHeader lund_header{
      .num_particles     = evt.size() - 1, // one less than `evt.size()`, since PDG == 90 (entry 0) represents the system
      .target_mass       = target_mass,
      .target_atomic_num = target_atomic_num,
      .target_spin       = spin_num[objTarget],
      .beam_spin         = spin_num[objBeam],
      .beam_type         = BEAM_PDG,
      .beam_energy       = beam_energy,
      .nucleon_pdg       = target_pdg,
      .process_id        = pyth.info.code(),
      .event_weight      = pyth.info.weight()
    };

    // true inclusive kinematics
    Pythia8::DISKinematics inc_kin(evt[1].p(), evt[5].p(), evt[2].p()); // TODO: write this to a separate file

    // stream to lund file
    lund_file.print("{} {:.{prec}} {:} {:} {:} {:} {:.{prec}} {:} {:} {:.{prec}}\n",
      lund_header.num_particles,
      lund_header.target_mass,
      lund_header.target_atomic_num,
      lund_header.target_spin,
      lund_header.beam_spin,
      lund_header.beam_type,
      lund_header.beam_energy,
      lund_header.nucleon_pdg,
      lund_header.process_id,
      lund_header.event_weight,
      fmt::arg("prec", float_precision)
      );
    for(auto const& lund_particle : lund_particles)
      lund_file.print("{:} {:.{prec}} {:} {:} {:} {:} {:.{prec}} {:.{prec}} {:.{prec}} {:.{prec}} {:.{prec}} {:.{prec}} {:.{prec}} {:.{prec}}\n",
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
          lund_particle.vz,
          fmt::arg("prec", float_precision)
          );

  } // end EVENT LOOP

  return 0;
}
