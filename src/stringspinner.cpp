#include <getopt.h>
#include <array>
#include <filesystem>
#include <fmt/format.h>
#include <Pythia8/Pythia.h>
#include <stringspinner/StringSpinner.h>

const int EXIT_ERROR = 1;
const int EXIT_SYNTAX = 2;

enum obj_enum { objBeam, objTarget, nObj };
const std::string obj_name[nObj] = { "beam", "target" };

static unsigned long num_events     = 10000;
static std::string out_file         = "out.lund";
static int verbose_mode             = 0;
static double beam_energy           = 10.60410;
static std::string target_type      = "proton";
static std::string pol_type         = "UU";
static std::string spin_type[nObj]  = {"", ""};
static double glgt_mag              = 0.2;
static double glgt_arg              = 0.0;
static int string_selection[2]      = {2, 2101};
static bool enable_string_selection = false;
static std::string config_file      = "clas12.cmnd";
static int seed                     = -1;

//////////////////////////////////////////////////////////////////////////////////

void Usage()
{
  fmt::print("USAGE: stringspinner [OPTIONS]...\n\n");
  fmt::print("  --numEvents NUM_EVENTS           number of events\n");
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
  fmt::print("  --selectString OBJ1,OBJ2         filter by strings, where OBJ1 and OBJ2 are PDG codes of quarks or diquarks;\n");
  fmt::print("                                   - PDG codes must be separated by a comma, with no spaces\n");
  fmt::print("                                   - examples:\n");
  fmt::print("                                       --selectString 2,2101  # selects 'u === (ud)_0' strings\n");
  fmt::print("                                       --selectString 0,0     # disable string selection\n");
  fmt::print("                                   default: {},{}\n\n", string_selection[0], string_selection[1]);
  fmt::print("  --config CONFIG_FILE             choose a configuration file from one of the following:\n");
  for(auto const& entry : std::filesystem::directory_iterator(STRINGSPINNER_ETC))
    fmt::print("                                       {}\n", entry.path().filename().string());
  fmt::print("                                   default: {:?}\n\n", config_file);
  fmt::print("  --seed SEED                      random number generator seed:\n");
  fmt::print("                                       default seed: -1\n");
  fmt::print("                                      based on time:  0\n");
  fmt::print("                                         fixed seed:  1 to 900_000_000\n\n");
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

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  // parse arguments
  struct option const opts[] = {
    {"numEvents",    required_argument, nullptr,       'n'},
    {"outFile",      required_argument, nullptr,       'o'},
    {"beamEnergy",   required_argument, nullptr,       'e'},
    {"targetType",   required_argument, nullptr,       'T'},
    {"polType",      required_argument, nullptr,       'p'},
    {"beamSpin",     required_argument, nullptr,       'b'},
    {"targetSpin",   required_argument, nullptr,       't'},
    {"glgtMag",      required_argument, nullptr,       'm'},
    {"glgtArg",      required_argument, nullptr,       'a'},
    {"selectString", required_argument, nullptr,       'q'},
    {"config",       required_argument, nullptr,       'c'},
    {"seed",         required_argument, nullptr,       's'},
    {"verbose",      no_argument,       &verbose_mode, 1},
    {"help",         no_argument,       nullptr,       'h'},
    {nullptr,        0,                 nullptr,       0}
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
      case 'q':
                {
                  std::istringstream token_stream(optarg);
                  std::string token;
                  int i=0;
                  while(getline(token_stream, token, ',') && i < 2)
                    string_selection[i++] = std::stoi(token);
                  if(i != 2)
                    return Error("value of option '--selectString' does not have 2 arguments");
                  break;
                }
      case 'c': config_file = std::string(optarg); break;
      case 's': seed = std::stoi(optarg); break;
      case 'h':
        Usage();
        return EXIT_SYNTAX;
      case '?':
        return EXIT_ERROR;
    }
  }

  enable_string_selection = ! (string_selection[0] == 0 && string_selection[1] == 0);

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
  Verbose(fmt::format("{:>30} = ({})===({})  [{}]", "selectString", string_selection[0], string_selection[1], enable_string_selection ? "enabled" : "disabled"));
  Verbose(fmt::format("{:>30} = {}", "seed", seed));
  Verbose(fmt::format("{:>30} = {}", "config", config_file));
  Verbose(fmt::format("{:=^82}", ""));

  // get path to configuration file
  // - must be installed in `STRINGSPINNER_ETC`
  // - take only `filename()` from user specified argument, to prevent them from using `../` to
  //   leave the `STRINGSPINNER_ETC` directory
  auto config_file_path = std::string(STRINGSPINNER_ETC) + "/" + std::filesystem::path{config_file}.filename().string();
  Verbose(fmt::format("config file path: {}", config_file_path));
  Verbose(fmt::format("{:=^82}", ""));

  // set target PDG and mass
  int target_pdg;
  double target_mass;
  if(target_type == "proton") {
    target_pdg = 2212;
    target_mass = 0.93827;
  }
  else if(target_type == "neutron") {
    target_pdg = 2112;
    target_mass = 0.93957;
  }
  else {
    return Error(fmt::format("unknown '--targetType' value {:?}", target_type));
  }

  // parse polarization type and spins -> set `spin_vec`, the spin vector for beam and target
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



  // start pythia with stringspinner hooks
  Pythia8::Pythia pythia;
  Pythia8::Event& EV = pythia.event;
  auto fhooks = std::make_shared<Pythia8::SimpleStringSpinner>();
  fhooks->plugInto(pythia);

  // configure pythia
  pythia.readFile(config_file_path);
  //// beam and target types
  pythia.readString(fmt::format("Beams:idA = 11"));
  pythia.readString(fmt::format("Beams:idB = {}", target_pdg));
  pythia.readString(fmt::format("Beams:eA = {}", beam_energy));
  pythia.readString(fmt::format("Beams:eB = {}", target_mass));
  //// seed
  pythia.readString("Random:setSeed = on");
  pythia.readString(fmt::format("Random:seed = {}", seed));
  //// beam polarization
  if(obj_is_polarized[objBeam]) {
    for(auto quark : std::vector<std::string>{"u", "d", "s", "ubar", "dbar", "sbar"})
      pythia.readString(fmt::format("StringSpinner:{}Polarisation = {}", quark, fmt::join(spin_vec[objBeam],",")));
  }
  //// target polarization
  if(obj_is_polarized[objTarget])
    pythia.readString(fmt::format("StringSpinner:targetPolarisation = {}", fmt::join(spin_vec[objTarget],",")));
  //// stringspinner free parameters
  pythia.readString(fmt::format("StringSpinner:GLGT = {}", glgt_arg));
  pythia.readString(fmt::format("StringSpinner:thetaLT = {}", glgt_mag));

  return 0;
}
