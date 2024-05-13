#include <getopt.h>
#include <fmt/format.h>
#include <Pythia8/Pythia.h>
#include <stringspinner/StringSpinner.h>

const int EXIT_ERROR = 1;
const int EXIT_SYNTAX = 2;

enum pol_type_enum { polU, polL, polT, nPol };
const std::string pol_type_name[nPol] = { "unpolarized", "longitudinal", "transverse" };
enum spin_enum { spinP, spinM };
enum obj_enum { objBeam, objTarget, nObj };
const std::string obj_name[nObj] = { "beam", "target" };

static unsigned long num_events    = 10000;
static std::string out_file        = "out.lund";
static int verbose_mode            = 0;
static std::string pol_type        = "UU";
static std::string spin_type[nObj] = {"", ""};

//////////////////////////////////////////////////////////////////////////////////

void Usage()
{
  fmt::print("USAGE: stringspinner [OPTIONS]...\n\n");
  fmt::print("--numEvents NUM_EVENTS           number of events\n");
  fmt::print("                                 default: {}\n\n", num_events);
  fmt::print("--outFile OUTPUT_FILE            output file name\n");
  fmt::print("                                 default: {:?}\n\n", out_file);
  fmt::print("--polType POLARIZATION_TYPE      beam and target polarization types\n");
  fmt::print("                                 - two characters: beam and target\n");
  fmt::print("                                 - types: 'U' = unpolarized\n");
  fmt::print("                                          'L' = longitudinally polarized\n");
  fmt::print("                                          'T' = transversely polarized\n");
  fmt::print("                                 default: {:?}\n\n", pol_type);
  fmt::print("--beamSpin BEAM_SPIN             the spin of the beam leptons\n");
  fmt::print("                                 - if longitudinally polarized ('L'):\n");
  fmt::print("                                   'p' = spin along +z axis\n");
  fmt::print("                                   'm' = spin along -z axis\n");
  fmt::print("                                 - if transversely polarized ('T'):\n");
  fmt::print("                                   'p' = spin along +y axis\n");
  fmt::print("                                   'm' = spin along -y axis\n");
  fmt::print("                                 - if unpolarized ('U'): no effect\n\n");
  fmt::print("--targetSpin TARGET_SPIN         the spin of the target nucleons\n");
  fmt::print("                                 - same usage as --beamSpin, applied to target\n\n");
  fmt::print("--verbose                        verbose printout\n\n");
  fmt::print("--help                           print this usage guide\n\n");
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
    {"numEvents",  required_argument, nullptr,       'n'},
    {"outFile",    required_argument, nullptr,       'o'},
    {"polType",    required_argument, nullptr,       'p'},
    {"beamSpin",   required_argument, nullptr,       'b'},
    {"targetSpin", required_argument, nullptr,       't'},
    {"verbose",    no_argument,       &verbose_mode, 1},
    {"help",       no_argument,       nullptr,       'h'},
    {nullptr,      0,                 nullptr,       0}
  };

  if(argc <= 1) {
    Usage();
    return EXIT_SYNTAX;
  };

  char opt;
  while((opt = getopt_long(argc, argv, "", opts, nullptr)) != -1) {
    switch(opt) {
      case 'n':
        num_events = std::stol(optarg);
        break;
      case 'o':
        out_file = std::string(optarg);
        break;
      case 'p':
        pol_type = std::string(optarg);
        break;
      case 'b':
        spin_type[objBeam] = std::string(optarg);
        break;
      case 't':
        spin_type[objTarget] = std::string(optarg);
        break;
      case 'h':
        Usage();
        return EXIT_SYNTAX;
      case '?':
        return EXIT_ERROR;
    }
  }

  Verbose(fmt::format("{:=^82}", " Arguments "));
  Verbose(fmt::format("{:>30} = {}", "numEvents", num_events));
  Verbose(fmt::format("{:>30} = {:?}", "outFile", out_file));
  Verbose(fmt::format("{:>30} = {:?}", "polType", pol_type));
  Verbose(fmt::format("{:>30} = {:?}", "beamSpin", spin_type[objBeam]));
  Verbose(fmt::format("{:>30} = {:?}", "targetSpin", spin_type[objTarget]));
  Verbose(fmt::format("{:=^82}", ""));

  // parse polarization type and spins
  std::array<double,3> spin_vec[nObj] = { {0, 0, 0}, {0, 0, 0} };
  if(pol_type.length() != 2)
    return Error(fmt::format("option '--polType' value {:?} is not 2 characters", pol_type));

  for(int i = 0; i < nObj; i++) {

    // parse the polarization type
    int pol_type_num;
    switch(std::toupper(pol_type.c_str()[i])) {
      case 'U': pol_type_num = polU; break;
      case 'L': pol_type_num = polL; break;
      case 'T': pol_type_num = polT; break;
      default:
        return Error(fmt::format("option '--polType' has unknown {} polarization type {:?}", obj_name[i], pol_type.c_str()[i]));
    }
    Verbose(fmt::format("{:>30} = {}", fmt::format("{} polarization type", obj_name[i]), pol_type_name[pol_type_num]));

    // if the polarization type is unpolarized, parse the spin type
    int spin_num = -1;
    if(pol_type_num != polU) {
      if(spin_type[i].empty())
        return Error(fmt::format("option '--{}Spin' must be set when {} polarization is {}", obj_name[i], obj_name[i], pol_type_name[pol_type_num]));
      if(spin_type[i].length() > 1)
        return Error(fmt::format("option '--{}Spin' value {:?} is not 1 character", obj_name[i], spin_type[i]));
      switch(std::tolower(spin_type[i].c_str()[0])) {
        case 'p': spin_num = spinP; break;
        case 'm': spin_num = spinM; break;
        default:
          return Error(fmt::format("option '--{}Spin' has unknown value {:?}", obj_name[i], spin_type[i]));
      }
    }

    // define the spin vector
    std::string spin_name = "0";
    if(pol_type_num == polU)
      spin_vec[i] = {0, 0, 0};
    else if(pol_type_num == polL) {
      switch(spin_num) {
        case spinP:
          spin_name = "+";
          spin_vec[i] = { 0, 0, i == objBeam ? -1.0 : 1.0 }; // beam has opposite sign, since quark momentum reversed after hard scattering
          break;
        case spinM:
          spin_name = "-";
          spin_vec[i] = { 0, 0, i == objBeam ? 1.0 : -1.0 };
          break;
      }
    }
    else if(pol_type_num == polT) {
      switch(spin_num) {
        case spinP:
          spin_name = "up";
          spin_vec[i] = { 0, i == objBeam ? -1.0 : 1.0, 0 }; // beam has opposite sign, since quark momentum reversed after hard scattering
          break;
        case spinM:
          spin_name = "down";
          spin_vec[i] = { 0, i == objBeam ? 1.0 : -1.0, 0 };
          break;
      }
    }
    Verbose(fmt::format("{:>30} = {}", fmt::format("{} spin", obj_name[i]), spin_name));
    Verbose(fmt::format("{:>30} = ({})", fmt::format("{} spin vector", i == objBeam ? "quark" : obj_name[i]), fmt::join(spin_vec[i], ", ")));
  }

  return 0;
}
