#include <getopt.h>
#include <fmt/format.h>
#include <Pythia8/Pythia.h>

#if defined (__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-compare"
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wmissing-braces"
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wmissing-braces"
#endif
#include <stringspinner/StringSpinner.h>
#if defined (__clang__)
#pragma clang diagnostic pop
#else
#pragma GCC diagnostic pop
#endif

static unsigned long num_events = 10000;
static std::string out_file     = "out.lund";
static int verbose              = 0;
static std::string pol_type     = "LU";
static std::string spin_beam    = "p";
static std::string spin_target  = "p";

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
  fmt::print("--spinBeam BEAM_SPIN             the spin of the beam leptons\n");
  fmt::print("                                 - if longitudinally polarized ('L'):\n");
  fmt::print("                                   'p' = spin along +z axis\n");
  fmt::print("                                   'm' = spin along -z axis\n");
  fmt::print("                                 - if transversely polarized ('T'):\n");
  fmt::print("                                   'p' = spin along +y axis\n");
  fmt::print("                                   'm' = spin along -y axis\n");
  fmt::print("                                 - if unpolarized ('U'): no effect\n");
  fmt::print("                                 default: {:?}\n\n", spin_beam);
  fmt::print("--spinTarget TARGET_SPIN         the spin of the target nucleons\n");
  fmt::print("                                 - same as --spinBeam, applied to target\n");
  fmt::print("                                 default: {:?}\n\n", spin_target);
  fmt::print("--verbose                        verbose printout\n\n");
  fmt::print("--help                           print this usage guide\n\n");
}

int main(int argc, char** argv)
{
  struct option const opts[] = {
    {"numEvents",  required_argument, nullptr,  'n'},
    {"outFile",    required_argument, nullptr,  'o'},
    {"polType",    required_argument, nullptr,  'p'},
    {"spinBeam",   required_argument, nullptr,  'b'},
    {"spinTarget", required_argument, nullptr,  't'},
    {"verbose",    no_argument,       &verbose, 1},
    {"help",       no_argument,       nullptr,  'h'},
    {nullptr,      0,                 nullptr,  0}
  };

  if(argc <= 1) {
    Usage();
    return 2;
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
        spin_beam = std::string(optarg);
        break;
      case 't':
        spin_target = std::string(optarg);
        break;
      case 'h':
        Usage();
        return 2;
      case '?':
        return 1;
    }
  }

  if(verbose==1) {
    fmt::print("{:=^82}\n", " Arguments ");
    fmt::print("{:>30} = {}\n", "numEvents", num_events);
    fmt::print("{:>30} = {:?}\n", "outFile", out_file);
    fmt::print("{:>30} = {:?}\n", "polType", pol_type);
    fmt::print("{:>30} = {:?}\n", "spinBeam", spin_beam);
    fmt::print("{:>30} = {:?}\n", "spinTarget", spin_target);
    fmt::print("{:=^82}\n", "");
  }

  return 0;
}
