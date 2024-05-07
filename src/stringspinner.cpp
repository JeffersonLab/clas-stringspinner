#include <getopt.h>
#include <fmt/format.h>
#include <Pythia8/Pythia.h>
#include <stringspinner/StringSpinner.h>

static unsigned long num_events = 10000;
static std::string out_file     = "out.lund";
static int verbose              = 0;

void Usage()
{
  fmt::print("USAGE: stringspinner [OPTIONS]...\n\n");
  fmt::print("--numEvents NUM_EVENTS        number of events\n");
  fmt::print("                              default: {}\n", num_events);
  fmt::print("--outFile OUTPUT_FILE         output file name\n");
  fmt::print("                              default: {}\n", out_file);
  fmt::print("--verbose                     verbose printout\n");
  fmt::print("--help                        print this usage guide\n");
}

int main(int argc, char** argv)
{
  struct option const opts[] = {
    {"numEvents", required_argument, nullptr,  'n'},
    {"outFile",   required_argument, nullptr,  'o'},
    {"verbose",   no_argument,       &verbose, 1},
    {"help",      no_argument,       nullptr,  'h'},
    {nullptr,     0,                 nullptr,  0}
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
    fmt::print("{:=^82}\n", "");
  }

  return 0;
}
