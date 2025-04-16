#include <fmt/format.h>
#include <functional>
#include <string>
#include <sstream>

int const EXIT_ERROR  = 1;
int const EXIT_SYNTAX = 2;

static bool enable_verbose_mode = false;
inline void Verbose(std::string msg)
{
  if(enable_verbose_mode)
    fmt::print(msg + "\n");
}

inline int Error(std::string msg)
{
  fmt::print(stderr, "[ERROR] " + msg + "\n");
  return EXIT_ERROR;
}

static int event_error_count = 0;
inline void EventError(std::string msg) {
  const int max_event_errors = 100;
  if(event_error_count <= max_event_errors) {
    Error(msg);
    if(event_error_count == max_event_errors)
      Error(fmt::format("More than {} event errors... suppressing the rest...", event_error_count));
    event_error_count++;
  }
}

inline void Tokenize(char const* str, std::function<void(std::string,int)> func)
{
  std::istringstream token_stream(str);
  std::string token;
  char const delim = ',';
  int i = 0;
  while(getline(token_stream, token, delim))
    func(token, i);
}
