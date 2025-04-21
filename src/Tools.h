#include <fmt/format.h>
#include <functional>
#include <string>
#include <sstream>

namespace clas {

  /// exit code for general error
  int const EXIT_ERROR  = 1;
  /// exit code for syntax error
  int const EXIT_SYNTAX = 2;

  /// if true, printouts are verbose
  static bool enable_verbose_mode = false;
  /// @param msg verbose printout of this message
  inline void Verbose(std::string msg)
  {
    if(enable_verbose_mode)
      fmt::println(msg);
  }

  /// @param msg error printout of this message
  /// @returns error exit code, so caller can use `return Error(...)`
  inline int Error(std::string msg)
  {
    fmt::println(stderr, "[ERROR] " + msg);
    return EXIT_ERROR;
  }

  /// number of event errors
  static int event_error_count = 0;
  /// @param msg print an error for this event; if there are too many, stop printing errors
  inline void EventError(std::string msg) {
    const int max_event_errors = 100;
    if(event_error_count <= max_event_errors) {
      Error(msg);
      if(event_error_count == max_event_errors)
        Error(fmt::format("More than {} event errors... suppressing the rest...", event_error_count));
      event_error_count++;
    }
  }

  /// @brief tokenize a string delimited by commas
  /// @param str the input string
  /// @param func what to do with each token
  inline void Tokenize(char const* str, std::function<void(std::string,int)> func)
  {
    std::istringstream token_stream(str);
    std::string token;
    char const delim = ',';
    int i = 0;
    while(getline(token_stream, token, delim))
      func(token, i++);
  }

}
