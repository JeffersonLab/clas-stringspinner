#include <fmt/format.h>
#include <functional>
#include <string>
#include <sstream>

namespace clas {

  /// exit code for general error
  int const EXIT_ERROR  = 1;
  /// exit code for syntax error
  int const EXIT_SYNTAX = 2;

  /// filter for verbose printouts
  static bool enable_verbose_mode = false;

  /// @param msg error printout of this message
  /// @returns error exit code, so caller can use `return Error(...)`
  template <typename... Args>
  inline int Error(fmt::format_string<Args...> fmt_str, Args&&... fmt_args)
  {
    fmt::println(stderr, "[ERROR] {}", fmt::format(fmt_str, std::forward<Args>(fmt_args)...));
    return EXIT_ERROR;
  }

  /// number of event errors
  static int event_error_count = 0;
  /// @param msg print an error for this event; if there are too many, stop printing errors
  template <typename... Args>
  inline void EventError(fmt::format_string<Args...> fmt_str, Args&&... fmt_args) {
    const int max_event_errors = 100;
    if(event_error_count <= max_event_errors) {
      Error(fmt_str, std::forward<Args>(fmt_args)...);
      if(event_error_count == max_event_errors)
        Error("More than {} event errors... suppressing the rest...", event_error_count);
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
