#pragma once

#include <fmt/format.h>
#include <functional>
#include <sstream>

/// namespace for `clas-stringspinner`
namespace string_spinner {

  /// event number type
  using evnum_t = unsigned long;

  /// filter for verbose printouts
  static bool enable_verbose_mode = false;

  /// @param msg error printout of this message
  template <typename... Args>
  inline void Error(fmt::format_string<Args...> fmt_str, Args&&... fmt_args)
  {
    fmt::println(stderr, "[ERROR] {}", fmt::format(fmt_str, std::forward<Args>(fmt_args)...));
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

  /// @brief expand `~` to the user's home directory
  /// @param path the file path
  /// @returns the file path with `~` expanded to `$HOME`
  inline std::string ExpandTilde(std::string const& path)
  {
    if(path.empty() || path[0] != '~')
      return path;
    if(char const* home_dir = std::getenv("HOME"); home_dir)
      return std::string(home_dir) + path.substr(1);
    throw std::runtime_error("cannot expand `~` since $HOME is not set");
  }


}
