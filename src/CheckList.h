#include <fmt/format.h>
#include <fmt/ranges.h>
#include <Pythia8/Event.h>
#include "Tools.h"

namespace clas {

  /// @brief checklist for particle cuts
  class CheckList {

    public:

      /// @param opt_name option name from command line
      CheckList(std::string const& opt_name) :
        m_opt_name(opt_name),
        m_enabled(false),
        m_enable_min_max(false) {}

      /// @brief setup the checklist by parsing command line arguments
      /// @param opt_arg argument string from command line
      /// @param parse_min_max if false, do not parse min,max from `opt_arg`
      void Setup(char const* opt_arg, bool const parse_min_max=true)
      {
        // parse command line arguments
        if(parse_min_max) {
          Tokenize(opt_arg, [&](auto token, auto i) {
            switch(i) {
              case 0: m_min = std::stod(token); break;
              case 1: m_max = std::stod(token); break;
              default: m_pdg_list.push_back(std::stoi(token));
            }
          });
          if(m_pdg_list.empty())
            throw std::runtime_error(fmt::format("value of option '--{}' does not have at least 3 arguments", m_opt_name));
          if(m_min >= m_max)
            throw std::runtime_error(fmt::format("option '--{}' has MIN >= MAX", m_opt_name));
        }
        else {
          Tokenize(opt_arg, [&](auto token, auto i) {
            m_pdg_list.push_back(std::stoi(token));
          });
          if(m_pdg_list.empty())
            throw std::runtime_error(fmt::format("value of option '--{}' does not have at least 1 argument", m_opt_name));
        }
        // initialize checklist
        m_checklist.clear();
        for(auto const& pdg : m_pdg_list)
          m_checklist.emplace_back(pdg, false);
        // local settings
        m_enabled = true;
        m_enable_min_max = parse_min_max;
      }

      // @returns true if this `CheckList` is enabled
      bool const& Enabled() { return m_enabled; }

      // @returns a string with information about this checklist
      std::string const GetInfoString() {
        if(!m_enabled)
          return "disabled";
        std::string result = fmt::format("for all PDGs ({})", fmt::join(m_pdg_list, ", "));
        if(m_enable_min_max)
          result = fmt::format("in ({}, {}) {}", m_min, m_max, result);
        return result;
      }

      /// @brief check the checklist
      /// @param evt the pythia event
      /// @param get_val a function to get the value to be checked; not needed if min,max is not used
      /// @param must_be_final if true, particle must be "final"
      /// @returns true if all particles pass the cut; if this checklist is not enabled, returns true
      bool const Check(
          Pythia8::Event const& evt,
          std::function<double(Pythia8::Particle const&)> get_val = [](Pythia8::Particle const& par){ return 0; },
          bool const& must_be_final = true)
      {
        // return true if not enabled
        if(!m_enabled)
          return true;
        // reset the checklist
        for(auto& [pdg, found] : m_checklist)
          found = false;
        // loop over event particles, and check the checkboxes
        std::size_t num_found = 0;
        for(auto const& par : evt) { // loop over event particles
          if(!must_be_final || (must_be_final && par.isFinal())) { // particle must be final, if `must_be_final==true`
            for(auto& [pdg, found] : m_checklist) { // loop over checklist
              if(!found && pdg == par.id()) { // if we haven't found this one yet, and the checklist PDG == particle PDG
                if(m_enable_min_max) { // check the box if val is in range (m_min, m_max)
                  auto val = get_val(par);
                  found = m_min < val && val < m_max;
                }
                else { // no value comparison, just check the box
                  found = true;
                }
                if(found) {
                  num_found++;
                  break; // break the loop over checklist
                }
              }
            }
          }
          // if all checkboxes are checked, stop looping over event particles and return true
          if(num_found == m_checklist.size()) {
            return true;
          }
        }
        // not all checkboxes were checked, return false
        return false;
      }

    private:
      /// the name of the checklist should match the CLI option name
      std::string m_opt_name;
      /// whether this checklist is enabled
      bool m_enabled;
      /// whether the min and max values are used for this checklist
      bool m_enable_min_max;
      /// the checklist: PDG -> checkbox
      std::vector<std::pair<int, bool>> m_checklist;
      /// list of PDGs in this checklist
      std::vector<int> m_pdg_list;
      /// cut minimum
      double m_min;
      /// cut maximum
      double m_max;
  };
}
