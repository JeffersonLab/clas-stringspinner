#include <fmt/format.h>
#include <fmt/ranges.h>
#include <Pythia8/Event.h>
#include "Tools.h"
#include "EventObjects.h"

namespace clas {

  /// @brief checklist for particle cuts
  class CheckList {

    public:

      /// modes for `CheckList` operations
      enum Mode {
        /// min and max values are not used
        kNoCuts,
        /// cuts applied to list of single hadrons
        k1hCuts,
        /// cuts applied to pair of hadrons
        k2hCuts
      };

      /// @param opt_name option name from command line
      /// @param the mode for this checklist
      CheckList(std::string const& opt_name, Mode mode) :
        m_opt_name(opt_name),
        m_mode(mode),
        m_enabled(false) {}

      /// @brief setup the checklist by parsing command line arguments
      /// @param opt_arg argument string from command line
      void Setup(char const* opt_arg)
      {
        // parse command line arguments
        switch(m_mode) {
          case kNoCuts:
            {
              Tokenize(opt_arg, [&](auto token, auto i) {
                m_pdg_list.push_back(std::stoi(token));
              });
              if(m_pdg_list.empty())
                throw std::runtime_error(fmt::format("value of option '--{}' does not have at least 1 argument", m_opt_name));
              break;
            }
          case k1hCuts:
          case k2hCuts:
            {
              Tokenize(opt_arg, [&](auto token, auto i) {
                switch(i) {
                  case 0: m_min = std::stod(token); break;
                  case 1: m_max = std::stod(token); break;
                  default: m_pdg_list.push_back(std::stoi(token));
                }
              });
              if(m_mode == k2hCuts && m_pdg_list.size() != 2)
                throw std::runtime_error(fmt::format("value of option '--{}' must have exactly 2 PDG values, which specify a dihadron", m_opt_name));
              if(m_pdg_list.empty())
                throw std::runtime_error(fmt::format("value of option '--{}' does not have at least 3 arguments", m_opt_name));
              if(m_min >= m_max)
                throw std::runtime_error(fmt::format("option '--{}' has MIN >= MAX", m_opt_name));
              break;
            }
        }
        m_enabled = true;
      }

      // @returns true if this `CheckList` is enabled
      bool const& Enabled() const { return m_enabled; }

      // @returns a string with information about this checklist
      std::string const GetInfoString() const {
        if(!m_enabled)
          return "disabled";
        switch(m_mode) {
          case kNoCuts:
            return fmt::format("for all PDGs ({})", fmt::join(m_pdg_list, ", "));
          case k1hCuts:
            return fmt::format("in [{:.3g}, {:.3g}], for all PDGs ({})", m_min, m_max, fmt::join(m_pdg_list, ", "));
          case k2hCuts:
            return fmt::format("in [{:.3g}, {:.3g}], for ({}) dihadrons", m_min, m_max, fmt::join(m_pdg_list, ", "));
        }
        throw std::runtime_error("GetInfoString failed");
      }

      /// @brief check the checklist for single particles
      /// @param evt the pythia event
      /// @param get_val a function to get the value to be checked, given a `Pythia8::Particle`; not used if `m_mode == kNoCuts`
      /// @param must_be_final if true, particle must be "final"
      /// @returns true if all particles pass the cut; if this checklist is not enabled, returns true
      bool const Check(
          Pythia8::Event const& evt,
          std::function<double(Pythia8::Particle const&)> const& get_val = [](Pythia8::Particle const& par){ return 0; },
          bool const& must_be_final = true) const
      {
        // return immediately if the checklist is disabled
        if(!m_enabled)
          return true;

        // “family-mode”: only require >=2 matches, not a full set
        bool const is_family_mode = (m_opt_name.find("family") != std::string::npos);
        std::size_t const required_matches = is_family_mode ? 2 : m_pdg_list.size();

        if(enable_verbose_mode)
          fmt::println("CHECKLIST for {} (need {} match{})",
                       m_opt_name, required_matches, required_matches==1?"":"es");

        // initialise checklist
        std::vector<std::pair<int, bool>> check_list;
        for(auto const& pdg : m_pdg_list)
          check_list.emplace_back(pdg, false);

        double num_found = 0; // running count of satisfied PDGs

        // loop over all particles in the event
        for(auto const& par : evt) {
          if(!must_be_final || par.isFinal() || par.id()==111) { // allow Pi0's
            for(auto& [pdg, found] : check_list) {
              if(!found && pdg == par.id()) {

                switch(m_mode) {
                  case kNoCuts:
                    found = true;
                    if(enable_verbose_mode) fmt::println("  [x] {} at idx {}", pdg, par.index());
                    break;

                  case k1hCuts:
                    {
                      double const val = get_val(par);
                      found = (m_min <= val && val <= m_max);
                      if(enable_verbose_mode && found)
                        fmt::println("  [x] {} at idx {}, value = {}", pdg, par.index(), val);
                    }
                    break;

                  case k2hCuts:
                    throw std::runtime_error("called single-particle `Check` but mode is dihadron");
                }

                if(found) {
                  // pi0's need two photons, so only increment found by 0.5 for photons
                  if(pdg==22) {num_found+=0.5;}
                  else {num_found++;}

                  if(num_found >= required_matches)
                    return true; // early exit once criterion is satisfied
                  break;         // stop scanning the checklist for this particle
                }
              }
            }
          }
        }

        // not enough PDGs satisfied
        return false;
      }

      /// @brief check if there's a dihadron which satisfies the cut
      /// @param evt the pythia event
      /// @param dih_kin a list of `DihadronKin` objects for this event
      /// @param get_val a function to get the value to be checked, given 2 `Pythia8::Particle`s
      /// @returns true if a dihadron which satisfies the cut is found; if this checklist is not enabled, returns true
      bool const Check(
          Pythia8::Event const& evt,
          std::vector<DihadronKin> const& dih_kin,
          std::function<double(Pythia8::Particle const&, Pythia8::Particle const&)> const& get_val
          ) const
      {
        // return true if not enabled
        if(!m_enabled)
          return true;
        // must be dihadron mode
        if(m_mode != k2hCuts)
          throw std::runtime_error("called dihadron `Check` but mode is not dihadron");
        // loop over dihadrons
        if(enable_verbose_mode) fmt::println("CHECKLIST for {}", m_opt_name);
        for(auto const& dih : dih_kin) {
          // if this is a dihadron the user wants
          if( (dih.pdgA == m_pdg_list.at(0) && dih.pdgB == m_pdg_list.at(1)) ||
              (dih.pdgA == m_pdg_list.at(1) && dih.pdgB == m_pdg_list.at(0)) )
          {
            // check the cuts
            auto const& parA = evt.at(dih.idxA);
            auto const& parB = evt.at(dih.idxB);
            auto val = get_val(parA, parB);
            if(m_min <= val && val <= m_max) {
              if(enable_verbose_mode) fmt::println("  [x] ({}, {}) at idxs ({}, {}), value = {}", parA.id(), parB.id(), dih.idxA, dih.idxB, val);
              return true; // satisfactory dihadron
            }
          }
        }
        return false; // no satisfactory dihadron
      }

    private:
      /// the name of the checklist should match the CLI option name
      std::string const m_opt_name;
      /// the mode for this checklist
      Mode const m_mode;
      /// whether this checklist is enabled
      bool m_enabled;
      /// list of PDGs in this checklist
      std::vector<int> m_pdg_list;
      /// cut minimum
      double m_min;
      /// cut maximum
      double m_max;
  };
}
