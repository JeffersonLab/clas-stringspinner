#include <fmt/os.h>
#include <fmt/ranges.h>
#include <optional>
#include <vector>
#include "Tools.h"

namespace clas {

  using evnum_t = unsigned long;

  /// Lund event header variables
  struct LundHeader {

    /// number of particles in the event
    int num_particles;
    /// target mass
    double target_mass;
    /// target atomic number
    int target_atomic_num;
    /// target spin
    double target_spin;
    /// beam spin (applies to the 1st `LundParticle`, actually)
    double beam_spin;
    /// beam PDG
    int beam_type;
    /// beam energy
    double beam_energy;
    /// target nucleon PDG
    int nucleon_pdg;
    /// pythia process code
    int process_id;
    /// event weight
    double event_weight;
    /// user values
    std::vector<int> user_values{};

    /// @brief stream to output file
    /// @param output the output file stream
    /// @param precision number of decimal places for floating point numbers
    void Stream(fmt::ostream& output, int const& precision) const {
      std::size_t const user_values_max_size = 90;
      std::optional<std::string> user_values_str = std::nullopt;
      if(!user_values.empty()) {
        if(user_values.size() <= user_values_max_size) {
          user_values_str = fmt::format(" {:d}", fmt::join(user_values, " "));
        }
        else {
          clas::EventError("LundHeader::user_values is too big, with size = {}; truncating this list", user_values.size());
          decltype(user_values) user_values_trun(user_values_max_size);
          std::copy(user_values.begin(), user_values.begin() + user_values_max_size, user_values_trun.begin());
          user_values_str = fmt::format(" {:d}", fmt::join(user_values_trun, " "));
        }
      }
      output.print("{num_particles:d} {target_mass:.{p}f} {target_atomic_num:d} {target_spin:.{p}f} {beam_spin:.{p}f} {beam_type:d} {beam_energy:.{p}f} {nucleon_pdg:d} {process_id:d} {event_weight:.{p}f}{user_values}\n",
          fmt::arg("num_particles", num_particles),
          fmt::arg("target_mass", target_mass),
          fmt::arg("target_atomic_num", target_atomic_num),
          fmt::arg("target_spin", target_spin),
          fmt::arg("beam_spin", beam_spin),
          fmt::arg("beam_type", beam_type),
          fmt::arg("beam_energy", beam_energy),
          fmt::arg("nucleon_pdg", nucleon_pdg),
          fmt::arg("process_id", process_id),
          fmt::arg("event_weight", event_weight),
          fmt::arg("user_values", user_values_str.value_or("")),
          fmt::arg("p", precision)
          );
    }
  };

  // ==================================================================================

  /// Lund particle variables
  struct LundParticle {

    /// particle index
    int index;
    /// particle lifetime
    double lifetime;
    /// particle type: 1 is propagated in Geant, 0 is not
    int type;
    /// particle PDG
    int pdg;
    /// first mother
    int mother1;
    /// first daughter
    int daughter1;
    /// particle momentum x-component
    double px;
    /// particle momentum y-component
    double py;
    /// particle momentum z-component
    double pz;
    /// particle energy
    double energy;
    /// particle mass
    double mass;
    /// particle vertex x-component
    double vx;
    /// particle vertex y-component
    double vy;
    /// particle vertex z-component
    double vz;

    /// @brief stream to output file
    /// @param output the output file stream
    /// @param precision number of decimal places for floating point numbers
    void Stream(fmt::ostream& output, int const& precision) const {
      output.print("{index:d} {lifetime:.{p}f} {type:d} {pdg:d} {mother1:d} {daughter1:d} {px:.{p}f} {py:.{p}f} {pz:.{p}f} {energy:.{p}f} {mass:.{p}f} {vx:.{p}f} {vy:.{p}f} {vz:.{p}f}\n",
          fmt::arg("index", index),
          fmt::arg("lifetime", lifetime),
          fmt::arg("type", type),
          fmt::arg("pdg", pdg),
          fmt::arg("mother1", mother1),
          fmt::arg("daughter1", daughter1),
          fmt::arg("px", px),
          fmt::arg("py", py),
          fmt::arg("pz", pz),
          fmt::arg("energy", energy),
          fmt::arg("mass", mass),
          fmt::arg("vx", vx),
          fmt::arg("vy", vy),
          fmt::arg("vz", vz),
          fmt::arg("p", precision)
          );
    }
  };

  // ==================================================================================

  /// Inclusive kinematics object
  struct InclusiveKin {

    /// event number
    evnum_t evnum;
    /// x
    double x;
    /// Q2
    double Q2;
    /// W
    double W;
    /// y
    double y;

    /// @brief print TTree branch-description
    static void Header(fmt::ostream& output) {
      output.print("{}\n", fmt::join(
            std::vector<std::string>{
            "evnum/I",
            "x/D",
            "Q2/D",
            "W/D",
            "y/D"
            }, ":"));
    }

    /// @brief stream to output file
    /// @param output the output file stream
    /// @param precision number of decimal places for floating point numbers
    void Stream(fmt::ostream& output, int const& precision) const {
      output.print("{evnum:d} {x:.{p}f} {Q2:.{p}f} {W:.{p}f} {y:.{p}f}\n",
          fmt::arg("evnum", evnum),
          fmt::arg("x", x),
          fmt::arg("Q2", Q2),
          fmt::arg("W", W),
          fmt::arg("y", y),
          fmt::arg("p", precision)
          );
    }

  };

  // ==================================================================================

  /// Single-hadron object
  struct SingleHadronKin {

    /// event number
    evnum_t evnum;
    /// index of the hadron
    int idx;
    /// PDG of the hadron
    int pdg;
    /// px
    double px;
    /// py
    double py;
    /// pz
    double pz;
    /// theta
    double theta;

    /// @brief print TTree branch-description
    static void Header(fmt::ostream& output) {
      output.print("{}\n", fmt::join(
            std::vector<std::string>{
            "evnum/I",
            "idx/I",
            "pdg/I",
            "px/D",
            "py/D",
            "pz/D",
            "theta/D"
            }, ":"));
    }

    /// @brief stream to output file
    /// @param output the output file stream
    /// @param precision number of decimal places for floating point numbers
    void Stream(fmt::ostream& output, int const& precision) const {
      output.print("{evnum:d} {idx:d} {pdg:d} {px:.{p}f} {py:.{p}f} {pz:.{p}f} {theta:.{p}f}\n",
          fmt::arg("evnum", evnum),
          fmt::arg("idx", idx),
          fmt::arg("pdg", pdg),
          fmt::arg("px", px),
          fmt::arg("py", py),
          fmt::arg("pz", pz),
          fmt::arg("theta", theta),
          fmt::arg("p", precision)
          );
    }

  };

  // ==================================================================================

  /// Dihadron object
  struct DihadronKin {

    /// event number
    evnum_t evnum;
    /// index of hadron A
    int idxA;
    /// index of hadron B
    int idxB;
    /// PDG of hadron A
    int pdgA;
    /// PDG of hadron B
    int pdgB;
    /// z of the pair
    double z;
    /// invariant mass
    double Mh;

    /// @brief print TTree branch-description
    static void Header(fmt::ostream& output) {
      output.print("{}\n", fmt::join(
            std::vector<std::string>{
            "evnum/I",
            "idxA/I",
            "idxB/I",
            "pdgA/I",
            "pdgB/I",
            "z/D",
            "Mh/D"
            }, ":"));
    }

    /// @brief stream to output file
    /// @param output the output file stream
    /// @param precision number of decimal places for floating point numbers
    void Stream(fmt::ostream& output, int const& precision) const {
      output.print("{evnum:d} {idxA:d} {idxB:d} {pdgA:d} {pdgB:d} {z:.{p}f} {Mh:.{p}f}\n",
          fmt::arg("evnum", evnum),
          fmt::arg("idxA", idxA),
          fmt::arg("idxB", idxB),
          fmt::arg("pdgA", pdgA),
          fmt::arg("pdgB", pdgB),
          fmt::arg("z", z),
          fmt::arg("Mh", Mh),
          fmt::arg("p", precision)
          );
    }

  };

}
