#include <fmt/os.h>
#include <fmt/ranges.h>
#include <optional>
#include <vector>
#include "Tools.h"

namespace clas {

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
    std::vector<double> user_values{};

    /// @brief stream to output file
    /// @param output the output file stream
    void Stream(fmt::ostream& output) const {
      std::size_t const user_values_max_size = 90;
      std::optional<std::string> user_values_str = std::nullopt;
      if(!user_values.empty()) {
        if(user_values.size() <= user_values_max_size) {
          user_values_str = fmt::format(" {:.5}", fmt::join(user_values, " "));
        }
        else {
          clas::EventError("LundHeader::user_values is too big, with size = {}; truncating this list", user_values.size());
          decltype(user_values) user_values_trun(user_values_max_size);
          std::copy(user_values.begin(), user_values.begin() + user_values_max_size, user_values_trun.begin());
          user_values_str = fmt::format(" {:.5}", fmt::join(user_values_trun, " "));
        }
      }
      output.print("{:} {:.5} {:} {:} {:} {:} {:.5} {:} {:} {:.5}{}\n",
          num_particles,
          target_mass,
          target_atomic_num,
          target_spin,
          beam_spin,
          beam_type,
          beam_energy,
          nucleon_pdg,
          process_id,
          event_weight,
          user_values_str.value_or("")
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
    void Stream(fmt::ostream& output) const {
      output.print("{:} {:.5} {:} {:} {:} {:} {:.5} {:.5} {:.5} {:.5} {:.5} {:.5} {:.5} {:.5}\n",
          index,
          lifetime,
          type,
          pdg,
          mother1,
          daughter1,
          px,
          py,
          pz,
          energy,
          mass,
          vx,
          vy,
          vz
          );
    }
  };

  // ==================================================================================

  /// Dihadron object
  struct DihadronKin {
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
  };

}
