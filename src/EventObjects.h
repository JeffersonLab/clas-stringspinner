#include <fmt/os.h>
#include <fmt/ranges.h>
#include <optional>
#include <vector>
#include <Pythia8/Event.h>
#include "Tools.h"

namespace clas {

  using evnum_t = unsigned long;

  /// Lund event header variables
  struct LundHeader {

    int num_particles; /// number of particles in the event
    double target_mass; /// target mass
    int target_atomic_num; /// target atomic number
    double target_spin; /// target spin
    double beam_spin; /// beam spin (applies to the 1st `LundParticle`, actually)
    int beam_type; /// beam PDG
    double beam_energy; /// beam energy
    int nucleon_pdg; /// target nucleon PDG
    int process_id; /// pythia process code
    double event_weight; /// event weight
    std::vector<int> user_values{}; /// user values

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

    int index; /// particle index
    double lifetime; /// particle lifetime
    int type; /// particle type: 1 is propagated in Geant, 0 is not
    int pdg; /// particle PDG
    int mother1; /// first mother
    int daughter1; /// first daughter
    double px; /// particle momentum x-component
    double py; /// particle momentum y-component
    double pz; /// particle momentum z-component
    double energy; /// particle energy
    double mass; /// particle mass
    double vx; /// particle vertex x-component
    double vy; /// particle vertex y-component
    double vz; /// particle vertex z-component

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

    Pythia8::Particle lep; /// scattered lepton
    evnum_t evnum; /// event number
    double x; /// x
    double Q2; /// Q2
    double W; /// W
    double y; /// y

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

    evnum_t evnum; /// event number
    int idx; /// index of the hadron
    int pdg; /// PDG of the hadron
    double px; /// px
    double py; /// py
    double pz; /// pz
    double theta; /// theta

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

    evnum_t evnum; /// event number
    int idxA; /// index of hadron A
    int idxB; /// index of hadron B
    int pdgA; /// PDG of hadron A
    int pdgB; /// PDG of hadron B
    double x{0}; /// DIS x
    double Q2{0}; /// DIS Q2
    double W{0}; /// DIS W
    double y{0}; /// DIS y
    double pxLep{0}; /// px of the scattered lepton
    double pyLep{0}; /// py of the scattered lepton
    double pzLep{0}; /// pz of the scattered lepton
    double pxA{0}; /// px of hadron A
    double pyA{0}; /// py of hadron A
    double pzA{0}; /// pz of hadron A
    double pxB{0}; /// px of hadron B
    double pyB{0}; /// py of hadron B
    double pzB{0}; /// pz of hadron B
    double z{0}; /// dihadron z
    double Mh{0}; /// dihadron invariant mass
    double MX{0}; /// dihadron missing mass
    double cosTheta{0}; /// dihadron partial wave cos(theta) (not scattering angle theta)

    static std::function<double(Pythia8::Particle const&, Pythia8::Particle const&)> GetZ(
        Pythia8::Vec4 const& vec_q,
        Pythia8::Vec4 const& vec_target
        )
    {
      return [&vec_target, &vec_q] (Pythia8::Particle const& parA, Pythia8::Particle const& parB) {
        return (vec_target * (parA.p()+parB.p())) / (vec_target * vec_q); // P.Ph / P.q
      };
    }

    void CalculateKinematics(
        InclusiveKin const& inc_kin,
        Pythia8::Particle const& hadA,
        Pythia8::Particle const& hadB,
        Pythia8::Vec4 const& vec_q,
        Pythia8::Vec4 const& vec_target
        )
    {
      auto const vec_Ph = hadA.p() + hadB.p();
      x  = inc_kin.x;
      Q2 = inc_kin.Q2;
      W  = inc_kin.W;
      y  = inc_kin.y;
      pxLep = inc_kin.lep.px();
      pyLep = inc_kin.lep.py();
      pzLep = inc_kin.lep.pz();
      pxA = hadA.px();
      pyA = hadA.py();
      pzA = hadA.pz();
      pxB = hadB.px();
      pyB = hadB.py();
      pzB = hadB.pz();
      z = GetZ(vec_q, vec_target)(hadA, hadB);
      Mh = vec_Ph.mCalc();
      MX = (vec_target + vec_q - vec_Ph).mCalc();
      // calculate cos(theta)
      auto vec_hadA__dih = hadA.p();
      vec_hadA__dih.rotbst(Pythia8::toCMframe(vec_Ph));  // boost to dihadron rest frame
      cosTheta = Pythia8::costheta(vec_hadA__dih, vec_Ph);
    }

    /// @brief print TTree branch-description
    static void Header(fmt::ostream& output) {
      output.print("{}\n", fmt::join(
            std::vector<std::string>{
            "evnum/I",
            "idxA/I",
            "idxB/I",
            "pdgA/I",
            "pdgB/I",
            "x/D",
            "Q2/D",
            "W/D",
            "y/D",
            "pxLep/D", "pyLep/D", "pzLep/D",
            "pxA/D", "pyA/D", "pzA/D",
            "pxB/D", "pyB/D", "pzB/D",
            "z/D",
            "Mh/D",
            "MX/D",
            "cosTheta/D"
            }, ":"));
    }

    /// @brief stream to output file
    /// @param output the output file stream
    /// @param precision number of decimal places for floating point numbers
    void Stream(fmt::ostream& output, int const& precision) const {
      output.print("{evnum:d} {ints:d} {doubles:.{p}f}\n",
          fmt::arg("evnum", evnum),
          fmt::arg("ints", fmt::join(
              std::vector<int>{
              idxA, idxB,
              pdgA, pdgB
              }, " ")),
          fmt::arg("doubles", fmt::join(
              std::vector<double>{
              x, Q2, W, y,
              pxLep, pyLep, pzLep,
              pxA, pyA, pzA,
              pxB, pyB, pzB,
              z,
              Mh,
              MX,
              cosTheta
              }, " ")),
          fmt::arg("p", precision)
          );
    }

  };

}
