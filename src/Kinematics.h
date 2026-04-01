#pragma once

#include <fmt/os.h>
#include <fmt/ranges.h>
#include <Pythia8/Event.h>

namespace css {

  using evnum_t = unsigned long;

  // ==================================================================================

  /// Inclusive kinematics object
  struct InclusiveKin {

    Pythia8::Particle lep; /// scattered lepton
    Pythia8::Vec4 vec_q; /// virtual photon 4-momentum
    Pythia8::Vec4 vec_target; /// target 4-momentum
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

    /// @brief generate a function to calculate dihadron z; this is separate from `CalculateKinematics` to enable its independent use
    /// @return a function to calculate dihadron z
    /// @param inc_kin the inclusive kinematics
    static std::function<double(Pythia8::Particle const&, Pythia8::Particle const&)> GetZfunction(InclusiveKin const& inc_kin)
    {
      return [&inc_kin] (Pythia8::Particle const& parA, Pythia8::Particle const& parB) {
        return (inc_kin.vec_target * (parA.p()+parB.p())) / (inc_kin.vec_target * inc_kin.vec_q); // P.Ph / P.q
      };
    }

    /// @brief calculate all dihadron kinematics
    /// @param inc_kin the inclusive kinematics
    /// @param hadA hadron A
    /// @param hadB hadron B
    void CalculateKinematics(
        InclusiveKin const& inc_kin,
        Pythia8::Particle const& hadA,
        Pythia8::Particle const& hadB
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
      z = GetZfunction(inc_kin)(hadA, hadB);
      Mh = vec_Ph.mCalc();
      MX = (inc_kin.vec_target + inc_kin.vec_q - vec_Ph).mCalc();
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
