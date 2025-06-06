// Kinematics.h
#ifndef KINEMATICS_H
#define KINEMATICS_H

#include <TLorentzVector.h>

class Kinematics {
public:
    double Q2(double beamE, double Eprime, double cosTh) const;
    double nu(double beamE, double Eprime) const;
    double y(double beamE, double Eprime) const;
    double x(double Q2, TLorentzVector q, TLorentzVector p) const;
    double W(double Q2, double Mp, double nu) const;
    double z(TLorentzVector, TLorentzVector, TLorentzVector) const;
    double xF(TLorentzVector, TLorentzVector, TLorentzVector, double) const;
    double com_th(TLorentzVector, TLorentzVector) const;
    double phi_R(TLorentzVector,TLorentzVector,TLorentzVector,TLorentzVector) const;
    double phi_h(TLorentzVector,TLorentzVector,TLorentzVector,TLorentzVector) const;
};

#endif // KINEMATICS_H