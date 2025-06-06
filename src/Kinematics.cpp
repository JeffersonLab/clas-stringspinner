// Kinematics.C
#include "Kinematics.h"
#include <cmath>


double Kinematics::Q2(double beamE, double Eprime, double cosTh) const {
    return 2.0 * beamE * Eprime * (1.0 - cosTh);
}

double Kinematics::nu(double beamE, double Eprime) const {
    return beamE - Eprime;
}

double Kinematics::y(double beamE, double Eprime) const {
    double n = nu(beamE, Eprime);
    return n > 0.0 ? n / beamE : 0.0;
}

double Kinematics::x(double Q2, TLorentzVector q, TLorentzVector p) const {
    return Q2/(2*p*q);
}

double Kinematics::W(double Q2, double Mp, double nu) const {
    double val = Mp*Mp + 2.0*Mp*nu - Q2;
    return val > 0.0 ? std::sqrt(val) : 0.0;
}

double Kinematics::z(TLorentzVector init_target, TLorentzVector part, TLorentzVector q) const {
  return (init_target*part)/(init_target*q);
}

double Kinematics::com_th(TLorentzVector P1, TLorentzVector P2) const{
  TLorentzVector Ptotal = P1+P2;
  TVector3 comBOOST = Ptotal.BoostVector();
  Ptotal.Boost(-comBOOST);
  P1.Boost(-comBOOST);
  return P1.Angle(comBOOST);
}

double Kinematics::phi_R(TLorentzVector Q, TLorentzVector L, TLorentzVector p1, TLorentzVector p2) const {
    TLorentzVector ph = p1 + p2;
    TLorentzVector r = 0.5*(p1-p2);
    
    TVector3 q(Q.Px(), Q.Py(), Q.Pz());
    TVector3 l(L.Px(), L.Py(), L.Pz());
    TVector3 R(r.Px(), r.Py(), r.Pz());
    
    TVector3 Rperp;
    
    // -- HERMES 0803.2367 angle, but used Matevosyan et al 1707.04999
    //    to obtain R_perp vector
    TLorentzVector init_target;
    init_target.SetPxPyPzE(0,0,0,0.938272);
    double z1 = (init_target*p1)/(init_target*Q);
    double z2 = (init_target*p2)/(init_target*Q);
    TVector3 P1(p1.Px(), p1.Py(), p1.Pz());
    TVector3 P2(p2.Px(), p2.Py(), p2.Pz());
    TVector3 P1perp = P1-(q*P1)/(q*q)*q;
    TVector3 P2perp = P2-(q*P2)/(q*q)*q;
    Rperp = (z2*P1perp-z1*P2perp)*((1)/(z1+z2));

    
    
    TVector3 qcrossl = q.Cross(l);
    TVector3 qcrossRperp = q.Cross(Rperp);
    
    double factor1 = (qcrossl*Rperp)/abs(qcrossl*Rperp);
    double factor2 = (qcrossl*qcrossRperp)/qcrossl.Mag()/qcrossRperp.Mag();
    
    return factor1*acos(factor2);
}


double Kinematics::phi_h(TLorentzVector Q, TLorentzVector L, TLorentzVector p1, TLorentzVector p2) const {
  TLorentzVector ph = p1 + p2;
  TLorentzVector r = 0.5*(p1-p2);

  TVector3 q(Q.Px(), Q.Py(), Q.Pz());
  TVector3 l(L.Px(), L.Py(), L.Pz());
  TVector3 Ph(ph.Px(), ph.Py(), ph.Pz());

  TVector3 qcrossl = q.Cross(l);
  TVector3 qcrossPh = q.Cross(Ph);

  double factor1 = (qcrossl*Ph)/abs(qcrossl*Ph);
  double factor2 = (qcrossl*qcrossPh)/qcrossl.Mag()/qcrossPh.Mag();
    
  return factor1*acos(factor2);
}

double Kinematics::xF(TLorentzVector q, TLorentzVector p, TLorentzVector init_target, double W) const {
  TLorentzVector com = q+init_target;
  TVector3 comBOOST = com.BoostVector();
  TLorentzVector qq = q;
  TLorentzVector pp = p;
  qq.Boost(-comBOOST);
  pp.Boost(-comBOOST);
  return 2*(qq.Vect().Dot(pp.Vect()))/(qq.Vect().Mag()*W);
}