struct LundHeader {
  int    num_particles;
  double target_mass;
  int    target_atomic_num;
  double target_spin;
  double beam_spin;
  int    beam_type;
  double beam_energy;
  int    nucleon_pdg;
  int    process_id;
  double event_weight;
};

struct LundParticle {
  int    index;
  double lifetime;
  int    status;
  int    pdg;
  int    mother1;
  int    daughter1;
  double px;
  double py;
  double pz;
  double energy;
  double mass;
  double vx;
  double vy;
  double vz;
};
