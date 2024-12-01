#ifndef SIM_H
#define SIM_H

#include <array>

struct Args {
  int num_particles = 0;
  int num_steps = 0;
  float box_dimension = 0.0;     // Angstroms
  float sigma = 3.3713;          // Angstroms
  float epsilon = 0.2261926386;  // kcal / mol
  float temperature = 120.0;     // Kelvin
  float mass = 39.948;           // g / mol
  std::string outdir = "ljsim_out";
};

struct Atom {
  float x;
  float y;
  float z;
  float vx;
  float vy;
  float vz;
};

float dist_atoms(Atom& a, Atom& b); 
float sqdist_atoms(Atom& a, Atom& b); 
std::array<float, 2> dist_sqdist_atoms(Atom& a, Atom& b); 

std::array<float, 2> lj_pot_force(Atom& a, Atom& b); // LJ 6-12 potential

template<int N>
class System {
 public:
  Atom atoms[N];
  Args opt;

  void init_pos();
  void time_step();
  void run_steps(int);
};

#endif