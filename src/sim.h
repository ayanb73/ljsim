#ifndef SIM_H
#define SIM_H

#include <array>
#include <string>
#include <tuple>
#include <vector>

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
  int id;
  float x = 0.0;
  float y = 0.0;
  float z = 0.0;
  float vx = 0.0;
  float vy = 0.0;
  float vz = 0.0;
  float Fx = 0.0;
  float Fy = 0.0;
  float Fz = 0.0;
};

float sqdist_atoms(Atom& a, Atom& b);

std::tuple<float, std::array<float, 3>> lj_pot_force(
    Atom& a, Atom& b);  // LJ 6-12 potential

class System {
 public:
  Args opt;
  std::vector<Atom> atoms;

  System(Args& o);

  void init_pos();
  void init_vel();
  float energy_forces();
};

#endif