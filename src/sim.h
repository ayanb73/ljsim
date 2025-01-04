#ifndef SIM_H
#define SIM_H

#include <array>
#include <string>
#include <tuple>
#include <vector>

struct Args {
  int num_particles = 0;
  int num_steps = 0;
  float box_dimension = 0.0;     // nm
  float sigma = 0.33713;          // nm
  float epsilon = 0.94639;  // kJ / mol
  float temperature = 100.0;     // Kelvin
  float mass = 39.948;           // amu
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

class System {
 public:
  Args opt;
  std::vector<Atom> atoms;
  float potential_energy;
  float kinetic_energy;
  float temperature;
  float half_box;

  float sigma_p12;
  float sigma_p6;

  System(Args& o);
  void init_pos();
  void init_vel();
  
  
  void update_pos(float dt);
  void update_vel(float half_dt);
  std::tuple<float, std::array<float, 3>> lj_pot_force(
    Atom& a, Atom& b);  // LJ 6-12 potential
  void compute_potential_energy_and_forces();
  void compute_kinetic_energy();
  void compute_temperature();
  void zero_com_velocity();

  void step_forward(float dt, float half_dt);
};

#endif