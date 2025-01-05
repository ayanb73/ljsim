#ifndef SIM_H
#define SIM_H

#include <array>
#include <string>
#include <tuple>
#include <vector>

struct Args {
  int num_particles = 0;
  int num_steps = 0;
  double box_dimension = 0.0;     // nm
  double sigma = 0.33713;          // nm
  double epsilon = 0.94639;  // kJ / mol
  double temperature = 100.0;     // Kelvin
  double mass = 39.948;           // amu
  std::string outdir = "ljsim_out";
};

struct Atom {
  int id;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double vz = 0.0;
  double Fx = 0.0;
  double Fy = 0.0;
  double Fz = 0.0;
};

class System {
 public:
  Args opt;
  std::vector<Atom> atoms;
  double potential_energy;
  double kinetic_energy;
  double temperature;
  double half_box;

  double sigma_p12;
  double sigma_p6;

  System(Args& o);
  void init_pos();
  void init_vel();
  
  
  void update_pos(double dt);
  void update_vel(double half_dt);
  std::tuple<double, std::array<double, 3>> lj_pot_force(
    Atom& a, Atom& b);  // LJ 6-12 potential
  void compute_potential_energy_and_forces();
  void compute_kinetic_energy();
  void compute_temperature();
  void zero_com_velocity();

  void step_forward(double dt);
  void report();
};

#endif