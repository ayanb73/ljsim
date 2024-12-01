#include <cmath>
#include <array>
#include "sim.h"
#include <tuple>
#include <iostream>

float sqdist_atoms(Atom& a, Atom& b) {
  float sqdist = pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2);
  return sqdist;
}

std::tuple<float, std::array<float, 3>> lj_pot_force(Atom& a, Atom& b) {
  // returns the force vector {Fx, Fy, Fz} acting on atom a
  // by newton's third law, atom b has the opposite force

  // if the atoms are the same, return zero
  if (&a == &b) {
    return std::make_tuple(0.0, std::array<float, 3>{0.0, 0.0, 0.0});
  }
  float sqd = sqdist_atoms(a, b);
  float r6 = 1 / pow(sqd, 3);
  float r12 = pow(r6, 2);
  float r14 = r12 / sqd;
  float r8 = r6 / sqd;
  float potE = 4 * (r12 - r6);
  float force_prefactor = 24 * r8 - 48 * r14;
  std::array<float, 3> force_vector = {force_prefactor * (a.x - b.x), force_prefactor * (a.y - b.y), force_prefactor * (a.z - b.z)};

  return std::make_tuple(potE, force_vector);
}

System::System(Args& o) {
  this->opt = o;
  for (int i = 0; i < o.num_particles; ++i) {
    Atom new_atom;
    new_atom.id = i;
    atoms.emplace_back(new_atom);
  }
}

void System::init_pos() {
  float min_d = - this->opt.box_dimension / 2;
  int fenceposts = round(cbrt(this->opt.num_particles) + 2);
  int limit = fenceposts - 2;
  float spacer = this->opt.box_dimension / (float)fenceposts;
  
  int k_min = pow(limit, 2);

  for (Atom& a : this->atoms) {
    // turn the id into a base (limit) number
    int i = 0;
    int j = 0;
    int k = 0;
    if (a.id >= k_min) {
      k = a.id / k_min;
    } 
    int rem = a.id % k_min;
    if (rem >= limit) {
      j = rem / limit;
    }
    i = rem % limit;

    a.x = min_d + i*spacer;
    a.y = min_d + j*spacer;
    a.z = min_d + k*spacer;
  }

}