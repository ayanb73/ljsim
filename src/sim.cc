#include <cmath>
#include <array>
#include "sim.h"
#include <tuple>

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
