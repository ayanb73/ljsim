#include <cmath>
#include <array>
#include "sim.h"


float dist_atoms(Atom& a, Atom& b) {
  float sqdist = pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2);
  float dist = sqrt(sqdist);
  return dist;
}

float sqdist_atoms(Atom& a, Atom& b) {
  float sqdist = pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2);
  return sqdist;
}

std::array<float, 2> dist_sqdist_atoms(Atom& a, Atom& b) {
  float sqdist = pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2);
  float dist = sqrt(sqdist);
  std::array<float, 2> res = {{dist, sqdist}};
  return res;
}

std::array<float, 2> lj_pot_force(Atom& a, Atom& b) {
  // if the atoms are the same, zero
  if (&a == &b) {
    return {{0.0, 0.0}};
  }
  std::array<float, 2> d_sqd = dist_sqdist_atoms(a, b);
  float c6 = 1 / pow(d_sqd[1], 3);
  float c12 = pow(c6, 2);
  float c13 = c12 / d_sqd[0];
  float c7 = c6 / d_sqd[0];
  float potE = 4 * (c12 - c6);
  float force = 24 * c7 - 48 * c13;

  return {{potE, force}};
}
