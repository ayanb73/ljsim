#include <cmath>
#include <array>
#include <tuple>
#include <iostream>

#include "sim.h"

/* Defining the System class methods */

System::System(Args& o) {
  opt = o;
  for (int i = 0; i < o.num_particles; ++i) {
    Atom new_atom;
    new_atom.id = i;
    atoms.emplace_back(new_atom);
  }
  n_atoms = atoms.size();
  half_box = opt.box_dimension / 2;
}

std::tuple<float, std::array<float, 3>> System::lj_pot_force(Atom& a, Atom& b) {
  // returns the force vector {Fx, Fy, Fz} acting on atom a
  // by newton's third law, atom b has the opposite force

  // if the atoms are the same, return zero
  if (&a == &b) {
    return std::make_tuple(0.0, std::array<float, 3>{0.0, 0.0, 0.0});
  }
  // calculate squared distance
  // checking for minimum image
  float dx = a.x - b.x;
  float dy = a.y - b.y;
  float dz = a.z - b.z;
  if (dx > half_box) {dx -= half_box;}
  if (dy > half_box) {dy -= half_box;}
  if (dz > half_box) {dz -= half_box;}
  float sqd = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
  
  // LJ but in reduced units of sigma = 1, epsilon = 1, mass = 1
  float r6 = 1 / pow(sqd, 3);
  float r12 = pow(r6, 2);
  float r14 = r12 / sqd;
  float r8 = r6 / sqd;
  float potE = 4  * (r12 - r6);
  float force_prefactor = (24 * r8 - 48 * r14);
  std::array<float, 3> force_vector = {force_prefactor * (a.x - b.x), force_prefactor * (a.y - b.y), force_prefactor * (a.z - b.z)};

  return std::make_tuple(potE, force_vector);
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

void System::init_vel() {
  // generating stuff randomly from the GROMACS algorithms
  float var_mwb = sqrt(this->opt.temperature * 0.00198720413 / this->opt.epsilon);
  for (Atom& a : this->atoms) {
    for (int i = 0; i < 3; ++i) {
      float random_sum = 0;
      for (int n = 0; n < 12; ++n) {
        random_sum += static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
      }
      random_sum -= 6.0;
      random_sum = random_sum * var_mwb;
      // so dumb
      if (i == 0) {a.vx = random_sum;}
      if (i == 1) {a.vy = random_sum;}
      if (i == 2) {a.vz = random_sum;}
    }

  }
  // fake time step w/ zero forces to remove com motion
  this->update_vel(0);

}

void System::compute_potential_energy_and_forces() {
  float potE = 0.0;
  // stupid double loop
  for (int i = 0; i < this->n_atoms; ++i) {
    for (int j = 0; j < this->n_atoms; ++j) {
      // update the velocity by half a time step

      // calculate pairwise energy and force
      std::tuple<float, std::array<float, 3>> pair_energy_force = this->lj_pot_force(this->atoms.at(i), this->atoms.at(j));
      // add to total Energy
      potE += std::get<0>(pair_energy_force);
      // update atom A forces
      this->atoms.at(i).Fx += std::get<1>(pair_energy_force).at(0);
      this->atoms.at(i).Fy += std::get<1>(pair_energy_force).at(1);
      this->atoms.at(i).Fz += std::get<1>(pair_energy_force).at(2);
      // update atom B forces
      this->atoms.at(j).Fx -= std::get<1>(pair_energy_force).at(0);
      this->atoms.at(j).Fy -= std::get<1>(pair_energy_force).at(1);
      this->atoms.at(j).Fz -= std::get<1>(pair_energy_force).at(2);

    }
  }
  this->potential_energy = potE;
}

void System::update_pos(float dt) {
  for (int i = 0; i < this->n_atoms; ++i) {
    this->atoms.at(i).x += this->atoms.at(i).vx * dt;
    this->atoms.at(i).y += this->atoms.at(i).vy * dt;
    this->atoms.at(i).z += this->atoms.at(i).vz * dt;
  }
}

void System::update_vel(float half_dt) {
  // our mass = 1 so we can
  // update velocities directly with 
  // forces and remove com velocity directly

  float com_vx = 0.0;
  float com_vy = 0.0;
  float com_vz = 0.0;
  for (int i = 0; i < this->n_atoms; ++i) {
    this->atoms.at(i).vx += this->atoms.at(i).Fx * half_dt;
    this->atoms.at(i).vy += this->atoms.at(i).Fy * half_dt;
    this->atoms.at(i).vz += this->atoms.at(i).Fz * half_dt;

    com_vx += this->atoms.at(i).vx;
    com_vy += this->atoms.at(i).vy;
    com_vz += this->atoms.at(i).vz;
  }

  float com_vx = com_vx / this->n_atoms;
  float com_vy = com_vy / this->n_atoms;
  float com_vz = com_vz / this->n_atoms;

  for (int i = 0; i < this->n_atoms; ++i) {
    this->atoms.at(i).vx -= com_vx;
    this->atoms.at(i).vy -= com_vy;
    this->atoms.at(i).vz -= com_vz;
  }
}


void System::compute_kinetic_energy() {
  float totE = 0;
  for (Atom& a : this->atoms) {
    totE += sqrt(pow(a.vx, 2) + pow(a.vy, 2) + pow(a.vz, 2));
  }
  totE = 0.5 * totE;
  this->kinetic_energy = totE;
}

void System::step_forward(float dt, float half_dt) {
  // update the velocities
  this->update_vel(half_dt);

  // update the positions
  this->update_pos(dt);
  
  // update forces (accelerations)
  this->compute_potential_energy_and_forces();

  // update the velocities again
  this->update_vel(half_dt);

  // compute kinetic energy
  this->compute_kinetic_energy();
}