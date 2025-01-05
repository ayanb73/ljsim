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
  half_box = opt.box_dimension / 2;
  sigma_p6 = pow(opt.sigma, 6);
}

std::tuple<double, std::array<double, 3>> System::lj_pot_force(Atom& a, Atom& b) {
  // returns the force vector {Fx, Fy, Fz} acting on atom a
  // by newton's third law, atom b has the opposite force

  // calculate squared distance checking for minimum image
  double dx = a.x - b.x;
  double dy = a.y - b.y;
  double dz = a.z - b.z;
  if (dx > this->half_box) {dx -= opt.box_dimension;} else if (dx <= -this->half_box) {dx += this->opt.box_dimension;}
  if (dy > this->half_box) {dy -= opt.box_dimension;} else if (dy <= -this->half_box) {dy += this->opt.box_dimension;}
  if (dz > this->half_box) {dz -= opt.box_dimension;} else if (dz <= -this->half_box) {dz += this->opt.box_dimension;}
  
  
  double sqd = dx*dx + dy*dy + dz*dz;
  double r6 = pow(sqd, 3);
  double c6 = this->sigma_p6 / r6;
  double potE = 4 * this->opt.epsilon * (c6 - 1)*(c6);
  double force_prefactor = 48 * this->opt.epsilon * c6 * (c6 - 0.5) / sqd;
  std::array<double, 3> force_vector = {force_prefactor * dx, force_prefactor * dy, force_prefactor * dz};

  return std::make_tuple(potE, force_vector);
}


void System::init_pos() {
  double min_d = 0;
  int fenceposts = round(cbrt(this->opt.num_particles) + 2);
  int limit = fenceposts - 2;
  double spacer = this->opt.box_dimension / (double)fenceposts;
  
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

    a.x = min_d + (i+1)*spacer;
    a.y = min_d + (j+1)*spacer;
    a.z = min_d + (k+1)*spacer;
  }

}

void System::init_vel() {
  // generating approx maxwell-boltzmann velocities
  double var_mwb = sqrt(this->opt.temperature * 8.3144621e-3 / this->opt.mass);
  for (Atom& a : this->atoms) {
    for (int i = 0; i < 3; ++i) {
      double random_sum = 0;
      for (int n = 0; n < 12; ++n) {
        random_sum += static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
      }
      random_sum -= 6.0;
      random_sum = random_sum * var_mwb;
      // so dumb
      if (i == 0) {a.vx = random_sum;}
      if (i == 1) {a.vy = random_sum;}
      if (i == 2) {a.vz = random_sum;}
    }

  }

  // remove com motion
  this->zero_com_velocity();

  // compute effective temperature
  this->compute_kinetic_energy();
  double target_E = 0.5 * 3 * (this->opt.num_particles - 1) * 8.3144621e-3 * this->opt.temperature;
  double correction_factor = sqrt(target_E / this->kinetic_energy);
  
  // rescale to target temperature
  for (Atom& a : this->atoms) {
    a.vx = a.vx * correction_factor;
    a.vy = a.vy * correction_factor;
    a.vz = a.vz * correction_factor;
  }

  // compute corrected temperature
  this->compute_kinetic_energy();
  this->compute_temperature();
  // remove com motion
  this->zero_com_velocity();

}

void System::compute_potential_energy_and_forces() {
  double potE = 0.0;
  for (Atom& a : this->atoms) {
    a.Fx = 0.0;
    a.Fy = 0.0;
    a.Fz = 0.0;
  }
  // double loop
  for (int i =0; i < this->opt.num_particles - 1; ++i) {
    Atom& a = this->atoms.at(i);
    for (int j = i + 1; j < this->opt.num_particles; ++j) {
      Atom& b = this->atoms.at(j);
      // calculate pairwise energy and force
      std::tuple<double, std::array<double, 3>> pair_energy_force = this->lj_pot_force(a, b);
      // add to total Energy
      potE += std::get<0>(pair_energy_force);
      // update atom A forces
      a.Fx += std::get<1>(pair_energy_force).at(0);
      a.Fy += std::get<1>(pair_energy_force).at(1);
      a.Fz += std::get<1>(pair_energy_force).at(2);
      // update atom B forces
      b.Fx -= std::get<1>(pair_energy_force).at(0);
      b.Fy -= std::get<1>(pair_energy_force).at(1);
      b.Fz -= std::get<1>(pair_energy_force).at(2);

    }
  }
  this->potential_energy = potE;
}

void System::update_pos(double dt) {
  for (Atom& a : this->atoms) {
    double new_x = a.x + a.vx * dt;
    while (new_x < 0) {
      new_x = new_x + this->opt.box_dimension;
    }
    while (new_x >= this->opt.box_dimension) {
      new_x = new_x - this->opt.box_dimension;
    }
    a.x = new_x;
    double new_y = a.y + a.vy * dt;
    while (new_y < 0) {
      new_y = new_y + this->opt.box_dimension;
    }
    while (new_y >= this->opt.box_dimension) {
      new_y = new_y - this->opt.box_dimension;
    }
    a.y = new_y;
    double new_z = a.z + a.vz * dt;
    while (new_z < 0) {
      new_z = new_z + this->opt.box_dimension;
    }
    while (new_z >= this->opt.box_dimension) {
      new_z = new_z - this->opt.box_dimension;
    }
    a.z = new_z;
  }
}

void System::update_vel(double half_dt) {
  for (Atom& a : this->atoms) {
    a.vx += (a.Fx / this->opt.mass) * half_dt;
    a.vy += (a.Fy / this->opt.mass) * half_dt;
    a.vz += (a.Fz / this->opt.mass) * half_dt;
  }
}

void System::zero_com_velocity() {
  double com_vx = 0.0;
  double com_vy = 0.0;
  double com_vz = 0.0;
  for (Atom&a : this->atoms) {
    com_vx += a.vx;
    com_vy += a.vy;
    com_vz += a.vz;
  }
  com_vx = com_vx / this->opt.num_particles;
  com_vy = com_vy / this->opt.num_particles;
  com_vz = com_vz / this->opt.num_particles;
  
  for (Atom&a : this->atoms) {
    a.vx -= com_vx;
    a.vy -= com_vy;
    a.vz -= com_vz;
  }
  
}


void System::compute_kinetic_energy() {
  double totE = 0;
  for (Atom& a : this->atoms) {
    totE += 0.5 * (pow(a.vx, 2) + pow(a.vy, 2) + pow(a.vz, 2)) * this->opt.mass;
  }
  this->kinetic_energy = totE;
}

void System::compute_temperature() {
  this->temperature = (2 * this->kinetic_energy) / (3 * (this->opt.num_particles - 1) * 8.3144621e-3);
}


void System::step_forward(double dt) {
  // update the velocities to v(t + dt/2)
  this->update_vel(0.5*dt);

  // update the positions to r(t + dt)
  this->update_pos(dt);
  
  // update forces and accelerations to be F(t + dt) / m = a(t + dt)
  this->compute_potential_energy_and_forces();

  // update the velocities again to v(t + dt)
  this->update_vel(0.5*dt);

  // remove com motion
  this->zero_com_velocity();

  // compute kinetic energy and temperature
  this->compute_kinetic_energy();
  this->compute_temperature();
}

void System::report() {
  printf("KE=%f PE=%f E=%f T=%f\n", this->kinetic_energy, this->potential_energy, this->kinetic_energy + this->potential_energy, this->temperature);
}