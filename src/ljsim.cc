#include <chrono>
#include <cmath>
#include <iostream>
#include <string>

#include "sim.h"

int main(int argc, char* argv[]) {
  char delimiter = '=';
  Args cli_args;
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    std::string begin = arg.substr(0, 2);
    if (begin != "--") {
      std::cout << "Expect -- for each command line arg." << std::endl;
      return 1;
    }
    size_t delim_index = arg.find(delimiter);
    if (delim_index == arg.npos) {
      std::cout << "Use --option=value for each command line arg." << std::endl;
      return 1;
    } else {
      std::string arg_key = arg.substr(2, delim_index - 2);
      std::string arg_value = arg.substr(delim_index + 1);
      if (arg_key == "num_particles") {
        cli_args.num_particles = std::stoi(arg_value, nullptr);
      } else if (arg_key == "num_steps") {
        cli_args.num_steps = std::stoi(arg_value, nullptr);
      } else if (arg_key == "box_dimension") {
        cli_args.box_dimension = std::stof(arg_value, nullptr);
      } else if (arg_key == "sigma") {
        cli_args.sigma = std::stof(arg_value, nullptr);
      } else if (arg_key == "epsilon") {
        cli_args.epsilon = std::stof(arg_value, nullptr);
      } else if (arg_key == "temperature") {
        cli_args.temperature = std::stof(arg_value, nullptr);
      } else if (arg_key == "mass") {
        cli_args.mass = std::stof(arg_value, nullptr);
      } else if (arg_key == "outdir") {
        cli_args.outdir = arg_value;
      } else {
        std::cout << "Unrecognized arg: " << arg_key << std::endl;
      }
    }
  }

  // if number of particles not provided
  if (cli_args.num_particles <= 0) {
    std::cout << "Must provide --num_particles=N option where N is an integer "
                 "greater than 0."
              << std::endl;
    return 1;
  }
  // if number of steps not provided
  if (cli_args.num_steps <= 0) {
    std::cout << "Must provide --num_steps=N option where N is an integer "
                 "greater than 0."
              << std::endl;
    return 1;
  }

  // if box dimension not provided
  if (cli_args.box_dimension == 0.0) {
    cli_args.box_dimension = std::cbrt(cli_args.num_particles * cli_args.sigma);
  }

  System simul_system(cli_args, 0);
  simul_system.init_pos();
  simul_system.init_vel();
  simul_system.compute_potential_energy_and_forces();
  simul_system.compute_kinetic_energy();

  double dt = 0.002;  // 2 fs
  double h_dt = 0.5 * dt;
  simul_system.report();
  for (int i = 1; i <= simul_system.opt.num_steps; ++i) {
    simul_system.step_forward(dt);
    simul_system.report();
  }
  return 0;
};
