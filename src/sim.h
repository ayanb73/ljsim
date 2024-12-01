#ifndef SIM_H
#define SIM_H

struct Args {
    int num_particles = 0;
    float box_dimension = 0.0;  // Angstroms
    float sigma = 3.4; // Angstroms
    float epsilon = 120.0; // Kelvin * k_B
    float temperature = 300.0; // Kelvin

    std::string outdir = "ljsim_output";
};

#endif