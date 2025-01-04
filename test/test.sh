#!/bin/bash
rm -rf ./ljsim
g++ src/*.cc -o ./ljsim -std=c++11
./ljsim --num_particles=2 --box_dimension=1.13526 --num_steps=100
