#!/bin/bash
rm -rf ./ljsim
g++ src/*.cc -o ./ljsim -std=c++11
./ljsim --num_particles=256 --box_dimension=2.25 --num_steps=1000
