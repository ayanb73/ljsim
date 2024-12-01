#!/bin/bash

g++ src/*.cc -o ljsim -std=c++11
./ljsim --num_particles=256 --box_dimension=21.6 --num_steps=10
