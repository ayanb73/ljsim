#!/bin/bash

g++ src/*.cc -o ljsim -std=c++11

./ljsim --num_particles=256 --num_steps=100
