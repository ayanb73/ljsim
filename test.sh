#!/bin/bash

g++ src/ljsim.cc -o ljsim -std=c++11

./ljsim --num_particles=1000 --sigma=3.4 --epsilon=120 --temperature=310 --outdir=out
