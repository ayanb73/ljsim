## ljsim 

`ljsim` is a pico-tool for running molecular dynamics (MD) simulations of homogenenous Lennard-Jones fluids in the NVE ensemble. You can do everything `ljsim` can do, and vastly more, in any of the common MD software packages; I wrote this mostly to learn things.

If you want to use `ljsim`, you can clone this repo, and compile the binary with `g++ src/*.cc -o ./ljsim -std=c++11`. Any `gcc` or `clang` that supports `c++11` should be fine. 

#### Roadmap

- [ ] Support for arbitrary mixtures
- [ ] Cleaner way to configure
- [ ] Data output file format
- [ ] Thermostat for NVT simulations
- [ ] Better testing harness
- [ ] Performance optimization
- [ ] Pivot into a library