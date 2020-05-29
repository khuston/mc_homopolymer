## Purpose
Sample configurations of a homopolymer, optionally at a sticky surface.

## Build and Run Tests
```bash
git clone https://github.com/khuston/mc_homopolymer.git
cd mc_homopolymer
python -m venv venv
venv/bin/pip install -e .[dev]
pytest
```

## Usage
This repo is a work in progress. The usage at the moment has just one option and no way to get data out.

```Python
import PyPolymers as Poly

polymer = Poly.MonteCarloChain(N=50, epsilon=1, sigma=1)
polymer.run(10000, tethered=False)
```

The target usage is:

```Python
import PyPolymers as Poly

N = number_of_particles = 50

particles = Poly.Particles.initialize_chain(N)

potentials = Poly.Potentials.add_harmonic_bonds(particles)
                            .add_harmonic_angles(particles)
                            .add_lennard_jones(particles, particles, epsilon=, sigma=)
                            .add_wall_10_4_3(particles, epsilon=, sigma=)

step1 = Poly.Steps.add_random_translations_to_all(particles)

step2 = Poly.Steps.add_crankshaft(particles, weight=1)
                  .add_pivot(particles, weight=1)
                  .add_pivot_z(particles, weight=1)
                  .add_slither(particles, weight=1)
                  .add_end_rotate(particles, weight=1)

steps = (step1, step2)

sampler = Poly.Samplers.UnbiasedSampler(steps, potentials)

num_steps = 1e5

position_logger = Poly.Loggers.PositionLogger(particles, )

sampler.add_logger()
sampler.run(num_steps)
```

## About
This library is implemented in C++ with Python bindings exposed through the [pybind11](https://github.com/pybind/pybind11) library.

## Requirements
- CMake 3.17
- Python 3.6
- C++14 compiler