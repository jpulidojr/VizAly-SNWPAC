# VizAly-SNWPAC: **S**pli**N**e **W**avelet **P**acking and **C**ompression 

**SN**o**WPAC** is a cubic B-Spline wavelet compressor for regular grid datasets (1D,2D,3D). Specializing in lossy compression for scientific data, the goal for this method is to provide a low cost, low overhead compression method for smart data reduction via a series of parameters achieving targetted levels of lossyness.

## Requirements
C++11 compatible compiler

CMake v3.6 Minimum

LZ4 1.8.2

## Usage
```
mkdir build
cd build
cmake ..

make -j
cd tests
./snwpac_test
```
## Build Status
[![Build Status](https://travis-ci.org/jpulidojr/VizAly-SNWPAC.svg?branch=master)](https://travis-ci.org/jpulidojr/VizAly-SNWPAC)

## Contributors

Lead Developer: Jesus Pulido

Supervisor: Bernd Hamann

# Copyright and License
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. The U.S. Government has rights to use, reproduce, and distribute this software. NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE. If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
