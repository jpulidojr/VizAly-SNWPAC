# VizAly-LossyWave: A lossy, wavelet-based compressor for scientific simulation data

LossyWave is a cubic B-Spline wavelet compressor for regular grid datasets (1D,2D,3D). The goal for this method is to provide a low cost, low overhead compression method for smart data reduction via a series of parameters achieving targetted levels of lossyness. Code is made in two flavors, C++ (multi-platform) and Matlab (Windows).

## Requirements

### C++
C++11 compatible compiler
CMake v3.6 Minimum

GSL 1.16 (with cubic B-Spline patches (included))

LZ4 1.8.2

### Matlab
Matlab R2012a+
- Wavelet Toolbox
- Image Processing Toolbox

7-zip (ZStandard) v15.xx or greater (Windows)

## Usage
```
mkdir build
cd build
cmake ..

make -j
cd tests
./lossywave_test
```
## Contributors

Lead Developer: Jesus Pulido

Supervisor: Bernd Hamann

## Copyright and license
LANS has asserted copyright on the software package C17078, entitled Framework for Analysis and Visualization of Simulation Data.
