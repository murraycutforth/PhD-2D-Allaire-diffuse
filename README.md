# Allaire-diffuse

This is a research code for two-material flow governed by the Euler equations, 
implementing the five-equation model of Allaire. Two methods are implemented, both
using the HLLC approximate Riemann solver:
1. Godunov's method, first order accurate
2. MUSCL-Hancock method (with Minbee slope limiting), second order accurate


## Setup instructions

1. (If necessary) Install Eigen and create a symlink to the Eigen headers on the default include path: 
```
sudo apt install eigen3-dev
sudo ln -s /usr/include/eigen3/Eigen /usr/local/include/Eigen
```
2. Clone these three repositories from the same base directory:
```
git clone https://github.com/murraycutforth/PhD-Common.git
git clone https://github.com/murraycutforth/PhD-2D-HCLsolver-framework.git
git clone https://github.com/murraycutforth/PhD-2D-Allaire-diffuse.git
```
3. Compile:
```
cd PhD-2D-Allaire-diffuse
make
```
4. Run code by editing the settingsfile.txt and then running ```./PhD-2D-Allaire-diffuse settingsfile.txt```


## Settings file options

TODO (you'll have to read the code for now- sorry!)
