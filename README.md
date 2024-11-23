# dhara_seq

Currently DHARA supports fully-compressible convection with ideal gas equation.

## Libraries required:
- Numpy
- CuPy (Version 12 recommended, `cupy-cuda12x`)
- h5py

For testing the code, we can switch off the saving fields and global .h5 files.

## Fully-compressible Convection Solver
- This is a sequential Object Oriented Python code which solves fully compressible convection in Cartesian box. 
- This is GPU enabled. 

### Steps to run the code:
1. Go to `para.py` file for changing the parameters. 
  - To make it run on GPU, put GPU for `device` and its rank. 
  - Make a output directory where you want to store output files and give its path to `output_dir`. 
  - You can change the time advance scheme from `Scheme`.
  - Put the grid parameters and control parameters.
2. To start the solver sequentially run `MainSeq.py` using Python `python3 MainSeq.py`.

### About the output files:
- The output field contains: rho, ux, uy, uz, theta and parameters. 
- theta is the perturbation from the adiabatic profile: T(x,y,z)-T_a(z).
- parameters contain: t, alpha, m, epsilon, Pr, Ra, beta in sequence.

### Notes:
- The solver needs the initial fields. I will recommend to leave them since they contain the equilribrium fields. 
