# README.md

This repository contains the data produced in the study of eccentric planet disc interactions. The underlying theory and methodology is detailed fully the paper. The data is divided into the circular and eccentric result directories. Each of the files contains the results across a variety of eccentricity and disc parameters.

- $p$ = radial surface density exponents
- $T$ = radial temperature exponents
- $h$ = disc aspect ratios at the planetary semi-major axis
- $e$ = eccentricity of the planetary orbit

## `circular/ `
- `parameters.out`: stores the eccentricity values e and disc parameters for q, p and h
- `t_net_m1.out`: net torque including the m = 1 mode
- `t_net_m2.out`: net torque excluding the m = 1 mode
- `t_l_m1.out`: Lindblad torque including the m = 1 mode
- `t_l_m2.out`: Lindblad torque excluding the m = 1 mode
- `t_c_m1.out`: Corotation torque including the m = 1 mode
- `t_c_m2.out`: Corotation torque including the m = 1 mode

## `eccentric/` 
- `parameters.out`: stores the eccentricity values e and disc parameters for q, p and h
- `t_c.out`: corotation torque
- `t_l.out`: Lindblad torque
- `t_net.out`: Net torque
- `tau_a_inv.out`: semi-major axis damping rates
- `tau_e_inv.out`: eccentricity damping rates
- `tau_L_inv.out`: angular momentum damping rates

In the root directory there is a sample Jupyter notebook script `read_data.ipynb` which shows how to read in the data and extract the corresponding parameters. There is also an example demonstrating how to interpolate the data to track the evolution of the orbital elements for a protoplanet in a typical disc, recreating Fig. 10 from the paper.



