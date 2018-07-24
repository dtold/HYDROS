# HYDROS
A dispersion relation solver for the hybrid-kinetic model of plasma physics.

===== HYDROS === HYbrid Dispersion RelatiOn Solver =====

The HYDROS code is intended to enable an easy-to-use interface for obtaining solutions of the hybrid-kinetic
dispersion relation. The derivation of the model is described in Told et al. (New J. Phys. 18 (2016) 075001), 
and is closely related to many space plasma physics simulation codes. The model includes a fully kinetic ion 
description, and a massless fluid electron species using a simple isothermal or adiabatic equation of state. 
The Ohm's law that determines the electric field can include a resistive term to model the wave physics seen 
in a turbulence simulation with finite resistivity.

The code is written in Python and runs on both Python 2.7/3, provided that the Numpy/Scipy modules are included. 
You are welcome to use this code, and to adapt it to your own purposes. If you produce enhancements which may be
useful to other members of the community, you are welcome to submit your changes to the Github repository at
https://github.com/dtold/HYDROS. Also, if you publish results based on this code, please cite the original paper
by Told et al. (New J. Phys., 2016). 

===== How to use the code =====

The HYDROS code is called by running "./hydros.py \<inputfile\>", where "\<inputfile\>" is an arbitrarily named
file modeled after the "1dscan" template shipped with the code. In the following, all parameters are described
in detail:

-- wavevector_mode: Determines whether kperp and kpar are used to define the wave vector (set the parameter to 1),
   or whether k (the magnitude of the wavenumber) and theta are used (set the parameter to 2).

-- kperp: The component of the wavenumber perpendicular to the background magnetic field.

-- kpar: The component of the wavenumber parallel to the background magnetic field.

-- k: The absolute value of the wavenumber.

-- theta: The wave propagation angle (relative to the background magnetic field) in degrees.

-- beta: The ratio of the parallel ion thermal pressure to the magnetic pressure. (i.e. beta = 8 * pi *n * T_i|| / B^2)

-- tau: The ratio of the parallel ion temperature to the electron temperature (assumed isotropic).

-- Tpar_Tperp: The ratio of parallel to perpendicular ion temperature.

-- gam: The polytropic coefficient used for the electron equation of state. gam=1 corresponds to isothermal electrons,
   and gam=5/3 corresponds to adiabatic electrons.

-- eta: A normalized resistivity that can be added to the Ohm's law determining the electric field, which is often
   used in hybrid-kinetic turbulence codes.

-- nb: Maximum number of terms to use in the Bessel sums. The code will automatically truncate the sums when a
   preset accuracy (1e-12 by default) is achieved. When introducing a too low limit, the code will report
   convergence problems in its output.log file. 

-- start: Initial guess for complex frequency. Should be set close to a known value of the mode that is to be
   tracked. 

-- normalization: The normalization to be used for in/output of wave numbers and frequencies. A setting of 1
   corresponds to wavenumbers normalized in d_i units and frequencies normalized to the ion cyclotron frequency.
   A setting of 2 corresponds to wavenumbers normalized in rho_i|| units (i.e. thermal gyroradius computed with
   the parallel ion temperature) and frequencies normalized to k_|| * rho_i||. This is the setting the code uses
   internally.

-- scan_type: Can be set to '1d', '2d', and 'DR-map'. The first two should be self-explanatory (with the additional
   defined in the following), and the last setting creates a complex frequency plot for the single parameter set
   specified at the top of the input file, with no additional scanning. The parameters of the frequency scan range
   presently have to be adjusted within the code itself. 

-- scanvar: Defined which variable is to be scanned over, and can be any of the physics parameters from the top
   part of the input file.
   
-- scan_range: Defines the range of the scan in the format (start, end). The scan range variable will supercede
   the value defined for this variable in the top part of the file.
   
-- nval: The number of scanned values within the scan range. Using a higher number can be helpful if the solver
   loses track of the mode.
   
-- log: Set to 1 to use logarithmic spacing, or 0 to use linear spacing of the scan nodes.

-- scan_list: Can be used in lieu of scan_range, nval, log to define a fixed list of values that are to be scanned.

The scanvar, scan_range, nval, log, scan_list variables for the second scan dimension operate in the same way.

-- radius: If the solver decides that the solution is too far away from the expectation in terms of frequency,
   damping rate, amplitudes or cross phases, it will scan around the expected position to find a better solution.
   "radius" determines the maximum radius (relative to the magnitude of the expected complex frequency).

-- dev_limit: Determines how many orders of magnitude the new solution may deviate from the old solution. Applies
   to field amplitudes and the frequency.

-- cp_limit: Determines how much the cross phases of dBy/dBz and dn/dBz may differ between the new and old
   solution. Units are in terms of the full 2*pi range.

-- method: Which method to use for predicting the position of the new root. Possible choices are
   -- "default": use old solution as starting point
   -- "predict": use parabolic extrapolation to predict new root position
   -- if "-ldr" is appended to the method name, the predicted damping rate will be modified to a lower value, in
      order to prefer less damped solutions

-- resultfile: Set the filename into which the results will be written. The output format is:
   var1 & var2 & frequency & damping rate & Re(dBy/dBz) & Im(dBy/dBz) & Re(dn/dBz) & Im(dn/dBz

-- logfile: Set the filename of the log file, into which the more verbose output of the root finding and scanning
   procedures are saved. It is recommended that this file be inspected in case of any unexpected problems. 

Finally, if there are any further questions about aspects not covered in this brief manual, do not hesitate to send
an e-mail to dtold &AT& ipp.mpg.de.
