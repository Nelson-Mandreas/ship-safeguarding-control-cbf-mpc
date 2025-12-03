# ship-safeguarding-control-cbf-mpc
This repository presents MATLAB code for simulating safeguarding control of a marine vessel. Control Barrier Functions (CBF) and Model Predictive Control (MPC) are implemented as distinct safety filters to a Dynamic Positioning (DP) nominal controller. The repository enables comparison of the performance of CBF and MPC as safety filters under environmental disturbances, including wind and ocean currents.

Scripts uses code from: 
  - T. I. Fossen and T. Perez (2004). Marine Systems Simulator (MSS)
    URL: https://github.com/cybergalactic/MSS
  - Mohamed W. Mehrez (2017). MPC-and-MHE-implementation-in-MATLAB-using-Casadi
    URL: https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi

Dependicies: 
  * MSS Toolbox
  * casADi Toolbox (version 3.7.2)

casADi can be downloaded from here: https://web.casadi.org/get/
