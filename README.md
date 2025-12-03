# ship-safeguarding-control-cbf-mpc
This repository presents Master Thesis MATLAB code for simulating safeguarding control of a marine vessel. Control Barrier Functions (CBF) and Model Predictive Control (MPC) are implemented as distinct safety filters to a Dynamic Positioning (DP) nominal controller. The repository enables comparison of the performance of CBF and MPC as safety filters under environmental disturbances, including wind and ocean currents. Master’s thesis conducted from August 2025 to January 2026.

Scripts in Master Thesis uses code from: 
  - T. I. Fossen and T. Perez (2004). Marine Systems Simulator (MSS)
    URL: https://github.com/cybergalactic/MSS
  - Mohamed W. Mehrez (2017). MPC-and-MHE-implementation-in-MATLAB-using-Casadi
    URL: https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi

Requirements: 
  * MATLAB
  * MSS Toolbox
  * casADi Toolbox (version 3.7.2)

casADi can be downloaded from here: https://web.casadi.org/get/

Previous work on safeguarding control of ship using CBF is given in the "Project Code" folder. Project conducted from January 2025 to May 2025.

Scripts in Project uses code from: 
  - T. I. Fossen and T. Perez (2004). Marine Systems Simulator (MSS)
    URL: https://github.com/cybergalactic/MSS
  - Jason J. Choi (2020). CBF-CLF-Helper 1.0: Library for Control Barrier Function (CBF) and Control Lyapunov Function (CLF) based control methods
    URL: https://github.com/HybridRobotics/CBF-CLF-Helper

Requirements: 
  * MATLAB
  * MSS Toolbox
  * Symbolic Math Toolbox

