function [tau, eta_d_out, z_int_out] = PIDmodified(eta,nu,eta_ref,M,wn,zeta,T_f,h,varargin) 
% Nonlinear MIMO PID regulator for dynamic positioning (DP). For 3-DOF
% models (surge, sway, and yaw) the output is: tau = [tau1, tau2, tau6]'. 
% In the general case, tau = [tau1, tau2, 0, 0, 0, tau6]', where
%    
%    z_int = eta - eta_d,     where eta_d = 1 / (T_f * s + 1) * eta_ref
%
%    tau = -R(psi)' * ( Kp * (eta - eta_d) + Ki * z_int ) - Kd * nu
%
%    Kp = M_diag * wn * wn,                  M_diag = diag(diag(M))
%    Kd = M_diag * 2 * zeta * wn
%    Ki = 1/10 * Kp * wn
%
% is based on Algorithm 15.2, MIMO nonlinear PID Pole-Placement Algorithm, 
% by Fossen (2021); see also Equation (15.82). 
%
% Persistent variables: 
% The setpoint eta_d and integral state z_int are persistent variables that 
% should be cleared by adding:
%
%    clear PIDmodified
%
% on the top in the script calling PIDmodified.m.
%
% Inputs:  
%   eta: generalized position vector, 3x1 (surge, sway, and yaw) or 6x1
%   nu:  generalized velocity vector, 3x1 (surge, sway, and yaw) or 6x1
%   eta_ref: vector [x_ref,y_ref,psi_ref] of setpoints in surge, sway, and yaw
%   M: system inertia matrix, 3x3 (surge, sway, and yaw) or 6x6
%   wn: closed-loop natural frequencies, scalar or diagonal matrix 3x3
%   zeta: closed-loop relative damping ratios, scalar or diagonal matrix 3x3
%   T_f: setpointlow-pass filter time constant (s)
%   h: sampling time (s)
%
% Outputs:  
%    tau: generalized control force, 3x1 (surge, sway, and yaw) or 6x1
%    eta_d_out: persistent variable eta_d value
%    z_int_out: persistent variable z_int value
%  
% Author:    Thor I. Fossen
% Date:      2 Sep 2023
%
% Code adapted/modified, and additional code written by Andreas Mork Nordsven, 03.12.2025
% Revisions: 
%   - Modified to extract persistent variables eta_d and z_int
%   - Modified to set persistent variables eta_d and z_int
%
% Original PIDnonlinearMIMO (T. I. Fossen and T. Perez (2004)): https://github.com/cybergalactic/MSS
%
% Dependicies: 
%   * MSS toolbox

persistent eta_d;  % LP-filtered commands
persistent z_int;  % integral states

% Initialization of desired state eta_d and integral state z_int 
if isempty(z_int)
    z_int = [0 0 0]';             
    eta_d = [0 0 0]';
end

eta_d_out = [];
z_int_out = [];
tau = [];

% Routines for getting values of persistent variables and set them
    if nargin >= 9 && strcmp(varargin{1}, 'get_persistents')
        eta_d_out = eta_d;
        z_int_out = z_int;
        return
    elseif nargin >= 11 && strcmp(varargin{1}, 'set_persistents')
        eta_d = varargin{2};
        z_int = varargin{3};
        return
    end

eta_ref = eta_ref(:);  % Make sure eta_ref is a column vector

% Reduce 6-DOF model to 3-DOF model
DOF = 3;
if length(nu) == 6
    eta = [eta(1) eta(2) eta(6)]';
    nu  = [nu(1) nu(2) nu(6)]';
    M = M([1 2 6],[1 2 6]);
    DOF = 6;
end

% Rotation matrix in yaw
R = [ cos(eta(3)) -sin(eta(3)) 0
      sin(eta(3))  cos(eta(3)) 0
      0            0           1 ];

% MIMO pole placement (Algorithm 15.2)
M_diag = diag(diag(M));             % diagonal M matrix
Kp = M_diag .* wn .* wn;
Kd = M_diag .* 2 .* zeta .* wn;
Ki = 1/10 * Kp .* wn;

% 3-DOF control law for surge, sway and yaw
e = eta - eta_d;
e(3) = ssa( e(3) );
tau_PID = -R' * ( Kp * e + Ki * z_int ) - Kd * nu;

if DOF == 6
    tau = [tau_PID(1), tau_PID(2), 0, 0, 0, tau_PID(3)]'; 
else  % 3-DOF
    tau = tau_PID;
end

% integral state: z_int[k+1]
z_int = z_int + h * (eta - eta_d);

% low-pass filter: eta_d[k+1]
eta_d = eta_d + h * (eta_ref - eta_d)/ T_f;


end