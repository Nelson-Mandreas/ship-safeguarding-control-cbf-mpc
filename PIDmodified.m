function [tau, eta_d_out, z_int_out] = PIDmodified(eta,nu,eta_ref,M,wn,zeta,T_f,h,varargin)

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

% Rutines for getting values of persistent variables and set them
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