% This script uses code from: 
%   - T. I. Fossen and T. Perez (2004). Marine Systems Simulator (MSS)
%     URL: https://github.com/cybergalactic/MSS
%
% Code adapted/modified, and additional code written by Andreas Mork Nordsven, 03.12.2025
%
% Dependicies: 
%   - MSS Toolbox
%   - CasADi Toolbox (version 3.7.2)
%   - Custom functions: 
%       * params_casADi.m
%       * PIDnonlinearMIMO (Original from MSS)
%       * navalvessel.m (Modified version of navalvessel.m from MSS)
%
% casADi can be downloaded from here: https://web.casadi.org/get/

clear PIDnonlinearMIMO
clearvars;
close all; 
clc;

addpath('C:\Users\andre\Downloads\casadi-3.7.2-windows64-matlab2018b')
import casadi.*

params = params_casADi();   % get parameters from params_casADi

%% ------------------------------------------------------------------------
sim_t = 30;               %        (*SET FINAL SIMULATION TIME HERE*)
dt = 0.05;                 % Sampling time

% Obstacle 1
xo_1 = -10;
yo_1 = -10;               %         (*SET OBSTACLE POSITIONS HERE*)
% Obstacle 2
xo_2 = -7;
yo_2 = -2;
% Obstacles radius        %         (*SET OBSTACLE RADIUS HERE*)
d = 3;

%------------------------       (*SET TUNE-ABLE CBF PARAMETERS HERE*)
cbf_gamma0 = 2;           % Parameter in CBF for obstacle 1
cbf_gamma0_2 = 2;         % Parameter in CBF for obstacle 2

cbf_rate = 0.2;           % CBF rate in CBF constraint for obstacle 1
cbf_rate2 = 0.15;         % CBF rate in CBF constraint for obstacle 2
%------------------------ 

% Use solver for obs 1 as long as time < time_interval1
% Use solver for obs 2 as long as time > time_interval1

time_interval1 = 15; %              (*SET TIME FOR SWITCH LOGIC HERE*)

Vc = 0;                   %           (*SET OCEAN CURRENT HERE*)
betaVc = deg2rad(-45);

%% ----------Building symbolic CBF-QP solver below this line-------------%%

% States
u = SX.sym('u');
v = SX.sym('v');
p = SX.sym('p');
r = SX.sym('r');
b = SX.sym('b');    % phi replaced with b
n = SX.sym('n');    % psi replaced with n
x_p = SX.sym('x_p');
y_p = SX.sym('y_p');

x1 = [u; v; p; r; b; n; x_p; y_p]; 

% Auxiliary variables
au = abs(u);
av = abs(v);
ar = abs(r);
ap = abs(p);

% Mass-Inertia Matrix 
M = [ (params.m-params.Xudot)      0                               0                               0                               0  0;
      0 (params.m-params.Yvdot)  -(params.m*params.zG+params.Ypdot) (params.m*params.xG-params.Yrdot) 0  0;
      0 -(params.m*params.zG+params.Kvdot) (params.Ixx-params.Kpdot) -params.Krdot                    0  0;
      0 (params.m*params.xG-params.Nvdot)  -params.Npdot (params.Izz-params.Nrdot)                   0  0;
      0 0 0 0 1 0;
      0 0 0 0 0 1 ];

% Hydrodynamic + Centripetal
Xh = params.Xuau*u*au + params.Xvr*v*r;

Yh = params.Yauv*au*v + params.Yur*u*r + params.Yvav*v*av + params.Yvar*v*ar + params.Yrav*r*av ...
   + params.Ybauv*b*abs(u*v) + params.Ybaur*b*abs(u*r) + params.Ybuu*b*u^2;

Kh = params.Kauv*au*v + params.Kur*u*r + params.Kvav*v*av + params.Kvar*v*ar + params.Krav*r*av ...
   + params.Kbauv*b*abs(u*v) + params.Kbaur*b*abs(u*r) + params.Kbuu*b*u^2 + params.Kaup*au*p ...
   + params.Kpap*p*ap + params.Kp*p + params.Kbbb*b^3 ...
   - (params.rho_water*params.g*params.gm*params.disp)*b;

Nh = params.Nauv*au*v + params.Naur*au*r + params.Nrar*r*ar + params.Nrav*r*av ...
   + params.Nbauv*b*abs(u*b) + params.Nbuar*b*u*ar + params.Nbuau*b*u*au;

Xc =  params.m*(r*v + params.xG*r^2 - params.zG*p*r);
Yc = -params.m*u*r;
Kc =  params.m*params.zG*u*r;
Nc = -params.m*params.xG*u*r;

F1 = Xh + Xc;
F2 = Yh + Yc;
F4 = Kh + Kc;
F6 = Nh + Nc;

% f in x_dot = f(x) + g(x)*tau
sola = M\[F1; F2; F4; F6; p; r];

% Position kinematics
x_dot_y_dot = [u*cos(n) - v*sin(n);
               u*sin(n) + v*cos(n)];

f = [sola; x_dot_y_dot];     % 8x1 drift vector

% g in x_dot = f(x) + g(x)*tau
Minv = inv(M);
g_body = Minv(:, [1 2 4]);  % we have forces/moment in Xe, Ye, Ne (surge, sway, yaw)

g = [g_body; zeros(2,3)];   % extend to 8x3

%% CBFs and Lie Derivatives 

e1_1 = x_p - xo_1;
e2_1 = y_p - yo_1;
e_1 = [e1_1, e2_1]';

B1_1 = -(e1_1)^2 - (e2_1)^2 + (d)^2;    % B1 for obstacle 1
B1_1_dot = -2 * e_1' * x_dot_y_dot;

cbf_1 = B1_1_dot + cbf_gamma0 * B1_1;   % B2 for obstacle 1

e1_2 = x_p - xo_2;
e2_2 = y_p - yo_2;
e_2 = [e1_2, e2_2]';

B1_2 = -(e1_2)^2 - (e2_2)^2 + (d)^2;  % B1 for obstacle 2
B1_2_dot = -2*e_2'*x_dot_y_dot; 

cbf_2 = B1_2_dot + cbf_gamma0_2 * B1_2;  % B2 for obstacle 2

% Lie Derivatives
dcbf_1 = jacobian(cbf_1, x1);
lf_cbf_1 = dcbf_1 * f;
lg_cbf_1 = dcbf_1 * g;

cbf_1_func = Function('cbf1', {x1}, {cbf_1});    % Call-able function B2 for obs 1
%Lf1_func = Function('Lf1_func', {x}, {lf_cbf_1});
%Lg1_func = Function('Lg1_func', {x}, {lg_cbf_1});
% To use function: cbf1_val = full(cbf_1_func(x_val));
B1_1_func = Function('B1_1', {x1}, {B1_1});

dcbf_2 = jacobian(cbf_2, x1);
lf_cbf_2 = dcbf_2 * f;
lg_cbf_2 = dcbf_2 * g;

cbf_2_func = Function('cbf2', {x1}, {cbf_2});    % Call-able function B2 for obs 2
%Lf2_func = Function('Lf2_func', {x}, {lf_cbf_2});
%Lg2_func = Function('Lg2_func', {x}, {lg_cbf_2});
B1_2_func = Function('B1_2', {x1}, {B1_2}); 

%% CBF-QP problem setup
tau = SX.sym('tau', 3);  % Optimization variables
tau_ref = SX.sym('tau_ref', 3);

% In order to pre-build the solver we need to put 
% the state vector and tau_ref into the parameters field
P = [x1; tau_ref];   

%---------------Create CBF solver for subproblem 1-------------------------
A = [lg_cbf_1];

b = [lf_cbf_1 + cbf_rate * cbf_1]; 

H = eye(3); 
obj = (tau - tau_ref)'*H*(tau - tau_ref);

prob_struct1 = struct('x', tau, 'f', obj, 'g', A*tau + b, 'p', P);  % x is opt variable (x = tau)

% Note: A*tau + b <= 0. lb = -inf, ub = 0

opts = struct();
opts.printLevel = 'none';

solver1 = qpsol('solver', 'qpoases', prob_struct1, opts);   % QP solver first subproblem

%---------------Create CBF solver for subproblem 2-------------------------
A = [lg_cbf_2];

b = [lf_cbf_2 + cbf_rate2 * cbf_2];

prob_struct2 = struct('x', tau, 'f', obj, 'g', A*tau + b, 'p', P);

solver1_2 = qpsol('solver', 'qpoases', prob_struct2, opts);   % QP solver second subproblem

%% --------------------------------------------------------------------- %%


%% -----------Build control allocation symbolic below this line ---------%% 
% The same control allocation algorithm as the one in MSS/VESSELS/SIMosv
% This code use casADi and solver IPOPT. MSS/VESSELS/SIMosv use SQP/fmincon.

alpha = SX.sym('alpha',2);
u = SX.sym('u', 4);
s = SX.sym('s', 3);
x2 = [alpha; u; s];    % Optimization variables

az_max = deg2rad(60);  % Max azimuth rotation angle (rad)
l_x = [37, 35, -51.5/2, -51.5/2];
l_y = [0, 0, 7, -7];
K_max = diag([300e3, 300e3, 655e3, 655e3]); % Max thrust for each propeller (N)
n_max = [140, 140, 150, 150]';  

args = struct;   % args contains ub, lb for x and g
args.lbx = [-az_max, -az_max, -1, -1, -1, -1, -inf, -inf, -inf];
args.ubx = [az_max, az_max, 1, 1, 1, 1, inf, inf, inf];

alpha_old = SX.sym('alpha_old',2);
u_old = SX.sym('u_old', 4);
tau_pid = SX.sym('tau_pid', 3);

% Need to put alpha_old, u_old, tau_pid in P matrix
P = [alpha_old; u_old; tau_pid];

w1 = 1;     % Weight for squared u
w2 = 100;   % Weight for squared slack variable s
w3 = 1;     % Weight for squared change in alpha
w4 = 0.1;   % Weight for squared change in u

% Cost function
obj = w1 * norm(u)^2 ...
    + w2 * norm(s)^2 ...
    + w3 * norm(alpha - alpha_old)^2 ...
    + w4 * norm(u - u_old)^2;

max_rate_alpha = 0.3;  % 0.3 rad/s = 17.2 deg/s
max_rate_u = 0.1;

% Inequality constraints
c1 = ( alpha(1)-alpha_old(1)) / dt - max_rate_alpha; 
c2 = -( alpha(1)-alpha_old(1)) / dt - max_rate_alpha;
c3 = ( alpha(2) - alpha_old(2) ) / dt - max_rate_alpha;
c4 = -( alpha(2) - alpha_old(2) ) / dt - max_rate_alpha;
c5 = ( u(1) - u_old(1) ) / dt - max_rate_u;
c6 = -( u(1) - u_old(1) ) / dt - max_rate_u;
c7 = ( u(2) - u_old(2) ) / dt - max_rate_u;
c8 = -( u(2) - u_old(2) ) / dt - max_rate_u;
c9 = ( u(3) - u_old(3) ) / dt - max_rate_u;
c10 = -( u(3) - u_old(3) ) / dt - max_rate_u;
c11 = ( u(4) - u_old(4) ) / dt - max_rate_u;
c12 = -( u(4) - u_old(4) ) / dt - max_rate_u;

% Equality constraint
T_thr = SX.zeros(3, 4);  % Thrust configuration matrix
T_thr(:,1) = [0 1 l_x(1)]';
T_thr(:,2) = [0 1 l_x(2)]';
T_thr(:,3) = [cos(alpha(1)) sin(alpha(1)) l_x(3)*sin(alpha(1))-l_y(3)*cos(alpha(1))]';
T_thr(:,4) = [cos(alpha(2)) sin(alpha(2)) l_x(4)*sin(alpha(2))-l_y(4)*cos(alpha(2))]';

ceq = T_thr * K_max * u - tau_pid + s; 

g = [c1; c2; c3; c4; c5; c6; c7; c8; c9; c10; c11; c12; ceq];
% lb for c1-c11: -inf, lb for ceq: 0, ub for c1-c11: 0, ub for ceq: 0

args.lbg = [-inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf,...
            0, 0, 0];   % last 3 is for ceq
args.ubg = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
            0, 0, 0];   % last 3 is for ceq

prob_struct = struct('x', x2, 'f', obj, 'g', g, 'p', P);

opts = struct;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.max_iter = 100;
opts.ipopt.tol = 1e-6;

solver2 = nlpsol('solver', 'ipopt', prob_struct, opts);  % Control alloc solver

% Initial guess: x0 = [az1 az2 u1 u2 u3 u4 s1 s2 s3]
args.x0 = [deg2rad(-28) deg2rad(28)  0 0 0 0  0 0 0];

%% ---Params for DP, PID controller, init for states, time, allocation---%%

% Define DP setpoints
x_ref = 0;                   % Reference North position in meters
y_ref = 0;                   % Reference East position in meters
psi_ref = deg2rad(0);        % Reference yaw angle in radians
eta_ref = [x_ref, y_ref, psi_ref]';  % Reference positions and heading

% Initialize the nonlinear MIMO PID controller          
wn =  0.1 * diag([1 1 3]);   % Natural frequencies for PID tuning
zeta = 1.0 * diag([1 1 1]);  % Damping ratios for PID tuning
T_f = 30;                    % Time constant for the setpoint low-pass filter (s)

total_k = ceil(sim_t / dt);     
t = 0;

% Initial states for ship 
x = [0; 0; 0; 0; 0; deg2rad(0)];       % x = [u; v; p; r; phi; psi]'
eta = [-15; -15; deg2rad(0)];            % north, east, psi (x, y, psi)

% Constant azimuth angles minimizing the condition number of T_thr                             
alpha0 = deg2rad([-28; 28]); 

alpha_old = alpha0;    % Initial values for dynamic optimization
u_old = [0, 0, 0, 0]'; % Initial propeller speeds 

ts   = zeros(total_k, 1);    % Store timesteps
xs   = zeros(total_k, 6);    % Store states for naval vessel
etas = zeros(total_k, 3);    % Store north, east, psi for naval vessel
us = zeros(total_k-1, 3);    % Store control inputs from CBF-QP controller
bs1 = zeros(total_k-1, 1);   % Store B1-CBF value for obs 1
bs2 = zeros(total_k-1, 1);   % Store B1-CBF value for obs 2
hs1 = zeros(total_k-1, 1);   % Store B2-CBF value for obs 1
hs2 = zeros(total_k-1, 1);   % Store B2-CBF value for obs 2

uis = zeros(total_k-1, 6);   % Store physical control inputs
uis(1, :) = [u_old', alpha_old'];   

xs(1, :)   = x';
etas(1, :) = eta;
ts(1)      = t;

h_waitbar = waitbar(0, 'Processing...');    % Display a wait bar 

% Store status from the two CBF solvers
CBF_solvers_info = struct( ...
    'success', cell(total_k-1, 1), ...         
    'return_status', cell(total_k-1, 1), ...
    'iter_count', cell(total_k-1, 1));
% Store status from the control alloc solver
controlalloc_solver_info = struct( ...
    'success', cell(total_k-1, 1), ...         
    'return_status', cell(total_k-1, 1), ...
    'iter_count', cell(total_k-1, 1));

%% --------------- Simulation loop (Main) starts here --------------------------%%
 
for k = 1:total_k-1

    % Update the progress bar every 10 iterations
    if mod(k, 10) == 0
       elapsedTime = ts(k) / sim_t;
       waitbar(elapsedTime, h_waitbar, ...
           sprintf('Progress: %3.0f%%', 100*elapsedTime));
    end

    % Control logic based on the elapsed simulation time
    if ts(k) > 80 
        eta_ref = [x_ref, y_ref, deg2rad(40)]'; % Change setpoint after 50 s
    end
    
    % 3 DOF model, 3 DOF tau
    nu = [x(1); x(2); x(4)];   % [u v r]
    tau_pid = PIDnonlinearMIMO(eta, nu, eta_ref, M([1 2 4],[1 2 4]), wn, zeta, T_f, dt); % ref controller
    
    
    % CBF-QP solver 2 subproblems. P = [x; tau_ref], x = [u v p r phi psi x_p y_p] 
    if ts(k) < time_interval1
        sol1 = solver1('p', [x; eta(1); eta(2); tau_pid], 'lbg', -inf(1,1), 'ubg', zeros(1,1));
        stats = solver1.stats();
    else
        sol1 = solver1_2('p', [x; eta(1); eta(2); tau_pid], 'lbg', -inf(1,1), 'ubg', zeros(1,1));
        stats = solver1_2.stats();
    end

    CBF_solvers_info(k).success = stats.success;
    CBF_solvers_info(k).return_status = stats.return_status;
    CBF_solvers_info(k).iter_count = stats.iter_count;

    tau_opt = full(sol1.x);  % surge, sway yaw
    us(k, :) = tau_opt';      % Could add wind forces (tau_w) here

    bs1(k) = full(B1_1_func([x; eta(1); eta(2)]));    % B1 obs1 values stored  
    bs2(k) = full(B1_2_func([x; eta(1); eta(2)]));    % B1 obs2 values stored
    hs1(k) = full(cbf_1_func([x; eta(1); eta(2)]));   % B2 obs1 values stored for each timestep
    hs2(k) = full(cbf_2_func([x; eta(1); eta(2)]));   % B2 obs2 values stored for each timestep 

    %% Control allocation
    
    sol2 = solver2('p', [alpha_old; u_old; tau_opt],'x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, ...
                   'lbg', args.lbg, 'ubg', args.ubg);

    stats2 = solver2.stats();
    controlalloc_solver_info(k).success = stats2.success;
    controlalloc_solver_info(k).return_status = stats2.return_status;
    controlalloc_solver_info(k).iter_count = stats2.iter_count;
    
    sol_opt = full(sol2.x);
    alpha_c = sol_opt(1:2);
    u_c = sol_opt(3:6);

    alpha_old = alpha_c;
    u_old = u_c;

    % Controls: ui = [ n_c(1) n_c(2) n_c(3) n_c(4) alpha_c(1) alpha_c(2) ]'
    u_c = n_max.^2 .* u_c;  % Scale control efforts to actual propeller speeds
    n_c = sign(u_c) .* sqrt(abs(u_c));  % Calculate each propeller's speed
    ui = [n_c; alpha_c];  

    %% Propagate naval vessel dynamics (using ode45).
    [ts_temp, xs_temp] = ode45(@(t, x) navalvessel(x, ui, Vc, betaVc), [t t+dt], x);
    x = xs_temp(end, :)';
    
    % Update ship position
    eta = eta + dt * Rzyx(0, 0, x(6)) * [x(1); x(2); x(4)];

    xs(k+1, :)   = x';
    etas(k+1, :) = eta;
    uis(k+1, :) = ui';
    ts(k+1)      = ts_temp(end);

    t = t + dt;
end

%% ---------------Print message in command window -----------------------%%
dist_obs1 = sqrt((etas(:,1) - xo_1).^2 + (etas(:,2) - yo_1).^2);
dist_obs2 = sqrt((etas(:,1) - xo_2).^2 + (etas(:,2) - yo_2).^2);

min_dist = d;    % constraint boundary
tol = 1e-8;             % tolerance

viol1 = dist_obs1 < min_dist - tol;  % - tol if we tolerate some error
viol2 = dist_obs2 < min_dist - tol;

% Check if any violations occurred and print messages
% Print violation times for Obs1
if any(viol1)
    fprintf('Ship violated safety distance for Obstacle 1 at:\n');
    fprintf('  t = %.3f s', ts(viol1));
end

% Print violation times for Obs2
if any(viol2)
    fprintf('\nShip violated safety distance for Obstacle 2 at:\n');
    fprintf('  t = %.3f s', ts(viol2));
end

% If no violations
if ~any(viol1) && ~any(viol2)
    fprintf('No safety distance violations detected.\n');
end

fprintf('\n');


%% -------Open CBF solvers and control alloc solver information----------%%
openvar('CBF_solvers_info');
openvar('controlalloc_solver_info');

%% -------------------------Plotting------------------------------------ %%        
close(h_waitbar);  % Close the progress indicator

legendLocation = 'best';

figure(1);
h1 = plot(etas(1,2), etas(1,1), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Start Ship');
hold on;
h2 = plot(yo_1, xo_1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Obstacle 1 Origin');
hold on;
h3 = plot(etas(:,2), etas(:,1), 'b-', 'LineWidth', 2, 'DisplayName', 'Ship Trajectory '); hold on;
hold on;
h4 = plot(yo_2, xo_2, 'r*', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Obstacle 2 Origin');
hold on

theta = linspace(0, 2*pi, 200);
circle_x = yo_1 + d * cos(theta);  % East
circle_y = xo_1 + d * sin(theta);  % North
plot(circle_x, circle_y, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Safety Radius (3 m)');

hold on; 

circle_x_2 = yo_2 + d *cos(theta);
circle_y_2 = xo_2 + d *sin(theta);
plot(circle_x_2, circle_y_2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Safety Radius (3 m)');

xlabel('East [m]');
ylabel('North [m]');
legend([h1, h2, h3 h4], 'Location', 'best');
title('Trajectory of Evasive Ship and Obstacle');
grid on;

figure(2);      % Plot control inputs for Xe, Ye, Ke, Ne
subplot(411)
plot(ts(1:end-1), us(:,1), 'b', 'LineWidth', 2);
title('Surge Control Force (Xe)')
grid on;
subplot(412)
plot(ts(1:end-1), us(:,2), 'g', 'LineWidth', 2);
title('Sway Control Force (Ye)')
grid on;
subplot(413)
plot(ts(1:end-1), 0*us(:,3), 'm', 'LineWidth', 2);  % No force in Ke
title('Roll Control Moment (Ke)')
grid on;
subplot(414)
plot(ts(1:end-1), us(:,3), 'r', 'LineWidth', 2);
title('Yaw Control Moment (Ne)')
grid on;

figure(3);  % Plot control barrier function
subplot(411)
plot(ts(1:end-1), bs1, 'r', 'LineWidth', 2);
title('Control Barrier Function (B1) for obstacle 1');
grid on;
subplot(412)
plot(ts(1:end-1), bs2, 'r', 'LineWidth', 2);
title('Control Barrier Function (B1) for obstacle 2');
grid on;
subplot(413)
plot(ts(1:end-1), hs1, 'b', 'LineWidth', 2);
title('Control Barrier Function (B2) for obstacle 1');
grid on;
subplot(414)
plot(ts(1:end-1), hs2, 'b', 'LineWidth', 2)
title('Control Barrier Function (B2) for obstacle 2')
grid on;

figure(4);
subplot(611) 
plot(ts, xs(:,1), 'b', 'LineWidth', 2); grid on;
title('Surge for evasive ship')
subplot(612) 
plot(ts, xs(:,2), 'g', 'LineWidth', 2);grid on;
title('Sway for evasive ship')
subplot(613) 
plot(ts, rad2deg(xs(:,3)), 'm', 'LineWidth', 2);grid on;
title('roll rate for evasive ship')
subplot(614) 
plot(ts, rad2deg(xs(:,4)), 'r', 'LineWidth', 2);grid on;
title('yaw rate for evasive ship')
subplot(615) 
plot(ts, rad2deg(xs(:,5)), 'k', 'LineWidth', 2);grid on;
title('roll angle for evasive ship')
subplot(616) 
plot(ts, rad2deg(xs(:,6)), 'k', 'LineWidth', 2);grid on;
title('yaw angle for evasive ship')

figure(5);
subplot(3,1,1) 
plot(ts(1:end), uis(:,1), 'LineWidth', 1.5); hold on;
plot(ts(1:end), uis(:,2), 'LineWidth', 1.5);

xlabel('Time [s]'), grid
legend('n_1','n_2','Location',legendLocation); title('Bow thrusters (RPM)'); grid on;

subplot(3,1,2)
plot(ts(1:end), uis(:,3), 'LineWidth', 1.5); hold on;
plot(ts(1:end), uis(:,4), 'LineWidth', 1.5);
xlabel('Time [s]'), title('Stern azimuth thrusters (RPM)')
legend('n_3','n_4','Location',legendLocation); grid on;

subplot(3,1,3)
plot(ts(1:end), rad2deg(uis(:,5)), 'LineWidth', 1.5); hold on;
plot(ts(1:end), rad2deg(uis(:,6)), 'LineWidth', 1.5);
xlabel('Time [s]'); 
legend('\alpha_1','\alpha_2','Location',legendLocation); title('Azimuth angles'); grid on;




