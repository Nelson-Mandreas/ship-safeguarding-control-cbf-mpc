% This script uses code from: 
%   - T. I. Fossen and T. Perez (2004). Marine Systems Simulator (MSS)
%     URL: https://github.com/cybergalactic/MSS
%
%   - Mohamed W. Mehrez (2017). MPC-and-MHE-implementation-in-MATLAB-using-Casadi
%     URL: https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi
%
% Code adapted/modified, and additional code written by Andreas Mork Nordsven, 03.12.2025
%
% Dependicies: 
%   - MSS Toolbox
%   - CasADi Toolbox (version 3.7.2)
%   - Custom functions: 
%       * params_casADi.m
%       * PIDmodified.m (Modified version of PIDnonlinearMIMO.m from MSS)
%       * navalvessel.m (Modified version of navalvessel.m from MSS)
%
% casADi can be downloaded from here: https://web.casadi.org/get/

clear PIDmodified
clearvars;
close all; 
clc;

addpath('C:\Users\andre\Downloads\casadi-3.7.2-windows64-matlab2018b')
import casadi.*

params = params_casADi();   % get parameters from params_casADi

%% --------------------------------------------------------------------- %%
% Set obstacle positions, size, prediction horizon (N), switch times,
% disturbances

T = 0.05;      % Sampling time
sim_t = 50;  %                         (*SET FINAL SIMULATION TIME HERE*)

N = 20;      % Prediction horizon                 (*SET N HERE*)

%                                     (*SET Initial states for ship HERE*) 
x0 = [0; 0; 0; 0; 0; deg2rad(0)];       % x = [u; v; p; r; phi; psi]'
eta = [-15; -15; deg2rad(0)];            % north, east, psi (x, y, psi)

% Set x-, y-position of obstacle 1 and 2
obs_x1 = -10; %                         (*SET POSITIONS OF OBSTACLES HERE*)      
obs_y1 = -10; 
obs_x2 = -7;
obs_y2 = -2;
obs_diam = 3; %                            (*SET OBSTACLE RADIUS HERE*)
% safety_margin = 0.5;                                    

% Use solver for obs 1 as long as time < time_interval1
% Use solver for obs 2 as long as time_interval1 < time < time_interval2
% When time > time_interval2 use only DP control law (PIDnonlinearMIMO)

time_interval1 = 20; %                  (*SET TIMES FOR SWITCH LOGIC HERE*)
time_interval2 = 50;

Vc = 0;          % Ocean current speed         (*SET OCEAN CURRENT HERE*)
betaVc = deg2rad(-45);      % Ocean current direction 


%% ----------Building symbolic MPC below this line-------------%%

% States
u = SX.sym('u');
v = SX.sym('v');
p = SX.sym('p');
r = SX.sym('r');
b = SX.sym('b');        % phi replaced with b
n = SX.sym('n');        % psi replaced with n
x_p = SX.sym('x_p');
y_p = SX.sym('y_p');

states = [u; v; p; r; b; n; x_p; y_p]; 
n_states = length(states);

tau = SX.sym('tau',3);
n_controls = length(tau);

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

fx = [sola; x_dot_y_dot];     % 8x1 drift vector

% g in x_dot = f(x) + g(x)*tau
Minv = inv(M);
g_body = Minv(:, [1 2 4]);  % we have forces/moment in Xe, Ye, Ne (surge, sway, yaw)

gg = [g_body; zeros(2,3)];   % extend to 8x3

rhs = fx + gg*tau;      % system r.h.s

f = Function('f',{states, tau}, {rhs}); 
U = SX.sym('U', n_controls, N);
P = SX.sym('P', n_states + n_controls*N);     % [x0; tau_pid_pred]

X = SX.sym('X', n_states, N+1);  % store prediction of states

X(:,1) = P(1:n_states); 
for k = 1:N
    st = X(:,k); con = U(:,k);
    f_value = f(st,con);
    st_next = st + (T*f_value);
    X(:,k+1) = st_next;
end

ff = Function('ff', {U,P}, {X});

obj = 0;
g = [];


R = zeros(3,3); R(1,1) = 0.1; R(2,2) = 0.1; R(3,3) = 0.5;
% compute objective
for k=1:N
    con = U(:,k);
    con_ref = P(n_states + (1:n_controls) + (k-1)*n_controls);
    obj = obj + (con-con_ref)' * R * (con-con_ref);
end

% Add constraints for collision avoidance

for k = 1:N+1
    g = [g; -sqrt((X(7,k)-obs_x1)^2+(X(8,k)-obs_y1)^2) + obs_diam];  % obs1
end

OPT_variables = reshape(U, 3*N,1); 
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);  

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level = 0;   % 0, 3
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob, opts);     % solver for 1st obstacle

args = struct;
args.lbg = [-inf*ones(1*(N+1),1)];  % inequality constraints    
args.ubg = [zeros(1*(N+1),1)];  % inequality constraints

g1 = [];
for k = 1:N+1
    g1 = [g1; -sqrt((X(7,k)-obs_x2)^2+(X(8,k)-obs_y2)^2) + obs_diam];  % obs2
end    

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g1, 'p', P);

solver1 = nlpsol('solver', 'ipopt', nlp_prob, opts);    % solver for 2nd obstacle


%% -----------Build control allocation symbolic below this line ---------%% 
% The same control allocation algorithm as the one in MSS/VESSELS/SIMosv
% This code use casADi and solver IPOPT. MSS/VESSELS/SIMos use SQP/fmincon.

alpha = SX.sym('alpha',2);
u = SX.sym('u', 4);
s = SX.sym('s', 3);
x2 = [alpha; u; s];    % Optimization variables

az_max = deg2rad(60);  % Max azimuth rotation angle (rad)
l_x = [37, 35, -51.5/2, -51.5/2];
l_y = [0, 0, 7, -7];
K_max = diag([300e3, 300e3, 655e3, 655e3]); % Max thrust for each propeller (N)
n_max = [140, 140, 150, 150]';  

args1 = struct;   % args contains ub, lb for x and g
args1.lbx = [-az_max, -az_max, -1, -1, -1, -1, -inf, -inf, -inf];
args1.ubx = [az_max, az_max, 1, 1, 1, 1, inf, inf, inf];

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

T = 0.05;      % Sampling time

% Inequality constraints
c1 = ( alpha(1)-alpha_old(1)) / T - max_rate_alpha; 
c2 = -( alpha(1)-alpha_old(1)) / T - max_rate_alpha;
c3 = ( alpha(2) - alpha_old(2) ) / T - max_rate_alpha;
c4 = -( alpha(2) - alpha_old(2) ) / T - max_rate_alpha;
c5 = ( u(1) - u_old(1) ) / T - max_rate_u;
c6 = -( u(1) - u_old(1) ) / T - max_rate_u;
c7 = ( u(2) - u_old(2) ) / T - max_rate_u;
c8 = -( u(2) - u_old(2) ) / T - max_rate_u;
c9 = ( u(3) - u_old(3) ) / T - max_rate_u;
c10 = -( u(3) - u_old(3) ) / T - max_rate_u;
c11 = ( u(4) - u_old(4) ) / T - max_rate_u;
c12 = -( u(4) - u_old(4) ) / T - max_rate_u;

% Equality constraint
T_thr = SX.zeros(3, 4);  % Thrust configuration matrix
T_thr(:,1) = [0 1 l_x(1)]';
T_thr(:,2) = [0 1 l_x(2)]';
T_thr(:,3) = [cos(alpha(1)) sin(alpha(1)) l_x(3)*sin(alpha(1))-l_y(3)*cos(alpha(1))]';
T_thr(:,4) = [cos(alpha(2)) sin(alpha(2)) l_x(4)*sin(alpha(2))-l_y(4)*cos(alpha(2))]';

ceq = T_thr * K_max * u - tau_pid + s; 

g = [c1; c2; c3; c4; c5; c6; c7; c8; c9; c10; c11; c12; ceq];
% lb for c1-c11: -inf, lb for ceq: 0, ub for c1-c11: 0, ub for ceq: 0

args1.lbg = [-inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf,...
            0, 0, 0];   % last 3 is for ceq
args1.ubg = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
            0, 0, 0];   % last 3 is for ceq

prob_struct = struct('x', x2, 'f', obj, 'g', g, 'p', P);

opts = struct;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.max_iter = 100;
opts.ipopt.tol = 1e-6;

solver2 = nlpsol('solver', 'ipopt', prob_struct, opts);

           
%% ----DP setpoints, inital conditions and precallocation--------- %%

% Define DP setpoints
x_ref = 0;                   % Reference North position in meters
y_ref = 0;                   % Reference East position in meters
psi_ref = deg2rad(0);        % Reference yaw angle in radians
eta_ref = [x_ref, y_ref, psi_ref]';  % Reference positions and heading

% Initialize the nonlinear MIMO PID controller          
wn =  0.1 * diag([1 1 3]);   % Natural frequencies for PID tuning
zeta = 1.0 * diag([1 1 1]);  % Damping ratios for PID tuning
T_f = 30;                    % Time constant for the setpoint low-pass filter (s)
 
total_k = ceil(sim_t / T);     
t = 0;

% Constant azimuth angles minimizing the condition number of T_thr                             
alpha0 = deg2rad([-28; 28]); 

alpha_old = alpha0;    % Initial values for dynamic optimization
u_old = [0, 0, 0, 0]'; % Initial propeller speeds 

u0 = zeros(N,3);    % 3 control inputs. Nx3 opt. variables

ts   = zeros(total_k, 1);    % Store timesteps
xs   = zeros(total_k, 6);    % Store states for naval vessel
etas = zeros(total_k, 3);    % Store north, east, psi for naval vessel
us = zeros(total_k-1, 3);

uis = zeros(total_k-1, 6);   % Store physical control inputs
uis(1, :) = [u_old', alpha_old'];  

xs(1, :)   = x0';
etas(1, :) = eta;
ts(1)      = t;

% Initial guess: x0 = [az1 az2 u1 u2 u3 u4 s1 s2 s3]
args1.x0 = [deg2rad(-28) deg2rad(28)  0 0 0 0  0 0 0];

tau_pid_pred = zeros(3,N);

h_waitbar = waitbar(0, 'Processing...');    % Display a wait bar 

% Store status from the two MPC solvers
MPC_solvers_info = struct( ...
    'success', cell(total_k-1, 1), ...         
    'return_status', cell(total_k-1, 1), ...
    'iter_count', cell(total_k-1, 1));
% Store status from the control alloc solver
controlalloc_solver_info = struct( ...
    'success', cell(total_k-1, 1), ...         
    'return_status', cell(total_k-1, 1), ...
    'iter_count', cell(total_k-1, 1));

%% --------------THE SIMULATION LOOP (MAIN) START FROM HERE-------------------- %%

for k = 1:total_k-1
  
    % Update the progress bar every 10 iterations
    if mod(k, 10) == 0
       elapsedTime = ts(k) / sim_t;
       waitbar(elapsedTime, h_waitbar, ...
           sprintf('Progress: %3.0f%%', 100*elapsedTime));
    end


    x_pred = [x0; eta(1); eta(2)];      % [u v p r phi psi x y]
    eta_pred = [x_pred(7); x_pred(8); x_pred(6)];   % [x y psi]
    nu_pred = [x_pred(1); x_pred(2); x_pred(4)];    % [u v r]
    
    if ts(k) < time_interval1 % Solver for the 1st obstacle
        
        for i = 1:N  % Predictive PID
            tau_pid_pred(:,i) = PIDmodified(eta_pred, nu_pred, eta_ref, ...
                                    M([1 2 4], [1 2 4]), wn, zeta, T_f, T);
        
            if i == 1    
                % Store persistent variables after first call
                [~, eta_d0, z_int0] = PIDmodified([], [], [], [], [], [], [], [], 'get_persistents'); 
            end

            f_dot = f(x_pred, tau_pid_pred(:,i));
            x_pred = x_pred + T*full(f_dot);    % propagate full state

            eta_pred = [x_pred(7); x_pred(8); x_pred(6)];
            nu_pred = [x_pred(1); x_pred(2); x_pred(4)];
        end

        % Set the persistent variables we got after one function call
        PIDmodified([], [], [], [], [], [], [], [], 'set_persistents', eta_d0, z_int0);

        args.p = [x0; eta(1); eta(2); reshape(tau_pid_pred, 3*N,1)];
        args.x0 = reshape(u0',3*N,1);

        sol = solver('x0', args.x0, ...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);

        stats = solver.stats();
        MPC_solvers_info(k).success = stats.success;
        MPC_solvers_info(k).return_status = stats.return_status;
        MPC_solvers_info(k).iter_count = stats.iter_count;

        u = reshape(full(sol.x)',3,N)';


    elseif ts(k) < time_interval2 && ts(k) > time_interval1  % solver for 2nd obstacle
        

        for i = 1:N  % Predictive PID
            tau_pid_pred(:,i) = PIDmodified(eta_pred, nu_pred, eta_ref, ...
                                    M([1 2 4], [1 2 4]), wn, zeta, T_f, T);
        
            if i == 1    
                % Store persistent variables after first call
                [~, eta_d0, z_int0] = PIDmodified([], [], [], [], [], [], [], [], 'get_persistents'); 
            end

            f_dot = f(x_pred, tau_pid_pred(:,i));
            x_pred = x_pred + T*full(f_dot);    % propagate full state

            eta_pred = [x_pred(7); x_pred(8); x_pred(6)];
            nu_pred = [x_pred(1); x_pred(2); x_pred(4)];
        end   
        
        % Set the persistent variables we got after one function call
        PIDmodified([], [], [], [], [], [], [], [], 'set_persistents', eta_d0, z_int0);

        args.p = [x0; eta(1); eta(2); reshape(tau_pid_pred, 3*N,1)];
        args.x0 = reshape(u0',3*N,1);

        sol = solver1('x0', args.x0, ...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);

        stats1 = solver1.stats();
        MPC_solvers_info(k).success = stats1.success;
        MPC_solvers_info(k).return_status = stats1.return_status;
        MPC_solvers_info(k).iter_count = stats1.iter_count;

        u = reshape(full(sol.x)',3,N)';
        
    else % use only PID
        
        u(1,:) = PIDmodified(eta_pred, nu_pred, eta_ref, ...
                                    M([1 2 4], [1 2 4]), wn, zeta, T_f, T);
    end
    
    us(k, :) = u(1,:);

    % Control allocation
    sol2 = solver2('p', [alpha_old; u_old; u(1,:)'],'x0', args1.x0, 'lbx', args1.lbx, 'ubx', args1.ubx, ...
                   'lbg', args1.lbg, 'ubg', args1.ubg);

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


    [ts_temp, xs_temp] = ode45(@(t, x0) navalvessel(x0, ui, Vc, betaVc), [t t+T], x0);
    x0 = xs_temp(end, :)';

    % Update ship position
    eta = eta + T * Rzyx(0, 0, x0(6)) * [x0(1); x0(2); x0(4)];

    xs(k+1, :)   = x0';
    etas(k+1, :) = eta;
    uis(k+1, :) = ui';
    ts(k+1)      = ts_temp(end);

    u0 = [u(2:size(u,1),:); u(size(u,1),:)];  % cut off the first u
    
    t = t + T;
end

%% ---------------Print message in command window -----------------------%%
dist_obs1 = sqrt((etas(:,1) - obs_x1).^2 + (etas(:,2) - obs_y1).^2);
dist_obs2 = sqrt((etas(:,1) - obs_x2).^2 + (etas(:,2) - obs_y2).^2);

min_dist = obs_diam;    % constraint boundary
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

%% -------Open MPC solvers and control alloc solver information----------%%
openvar('MPC_solvers_info');
openvar('controlalloc_solver_info');

%% -------------------------Plotting------------------------------------ %%        
close(h_waitbar);  % Close the progress indicator

legendLocation = 'best';

figure(1);
h1 = plot(etas(1,2), etas(1,1), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Start Ship');
hold on;
h2 = plot(etas(:,2), etas(:,1), 'b-', 'LineWidth', 2, 'DisplayName', 'Ship Trajectory ');
hold on
h3 = plot(obs_y1, obs_x1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Obstacle 1 Origin');
hold on;
h4 = plot(obs_y2, obs_x2, 'r*', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Obstacle 2 Origin');
hold on

theta = linspace(0, 2*pi, 200);
circle_x = obs_y1 + obs_diam * cos(theta);  % East
circle_y = obs_x1 + obs_diam * sin(theta);  % North
plot(circle_x, circle_y, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Safety Radius (3 m)');

hold on; 

circle_x_2 = obs_y2 + obs_diam *cos(theta);
circle_y_2 = obs_x2 + obs_diam *sin(theta);
plot(circle_x_2, circle_y_2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Safety Radius (3 m)');

grid on;
xlabel('East [m]');
ylabel('North [m]');
legend([h1, h2, h3, h4], 'Location', 'best');

figure(2);     % Plot control inputs for Xe, Ye, Ke, Ne
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

figure(3);
subplot(3,1,1); grid
plot(ts(1:end), uis(:,1), 'LineWidth', 1.5); hold on;
plot(ts(1:end), uis(:,2), 'LineWidth', 1.5);

xlabel('Time [s]'), grid
legend('n_1','n_2','Location',legendLocation); title('Bow thrusters (RPM)'); grid on;

subplot(3,1,2)
plot(ts(1:end), uis(:,3), 'LineWidth', 1.5); hold on;
plot(ts(1:end), uis(:,4), 'LineWidth', 1.5);
xlabel('Time [s]'), title('Stern azimuth thrusters (RPM)'), grid
legend('n_3','n_4','Location',legendLocation)

subplot(3,1,3), grid
plot(ts(1:end), rad2deg(uis(:,5)), 'LineWidth', 1.5); hold on;
plot(ts(1:end), rad2deg(uis(:,6)), 'LineWidth', 1.5);
xlabel('Time [s]');
legend('\alpha_1','\alpha_2','Location',legendLocation); title('Azimuth angles'); grid on;

% figure(4)
% stairs(ts(1:length(u_cl)), u_cl(:,1), 'LineWidth', 1.5); hold on;
% stairs(ts(1:length(u_cl)), u_cl(:,2), 'LineWidth', 1.5);
% stairs(ts(1:length(u_cl)), u_cl(:,3), 'LineWidth', 1.5);
% grid on;
% xlabel('Time [s]');
% ylabel('Control inputs');
% title('MPC control inputs (u\_cl)');
% legend('X_e','Y_e','N_e','Location','best');


figure(5); clf
hold on; grid on;
% Plot distances
plot(ts(1:end), dist_obs1, 'r', 'LineWidth', 2);
plot(ts(1:end), dist_obs2, 'b', 'LineWidth', 2);
% Add safety line
yline(min_dist, '--k', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Distance [m]');
legend('Distance to Obs 1', 'Distance to Obs 2', 'Safety Limit','Location','best');
title('Distance Between Ship and Obstacles');
ylim([min_dist - 1, max([dist_obs1; dist_obs2]) + 1]); % zoom to relevant range

figure(6);
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



