clearvars;
close all; 
clc;

params = getParams_final();  % Retrieves ship, wind, current, obstacle parameters
%%
dt = 0.02;
sim_t = 20;
total_k = ceil(sim_t / dt);
t = 0;

%% Initial states for ship (w/ cbf)
x0 = [5; 0; 0; 0; 0; deg2rad(45)];       % x = [u; v; p; r; phi; psi]'
eta0 = [15; -10; deg2rad(45)];            % north, east, psi

u_ref = zeros(8,1);               

x = x0;            
eta = eta0;

psi_ref1 = deg2rad(45);        % Want ship to reach this heading eventually
psi_d1 = eta(3);

psi_int1 = 0;                   % Integral state for yaw control of naval vessel 1

%% Initial states for second naval vessel (w/o cbf)
x20 = [2; 0; 0; 0; 0; deg2rad(-45)];       
eta20 = [10; 60; deg2rad(-45)];          

x2 = x20;
eta2 = eta20;

% PID controller (Naval vessel 1 and 2 use the same gains)
Kp = 1e6;                    % Controller proportional gain
Ki = 1e4;                    % Controller integral gain
Td = 1;                      % Controller derivative time

% Reference signal
w_n = 0.1;                  % Low-pass filter cut-off frequency (rad/s)
psi_ref = deg2rad(-45);        % Desired heading
psi_d = eta2(3);                 

%% Initialize traces
ts   = zeros(total_k, 1);    % Store timesteps

xs   = zeros(total_k, 6);    % Store states for evasive naval vessel
etas = zeros(total_k, 3);    % Store north, east, psi for evasive naval vessel

xs2  = zeros(total_k, 6);    
etas2 = zeros(total_k, 3);

hs = zeros(total_k-1, 1);    % Store control barrier function value
us = zeros(total_k-1, 8);    % Store control inputs from CBF-QP controller

psi_ds = zeros(total_k-1, 1);   % Store psi_d for naval vessel 1

distances = zeros(total_k-1, 1);    % Store distance between the two ships

%% Create an instance of EvasiveManeuvering
maneuver = EvasiveManeuvering(params);

% Starting point for traces
xs(1, :)   = x0';
etas(1, :) = eta;
ts(1)      = t;
xs2(1, :)  = x20';
etas2(1, :) = eta2;
psi_ds(1) = psi_d1;
distances(1) = norm(eta(1:2)-eta2(1:2));

%% Simulation loop 
for k = 1:total_k-1
    t        % Show how far the simulation has come
    
    % Calculate wind forces/moments.
    tau_wind = Tau_Wind(x, params.Vw, params.betaVw);
    
    %% Ship measurements of naval vessel w/ cbf
    u_f    = x(1) + 0.001*randn;     % surge           % NOTE: Uten forstyrrelser på målingene
    v_f    = x(2) + 0.001*randn;     % sway              så kan jeg ha stor CBF-rate. 
    p_f    = x(3) + 0.001*randn;     % roll rate
    r_f    = x(4) + 0.001*randn;     % yaw rate
    phi_f  = x(5) + 0.001*randn;     % roll angle
    psi_f  = x(6) + 0.001*randn;     % yaw angle
    x_pos_f = eta(1) + 0.001*randn;  % North position (ship)
    y_pos_f = eta(2) + 0.001*randn;  % East position (ship)
    
    x_f = [u_f; v_f; p_f; r_f; phi_f; psi_f; x_pos_f; y_pos_f];
    
    %% Second naval vessel states (assumed perfect measurements)
    u_2     = x2(1);
    v_2     = x2(2);
    p_2     = x2(3);
    r_2     = x2(4);
    phi_2   = x2(5);
    psi_2   = x2(6);
    x_pos_2 = eta2(1);
    y_pos_2 = eta2(2);
    
    %% Second naval vessel control
    % Control system
    tauX = 1e5;                                    % Thrust
    tauN = -Kp * ( ssa(psi_2 - psi_d) + Td * r_2);  % PD heading controller
    
    tau_2 = [tauX+tau_wind(1) tau_wind(2) 0 tauN+tau_wind(4)]';
   
    %% Call the CBF-QP controller.
    tauX1 = 1e5;
    tauN1 = -Kp * ( ssa(psi_f - psi_d1) + Td * r_f + (Ki/Kp) * psi_int1);  % PID heading controller
    
    u_ref = [tauX1, 0, 0, tauN1, 0, 0, 0, 0]';
    % The second ship position stored in x_naval2, used in CBF-QP controller
    x_naval2 = [eta2(1)+0.001*randn; eta2(2)+0.001*randn]; 
    [Inputs, B, feas, comp_time] = maneuver.ctrlCbfQp(x_f, u_ref, 0, x_naval2);  % CBF-QP cpntroller
    us(k, :) = Inputs';        % Control inputs given by the CBF-QP controller
    hs(k)   = B;               % CBF values stored for each timestep
    
    % For the ship dynamics, tau = Inputs + wind.
    tau = Inputs(1:4) + tau_wind;
    
    %% Propagate naval vessel dynamics (using ode45).
    [ts_temp, xs_temp] = ode45(@(t, x) navalvessel(x, tau), [t t+dt], x);
    x = xs_temp(end, :)';
    
    % Update ship position
    eta = eta + dt * Rzyx(0, 0, x(6)) * [x(1); x(2); x(4)];
    psi_d1 = lowPassFilter(psi_d1, psi_ref1, w_n, dt);             % psi_d1 eventually reaches psi_ref1
    
    xs(k+1, :)   = x';
    ts(k+1)      = ts_temp(end);
    etas(k+1, :) = eta;
    psi_ds(k+1) = psi_d1;
    
    %% Propagate naval vessel 2 dynamics
    [tss_temp, xss_temp] = ode45(@(t,x2) navalvessel(x2, tau_2), [t t+dt], x2);
    x2 = xss_temp(end, :)';
    
    eta2 = eta2 + dt * Rzyx(0, 0, x2(6)) * [x2(1); x2(2); x2(4)];
    psi_d = lowPassFilter(psi_d, psi_ref, w_n, dt);        
    
    xs2(k+1, :) = x2';
    etas2(k+1, :) = eta2;
    
    distances(k+1) = norm(eta(1:2) - eta2(1:2));
   
    t = t + dt;
    
    psi_int1 = psi_int1 + dt * ssa(psi_f - psi_d1);   % Yaw integral state for naval vessel 1
    
end

figure(1);   % Plot north-east position of the two ships
% Plot the starting points as dots:
h1 = plot(etas(1,2), etas(1,1), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Start Vessel 1');
hold on;
h2 = plot(etas2(1,2), etas2(1,1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Start Vessel 2');
hold on;
h3 = plot(etas(:,2), etas(:,1), 'b-', 'LineWidth', 2, 'DisplayName', 'Vessel 1 (Evasive) Trajectory '); hold on;
h4 = plot(etas2(:,2), etas2(:,1), 'r--', 'LineWidth', 2, 'DisplayName', 'Vessel 2 (Intruder) Trajectory');
xlabel('East [m]');
ylabel('North [m]');
legend([h1, h2, h3, h4], 'Location', 'best');
title('Trajectories of the Two Naval Vessels');
grid on;

figure(2);  % Plot control barrier function
plot(ts(1:end-1), hs, 'b', 'LineWidth', 2);
xlabel('t [s]');
ylabel('B(x)');
title('Control Barrier Function');
grid on;

figure(3);  % Plot control inputs for Xe, Ye, Ke, Ne
subplot(411)
plot(ts(1:end-1), us(:,1), 'b', 'LineWidth', 2);
title('Surge Control Force (Xe)')
subplot(412)
plot(ts(1:end-1), us(:,2), 'g', 'LineWidth', 2);
title('Sway Control Force (Ye)')
subplot(413)
plot(ts(1:end-1), us(:,3), 'm', 'LineWidth', 2);
title('Roll Control Moment (Ke)')
subplot(414)
plot(ts(1:end-1), us(:,4), 'r', 'LineWidth', 2);
title('Yaw Control Moment (Ne)')

figure(4); % Plot states for Evasive Ship
subplot(611)
plot(ts, xs(:,1), 'b', 'LineWidth', 2);
title('Surge for evasive ship')
subplot(612)
plot(ts, xs(:,2), 'g', 'LineWidth', 2);
title('Sway for evasive ship')
subplot(613)
plot(ts, rad2deg(xs(:,3)), 'm', 'LineWidth', 2);
title('roll rate for evasive ship')
subplot(614)
plot(ts, rad2deg(xs(:,4)), 'r', 'LineWidth', 2)
title('yaw rate for evasive ship')
subplot(615)
plot(ts, rad2deg(xs(:,5)), 'k', 'LineWidth', 2)
title('roll angle for evasive ship')
subplot(616)
plot(ts, rad2deg(xs(:,6)), 'y', 'LineWidth', 2)
title('yaw angle for evasive ship')

% figure(5); % Plot states for Intruder Ship
% subplot(611)
% plot(ts, xs2(:,1), 'b', 'LineWidth', 2);
% subplot(612)
% plot(ts, xs2(:,2), 'g', 'LineWidth', 2);
% subplot(613)
% plot(ts, xs2(:,3), 'm', 'LineWidth', 2);
% subplot(614)
% plot(ts, xs2(:,4), 'r', 'LineWidth', 2)
% subplot(615)
% plot(ts, xs2(:,5), 'k', 'LineWidth', 2)
% subplot(616)
% plot(ts, xs2(:,6), 'y', 'LineWidth', 2)

figure(6); % Plot distance between the ships
plot(ts, distances, 'b', 'LineWidth', 2);
hold on;
plot(ts, ones(size(ts)) * params.d, 'r', 'LineWidth', 2);
hold on;
plot(ts, ones(size(ts))* 35, 'k', 'LineWidth', 2);
grid on; 
title("Distance between the ships")
ylabel("Distance [m]");
xlabel("Time [s]");
legend('Actual Distance', 'Safety Distance (d) used in CBF', 'Actual Safety Distance');


figure(7); % Plot psi_d for intruding naval vessel
plot(ts, rad2deg(psi_ds), 'b', 'LineWidth', 2);
hold on; 
plot(ts, rad2deg(xs(:,6)), 'r', 'LineWidth', 2);
grid on;
title("Yaw angle (deg)");
ylabel("Angle [deg]");
xlabel("Time [s]")
legend('Desired heading angle psi_d', 'Actual heading angle psi')