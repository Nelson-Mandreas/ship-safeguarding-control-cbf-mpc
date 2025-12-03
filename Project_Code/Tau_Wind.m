function tau_wind = Tau_Wind(x, Vw, betaVw)

    % Relative wind velocities
    u_rw = x(1) - Vw * cos(betaVw - x(6));   % x(6) = psi
    v_rw = x(2) - Vw * sin(betaVw - x(6));   % x(1) = u,  x(2) = v
    V_rw = sqrt( u_rw^2 + v_rw^2);
    gamma_rw = -atan2(v_rw, u_rw);

    AFw = 160.70;
    ALw = 434.80;
    sH = 5.10;
    sL = 1.48;
    L = 48.0;     % Lpp used

    % Wind forces
    tau_w = blendermann94(gamma_rw, V_rw, AFw, ALw, sH, sL, L, 13);  % Research vessel is no. 13
    
    tau_wind = [tau_w(1), tau_w(2), 0, tau_w(3)]';  % X, Y, N
    
end