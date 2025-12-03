classdef EvasiveManeuvering < ctrlAffineSys
    methods       
        function [x, f, g] = defineSystem(obj, params)
            % Outputs: x: symbolic state vector
            %          f: drift term, expressed symbolically wrt x.
            %          g: control vector fields, expressed symbolically wrt x.
            
        syms u v p r b n x_pos y_pos        % phi is replaced with b, psi is replaced with n
        x = [u; v; p; r; b; n; x_pos; y_pos];   % ADD: u_t, v_t, x_t, y_t
        
        % Auxiliary variables
        au = abs(u);     
        av = abs(v);     
        ar = abs(r); 
        ap = abs(p); 
        ab = abs(b);
        L2 = params.Lpp^2; 
        
        M =[ (params.m-params.Xudot)  0   0   0   0   0;
           0 (params.m-params.Yvdot) -(params.m*params.zG+params.Ypdot) (params.m*params.xG-params.Yrdot) 0 0;
           0 -(params.m*params.zG+params.Kvdot) (params.Ixx-params.Kpdot) -params.Krdot 0 0;
           0 (params.m*params.xG-params.Nvdot) -params.Npdot (params.Izz-params.Nrdot) 0 0;
           0 0 0 0 1 0; 
           0 0 0 0 0 1] ;
        
        Xh  = params.Xuau*u*au+params.Xvr*v*r;

        Yh = params.Yauv*au*v + params.Yur*u*r + params.Yvav*v*av + params.Yvar*v*ar + params.Yrav*r*av ...
          + params.Ybauv*b*abs(u*v) + params.Ybaur*b*abs(u*r) + params.Ybuu*b*u^2;

        Kh = params.Kauv*au*v +params.Kur*u*r + params.Kvav*v*av + params.Kvar*v*ar + params.Krav*r*av ...
           + params.Kbauv*b*abs(u*v) + params.Kbaur*b*abs(u*r) + params.Kbuu*b*u^2 + params.Kaup*au*p...
           + params.Kpap*p*ap +params.Kp*p +params.Kbbb*b^3-(params.rho_water*params.g*params.gm*params.disp)*b;

        Nh = params.Nauv*au*v + params.Naur*au*r + params.Nrar*r*ar + params.Nrav*r*av...
           +params.Nbauv*b*abs(u*b) + params.Nbuar*b*u*ar + params.Nbuau*b*u*au;
       
        Xc =   params.m*(r*v+params.xG*r^2-params.zG*p*r);  
        Yc = - params.m*u*r;
        Kc =   params.m*params.zG*u*r;
        Nc = - params.m*params.xG*u*r;
       
        F1 = Xh+Xc;
        F2 = Yh+Yc;
        F4 = Kh+Kc;
        F6 = Nh+Nc;                                     
        
        % Dynamics
        sol = M\[F1; F2; F4; F6; p; r];
        extra_rows = [u*cos(n)-v*sin(n); u*sin(n)+v*cos(n)]; % x_pos_dot and y_pos_dot in NED
        f = [sol; extra_rows];
        
        sol_g = M\[1, 0, 0, 0, 0, 0;
                   0, 1, 0, 0, 0, 0;
                   0, 0, 1, 0, 0, 0;             % If you only want control forces/moments
                   0, 0, 0, 1, 0, 0;             % in e.g surge and yaw you put zeros in (2,2)
                   0, 0, 0, 0, 0, 0;             % and (3,3) and leave ones in (1,1), (4,4)
                   0, 0, 0, 0, 0, 0];
               
        g = [sol_g, zeros(size(M,1), 2); zeros(2, size(M,2) + 2)];
        % u = [Xe,Ye,Ke,Ne,0,0,0,0]
        
        end
        
        function [u, B, feas, comp_time] = ctrlCbfQp(obj, x, u_ref, verbose, x_obstacle)
            % Implementation of a basic CBF-QP.
            % Inputs:
            %   x: evasive ship state vector
            %   u_ref: reference control input
            %   verbose: flag for printing info (1) or silent (0)
            %   x_obstacle: obstacle position [xo; yo]
            %
            % Outputs:
            %   u: control input from the QP
            %   B: CBF value at state x 
            %   feas: feasibility flag (1 if feasible)
            %   comp_time: computation time for the QP.
            
            if isempty(obj.cbf)
                error('CBF is not defined. Provide a symbolic CBF in defineCbf.');
            end
            if nargin < 3 || isempty(u_ref)
                u_ref = zeros(obj.udim, 1);
            end
            if nargin < 4 || isempty(verbose)
                verbose = 0;
            end
            if nargin < 5
                error('Please provide the obstacle positions as a 2×1 vector [xo;yo].');
            end
            
            tstart = tic;
            
            % CBF and its Lie derivatives
            B = obj.cbf(x, x_obstacle);
            LfB = obj.lf_cbf(x, x_obstacle);
            LgB = obj.lg_cbf(x, x_obstacle);
            
            %% QP Constraint: -LgB * u <= LfB + rate * B
            A = -LgB;
            b = LfB + obj.params.cbf.rate * B;
            
            % Add input constraints if specified.
            if isfield(obj.params, 'u_max')
                A = [A; eye(obj.udim)];
                if isscalar(obj.params.u_max)
                    b = [b; obj.params.u_max * ones(obj.udim, 1)];
                elseif numel(obj.params.u_max) == obj.udim
                    b = [b; obj.params.u_max];
                else
                    error('params.u_max must be a scalar or a (udim, 1) array.');
                end
            end
            if isfield(obj.params, 'u_min')
                A = [A; -eye(obj.udim)];
                if isscalar(obj.params.u_min)
                    b = [b; -obj.params.u_min * ones(obj.udim, 1)];
                elseif numel(obj.params.u_min) == obj.udim
                    b = [b; -obj.params.u_min];
                else
                    error('params.u_min must be a scalar or a (udim, 1) array.');
                end
            end
            
            %% Cost function for the QP: minimize 0.5 * u' * H * u + f_cost' * u.
            if isfield(obj.params.weight, 'input')
                if isscalar(obj.params.weight.input)
                    weight_input = obj.params.weight.input * eye(obj.udim);
                elseif all(size(obj.params.weight.input)==[obj.udim obj.udim])
                    weight_input = obj.params.weight.input;
                else
                    error('params.weight.input must be a scalar or an (udim,udim) array.');
                end
            else
                weight_input = eye(obj.udim);
            end
            
            if verbose
                options = optimset('Display', 'notify');
            else
                options = optimset('Display', 'off');
            end
            
            H = weight_input;
            f_cost = -weight_input * u_ref;
            [u, ~, exitflag, ~] = quadprog(H, f_cost, A, b, [], [], [], [], [], options);
            if exitflag == -2
                feas = 0;
                disp('Infeasible QP. The CBF constraint conflicts with input constraints.');
            else
                feas = 1;
            end
            comp_time = toc(tstart);
        end  
    end
end