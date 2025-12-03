classdef ctrlAffineSys < handle    
    %% Control-Affine Dynamic System Class.
    properties
        % State dimension
        xdim 
        % Control input dimension
        udim
        % Model parameters as a structure
        params 
        
        % Dynamics: xdot = f(x) + g(x) u
        f 
        g 
        
        % Function handles for CBF/CLF and its Lie derivatives.
        cbf 
        lf_cbf 
        lg_cbf 
        clf
        lf_clf
        lg_clf        
    end
    
    methods
        function obj = ctrlAffineSys(params)
            if nargin < 1
                obj.params = [];
                disp("Warning: params argument is missing.");
            else
                obj.params = params;
            end
            
            [x, f, g] = obj.defineSystem(params);
            
            % symbolic variables for the obstacle's North and East position.
            syms xo yo real
            symbolic_x_obstacle = [xo; yo];
            
            clf = obj.defineClf(params, x);
            cbf = obj.defineCbf(params, x, symbolic_x_obstacle);
            
            % Initialize the system 
            obj.initSys(x, f, g, cbf, clf, symbolic_x_obstacle);
        end
        
        
        function [x, f, g] = defineSystem(~, params)
            % Outputs:
            %   x: symbolic state vector
            %   f: drift term (symbolic)
            %   g: control vector fields (symbolic)
            x = [];
            f = [];
            g = [];
        end
        
        function clf = defineClf(obj, params, symbolic_state)
            % Returns a symbolic expression for the CLF.
            clf = [];
        end
        
        function cbf = defineCbf(obj, params, symbolic_state, symbolic_x_obstacle)
            % Returns a symbolic expression for the CBF.
            %
            % Inputs:
            %   symbolic_state: the evasive ship’s symbolic state vector
            %   symbolic_x_obstacle: obstacle position in the x-y plane
    
            x = symbolic_state;
            
            x_pos = x(7);
            y_pos = x(8);
            
            xo = symbolic_x_obstacle(1);
            yo = symbolic_x_obstacle(2);
            
            % Extract parameters
            d = params.d;                           % desired distance
            cbf_gamma0 = params.cbf_gamma0;
            
            
            distance = (x_pos - xo)^2 + (y_pos - yo)^2 - d^2; 
            
            n = x(6);  
            u = x(1);
            v = x(2);
            derivDistance = 2*(x_pos - xo)*(u*cos(n) - v*sin(n)) + 2*(y_pos - yo)*(u*sin(n) + v*cos(n));
            
            cbf = derivDistance + cbf_gamma0 * distance;
        end
        
        function initSys(obj, symbolic_x, symbolic_f, symbolic_g, symbolic_cbf, symbolic_clf, symbolic_x_obstacle)
            
            if isempty(symbolic_x) || isempty(symbolic_f) || isempty(symbolic_g)
                error('x, f, or g are empty. Define your dynamics using symbolic expressions.');
            end
            
            if ~isa(symbolic_f, 'sym')
                f_sym = sym(symbolic_f);
            else
                f_sym = symbolic_f;
            end
            if ~isa(symbolic_g, 'sym')
                g_sym = sym(symbolic_g);
            else
                g_sym = symbolic_g;
            end
            
            x = symbolic_x;
            obj.xdim = size(x, 1);
            obj.udim = size(g_sym, 2);
            
            % Dynamics.
            obj.f = matlabFunction(f_sym, 'vars', {x});
            obj.g = matlabFunction(g_sym, 'vars', {x});
            
            % CBF and its Lie derivatives
            if ~isempty(symbolic_cbf)
                dcbf = simplify(jacobian(symbolic_cbf, x));
                lf_cbf_sym = dcbf * f_sym;
                lg_cbf_sym = dcbf * g_sym;
                
                save('cbf_terms.mat', 'lf_cbf_sym', 'lg_cbf_sym');
                
                obj.cbf = matlabFunction(symbolic_cbf, 'vars', {x, symbolic_x_obstacle});
                obj.lf_cbf = matlabFunction(lf_cbf_sym, 'vars', {x, symbolic_x_obstacle});
                obj.lg_cbf = matlabFunction(lg_cbf_sym, 'vars', {x, symbolic_x_obstacle});
            end
            
            % CLF and its Lie derivatives
            if ~isempty(symbolic_clf)
                dclf = simplify(jacobian(symbolic_clf, x));
                lf_clf_sym = dclf * f_sym;
                lg_clf_sym = dclf * g_sym;
                obj.clf = matlabFunction(symbolic_clf, 'vars', {x});
                obj.lf_clf = matlabFunction(lf_clf_sym, 'vars', {x});
                obj.lg_clf = matlabFunction(lg_clf_sym, 'vars', {x});
            end
        end
    end
end

