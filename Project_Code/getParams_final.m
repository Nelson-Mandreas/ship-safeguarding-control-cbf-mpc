function params = getParams_final()

%% Parameters
% Ocean current
params.Vc = 2;                    % Ocean current speed (m/s) (Horizontal)
params.betaVc = deg2rad(90);       % Ocean current direction (rad) (Horizontal direction)
params.wc = 1;                             % Vertical speed (m/s)

% Wind
params.Vw = 10;                      % Wind speed (m/s)
params.betaVw = deg2rad(0);       % Wind direction (rad)

%%

% Parameters copied from navalvessel(.)
params.rho_water     =	1014.0;	        %	water density	[kg/m^3]	
params.rho_air		=	1.225	;	    %	air density		[kg/m^3]	
params.g				=	9.81;	        %	gravity constant	[m/s^2]	
params.deg2rad 		=	pi/180;	        %	degrees to radians	
params.rad2deg 		=	180/pi;	        %	rad to degrees		
params.ms2kt			=	3600/1852;	    % 	m/s to kt 			
params.kt2ms 		=	1852/3600;	    %	kt to m/s			
params.RPM2rads		=	2*pi/60;	    %	RPM to rad/s		
params.rads2RPM		=   60/(2*pi);	    %	rad/s to RPM		
params.HP2W			=	745.700;	    %	HP to Watt
        
params.sp    =1.5;                   % span
params.A     =1.5;                   % Area
params.ar    =3;                     % aspect ratio
params.dCL   =0.054; % 1/deg         % dCL/d a_e
params.stall =23;                    % a_stall 
        
params.Lpp    =  51.5 ;                  % Length between perpendiculars [m]
params.B      =  8.6  ;                  % Beam over all  [m]
params.D	     =  2.3  ;                  % Draught [m]
        
params.disp   =  357.0;                   % Displacement  [m^3]
params.m      =  params.disp*params.rho_water;  % Mass [Kg]
params.Izz    =  47.934*10^6 ;            % Yaw Inertia
params.Ixx    =  2.3763*10^6 ;            % Roll Inertia
params.U_nom  =  8.0   ;	                 % Speed nominal [m/sec] (app 15kts) 
params.KM		=  4.47;	             %  [m] Transverse metacentre above keel
params.KB		=  1.53;	             %  [m] Transverse centre of bouancy
params.gm 		=  1.1;	                 %  [m]	Transverse Metacenter
params.bm 		=  params.KM - params.KB;
params.LCG       = 20.41 ;                % [m]	Longitudinal CG (from AP considered at the rudder stock)
params.VCG       = 3.36  ;                % [m]	Vertical  CG  above baseline
params.xG        = -3.38  ;               % Coordinate of CG from the body fixed frame adopted for the PMM test  
params.zG  	    = -(params.VCG-params.D);          % Coordinate of CG from the body fixed frame adopted for the PMM test  
params.m_xg	    = params.m * params.xG;
params.m_zg	    = params.m * params.zG;
params.Dp        = 1.6 ;                  % Propeller diameter [m]
        
params.Xudot  	= -17400.0 ;
params.Xuau     	= -1.96e+003 ;
params.Xvr    	=  0.33 * params.m ;
        
params.Yvdot = -393000 ; 
params.Ypdot = -296000 ; 
params.Yrdot = -1400000 ; 
params.Yauv  = -11800 ; 
params.Yur   =  131000 ; 
params.Yvav  = -3700 ; 
params.Yrar   =  0 ;
params.Yvar  = -794000 ; 
params.Yrav  = -182000 ; 
params.Ybauv =  10800 ; % Y_{\phi |v u|}
params.Ybaur =  251000 ; 
params.Ybuu  = -74 ; 
        
params.Kvdot =  296000 ;
params.Kpdot = -774000 ;
params.Krdot =  0 ;
params.Kauv  =  9260 ;
params.Kur   = -102000 ;
params.Kvav  =  29300 ;
params.Krar  =  0 ;
params.Kvar  =  621000 ;
params.Krav  =  142000 ;
params.Kbauv =  -8400 ;
params.Kbaur =  -196000 ;
params.Kbuu  =  -1180 ;
params.Kaup  =  -15500 ;
params.Kpap  =  -416000 ;
params.Kp    =  -500000 ;
params.Kb    =  0.776*params.m*params.g;
params.Kbbb  =  -0.325*params.m*params.g ;
        
params.Nvdot =  538000 ;
params.Npdot =  0 ;
params.Nrdot = -38.7e6;
params.Nauv  = -92000 ;
params.Naur  = -4710000 ;
params.Nvav  =  0 ;
params.Nrar  = -202000000 ;
params.Nvar  =  0 ;
params.Nrav  = -15600000 ;
params.Nbauv = -214000 ;
params.Nbuar = -4980000 ;
params.Nbuau = -8000 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obstacle radius
params.d = 40;            % Want to be atleast ... m away from the remus100
params.cbf_gamma0 = 8;

params.cbf.rate = 5;      
params.weight.slack = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

