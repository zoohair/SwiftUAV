%%
%Zouhair Mahboubi, 11-11-2009
%Stability derivatives approximations for Swift
%Reference Etkins, Dynamics of Flight: stability and control

%%
%Defining flight conditions
%All units are SI unless otherwise specified
%angles in radians and frequencies in radians/sec

g = 9.81; %m/s^2
rho = 1.22; %density at sea-level
mu_d = 2e-5; %dynamic viscocity at sea-level

M = 190; %Mass
W = M*g; %Weight ~418lbs
%Moments of Inertia use empty weight and assume that the major contribution
%comes from the wing. These are likely best computed from the CAD model...
Ix=882.;Iy=97.4; Iz= 990.; 
Iyz = 0; Ixy = 0; Ixz = 0; %Ixz should not be zero, but not needed for now


Sref = 12.5; %reference area ~ 134.5 ft^2
bref = 12.8; %reference span ~ 42 ft
cavg = Sref/bref; %average chord ~ 3.2 ft
cbar = 1.07; %mean aerodynamic chord, 3.5 ft (given by John)
AR = bref^2/Sref; %Aspect ratio

U = 50*1e3/3600; %Cruising at 50km/h
L_D = 20; %L over D


%Non-dimensional quantities
Re = rho*U*cbar/mu_d; %Reynolds number
q = .5*rho*U^2;
CL = W / (q*Sref);
CD = CL/L_D;

mu = M/(rho*Sref*cbar/2); %non-dimensionalized mass
t_star = cbar/(2*U);
Iy_bar = Iy / (rho*Sref*(cbar/2)^3);

%Estimated stability derivatives 
CL_alpha = 5.09; %Given by John Melton, Linair Computes 5.98
Cz_alpha = - (CL_alpha + CD); % Eq 5.2,4

Cm_alpha = -.406; %Given by John Melton, LinAir computes -1.347
Cm_alphadot = 0; %Given by John Melton
Cm_q = -0.0649 ; %Given by John Melton, LinAir computes -5.65


%%
%Lancaster approximation of phugoid period
%equations 6.3,5 and 6.3,12 of Etkins
T_phugoid = sqrt(2)*pi*U/g;
w_phugoid = 2*pi/T_phugoid;
zeta_phugoid = CD/CL;

%short period approximation
C = -(1/t_star^2/Iy_bar) * (Cm_alpha - Cm_q*Cz_alpha/2/mu );
B = -(1/t_star)*(Cz_alpha/2/mu+ (Cm_q + Cm_alphadot)/Iy_bar);
w_sp = sqrt(C);
zeta_sp = B/(2*w_sp);

fprintf('Phugoid frequency %f rad/s , damping %f \n', w_phugoid,zeta_phugoid);
fprintf('SP frequency %f rad/s , damping %f \n', w_sp,zeta_sp);

%results 
%Phugoid frequency 0.998887 rad/s , damping 0.050000 
%SP frequency 2.58 rad/s , damping 0.56

%Linair computes:
% Phugoid     : -0.0836692	+/- 0.5058501j => wn = 0.51 rad/s zeta = 0.16
% Short period: -10.46046 +/- 6.9383464j   => wn = 12.55 rad/s zeta = 0.83
%It seems like things are more damped according to LinAir. Not to be
%trusted for the moment being...
%%
s = 4; %zeta*wn
zw = 4.6; %sqrt(1-zeta^2)*wn

lin_wn = sqrt(s^2+zw^2)
lin_zeta = s/lin_wn

% Mode    	Real(s) 	Imag(s) 
% Short Period	-04.073607	4.6035423
% Short Period	-04.073607	-04.603542
% Phugoid 	-0.0031249	0.6405237
% Phugoid 	-0.0031249	-0.6405237
% Dutch Roll	-0.3504365	1.4492235
% Dutch Roll	-0.3504365	-01.449223
% Spiral  	0.0521827	0.0000000
% Roll    	-08.931337	0.0000000


