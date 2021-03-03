fn = 50;
load_damping = 400; % MW/Hz
S_base = 100;  % Not used?

time_activate_wind = 40; % Let wind power models stabilize at Popt before activating them in the grid model
time_event = 45;

%% Wind
% v_wind = 8;
data_NREL = load_wind_para('model_data.mat',v_wind,kstab_multiplier);

% ww = 1.3091;
% data_NREL.wt.rotor.ratedspeed = ww;
% data_NREL.wt.gen.ratedspeed = ww*data_NREL.wt.gen.N;
% data_NREL.wt.ctrl.gen.rated = ww*data_NREL.wt.gen.N;

Pwind_rated = data_NREL.Pb*1e-6; % Rated power per wind turbine (default 5 MW) 
Pe_opt = data_NREL.ctrl.Pe_opt*1e-6; %Optimal power at the current wind speed.
w_opt = data_NREL.w_initial/data_NREL.wt.gen.N;
w_rated = data_NREL.wt.rotor.ratedspeed;

%    BUS   BSKV   MBASE  V_WIND ROT_OPT ROT_RATED
M = [6     400    25    v_wind  w_opt   w_rated]; % OBS MAKE SURE THAT LOAD FLOW IS UPDATED AS WELL...
BUS = M(:,1);
BSKV = M(:,2);
MBASE = M(:,3);  POPT = Pe_opt*MBASE/Pwind_rated; % [MW]
V_WIND = M(:,4); ROT_OPT = M(:,5); ROT_RATED = M(:,6);

POPT = min(POPT,MBASE);
wind_turbines = table(BUS,BSKV,MBASE,POPT,V_WIND,ROT_OPT, ROT_RATED);

% Linearized wind 
s = tf('s');
% z_wind = 0.05;
z_wind = data_NREL.ctrl.z_wind;
a =v_wind * data_NREL.ctrl.C*(data_NREL.ctrl.k_stab) / data_NREL.ctrl.x_min;
p_wind = a-z_wind;
Gwind = (s-z_wind)/(s+p_wind);
Gwind = tf(Gwind);

J = data_NREL.ctrl.J;
Gwind_speed = -Pe_opt*1e6/(J*w_opt*w_opt) * 1/(s+p_wind);

FCR_wind = tf(1);

% v_lin = 8.5 ; % approx = z_wind/0.006
% r = data_NREL.wt.rotor.radius;
% P_lin = v_lin^3*(data_NREL.env.rho*pi*r^2)/2*data_NREL.ctrl.c_opt*1e-6;
% w_lin = v_lin*data_NREL.ctrl.lambda_opt/r;




%% Lines
from = [1];
to   = [1];
V = ones(size(from))*400;
km = [20];
lines = table(from,to,V,km);

