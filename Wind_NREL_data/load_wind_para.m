function [data_NREL] = load_wind_para(model_name,v_wind,varargin)
% model_name = 'model_data.mat';
% v_wind = 10; % wind speed m/s

if nargin == 2
    kstab_multiplier = 2;
else
    kstab_multiplier = varargin{1};
end

data_NREL=load(model_name);
data_NREL = data_NREL.l;
env=data_NREL.env;
wt=data_NREL.wt;

% Tgen = wt.gen.timeconstant;
n = wt.ctrl.gen.effeciency;
Pb = data_NREL.public.rated; % 5 MW
Pb_m = Pb/n;

data_NREL.env.v_wind = v_wind;
data_NREL.Pb = Pb;

J = data_NREL.wt.gen.N^2*data_NREL.wt.gen.inertia + data_NREL.wt.rotor.inertia;
data_NREL.ctrl.J = J;

rho =env.rho;
r = wt.rotor.radius;

data_NREL.ctrl.v_wind = v_wind;

%% Stabilizing Controller
cp = wt.cp.table;
lambda = wt.cp.tsr;
cp = cp(1,:);
[c_opt,idx] = max(cp);
lambda_opt = lambda(idx);

x = lambda/lambda_opt;
N = length(x);
x_min = 0.8; % Allowed deviation from optimal speed
for i = 2:N-1
    i1 = i-1;
    i2 = i+1;
    dx = x(i2)-x(i1);
    d_cp = cp(i2)-cp(i1);
    d_cp_x(i-1) = d_cp/dx;
end
[~,idx] = min(abs(x(2:N-1)-x_min)); % Range 2:N-1 since we take values before and after (could be changed)
d_c_max =  d_cp_x(idx); % cp derivative at x = x_min
c_x_min =  cp(idx+1); % cp value at x = x_min

k_stab = d_c_max*kstab_multiplier;
w_opt = min(v_wind*lambda_opt/r, data_NREL.wt.rotor.ratedspeed);
data_NREL.ctrl.k_stab = k_stab;


data_NREL.ctrl.lambda_opt = lambda_opt;
data_NREL.ctrl.c_opt = c_opt;


% Linearized model
C = rho/2 * pi*r^2/J * r^2/lambda_opt^2;
data_NREL.ctrl.C = C;
z_wind = v_wind * C * d_c_max / x_min;
data_NREL.ctrl.z_wind = z_wind;

a = v_wind * C * k_stab / x_min;
data_NREL.ctrl.p_wind = a-z_wind;

% Speed
data_NREL.ctrl.w_opt = w_opt;
data_NREL.w_initial  = w_opt*data_NREL.wt.gen.N; % Generator initial speed

% Power
P_wind = v_wind^3*(rho*pi*r^2)/2;
P_opt = min(P_wind*c_opt,Pb_m);
data_NREL.torque_initial  = P_opt/data_NREL.w_initial;
data_NREL.ctrl.Pm_opt = P_opt;
data_NREL.ctrl.Pe_opt = P_opt*n;

% Low speed safety controller
data_NREL.ctrl.x_min = x_min;
data_NREL.ctrl.low_speed_safety = 1;
data_NREL.ctrl.c_x_min = c_x_min;
