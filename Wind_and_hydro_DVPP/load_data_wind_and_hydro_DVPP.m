fn = 50;
load_damping = 400; % MW/Hz
S_base = 100;  % Not used?

time_activate_wind = 40; % Let wind power models stabilize at Popt before activating them in the grid model
time_event = 45;

%% Wind
v_wind = 8;
data_NREL = load_wind_para('model_data.mat',v_wind);

Pwind_rated = data_NREL.Pb*1e-6; % Rated power per wind turbine (default 5 MW) 
Pe_opt = data_NREL.ctrl.Pe_opt*1e-6; %Optimal power at the current wind speed.
w_opt = data_NREL.w_initial/data_NREL.wt.gen.N;
w_rated = data_NREL.wt.rotor.ratedspeed;

%    BUS   BSKV   MBASE  V_WIND ROT_OPT ROT_RATED
M = [6     400    20    v_wind  w_opt   w_rated]; % OBS MAKE SURE THAT LOAD FLOW IS UPDATED AS WELL...
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


%% Power system
%% Production
%             
Pgen = [80; 80*100]; 
Sgen = Pgen./[0.8;0.9];% Assume hydro plant operate at 80% of rated power
Pgen(1) = 80;
Wkin = Sgen.*[3;6]; 
Wkin_sum = sum(Wkin)

%% Loads
n = 2;
BUS = (1:n)';
V =[400;400;]; 

PL =        [0; Pgen(2)];
QL =         PL*0.0; 
QSHUNT =    [0; 0];  % Adjust these to fix voltage/reactive power issues, if necessary
PSHUNT = ones(n,1)*0;
PL_freq = PL/sum(PL)*load_damping; % Freqeuncy dependent load proportionally distributed

np =  ones(n,1)*2; % Dynamic three phase load should be able to deal with non-impedance loads. But, this does not seem to work that well for larger simulations.
nq =  ones(n,1)*2;
     
loads = table(BUS,V,PL,QL,PSHUNT,QSHUNT,PL_freq,np,nq);

%% Lines
from = [1];
to   = [1];
V = ones(size(from))*400;
km = [20];
lines = table(from,to,V,km);

%% Exciter
n = 2; 
BUS = (1:n)';
K =         [300; 300];
MIN =       [0;  0 ];
MAX =       [5;  5 ];

ExcData = table(BUS,K,MIN,MAX);

%% PSS 
Kpss = [1; 1];

% T1 = [0.07358; 0.04903; 0.07362; 0.1043; 0.0741];
% T2 = [0.9843; 0.6239; 0.2461 ; 0.4696 ; 0.9773];
% 
% T3 = [0; 0; 0; 0; 0];
% T4 = [0; 0; 0; 0; 0];
T1 = [0.1043; 0.1043];
T2 = [0.4696; .4696];
T3 = [0; 0];
T4 = [0; 0];

Tw = ones(n,1)*4.5; % Washout
Tf = ones(n,1)*0.01; % Low-pass

VMIN = ones(n,1)*0.05;
VMAX = ones(n,1)*0.05;

PSSData = table(BUS,Tw,Tf,Kpss,T1,T2,T3,T4,VMIN,VMAX);

%% Generator
% type: 1 -> Salient pole rotor (hydro), 2 -> Round rotor (thermal)
%    BUS type Td0p  Td0pp  Tq0p Tq0pp  H  Xd    Xq    Xdp   Xqp   Xdpp  Xqpp  XL     Xt      
M = [1   1     5    0.05   NaN  0.1    3  1.25  0.85  0.4   NaN   0.35  0.35  0.15   0.15;
     2   2     7    0.05   1.5  0.05   6  2.35  2.15  0.45  0.55  0.35  0.35  0.15   0.15 ];  
% BUS = M(:,1); 
type = M(:,2); 
Td0p = M(:,3); Td0pp = M(:,4); Tq0p = M(:,5); Tq0pp = M(:,6);
H = M(:,7);
Xd = M(:,8); Xq = M(:,9); Xdp = M(:,10); Xqp = M(:,11); Xdpp = M(:,12); Xqpp = M(:,13);
XL = M(:,14); Xt = M(:,14);

%    BSKV  V0  PGEN    QGEN       QMAX        QMIN          MBASE       
M = [400  1    Pgen(1) Pgen(1)/4  Sgen(1)/3  -Sgen(1)/6   Sgen(1);
     400  1    Pgen(2) Pgen(2)/4  Sgen(2)/3  -Sgen(2)/6   Sgen(2);];
  
BSKV = M(:,1); V0 = M(:,2) ;
PGEN = M(:,3); QGEN = M(:,4); QMAX = M(:,5); QMIN = M(:,6) ;
MBASE = M(:,7);
%      Ta     vgMAX    gMAX    gMIN   TW   DREF   GenType  GovType           
M = [ 0.2     0.1     0.95      0     2    1      1        1; % <------------- Change TW for more illustrative simulation
      NaN     NaN      NaN    NaN    NaN   NaN    2        0];
Ta = M(:,1);
vgMAX = M(:,2); gMAX = M(:,3);  gMIN = M(:,4);
TW = M(:,5);
DREF = M(:,6); % unused
GenType = M(:,7);
GovType = M(:,8);

GenData = table(BUS, type, Td0p,  Td0pp,  Tq0p,   Tq0pp,...
                    H, Xd, Xq, Xdp, Xqp, Xdpp, Xqpp, XL, Xt,...
                    BSKV, V0, PGEN, QGEN, QMAX, QMIN, MBASE,...
                    Ta, vgMAX, gMAX, gMIN, TW, DREF, GenType, GovType); 

%% FCR controller design
s = tf('s');

% fdcmax = 0.4; % Max allowed steady state freq deviation
% d_dim = 1400;
% R = (d_dim/fdcmax - load_damping); % MW/Hz to provide with FCR
% R = R*fn;

R = 20*fn;

% K_ideal = R/(s*7.2+1); % Ideal FCR
K_ideal = R*(1+6.5*s)/((1+2*s)*(1+17*s)); % Second order FCR (avoids second freqeuncy dip)

% c_slow = [4230;2500;800]; % Distribution of FCR-N in Saarinen's paper
%    hydro, 
c_slow = [1]; % Distribution of steady-state FCR
c_tot = sum(c_slow);
c_slow = c_slow/c_tot;

c_fast = [1]; % Distribution of fast FCR
c_tot = sum(c_fast);
c_fast = c_fast/c_tot;

n_hydro = length(c_slow);
n_wind = length(c_fast);

Hhydro = []; c_hydro = [];
Hwind = []; c_wind = [];
for i = 1:(n_hydro+n_wind)
    if i == 1 % Hydro turbine
        % Assumes that hydro turbines are listed first
        Ty = GenData.Ta(i);
        Tw = GenData.TW(i);
        g0 = GenData.PGEN(i)/GenData.MBASE(i);
        Tz_hydro = Tw*g0; % 1/z_hydro

        Gservo = 1/(s*Ty+1); 
        Ghydro = (1 - s*Tz_hydro)/(1 + s*Tz_hydro/2);
        Hhydro{i} = Ghydro*Gservo*Sgen(i);
        
        % Participation factor need to include the RHP zero
        w_x = 1; % Add a pole at w_x times the freqeuncy as the RHP zero. 
        c_hydro{i} = c_slow(i)*(1 - s*Tz_hydro)/(1 + s*Tz_hydro/w_x); 
        
        if i == n_hydro
            c_slow_sum = tf(0);
            for ii = 1:n_hydro
                c_slow_sum = c_slow_sum + c_hydro{ii};
            end
            c_fast_sum = 1- c_slow_sum;
        end
       
    elseif i == 2 % Wind turbine
        ii = i - n_hydro;
        z_wind = data_NREL.ctrl.z_wind;
        a =v_wind * data_NREL.ctrl.C*(data_NREL.ctrl.k_stab) / data_NREL.ctrl.x_min;
        p_wind = a-z_wind;
        Gwind = (s-z_wind)/(s+p_wind);
        Hwind{ii} = Gwind*wind_turbines.POPT(ii);
        
        % Participation factor need to include the RHP zero
        c_wind{ii} = c_fast_sum*c_fast(ii)*Gwind;
    end       
end

c_fast_sum = tf(0);
for i = 1:n_wind
    c_fast_sum = c_fast_sum + c_wind{i};
end

% Normalize so that c_hydro+c_wind = 1
c_sum = minreal(c_slow_sum+c_fast_sum);
c_result = tf(0);
for i = 1:n_hydro
    c_hydro{i} = c_hydro{i}/c_sum;
    c_result = c_result + c_hydro{i};
end
for i = 1:n_wind
    c_wind{i} = c_wind{i}/c_sum;
    c_result = c_result + c_wind{i};
end

c_result = minreal(c_result);

FCR = [];
for i = 1:(n_hydro+n_wind)
    if i == 1 % Hydro turbine
        ii = i;
        FCR{i} = minreal(K_ideal*c_hydro{ii}/Hhydro{ii});    
    elseif i == 2 % Wind turbine
        ii = i - n_hydro;
        FCR{i} = minreal(K_ideal*c_wind{ii}/Hwind{ii});  
    end
end

%% Linear sim
% F1 = FCR{1};
% H1 = Hhydro{1};
% 
% mult1 = 1/(Sgen(1)*fn);
% 
% step(F1*H1*mult1*10, K_ideal/50/Sgen(1)*10);

