
fn = 50;
load_damping = 400; % MW/Hz
S_base = 100;  % Not used?

time_event = 5;
d_fault = 1400;

ideal_FCR = false;
unlimited_Q = false; % Useful for finding a stable load flow

%% Production
%       1         2         3         4              5        
%       NO-hydro, SE-hydro, FI-hydro, SE-thermal, FI-thermal,
% Pgen = [9000;     6000;     2000;     5000;       2000]; % 110GWS case
Pgen = [18000; 12000; 2000; 11000; 7000]; % 240GWS case
Sgen = Pgen./[0.8;0.8;0.8;0.9;0.9];
Wkin = Sgen.*[3;3;3;6;6]; % Gives kinetic energy =  110 GWs 
sum(Wkin)

%% Loads
n = 5;
BUS = (1:n)';
V =[400;400;400;400;400]; 

%            1          2          3         4           5        
%            NO-hydro,  SE-hydro,  FI-hydro, SE-thermal, FI-thermal
PL = Pgen + [0;         -5500;         0;      5000;          500];
QL = PL*0.0; 
QSHUNT = PL*0.03+ [0;         1000;         0;        2000;   500];  % Adjust these to fix voltage/reactive power issues, if necessary
PSHUNT = ones(n,1)*0;
PL_freq = PL/sum(PL)*load_damping; % Freqeuncy dependent load proportionally distributed

np =  ones(n,1)*2; % Dynamic three phase load should be able to deal with non-impedance loads. But, this does not seem to work that well for larger simulations.
nq =  ones(n,1)*2;
     
loads = table(BUS,V,PL,QL,PSHUNT,QSHUNT,PL_freq,np,nq);

%% Lines
from = [1;2;2;3;1];
to   = [4;4;3;5;2];
V = ones(size(from))*400;
km = [250;100;100;150;400];
% km = [220; 100; 100; 150; 300];
lines = table(from,to,V,km);

%% Transformers (not used)
% from = [1;4;5];
% to = [11;41;51];
% V1 = [400;400;400];
% V2 = [130;130;130];
% L2 = [0.01; 0.01; 0.01]; 
% ratio = [1; 1; 1];
% Sbase = [PL(1);PL(4);PL(5)]*1;
% 
% trafos = table(from,to,V1,V2,L2,ratio,Sbase);

%% Exciter
n = 5; 
BUS = (1:n)';
K =         [300;  300;  300;  300; 300];
MIN =       [0;   0;   0;   0;   0];
MAX =       [5;   5;   5;   5;   5];

ExcData = table(BUS,K,MIN,MAX);

%% PSS (prio 2, 3, 5, 4, 1)
Kpss = [0.15; 0.6; 0.15; 0.15; 0.15];

T1 = [sqrt(0.01751); sqrt(0.01159); 0.05269; sqrt(0.01729); sqrt(0.03346)];
T2 = [sqrt(0.4547); sqrt(0.1986); 0.3223 ; sqrt( 0.4603);  sqrt(0.9547)];

T3 = [sqrt(0.01751); sqrt(0.01159);      0;       sqrt(0.01729); sqrt(0.03346)];
T4 = [sqrt(0.4547); sqrt(0.1986);      0;       sqrt( 0.4603);  sqrt(0.9547)];

Tw = ones(n,1)*4.5; % Washout
Tf = ones(n,1)*0.01; % Low-pass

VMIN = ones(n,1)*0.05;
VMAX = ones(n,1)*0.05;

PSSData = table(BUS,Tw,Tf,Kpss,T1,T2,T3,T4,VMIN,VMAX);

%% Generator
% type: 1 -> Salient pole rotor (hydro), 2 -> Round rotor (thermal)
%    BUS type Td0p  Td0pp  Tq0p Tq0pp  H  Xd    Xq    Xdp   Xqp   Xdpp  Xqpp  XL     Xt      
M = [1   1     5    0.05   NaN  0.1    3  1.25  0.85  0.4   NaN   0.35  0.35  0.15   0.15; 
     2   1     5    0.05   NaN  0.1    3  1.25  0.85  0.4   NaN   0.35  0.35  0.15   0.15;  
     3   1     5    0.05   NaN  0.1    3  1.25  0.85  0.4   NaN   0.35  0.35  0.15   0.15;  
     4   2     7    0.05   1.5  0.05   6  2.35  2.15  0.45  0.55  0.35  0.35  0.15   0.15;
     5   2     7    0.05   1.5  0.05   6  2.35  2.15  0.45  0.55  0.35  0.35  0.15   0.15];  
% BUS = M(:,1); 
type = M(:,2); 
Td0p = M(:,3); Td0pp = M(:,4); Tq0p = M(:,5); Tq0pp = M(:,6);
H = M(:,7);
Xd = M(:,8); Xq = M(:,9); Xdp = M(:,10); Xqp = M(:,11); Xdpp = M(:,12); Xqpp = M(:,13);
XL = M(:,14); Xt = M(:,14);

%    BSKV  V0  PGEN    QGEN       QMAX        QMIN          MBASE       
M = [400  1    Pgen(1) Pgen(1)/4  Sgen(1)/3  -Sgen(1)/6   Sgen(1);
     400  1    Pgen(2) Pgen(2)/4  Sgen(2)/3  -Sgen(2)/6   Sgen(2) ;
     400  1    Pgen(3) Pgen(3)/4  Sgen(3)/3  -Sgen(3)/6   Sgen(3) ;
     400  1    Pgen(4) Pgen(4)/4  Sgen(4)/3  -Sgen(4)/6   Sgen(4);
     400  1    Pgen(5) Pgen(5)/4  Sgen(5)/3  -Sgen(5)/6   Sgen(5)];
  
BSKV = M(:,1); V0 = M(:,2) ;
PGEN = M(:,3); QGEN = M(:,4); QMAX = M(:,5); QMIN = M(:,6) ;
MBASE = M(:,7);
%      Ta     vgMAX    gMAX    gMIN   TW     DREF   GenType  GovType           
M = [ 0.2     0.1     0.95      0     0.7    1      1        1 ;   
      0.2     0.1     0.95      0     1.4    1      1        1 ; 
      0.2     0.1     0.95      0     1.4    1      1        1 ; 
      NaN     NaN      NaN    NaN     NaN    NaN    2        0;
      NaN     NaN      NaN    NaN     NaN    NaN    2        0] ; 
Ta = M(:,1);
vgMAX = M(:,2); gMAX = M(:,3);  gMIN = M(:,4);
TW = M(:,5);
DREF = M(:,6); % unused
GenType = M(:,7);
GovType = M(:,8);

if unlimited_Q % Useful for finding a stable load flow
    QMAX = QMAX*1e6;
    QMIN = QMIN*1e6;
end

GenData = table(BUS, type, Td0p,  Td0pp,  Tq0p,   Tq0pp,...
                    H, Xd, Xq, Xdp, Xqp, Xdpp, Xqpp, XL, Xt,...
                    BSKV, V0, PGEN, QGEN, QMAX, QMIN, MBASE,...
                    Ta, vgMAX, gMAX, gMIN, TW, DREF, GenType, GovType); 
                
%%
fcr_controller_design
                