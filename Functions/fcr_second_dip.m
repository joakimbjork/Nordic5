% Script to tune FCR design target
Wsum = 1.1042*10^(5);
M = 2*Wsum; % 110 MWs
% M = 240e3;

D = 400; % MW/Hz

s = tf('s');
fdcmax = 0.4; % Max allowed steady state freq deviation
d_dim = 1400;
fn = 50;
R = (d_dim/fdcmax - D); % MW/Hz to provide with FCR
K = R/(s*7.2+1); % ideal FCR MW/Hz

G = fn/(s*M+D*fn);

L = G*K;
S = inv(1+L);
SG = S*G;

Gcl = feedback(G,K);

% step(-Gcl*d_dim-0.1)

% Lead filter'
K2 = R/(s*17+1); % ideal FCR MW/Hz
L2 = G*K2;
[~,Pm,~,wc] = margin(L2);
d_phase = (75-Pm)*pi/180;
% d_phase = (50-Pm)*pi/180;
ang = d_phase*180/pi
F = leadFilter(d_phase,wc*1.3);

% b = wc/10;
% a = b/k;
% Fcomp = (s+a)/(s+b)*k;

K2 = K2*F;
Kmod = R*(1+6.5*s)/((1+2*s)*(1+17*s));
L2 = G*K2;
Lmod = G*Kmod;

% margin(L2)

Gcl2 = feedback(G,K2);
Gcl3 = feedback(G,Kmod);
figure()
step(-Gcl*d_dim-0.1,-Gcl2*d_dim-0.1,-Gcl3*d_dim-0.1)
legend('First order FCR', 'First order FCR + Lead filter', 'Modified')
% step(-Gcl*d_dim-0.1,-Gcl2*d_dim-0.1,60)
% legend('First order FCR', 'First order FCR + Lead filter')


%%
K2 = zpk(K2);
K2_=zpk(K2.Z{1}, K2.P{1},K2.K,'DisplayFormat','time constant')