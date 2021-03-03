% Controllers in the DVPP example

%        0.13889 (s+5) (s+1.25) (s+0.04786)
% K1 =  ------------------------------------   [pu power/ pu speed]
%        (s+0.4081) (s+0.1389) (s+0.0733)

T = 7.2;
R = 20*50; % 20MW/Hz

Ty = 0.2;

Sbase = 100; % Base power of gen

k1 = 0.13889 * 100 * 50


%            39.007 s (s+0.04786)
% K1 =  --------------------------------
%       (s+0.4081) (s+0.1389) (s+0.0733)

Popt = 7.1212; % Base power of wind turbine

k2 = 39.007 * Popt * 50