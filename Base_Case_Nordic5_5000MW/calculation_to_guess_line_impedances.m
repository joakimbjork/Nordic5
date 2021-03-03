%% Estimate of line impedance (OBS 240 GWS case)
% All lines and bachines are at 400 kV base power
Sbase = 100*1e6; % MW
Vbase = 400*1e3;
Sr = [22500 15000 2500 12222 7778]*1e6; % MW
Xg = 0.55; % Xdp + XL assumed to be the same for all generators
H = [3, 3, 3, 6, 6];
w0 = 50*2*pi;

L_per_km = 0.9337e-3; % Henry per km

% 0.33Hz Mode (NO_1+SE_4) -- (FI_3+FI_5)
mode = 0.33*2*pi;
idx1 = [1,4]; idx2 = [3,5];

S1 = sum(Sr(idx1)); S2 = sum(Sr(idx2));

M1 = 2*sum(H(idx1).*Sr(idx1))/w0;
M2 = 2*sum(H(idx2).*Sr(idx2))/w0;

M_1 = 4.5*S1/w0;

X1 = Xg/S1 ;X2 = Xg/S2;

Xline = ( (M1+M2)/(mode^2*M1*M2) - X1 - X2) * Vbase^2; % Ohm
Lline = Xline/w0;

dist_1 = Lline/L_per_km

% 0.48Hz Mode (NO_1) -- (SE_4)
mode = 0.48*2*pi;
idx1 = [1]; idx2 = [4];

S1 = sum(Sr(idx1)); S2 = sum(Sr(idx2));

M1 = 2*sum(H(idx1).*Sr(idx1))/w0;
M2 = 2*sum(H(idx2).*Sr(idx2))/w0;

X1 = Xg/S1 ;X2 = Xg/S2;

Xline = ( (M1+M2)/(mode^2*M1*M2) - X1 - X2) * Vbase^2; % Ohm
Lline = Xline/w0;

dist_2 = Lline/L_per_km