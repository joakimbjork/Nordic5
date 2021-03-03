% PSS parameter settings.
%
% Filter time constants are set using main_PSS_controller_design.m. 
% 
% The consired case are the 240 GWs (main case) and the 110 Gws load 5000
% MW load flows.
%
% PSS's are tuned for linearized model (without PSS) using hydro FCR.

%% PSS tuned for load_data_Nordic5_5000MW_240 (Main case)
Kpss = [0.15; 0.6; 0.15; 0.15; 0.15];

T1 = [sqrt(0.01751); sqrt(0.01159); 0.05269; sqrt(0.01729); sqrt(0.03346)];
T2 = [sqrt(0.4547); sqrt(0.1986); 0.3223 ; sqrt( 0.4603);  sqrt(0.9547)];

T3 = [sqrt(0.01751); sqrt(0.01159);      0;       sqrt(0.01729); sqrt(0.03346)];
T4 = [sqrt(0.4547); sqrt(0.1986);      0;       sqrt( 0.4603);  sqrt(0.9547)];


%% PSS tuned for load_data_Nordic5_5000MW_110
Kpss = [0.15; 0.6; 0.15; 0.15; 0.15];

T1 = [sqrt(0.01086); 0.0519 ; 0.05397; sqrt(0.01146) ; sqrt(0.01663)];
T2 = [sqrt(0.1929); 0.5546; 0.3156; sqrt(0.1827); sqrt(0.2941)];

T3 = [sqrt(0.01086); 0 ; 0; sqrt(0.01146) ; sqrt(0.01663)];
T4 = [sqrt(0.1929); 0; 0; sqrt(0.1827); sqrt(0.2941)];

