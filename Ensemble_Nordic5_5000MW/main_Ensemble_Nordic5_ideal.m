% Copyright (C) 2021  Joakim Björk
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

%% Simulation
clear all
path = '';

%% Run Simulation
simulate_model = true;
run_load_flow = false;

s = tf('s');
fdcmax = 0.4; % Max allowed steady state freq deviation
d_dim = 1400;
load_damping = 400;
fn = 50;
R = (d_dim/fdcmax - load_damping); % MW/Hz to provide with FCR
R = R*fn;

% K_ideal = R/(s*7.2+1); % First order FCR
K_ideal = R*(1+6.5*s)/((1+2*s)*(1+17*s)); % Second order FCR (avoids second freqeuncy dip)
runs = [];
for i_run = 1
    if simulate_model
        if i_run == 1
                model = 'Ensemble_Nordic5_5000MW'; 
                load_data_Ensemble_Nordic5_5000MW_110
                fcr_mode = 0;        
        elseif i_run == 2
                model = 'Ensemble_Nordic5_240_5000MW';
                load_data_Ensemble_Nordic5_5000MW_240
                fcr_mode = 0;                              
        end
        
        c_slow = [60;30;10]; % Distribution of FCR-D
        c_fast = [1,2]; % Distribution of fast FCR

        c{1} = c_slow;
        c{2} = c_fast;
        
        [FCR, FCR_ideal] = ensemble_fcr_controller(...
                       K_ideal,c,GenData,fcr_mode,data_NREL,wind_turbines);        
     
        % For hydro + wind, fcr_mode = 1,   
        if fcr_mode == 1
            i_fcr=1; FCR{i_fcr} = modelred_hsv(FCR{i_fcr});
            i_fcr=2; FCR{i_fcr} = modelred_hsv(FCR{i_fcr});
            i_fcr=3; FCR{i_fcr} = modelred_hsv(FCR{i_fcr});
            i_fcr=4; FCR{i_fcr} = modelred_hsv(FCR{i_fcr});
            i_fcr=5; FCR{i_fcr} = modelred_hsv(FCR{i_fcr});           
        end

    %         [FCR, FCR_ideal] = create_fcr_controller(...
    %                        K_ideal,c,GenData,fcr_mode,data_NREL,wind_turbines);

        time_activate_wind = 0; % Let wind power models stabilize at Popt before activating them in the grid model
        time_event = time_activate_wind + 1;
        
        if run_load_flow
            LF = power_loadflow('-v2',model,'solve');                           
        end

        sim(model);   

        % Save Parameters
        runs{i_run}.t = t-time_activate_wind ; 
        % Freqeuncy
        runs{i_run}.w = w; 
        runs{i_run}.w_coi = w_coi;
        % FCR [MW]
        runs{i_run}.P_FCR = P_FCR; 
        runs{i_run}.P_FCR_tot = P_FCR_tot;
        % Freqeuncy dependent load
        runs{i_run}.PL_freq_tot = PL_freq_tot;
        if fcr_mode == 1
            % Wind data
            runs{i_run}.P_wind = P_wind;
            runs{i_run}.w_wind = w_wind;
        end
        
        M = 2*sum(Wkin);
        G = fn/(s*M+load_damping*fn);
        L = G*K_ideal/fn;
        S = inv(1+L);
        % Gcl = feedback(G,K_ideal/fn);
        % step(-S*G*d_dim-0.1)
        
        runs{i_run}.L = L;
        runs{i_run}.SG = S*G;
    end      
end

%

%% Load From File
load_model = true;
if load_model
%     load('run_01_07.mat')
    load('saved_example.mat')
    disp('Load simulation result from file')
    load_data_Ensemble_Nordic5_5000MW_110
end

%% Step
xlims=[0,40];

ylim1=[48.8,50];
ylim2=[-100,1600];

xlims2=[0,60];
ylim3=[200,850];
ylim4=[0.8,1.05];

n_down = 5; % For down sampling signals

% t = downsample(t,n_down);
% y = downsample(y,n_down);
% Ideal and Hydro
h=1/sqrt(2);
figureLatex
set(hfig,'position',[pos(1:2),3,3])

co = ([0.6,0.6,0.6;0,0,0;0,0.5,0.7;1,0,0;]);
co_gray = ([0.6,0.6,0.6]);
co_light = [1,0.8,0];

subplot(2,1,1)
for i = [2,1]
    run = runs{1,i};
    t = run.t;    
    y = run.w_coi;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    plot(t,y(:,1),'color',co(i,:)); hold on
end

for i = [2,1]
    run = runs{1,i};
    t = run.t;    
    y = run.w;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    plot(t,y(:,[5]),'color',co_light); hold on
end

for i = [2,1]
    run = runs{1,i};
    t = run.t;    
    y = run.w_coi;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    plot(t,y(:,1),'color',co(i,:)); hold on
end

yline(49,'k')
l = legend('$\omega_\mathrm{COI}$ (240~GWs)',...
           '$\omega_\mathrm{COI}$ (110~GWs)',...
           '$\omega_5$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','southeast')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left

subplot(2,1,2)
for i = [2,1]
    run = runs{1,i};
    t = run.t;    
    y = run.P_FCR_tot;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    plot(t,y,'color',co(i,:)); hold on
%     if i == 2
%         y2 = run.P_FCR(:,1:3);    
%         y2 = downsample(y2,n_down);
%         plot(t,y2,'color',co(3,:));
%         plot(t,y,'color',co(i,:));
%     end
end


yline(0,'k')
l = legend('$P_\mathrm{FCR}$ (240~GWs)',...
           '$P_\mathrm{FCR}$ (110~GWs)');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','southeast')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left

subplot(2,1,1)
ylabel('Speed [Hz]')
xlim(xlims)
ylim(ylim1)
subplot(2,1,2)
ylabel('Active power [MW]')
xlabel('Time [s]')
xlim(xlims)
ylim(ylim2)

xticks([0:10:60])

saveFigure = false;
if saveFigure
    saveas(hfig,fullfile(path, 'Nordic5_fcr'),'epsc');
end

%% Bode
omega = logspace(-3,1,100);
h = 1/sqrt(2);
figureLatex
set(hfig,'position',[pos(1:2),3,3])
co = ([0.6,0.6,0.6;0,0,0;0,0.5,0.7;1,0,0;]);

subplot(2,1,1);
wc=[];
for i = [2,1]
    L = runs{1,i}.L;
    [MAG,PHASE] = bode(L,omega);
    [~,Pm,~,wc(i)] = margin(L);
    MAG = squeeze(MAG);
    loglog(omega,MAG,'color',co(i,:)); hold on    
end
for i = [2,1]
    xline(wc(i),'color',co(i,:))
end
yline(1);

l = legend('240~GWs',...
           '110~GWs',...
           num2str(round(wc(2),2)),...
           num2str(round(wc(1),2)));
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','southwest')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
ylim([0.001,20])
xlim([1e-3,10])

ylabel('Amplitude')
xlabel('Angular frequency [rad/s]')

saveFigure = false;
if saveFigure
    saveas(hfig,fullfile(path, 'Nordic5_bode'),'epsc');
end


