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
simulate_model = false;

s = tf('s');
fdcmax = 0.4; % Max allowed steady state freq deviation
d_dim = 1400;
load_damping = 400;
fn = 50;
R = (d_dim/fdcmax - load_damping); % MW/Hz to provide with FCR
R = R*fn;

% K_ideal = R/(s*7.2+1); % First order FCR
K_ideal = R*(1+6.5*s)/((1+2*s)*(1+17*s)); % Second order FCR (avoids second freqeuncy dip)

if simulate_model
%     runs = [];
    load_data_Wind_Nordic5_5000MW_110
        
    for i_run = 4
        if i_run == 1
            model = 'Nordic5_5000MW'; 
            fcr_mode = 0;    
            
        elseif i_run == 2
            model = 'Nordic5_5000MW'; 
            fcr_mode = 2;                               
        elseif i_run == 3 
            model = 'Wind_Nordic5_5000MW';
            fcr_mode = 1;     
        elseif i_run == 4
            model = 'Wind_Nordic5_5000MW';
            load_data_Wind_Nordic5_5000MW_110_sens
            fcr_mode = 1;   
        end   
        c_slow = [60;30;10]; % Distribution of FCR-D
        c_fast = [1,2]; % Distribution of fast FCR

        c{1} = c_slow;
        c{2} = c_fast;
        
        [FCR, FCR_ideal] = ensemble_fcr_controller(...
                       K_ideal,c,GenData,fcr_mode,data_NREL,wind_turbines);
        
        % For hydro + wind, fcr_mode = 1,   
        if fcr_mode == 1;
            i_fcr=1;
            FCR{i_fcr} = modelred_hsv(FCR{i_fcr},4);
            i_fcr=2;
            FCR{i_fcr} = modelred_hsv(FCR{i_fcr},4);
            i_fcr=3;
            FCR{i_fcr} = modelred_hsv(FCR{i_fcr},4);
            i_fcr=4;
            FCR{i_fcr} = modelred_hsv(FCR{i_fcr},4);
            i_fcr=5;
            FCR{i_fcr} = modelred_hsv(FCR{i_fcr},4);           
        end
             
        sim(model);   
        
        % Save Parameters
        runs{i_run}.t = t-time_activate_wind-4 ; 
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
            runs{i_run}.P_wind_set = P_wind_set;
        end
    end      
end


%% Load From File
load_model = true;
if load_model
%     load('run_01_07.mat')
    load('saved_example.mat')
    disp('Load simulation result from file')
    load_data_Wind_Nordic5_5000MW_110
end


%%

xlims=[0,40];

ylim1=[48.5,50];
ylim2=[-300,2000];

xlims2=[0,60];
ylim3=[200,850];
ylim4=[0.8,1.05];

n_down = 5; % For down sampling signals
% t = downsample(t,n_down);
% y = downsample(y,n_down);
%% Ideal and Hydro
h=1;
figureLatex

co = ([0.6,0.6,0.6;0,0,0;0,0.5,0.7;1,0,0;]);
co_gray = ([0.6,0.6,0.6]);
co_light = [1,0.8,0];

subplot(2,1,1)
for i = 1:2
    run = runs{1,i};
    t = run.t;    
    y = run.w_coi;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    plot(t,y(:,1),'color',co(i,:)); hold on
end

for i = 1:2
    run = runs{1,i};
    t = run.t;    
    y = run.w;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    plot(t,y(:,[5]),'color',co_light); hold on
end



for i = 1:2
    run = runs{1,i};
    t = run.t;    
    y = run.w_coi;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    plot(t,y(:,1),'color',co(i,:)); hold on
end
yline(49,'k')
l = legend('$\omega_\mathrm{COI}$ (ideal)',...
           '$\omega_\mathrm{COI}$ (hydro)',...
           '$\omega_5$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','southeast')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left

subplot(2,1,2)
for i = 1:2
    run = runs{1,i};
    t = run.t;    
    y = run.P_FCR_tot;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    plot(t,y,'color',co(i,:)); hold on
    if i == 2
        y2 = run.P_FCR(:,1:3);    
        y2 = downsample(y2,n_down);
        plot(t,y2,'color',co(3,:));
        plot(t,y,'color',co(i,:));
    end
end


yline(0,'k')
l = legend('$P_\mathrm{ideal}$ ',...
           '$P_\mathrm{hydro}$',...
           'Bus 1-3');
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

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'Nordic5_hydro'),'epsc');
end

%% Ideal and Hydro ZOOM
h=0.99;
figureLatex

co = ([0.6,0.6,0.6;0,0,0;0,0.5,0.7;1,0,0;]);
co_gray = ([0.6,0.6,0.6]);

run = runs{1,2};
for i = 1:3    
    t = run.t;    
    y = run.P_FCR;    
%     t = downsample(t,n_down);
%     y = downsample(y,n_down);
    plot(t,y(:,i),'color',co(i,:)); hold on
end


yline(0,'k')
l = legend('Bus 1',...
           'Bus 2',...
           'Bus 3');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','southeast')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left


ylabel('Active power [MW]')
xlabel('Time [s]')
xlim([0,10])
ylim([-50,100])

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'Nordic5_hydro_zoom'),'epsc');
end

%% Ideal and Hydro + Wind
h=1;
figureLatex

co = ([0.6,0.6,0.6;0,0.5,0.7;0,0,0;1,0,0;]);
co_gray = ([0.6,0.6,0.6]);
co_light = [1,0.8,0];

subplot(2,1,1)
for i = [1,3]
    run = runs{1,i};
    t = run.t;    
    y = run.w_coi;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    if i == 3
        plot(t,y,'--','color',co(i,:)); hold on
    else
        plot(t,y,'color',co(i,:)); hold on
    end 
end

for i = [1,3]
    run = runs{1,i};
    t = run.t;    
    y = run.w;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    plot(t,y(:,[5]),'color',co_light); hold on
end

for i = [1,3]
    run = runs{1,i};
    t = run.t;    
    y = run.w_coi;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    if i == 3
        plot(t,y,'--','color',co(i,:)); hold on
    else
        plot(t,y,'color',co(i,:)); hold on
    end 
end

yline(49,'k')
l = legend('$\omega_\mathrm{COI}$ (ideal)',...
           '$\omega_\mathrm{COI}$ (hydro+wind)',...
           '$\omega_5$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','southeast')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left


subplot(2,1,2)
for i = [1,3]
    run = runs{1,i};
    t = run.t;    
    y = run.P_FCR_tot;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    if i == 3
        plot(t,y,'--','color',co(i,:)); hold on
    else
        plot(t,y,'color',co(i,:)); hold on
    end 
    if i == 3
        y2 = sum(run.P_FCR(:,1:3),2);    
        y2 = downsample(y2,n_down);
        y3 = sum(run.P_FCR(:,4:5),2);    
        y3 = downsample(y3,n_down);
        plot(t,y2,'color',[0,0.5,0.7]);
        plot(t,y3,'color',[0,0,0]);
    end
end

yline(0,'k')
l = legend('$P_\mathrm{ideal}$ ',...
           '$P_\mathrm{hydro+wind}$',...
           '$P_\mathrm{hydro}$',...
           '$P_\mathrm{wind}$');
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

if saveFigure
    saveas(hfig,fullfile(path, 'Nordic5_hydro_and_wind'),'epsc');
end

%% Plot wind response
h=1;
figureLatex
       
co = ([;0,0,0;0,0.5,0.7]);

subplot(2,1,1)
for i = 3
    run = runs{1,i};
    t = run.t;    
    y = run.P_wind;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    plot(t,y(:,1),'color',co(1,:)); hold on
    plot(t,y(:,2),'color',co(2,:));
end


l = legend('Bus 2' , 'Bus 4');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','northeast')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left

subplot(2,1,2)
for i = 3
    run = runs{1,i};
    t = run.t;
    y = run.w_wind;
    n = 10;
    t = downsample(t,n);
    y = downsample(y,n);
    plot(t,y(:,1),'color',co(1,:)); hold on
    plot(t,y(:,2),'color',co(2,:));
end
% yyaxis left

subplot(2,1,1)
ylabel('Active power [MW]')
xlim(xlims2)
ylim(ylim3)
subplot(2,1,2)
ylabel('Normalized speed')
xlabel('Time [s]')
xlim(xlims2)
ylim(ylim4)
%
saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'Nordic5_wind'),'epsc');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ideal and Hydro + Wind (less wind resources)
h=1;
figureLatex
n_down = 5;
co = ([0.6,0.6,0.6;0,0.5,0.7;0,0,0;0,0,0;]);
co_gray = ([0.6,0.6,0.6]);
co_light = [1,0.8,0];

subplot(2,1,1)
for i = [1,4]
    run = runs{1,i};
    t = run.t;    
    y = run.w_coi;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    if i == 4
        plot(t,y,'--','color',co(i,:)); hold on
    else
        plot(t,y,'color',co(i,:)); hold on
    end 
end

for i = [1,4]
    run = runs{1,i};
    t = run.t;    
    y = run.w;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    plot(t,y(:,[5]),'color',co_light); hold on
end

for i = [1,4]
    run = runs{1,i};
    t = run.t;    
    y = run.w_coi;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    if i == 4
        plot(t,y,'--','color',co(i,:)); hold on
    else
        plot(t,y,'color',co(i,:)); hold on
    end 
end

yline(49,'k')
l = legend('$\omega_\mathrm{COI}$ (ideal)',...
           '$\omega_\mathrm{COI}$ (hydro+wind)',...
           '$\omega_5$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','southeast')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left


subplot(2,1,2)
for i = [1,4]
    run = runs{1,i};
    t = run.t;    
    y = run.P_FCR_tot;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    if i == 4
        plot(t,y,'--','color',co(i,:)); hold on
    else
        plot(t,y,'color',co(i,:)); hold on
    end 
    if i == 4
        y2 = sum(run.P_FCR(:,1:3),2);    
        y2 = downsample(y2,n_down);
        y3 = sum(run.P_FCR(:,4:5),2);    
        y3 = downsample(y3,n_down);
        plot(t,y2,'color',[0,0.5,0.7]);
        plot(t,y3,'color',[0,0,0]);
    end
end

yline(0,'k')
l = legend('$P_\mathrm{ideal}$ ',...
           '$P_\mathrm{hydro+wind}$',...
           '$P_\mathrm{hydro}$',...
           '$P_\mathrm{wind}$');
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

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'Nordic5_hydro_and_wind_x2'),'epsc');
end

%% Plot wind response (less wind resources)
h=1;
figureLatex
       
co = ([;0,0,0;0,0.5,0.7]);

subplot(2,1,1)
for i = 4
    run = runs{1,i};
    t = run.t;    
    y = run.P_wind; 
    y2 = run.P_wind_set*0.96*1e-5;  
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    y2 = downsample(y2,n_down);
    y2 = y2-y2(1,:)+y(1,:);
    plot(t,y(:,1),'color',co(1,:)); hold on
    plot(t,y(:,2),'color',co(2,:));
    tspan = 321:374;
    plot(t(tspan),y2(tspan-1,1),'-','color',[0.6,0.6,0.6]);
    tspan = 377:423;
    plot(t(tspan),y2(tspan-1,2),'-','color',[0.6,0.6,0.6]);
%     plot(t,y(:,2),'color',co(2,:));
end


l = legend('Bus 2' , 'Bus 4');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','northeast')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left

subplot(2,1,2)
for i = 4
    run = runs{1,i};
    t = run.t;
    y = run.w_wind;
    n = 10;
    t = downsample(t,n);
    y = downsample(y,n);
    plot(t,y(:,1),'color',co(1,:)); hold on
    plot(t,y(:,2),'color',co(2,:));
end
% yyaxis left

subplot(2,1,1)
ylabel('Active power [MW]')
xlim(xlims2)
ylim([100,550])
subplot(2,1,2)
ylabel('Normalized speed')
xlabel('Time [s]')
xlim(xlims2)
ylim(ylim4-[0.1,0])
%
saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'Nordic5_wind_x2'),'epsc');
end