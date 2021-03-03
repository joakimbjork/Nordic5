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
if simulate_model
    run = [];

    v_wind = 8;
    load_data_wind_and_hydro_DVPP

    model = 'model_wind_and_hydro_DVPP';
    sim(model);

    run.t = t-time_activate_wind ; 
    % wind x/x_opt
    run.Pe = Pe; run.w_rot = w_rot;
        % line measurements [MW]
    run.Pe_wind_hydro = Pe_wind_hydro; run.Pe_tot = Pe_tot;

end

%% Load From File
load_model = false;
if load_model
    v_wind = 8;
    load_data_wind_and_hydro_DVPP
    load('saved_example.mat')
    disp('Load simulation result from file')
end

%% Plot wind response
ylim1 = [0.70,1.55]*POPT(1);
ylim2 = [0.8, 1.05];
xlim1 = [0,100];
h=1;
figureLatex
       
co = ([0,0,0;0,0.5,0.7]);

subplot(2,1,1)

t = run.t;
y = run.Pe*POPT(1);
n = 10;
t = downsample(t,n);
y = downsample(y,n);
plot(t,y(:,1),'color',co(1,:)); hold on


y = run.Pe*POPT(1);
n = 10;
y = downsample(y,n);
plot(t,y(:,2),'--','color',co(1,:))
plot(t,y(:,3),'color',[0.6,0.6,0.6])

l = legend('Nonlinear', 'Linear','$P_\mathrm{ref}$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','northeast')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left
ylabel('Active power [MW]')

xlim(xlim1)
ylim(ylim1)
%
co = ([0,0,0;0,0.5,0.7]);

subplot(2,1,2)

y = run.w_rot;
n = 10;
y = downsample(y,n);
plot(t,y(:,1),'color',co(1,:)); hold on
plot(t,y(:,2),'--','color',co(1,:));

% yyaxis left
ylabel('Normalized speed')
xlabel('Time [s]')

xlim(xlim1)
ylim(ylim2)
%
saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'DVPP_wind'),'epsc');
end

%% Plot response

ylim1 = [86,98];
ylim2 = [-2, 12];
xlim1 = [0,100];
h=1;
figureLatex

Phydro = 80;
Pwind = 7.1212;
       
co = ([0,0,0;0,0.5,0.7]);

subplot(2,1,1)

y = run.Pe_tot;
n = 10;
t = run.t;
t = downsample(t,n);
y = downsample(y,n);
plot(t,y(:,1),'color',co(1,:)); hold on
plot(t,y(:,2),'color',[0.6,0.6,0.6])

l = legend('$P_\mathrm{DVPP}$', '$P_\mathrm{ideal}$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','east')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left
ylabel('Active power [MW]')
yticks([82:5:120])
xlim(xlim1)
ylim(ylim1)
%

subplot(2,1,2)

y = run.Pe_wind_hydro;
n = 10;
y = downsample(y,n);
yyaxis left
plot(t,y(:,2)+Phydro,'-','color',co(2,:)); hold on
plot(t,y(:,1)+Pwind,'-','color',co(1,:)); 
yyaxis right
plot(t,y(:,2)+Phydro,'-','color',[0,0.5,0.7]); 

yyaxis left
ylabel('Active power [MW]')
xlabel('Time [s]')

l = legend( '$P_1$ (hydro)','$P_2$ (wind)');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','east')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')

xlim(xlim1)
yyaxis left
yticks([7:5:20])
ylim(ylim2+Pwind)
yyaxis right
ylim(ylim2+Phydro)

ax = gca;
 ax.YAxis(2).Color = co(2,:);
%
saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'DVPP_line'),'epsc');
end
