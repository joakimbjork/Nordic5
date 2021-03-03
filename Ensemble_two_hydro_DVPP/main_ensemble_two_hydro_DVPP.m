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

%%
v_wind = 8;
load_data_ensemble_two_hydro_DVPP

% FCR
s = tf('s');
R = 20; % MW/Hz
R = R*fn;
% K_ideal = R/(s*7.2+1); % First order FCR
K_ideal = R*(1+6.5*s)/((1+2*s)*(1+17*s)); % Second order FCR (avoids second freqeuncy dip)

c_slow = [50;50]; % Distribution of FCR-D
c_fast = [100]; % Distribution of fast FCR
    
c_stat{1} = c_slow;
c_stat{2} = c_fast;



%% Hydro and Battery
fcr_mode = 0;

[FCR, FCR_ideal, c] = ensemble_fcr_controller(...
                       K_ideal,c_stat,GenData,fcr_mode,data_NREL,wind_turbines);

% For hydro + wind, fcr_mode = 1,   
if fcr_mode == 1;
    for i = 1:3
        FCR{i} = modelred_hsv(FCR{i},4);
    end
end

Hhydro = [];
s = tf('s');
for i = 1:2
   Ty = GenData.Ta(i);
   Tw = GenData.TW(i);
   g0 = GenData.PGEN(i)/GenData.MBASE(i);
   Tz_hydro = Tw*g0; % 1/z_hydro

   Gservo = 1/(s*Ty+1); 
   Ghydro = (1 - s*Tz_hydro)/(1 + s*Tz_hydro/2);
   Hhydro{i} = Ghydro*Gservo*GenData.MBASE(i); 
end

c1 = zpk(c.c_hydro{1});
c2 = zpk(c.c_hydro{2});
c3 = zpk(c.c_wind{1});
F1 = Hhydro{1}*FCR{1};
F2 = Hhydro{2}*FCR{2};
F3 = Gwind*FCR{3};

if fcr_mode == 1;
    c1 = modelred_hsv(c1,3);
    c2 = modelred_hsv(c2,3);
    c3 = modelred_hsv(c3,3);
end


%% Bode hydro only 
ylim1 = [0.40,1.1];
ylim2 = [-180-20, 20];
xlim1 = [1e-2,1e2];

omega = logspace(-2,2,100);
[MAG,PHASE] = bode([c1;c2;c1+c2],omega); 

MAG = squeeze(MAG); PHASE = squeeze(PHASE)-360;

h=1;
figureLatex
       
co = ([0,0,0;0.6,0.6,0.6;0,0.5,0.7]);
set(groot,'defaultAxesColorOrder',co)

subplot(2,1,1)
semilogx(omega,MAG(1:2,:)), hold on
semilogx(omega,MAG(3,:),'k--')

ylim(ylim1)

subplot(2,1,2)
semilogx(omega,PHASE(1:2,:)), hold on;
semilogx(omega,PHASE(3,:),'k--')
ylim(ylim2)
yticks([-3*90,-2*90,-90,0,90])

xline(1,'k');
xline(0.5,'color',co(2,:));
yline(-90,'k');

subplot(2,1,1)
ylabel('Amplitude')

l = legend('$c_1$', '$c_2$','$c_1+c_2$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','east')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')

subplot(2,1,2);
ylabel('Phase [$^\circ$]')
xlabel('Angular frequency [rad/s]')

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'c_factors_two_hydro'),'epsc');
end

%% Bode hydro+battery
ylim1 = [0,2];
ylim2 = [-190, 100];
xlim1 = [1e-2,1e2];

omega = logspace(-2,2,100);
c3 = 1-(c1+c2);
[MAG,PHASE] = bode([c1;c2;c3;c1+c2+c3],omega); 

MAG = squeeze(MAG); PHASE = squeeze(PHASE);

PHASE([1,2],:) = PHASE([1,2],:)-360;

h=1;
figureLatex
       
co = ([0,0,0;0.6,0.6,0.6;0,0.5,0.7]);
set(groot,'defaultAxesColorOrder',co)

subplot(2,1,1)
semilogx(omega,MAG(1:3,:)), hold on
semilogx(omega,MAG(4,:),'k--')

 ylim(ylim1)

subplot(2,1,2)
semilogx(omega,PHASE(1:3,:)), hold on;
semilogx(omega,PHASE(4,:),'k--')
ylim(ylim2)
yticks([-3*90,-2*90,-90,0,90])

xline(1,'k');
xline(0.5,'color',co(2,:));
yline(-90,'k');

subplot(2,1,1)
ylabel('Amplitude')

l = legend('$c_1$', '$c_2$','$c_3$','$c$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','northwest')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')

subplot(2,1,2);
ylabel('Phase [$^\circ$]')
xlabel('Angular frequency [rad/s]')

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'c_factors_hydro_battery'),'epsc');
end

%% Hydro and Wind
fcr_mode = 1;

[FCR, FCR_ideal, c] = ensemble_fcr_controller(...
                       K_ideal,c_stat,GenData,fcr_mode,data_NREL,wind_turbines);

% For hydro + wind, fcr_mode = 1,   
if fcr_mode == 1;
    for i = 1:3
        FCR{i} = modelred_hsv(FCR{i},4);
    end
end

Hhydro = [];
s = tf('s');
for i = 1:2
   Ty = GenData.Ta(i);
   Tw = GenData.TW(i);
   g0 = GenData.PGEN(i)/GenData.MBASE(i);
   Tz_hydro = Tw*g0; % 1/z_hydro

   Gservo = 1/(s*Ty+1); 
   Ghydro = (1 - s*Tz_hydro)/(1 + s*Tz_hydro/2);
   Hhydro{i} = Ghydro*Gservo*GenData.MBASE(i); 
end

c1 = zpk(c.c_hydro{1});
c2 = zpk(c.c_hydro{2});
c3 = zpk(c.c_wind{1});
F1 = Hhydro{1}*FCR{1};
F2 = Hhydro{2}*FCR{2};
F3 = Gwind*FCR{3};

if fcr_mode == 1;
    c1 = modelred_hsv(c1,3);
    c2 = modelred_hsv(c2,3);
    c3 = modelred_hsv(c3,3);
end

%% Bode hydro+wind
ylim1 = [0,2];
ylim2 = [-190, 280];
xlim1 = [1e-3,1e1];

omega = logspace(-3,1,100);
c3 = 1-(c1+c2);
[MAG,PHASE] = bode([c1;c2;c3;c1+c2+c3],omega); 

MAG = squeeze(MAG); PHASE = squeeze(PHASE);

PHASE([1,2],:) = PHASE([1,2],:)-360;

h=1;
figureLatex
       
co = ([0,0,0;0.6,0.6,0.6;0,0.5,0.7]);
set(groot,'defaultAxesColorOrder',co)

subplot(2,1,1)
semilogx(omega,MAG(1:3,:)), hold on
semilogx(omega,MAG(4,:),'k--')

ylim(ylim1)

subplot(2,1,2)
semilogx(omega,PHASE(1:3,:)), hold on;
semilogx(omega,PHASE(4,:),'k--')
ylim(ylim2)
yticks([-3*90,-2*90,-90,0,90,2*90,3*90,4*90])

xline(1,'k');
xline(0.5,'color',co(2,:));
xline(0.04786,'color',co(3,:))
yline(-90,'k');

subplot(2,1,1)
ylabel('Amplitude')

l = legend('$c_1''$', '$c_2''$','$c_3''$','$c''$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','northwest')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')

subplot(2,1,2);
ylabel('Phase [$^\circ$]')
xlabel('Angular frequency [rad/s]')

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'c_factors_hydro_wind'),'epsc');
end

%% Run Simulation
simulate_model = false;
if simulate_model
    v_wind = 8;
    load_data_ensemble_two_hydro_DVPP

    % FCR
    s = tf('s');
    R = 20; % MW/Hz
    R = R*fn;
    % K_ideal = R/(s*7.2+1); % First order FCR
    K_ideal = R*(1+6.5*s)/((1+2*s)*(1+17*s)); % Second order FCR (avoids second freqeuncy dip)

    c_slow = [50;50]; % Distribution of FCR-D
    c_fast = [100]; % Distribution of fast FCR

    c_stat{1} = c_slow;
    c_stat{2} = c_fast;

    for i = 1:3
        if i == 1
            fcr_mode = 2;
            [FCR, FCR_ideal, c] = ensemble_fcr_controller(...
                       K_ideal,c_stat,GenData,fcr_mode,data_NREL,wind_turbines);
            Pwind0 = 0;
            Fbat = tf(0);
        elseif i == 2
            fcr_mode = 2;
            [FCR, FCR_ideal, c] = ensemble_fcr_controller(...
                       K_ideal,c_stat,GenData,fcr_mode,data_NREL,wind_turbines);
            Pwind0 = 0;
            Fbat = (1-(c.c_hydro{1}+c.c_hydro{2}))*K_ideal;
        else
            fcr_mode = 1;
            [FCR, FCR_ideal, c] = ensemble_fcr_controller(...
                       K_ideal,c_stat,GenData,fcr_mode,data_NREL,wind_turbines);
            for ii = 1:3
                 FCR{ii} = modelred_hsv(FCR{ii},4);
            end
            Pwind0 = POPT(1);
            Fbat = tf(0);
        end
        model = 'model_ensemble_two_hydro_DVPP';
        sim(model);

        runs{i}.t = t-time_activate_wind ; 
        % wind x/x_opt
        runs{i}.Pe = Pe; runs{i}.w_rot = w_rot;
            % line measurements [MW]
        runs{i}.Pe_wind_hydro = Pe_wind_hydro; 
        runs{i}.Pe_tot = Pe_tot;
        runs{i}.Pe_bat_storage = Pe_bat_storage;
    end

end

%% Load From File
load_model = true;
if load_model
    load('saved_example.mat')
    disp('Load simulation result from file')
    
    v_wind = 8;
    load_data_ensemble_two_hydro_DVPP
end

%% Plot hydro only response
run = runs{1};
ylim1 = [75,105];
ylim2 = [-2, 12];
xlim1 = [0,100];
xlim2 = [0,100];
h=1;
figureLatex

Phydro = Pgen(1);
Pwind = POPT(1);
       
co = ([0,0,0;0,0.5,0.7]);

subplot(2,1,1)

y = run.Pe_tot;
n = 10;
t = run.t;
t = downsample(t,n);
y = downsample(y,n);
plot(t,y(:,1),'color',co(1,:)); hold on
plot(t,y(:,2),'--','color',[0.6,0.6,0.6])

l = legend('$P_1+P_2$', '$P_\mathrm{des}$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','east')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left
ylabel('Active power [MW]')
yticks([60:10:120])

%

subplot(2,1,2)

y = run.Pe_wind_hydro;
n = 10;
y = downsample(y,n);
plot(t,y(:,2)+Phydro,'-','color',co(1,:)); hold on
plot(t,y(:,3)+Phydro,'-','color',[0.6,0.6,0.6]);

ylabel('Active power [MW]')
xlabel('Time [s]')

l = legend( '$P_1$','$P_2$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','east')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')

subplot(2,1,1)
xlim(xlim1)
ylim(ylim1)
subplot(2,1,2)
xlim(xlim2)
ylim(ylim2+Phydro)

%
saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'DVPP_hydro_only_line'),'epsc');
end

%% Plot hydro + storage response
run = runs{2};
ylim1 = [75,105];
ylim2 = [-2, 12];
xlim1 = [0,100];
xlim2 = [0,100];
% Plot line response

h=1;
figureLatex

Phydro = Pgen(1);
Pwind = POPT(1);
       
co = ([0,0,0;0,0.5,0.7]);

subplot(2,1,1)

y = run.Pe_tot;
n = 10;
t = run.t;
t = downsample(t,n);
y = downsample(y,n);
plot(t,y(:,1),'color',co(1,:)); hold on
plot(t,y(:,2),'--','color',[0.6,0.6,0.6])

l = legend('$P_\mathrm{DVPP}$', '$P_\mathrm{des}$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','east')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left
ylabel('Active power [MW]')
yticks([60:10:120])
xlim(xlim1)
ylim(ylim1)
%

subplot(2,1,2)

y = run.Pe_wind_hydro;
n = 10;
y = downsample(y,n);
yyaxis left
plot(t,y(:,2)+Phydro,'-','color',co(1,:)); hold on
plot(t,y(:,3)+Phydro,'-','color',[0.6,0.6,0.6]);

yyaxis right
plot(t,y(:,4),'-','color',co(2,:)); 
% plot(t,y(:,2:3)+Phydro,'-','color',[0,0.5,0.7]); 

yyaxis left
ylabel('Active power [MW]')
xlabel('Time [s]')

l = legend( '$P_1$','$P_2$','$P_3$ (battery)');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','east')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')

xlim(xlim1)
yyaxis right
yticks([0:5:10])
ylim(ylim2)
yyaxis left
ylim(ylim2+Phydro)

ax = gca;
 ax.YAxis(2).Color = co(2,:);
%
saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'DVPP_storage_line'),'epsc');
end

%% Plot Storage
ylim1 = [-1,6];
ylim2 = [-18, 2];
xlim1 = [0,100];
xlim2 = [0,100];
% Plot line response

h=1;
figureLatex
       
co = ([0,0,0;0,0.5,0.7]);

subplot(2,1,1)

y = run.Pe_wind_hydro;
n = 10;
t = run.t;
t = downsample(t,n);
y = downsample(y,n);
plot(t,y(:,4),'color',co(2,:)); hold on
% plot(t,y(:,2)+y(:,3),'-','color',co(1,:)); hold on


% yyaxis left
ylabel('Active power [MW]')
yticks([0:2:10])
xlim(xlim1)
ylim(ylim1)
%

subplot(2,1,2)

y = -run.Pe_bat_storage*1e-3/(60*60);
n = 10;
y = downsample(y,n);
plot(t,y(:,1),'-','color',co(2,:)); hold on

ylabel('Energy [kWh]')
xlabel('Time [s]')


xlim(xlim1)
ylim(ylim2)
%
saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'DVPP_storage'),'epsc');
end

%% Plot wind response
run = runs{3};
ylim1 = [0.70,1.7]*POPT(1);
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
% plot(t,y(:,3),'color',[0.6,0.6,0.6]); hold on
plot(t,y(:,1),'color',co(2,:)); ; hold on


y = run.Pe*POPT(1);
n = 10;
y = downsample(y,n);
plot(t,y(:,2),'--','color',co(2,:))


l = legend('Nonlinear', 'Linear');
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
plot(t,y(:,1),'color',co(2,:)); hold on
plot(t,y(:,2),'--','color',co(2,:));

% yyaxis left
ylabel('Normalized rotor speed')
xlabel('Time [s]')

xlim(xlim1)
ylim(ylim2)
%
saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'DVPP_wind'),'epsc');
end

%% Plot wind + hydro

ylim1 = [87,113];
ylim2 = [-2, 12];
xlim1 = [0,100];
h=1;
figureLatex

Phydro = Pgen(1);
Pwind = POPT(1);
       
co = ([0,0,0;0,0.5,0.7]);

subplot(2,1,1)

y = run.Pe_tot;
n = 10;
t = run.t;
t = downsample(t,n);
y = downsample(y,n);
plot(t,y(:,1),'color',co(1,:)); hold on
plot(t,y(:,2),'--','color',[0.6,0.6,0.6])

l = legend('$P_\mathrm{DVPP}$', '$P_\mathrm{des}$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','east')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left
ylabel('Active power [MW]')
yticks([60:10:120])
xlim(xlim1)
ylim(ylim1)
%

subplot(2,1,2)

y = run.Pe_wind_hydro;
n = 10;
y = downsample(y,n);
yyaxis left
plot(t,y(:,2)+Phydro,'-','color',co(1,:)); hold on
plot(t,y(:,3)+Phydro,'-','color',[0.6,0.6,0.6]);

yyaxis right
plot(t,y(:,1)+Pwind,'-','color',co(2,:)); 
% plot(t,y(:,2:3)+Phydro,'-','color',[0,0.5,0.7]); 

yyaxis left
ylabel('Active power [MW]')
xlabel('Time [s]')

l = legend( '$P_1$ ','$P_2$ ','$P_3$ (wind)');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','east')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')

xlim(xlim1)
yyaxis right
yticks([10:5:20])
ylim(ylim2+Pwind)
yyaxis left
ylim(ylim2+Phydro)

ax = gca;
 ax.YAxis(2).Color = co(2,:);
%
saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'DVPP_wind_line'),'epsc');
end
