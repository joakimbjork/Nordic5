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
%% Diferent wind speeds
v = [8,10];
kstab_multiplier = 2;
pstep = 0.2;

simulate_model = true;
if simulate_model
    for i = 1:2
        v_wind = v(i);
        load_data_wind_open_loop
        
        model = 'model_wind_open_loop';
        sim(model);
        runs{i}.t = t-time_activate_wind ; 
        % x/x_opt
        runs{i}.Pe = Pe; runs{i}.w_rot = w_rot;
        % x
        runs{i}.Pe_abs = Pe_abs; runs{i}.w_rot_abs = w_rot_abs;
    end
end

%% Plot response
ylim1 = [1.0,4.5];
ylim2 = [0.6, 1.25];
xlim1 = [0,100];
h=1;
figureLatex
       
co = ([0,0,0;0,0.5,0.7;1,0,0]);

subplot(2,1,1)
for i = 1:2
    run = runs{1,i};
    t = run.t;
    y = run.Pe_abs*5;
    n = 10;
    t = downsample(t,n);
    y = downsample(y,n);
    plot(t,y(:,1),'color',co(i,:)); hold on
end
for i = 1:2
    run = runs{1,i};
    t = run.t;
    y = run.Pe_abs*5;
    n = 10;
    t = downsample(t,n);
    y = downsample(y,n);
    plot(t,y(:,2),'--','color',co(i,:))
    plot(t,y(:,3),'color',[0.6,0.6,0.6])
%     yline(y(1,1),'k')
end

l = legend('$v = 8$ m/s', '$v = 10 $ m/s','Linear','$P_\mathrm{ref}$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','northeast')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left
ylabel('Active power, $P_e$ [MW]')

xlim(xlim1)
ylim(ylim1)
%

subplot(2,1,2)
for i = 1:2
    run = runs{1,i};
    t = run.t;
    y = run.w_rot_abs*data_NREL.wt.rotor.ratedspeed;
    n = 10;
    t = downsample(t,n);
    y = downsample(y,n);
    plot(t,y(:,1),'color',co(i,:)); hold on
    plot(t,y(:,2),'--','color',co(i,:));
end
% yyaxis left
ylabel('Turbine speed, $\Omega$ [rad/s]')
xlabel('Time [s]')

xlim(xlim1)
ylim(ylim2)


%
saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'speed_8_and_10_abs'),'epsc');
end

%% Plot response

ylim1 = [0.75,1.25];
ylim2 = [0.7, 1.05];
xlim1 = [0,100];
h=1;
figureLatex
       
co = ([0,0,0;0,0.5,0.7;1,0,0;]);

subplot(2,1,1)
for i = 1:2
    run = runs{1,i};
    t = run.t;
    y = run.Pe;
    n = 10;
    t = downsample(t,n);
    y = downsample(y,n);
    plot(t,y(:,1),'color',co(i,:)); hold on
end
for i = 1:2
    run = runs{1,i};
    t = run.t;
    y = run.Pe;
    n = 10;
    t = downsample(t,n);
    y = downsample(y,n);
    plot(t,y(:,2),'--','color',co(i,:))
    plot(t,y(:,3),'color',[0.6,0.6,0.6])
end

l = legend('$v = 8$ m/s', '$v = 10 $ m/s','Linear','$P_\mathrm{ref}$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','northeast')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left
ylabel('$P_e/P_\mathrm{opt}$')

xlim(xlim1)
ylim(ylim1)
%

subplot(2,1,2)
for i = 1:2
    run = runs{1,i};
    t = run.t;
    y = run.w_rot;
    n = 10;
    t = downsample(t,n);
    y = downsample(y,n);
    plot(t,y(:,1),'color',co(i,:)); hold on
    plot(t,y(:,2),'--','color',co(i,:));
end
% yyaxis left
ylabel('$\Omega/\Omega_\mathrm{opt}$')
xlabel('Time [s]')

xlim(xlim1)
ylim(ylim2)
%
saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'speed_8_and_10'),'epsc');
end

%% Different gain
k0 = 2*[1,1.5];
v_wind = 8;
pstep = 0.2;

simulate_model = true;
if simulate_model
    for i = 1:2
        kstab_multiplier = k0(i);
        load_data_wind_open_loop
               
        model = 'model_wind_open_loop';        
        sim(model);
        runs{i}.t = t-time_activate_wind ; 
        % x/x_opt
        runs{i}.Pe = Pe; runs{i}.w_rot = w_rot;
        % x
        runs{i}.Pe_abs = Pe_abs; runs{i}.w_rot_abs = w_rot_abs;
    end
end

%% Plot response
ylim1 = [ 0.7, 1.4 ];
ylim2 = [ 0.7, 1.05 ];

% ylim1 = [ 0.7, 1.4 ];
% ylim2 = [ 0.7, 1.05 ];
xlim1 = [0,100];
h=1;
figureLatex
       
co = ([0,0,0;0,0.5,0.7;1,0,0;]);

subplot(2,1,1)
for i = 1:2
    run = runs{1,i};
    t = run.t;
    y = run.Pe;
    n = 10;
    t = downsample(t,n);
    y = downsample(y,n);
    plot(t,y(:,1),'color',co(i,:)); hold on
end
for i = 1:2
    run = runs{1,i};
    t = run.t;
    y = run.Pe;
    n = 10;
    t = downsample(t,n);
    y = downsample(y,n);
    plot(t,y(:,2),'--','color',co(i,:)); 
    if i == 1
        plot(t,y(:,3),'color',[0.6,0.6,0.6])
    end
end


l = legend('$k = 0.72 $', '$k = 1.08$','Linear','$P_\mathrm{ref}$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','best')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
ylabel('$P_e/P_\mathrm{opt}$')
% plot(t,y(:,3),'k--')
% plot(t,y(:,2),'k--')
xlim(xlim1)
ylim(ylim1)
%

subplot(2,1,2)
for i = 1:2
    run = runs{1,i};
    t = run.t;
    y = run.w_rot;
    n = 10;
    t = downsample(t,n);
    y = downsample(y,n);
    plot(t,y(:,1),'color',co(i,:)); hold on
    plot(t,y(:,2),'--','color',co(i,:));
end

ylabel('$\Omega/\Omega_\mathrm{opt}$')
xlabel('Time [s]')

xlim(xlim1)
ylim(ylim2)
%
saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'gain_1_5'),'epsc');
end

%% Low Speed Protection test
v_wind = 8;
kstab_multiplier = 2;
pstep_ = [0.3,0.3];
safety_on = [0,1];

simulate_model = true;
if simulate_model
    for i = 1:2      
        load_data_wind_open_loop
        data_NREL.ctrl.low_speed_safety=safety_on(i);
        pstep = pstep_(i);
        
        model = 'model_wind_open_loop';
        sim(model);
        runs{i}.t = t-time_activate_wind ; 
        % x/x_opt
        runs{i}.Pe = Pe; runs{i}.w_rot = w_rot;
        % x
        runs{i}.Pe_abs = Pe_abs; runs{i}.w_rot_abs = w_rot_abs;
    end
end

%% Plot
ylim1 = [ 0.7, 1.4 ];
ylim2 = [ 0.7, 1.05 ];

load_model = false;
if load_model
    load('rund_.mat')
    disp('Load simulation result from file')
end

% ylim1 = [ 0.7, 1.4 ];
% ylim2 = [ 0.7, 1.05 ];
xlim1 = [0,100];
h=1;
figureLatex
       
co = ([0,0,0;0,0.5,0.7;1,0,0;]);

subplot(2,1,1)
for i = 1:2
    run = runs{1,i};
    t = run.t;
    y = run.Pe;
    n = 10;
    t = downsample(t,n);
    y = downsample(y,n);
    plot(t,y(:,1),'color',co(i,:)); hold on   
    if i == 1
        plot(t,y(:,2),'--','color',co(i,:));
        plot(t,y(:,3),'color',[0.6,0.6,0.6])
    end
end

l = legend('k=0.72', 'Linear', '$P_\mathrm{ref}$', 'Low speed protection');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','best')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
ylabel('$P_e/P_\mathrm{opt}$')

xlim(xlim1)
ylim(ylim1)
%
subplot(2,1,2)
for i = [1,2]
    run = runs{1,i};
    t = run.t;
    y = run.w_rot;
    n = 10;
    t = downsample(t,n);
    y = downsample(y,n);
    plot(t,y(:,1),'color',co(i,:)); hold on
    if i ==1
     plot(t,y(:,2),'--','color',co(i,:));
    end
end




ylabel('$\Omega/\Omega_\mathrm{opt}$')
xlabel('Time [s]')

xlim(xlim1)
ylim(ylim2)

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'step_20_30'),'epsc');
end

