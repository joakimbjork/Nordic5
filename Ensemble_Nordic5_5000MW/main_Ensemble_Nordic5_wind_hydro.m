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

%% Load From File
load_model = false;
if load_model
    load('saved_example_wind_hydro.mat')
    disp('Load simulation result from file')
end

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
K_fcr = R*(1+6.5*s)/((1+2*s)*(1+17*s)); % Second order FCR (avoids second freqeuncy dip)

for i_run = 1:3
    
    if i_run == 1
        model = 'Ensemble_Nordic5_5000MW'; 
        load_data_Ensemble_Nordic5_5000MW_110
        fcr_mode = 0;        
    elseif i_run == 2
        model = 'Ensemble_Nordic5_5000MW'; 
        load_data_Ensemble_Nordic5_5000MW_110
        fcr_mode = 2;  
    elseif i_run == 3
        model = 'Ensemble_Nordic5_5000MW_wind_hydro'; 
        load_data_Ensemble_Nordic5_5000MW_110
        fcr_mode = 1;                                         
    end

    c_slow = [60;30;10]; % Distribution of FCR-D
    c_fast = [1,2]; % Distribution of fast FCR
    c = {};
    c{1} = c_slow;
    c{2} = c_fast;

    [FCR, FCR_ideal,~,H] = ensemble_fcr_controller(...
                   K_fcr,c,GenData,fcr_mode,data_NREL,wind_turbines);        

    % For hydro + wind, fcr_mode = 1,   
    if fcr_mode == 1
        i_fcr=1; FCR{i_fcr} = modelred_hsv(FCR{i_fcr},4);
        i_fcr=2; FCR{i_fcr} = modelred_hsv(FCR{i_fcr},4);
        i_fcr=3; FCR{i_fcr} = modelred_hsv(FCR{i_fcr},4);
        i_fcr=4; FCR{i_fcr} = modelred_hsv(FCR{i_fcr},4);
        i_fcr=5; FCR{i_fcr} = modelred_hsv(FCR{i_fcr},4);           
    end

    if run_load_flow
        LF = power_loadflow('-v2',model,'solve');                           
    end
%         if linearize
%             time_event = 120;            
%         end

%%

    %%
    M = 2*sum(Wkin);
    G = fn/(s*M+load_damping*fn);
    Ldes = G*K_fcr/fn;
    L = Ldes;  
    F = tf(0);
    if fcr_mode>=1             
        for i_plant = 1:length(H)
            F = F + H{i_plant}*FCR{i_plant};
        end
        L = G*F/fn;            
    end 
    S = inv(1+L);
    % Gcl = feedback(G,K_ideal/fn);
    % step(-S*G*d_dim-0.1)

    runs{i_run}.L = L;
    runs{i_run}.Ldes = Ldes;
    runs{i_run}.SG = S*G;
    
        
    F = tf(0);
    if fcr_mode==1             
        for i_plant = 1:3
            F = F + H{i_plant}*FCR{i_plant};
        end
        L = G*F/fn;            
    end 
    runs{i_run}.Lhydro = L;
    
    F = tf(0);
    if fcr_mode==1             
        for i_plant = 4:5
            F = F + H{i_plant}*FCR{i_plant};
        end
        L = G*F/fn;            
    end 
    runs{i_run}.Lwind = L;
    
    %% alt open-loop
    M = 2*sum(Wkin);
    G = fn/(s*M);
    Ldes = G*(K_fcr+load_damping)/fn;
    L = Ldes;  
    F = tf(load_damping);
    if fcr_mode>=1             
        for i_plant = 1:length(H)
            F = F + H{i_plant}*FCR{i_plant};
        end
        L = G*F/fn;            
    end
    S = inv(1+L);
    % Gcl = feedback(G,K_ideal/fn);
    % step(-S*G*d_dim-0.1)

    runs{i_run}.L_alt = L;
    runs{i_run}.Ldes_alt = Ldes;
    runs{i_run}.SG_alt = S*G;

    %%
    if simulate_model
        time_event = time_activate_wind+5;

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
        end
        
        
    end      
end

%% Time hydro
xlims=[0,40];

ylim1=[48.5,50];
ylim2=[-300,1800];

xlims2=[0,60];
ylim3=[200,850];
ylim4=[0.8,1.05];

n_down = 5; % For down sampling signals

% t = downsample(t,n_down);
% y = downsample(y,n_down);
% Ideal and Hydro
h=1;
figureLatex

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
    if i == 1
        plot(t,y(:,1),'--','color',co(i,:)); hold on
    else
        plot(t,y(:,1),'color',co(i,:)); hold on
    end
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
    if i == 1
        plot(t,y(:,1),'--','color',co(i,:)); hold on
    else
        plot(t,y(:,1),'color',co(i,:)); hold on
    end
end

yline(49,'k')
l = legend('$\omega_\mathrm{COI}$ (hydro)',...
           '$\omega_\mathrm{COI}$ (desired)',...
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
    y2 = run.P_FCR(:,1);    
    y2 = downsample(y2,n_down);
    if i == 1
        plot(t,y(:,1),'--','color',co(i,:)); hold on        
    else
        plot(t,y(:,1),'color',co(i,:)); hold on
        plot(t,y2,'color',co(3,:));
    end
end
for i = 2
    run = runs{1,i};
    t = run.t;    
    t = downsample(t,n_down);
    y2 = run.P_FCR(:,2:3);    
    y2 = downsample(y2,n_down);
    plot(t,y2,'color',co(3,:));
end


yline(0,'k')
l = legend('$P_\mathrm{hydro}$',...
           'Bus 1-3',...
           '$P_\mathrm{des}$');
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

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'Nordic5_fcr_hydro'),'epsc');
end

%% Ideal and Hydro bode
omega = logspace(-3,1,100);
h = 1;
figureLatex
co = ([0.6,0.6,0.6;0,0.5,0.7;0,0,0;0,0.5,0.7;1,0,0;]);

subplot(2,1,1);
wc=[]; Pm = [];
for i = [2,1]     
    L = runs{i}.L;
%     z(i)=100000;
%     if i ==2
%         [p,z_] = pzmap(L);
%         z(i) = z_(z_>0);
%     end
        
    [MAG,PHASE] = bode(L,omega);
    [~,Pm(i),~,wc(i)] = margin(L);
    MAG = squeeze(MAG);
    PHASE = squeeze(PHASE);
    adj = PHASE(1,1)/360;
    adj = round(adj,0)*360;
    PHASE = PHASE-adj;
    subplot(2,1,1);
    if i == 1
        loglog(omega,MAG,'--','color',co(i,:)); hold on     
    else
        loglog(omega,MAG,'color',co(i,:)); hold on 
    end
    subplot(2,1,2);
    if i == 1
        semilogx(omega,PHASE,'--','color',co(i,:)); hold on  
    else
        semilogx(omega,PHASE,'color',co(i,:)); hold on 
    end
end
for i = [1,2]
    semilogx([wc(i),wc(i)],[0,Pm(i)]-180,'color',co(i,:)); hold on 
end
for i = [2]
    subplot(2,1,1);
    xline(wc(i))
%     xline(z(i),'color',co(3,:))
    yline(1)
    subplot(2,1,2);
    yline(-180)
%     xline(z(i),'color',co(3,:))
end

subplot(2,1,1);
l = legend('$L_\mathrm{hydro}$',...
           '$L_\mathrm{des}$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','southwest')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
ylim([0.001,20])
xlim([1e-3,10])
ylabel('Amplitude')

subplot(2,1,2);
l = legend([num2str(round(Pm(2),0)),'$^\circ$'],...
           [num2str(round(Pm(1),0)),'$^\circ$']);
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','southwest')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
ylim([-90*3,10])
yticks(-90*6:90:90)
xlim([1e-3,10])
ylabel('Phase [$^\circ$]')
xlabel('Angular frequency [rad/s]')

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'Nordic5_bode_hydro'),'epsc');
end

%% Time hydro+wind
xlims=[0,40];

ylim1=[48.5,50];
ylim2=[-300,1800];

xlims2=[0,60];
ylim3=[200,850];
ylim4=[0.8,1.05];

n_down = 5; % For down sampling signals

% t = downsample(t,n_down);
% y = downsample(y,n_down);
% Ideal and Hydro
h=1;
figureLatex

co = ([0.6,0.6,0.6;0,0,0;0,0.5,0.7;1,0,0;]);
co_gray = ([0.6,0.6,0.6]);
co_light = [1,0.8,0];

subplot(2,1,1)
for i = [3,1]
    run = runs{1,i};
    t = run.t;    
    y = run.w_coi;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    if i == 1
        plot(t,y(:,1),'--','color',co(1,:)); hold on
    else
        plot(t,y(:,1),'color',co(2,:)); hold on
    end
end

for i = [3,1]
    run = runs{1,i};
    t = run.t;    
    y = run.w;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    plot(t,y(:,[5]),'color',co_light); hold on
end

for i = [3,1]
    run = runs{1,i};
    t = run.t;    
    y = run.w_coi;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);
    if i == 1
        plot(t,y(:,1),'--','color',co(1,:)); hold on
    else
        plot(t,y(:,1),'color',co(2,:)); hold on
    end
end

yline(49,'k')
l = legend('$\omega_\mathrm{COI}$ (hydro+wind)',...
           '$\omega_\mathrm{COI}$ (des)',...
           '$\omega_5$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','southeast')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')
% yyaxis left

subplot(2,1,2)
for i = [3,1]
    run = runs{1,i};
    t = run.t;    
    y = run.P_FCR_tot;    
    t = downsample(t,n_down);
    y = downsample(y,n_down);    
    if i == 1
        plot(t,y(:,1),'--','color',co(1,:)); hold on        
    else
        y2 = run.P_FCR(:,1:3);    
        y2 = downsample(y2,n_down);
        y3 = run.P_FCR(:,4:5);    
        y3 = downsample(y3,n_down);
        plot(t,y(:,1),'color',co(2,:)); hold on
        plot(t,sum(y2,2),'color',co(3,:));
        plot(t,sum(y3,2),'color',co(1,:));
    end
end
% for i = 3
%     run = runs{1,i};
%     t = run.t;    
%     t = downsample(t,n_down);
%     y2 = run.P_FCR(:,1:3);    
%     y2 = downsample(y2,n_down);
%     plot(t,y2,'-','color',co(3,:));
%     y2 = run.P_FCR(:,4:5);    
%     y2 = downsample(y2,n_down);
%     plot(t,y2,'-','color',co(1,:));
% end


yline(0,'k')
l = legend('$P_\mathrm{hydro}+P_\mathrm{wind}$',...
           '$P_\mathrm{hydro}$',...
           '$P_\mathrm{wind}$',...
           '$P_\mathrm{des}$');
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

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'Nordic5_fcr_hydro_wind'),'epsc');
end

%% Ideal and Hydro-Wind bode
omega = logspace(-3,1,100);
h = 1;
figureLatex
co = ([0.6,0.6,0.6;...
    0,0,0;...
    0,0.5,0.7;1,0,0;]);

subplot(2,1,1);
wc=[]; Pm = [];
for i = [3,1]
    if i==3
        L1 = runs{i}.Lhydro;
        [MAG1,PHASE1] = bode(L1,omega);
        MAG1 = squeeze(MAG1);
        PHASE1 = squeeze(PHASE1);
        adj = PHASE1(1,1)/360;
        adj = round(adj,0)*360;
        PHASE1 = PHASE1-adj;
        L2 = runs{i}.Lwind;
        [MAG2,PHASE2] = bode(L2,omega);
        MAG2 = squeeze(MAG2);
        PHASE2 = squeeze(PHASE2);
        adj = PHASE2(1,1)/360;
        adj = round(adj,0)*360;
        PHASE2 = PHASE2-adj+360;
        L = L1+L2;
    else
        L = runs{i}.L;
    end
    [MAG,PHASE] = bode(L,omega);
    [~,Pm(i),~,wc(i)] = margin(L);
    MAG = squeeze(MAG);
    PHASE = squeeze(PHASE);
    adj = PHASE(1,1)/360;
    adj = round(adj,0)*360;
    PHASE = PHASE-adj;
    
    
    
    subplot(2,1,1);
    if i == 1
        loglog(omega,MAG,'--','color',co(1,:)); hold on     
    else
        loglog(omega,MAG,'color',co(2,:)); hold on 
        loglog(omega,MAG1,'color',co(3,:));
        loglog(omega,MAG2,'color',co(1,:));
    end
    
    subplot(2,1,2);
    if i == 1
        semilogx(omega,PHASE,'--','color',co(1,:)); hold on  
    else
        semilogx(omega,PHASE,'color',co(2,:)); hold on 
        loglog(omega,PHASE1,'color',co(3,:));
        loglog(omega,PHASE2,'color',co(1,:));
    end
    
    
end

subplot(2,1,2);
semilogx([wc(3),wc(3)],[0,Pm(3)]-180,'color',co(2,:)); hold on 
semilogx([wc(3),wc(3)],[0,Pm(3)]-180,'--','color',co(1,:));


for i = [3]
    subplot(2,1,1);
    xline(wc(i),'color',co(2,:))
    yline(1)
    subplot(2,1,2);
    yline(-180)

end

subplot(2,1,1);

ylim([0.001,20])
xlim([1e-3,10])
ylabel('Amplitude')

subplot(2,1,2);
l = legend('$L_\mathrm{hydro}+L_\mathrm{wind}$',...
            '$L_\mathrm{hydro}$',...
             '$L_\mathrm{wind}$',...
            '$L_\mathrm{des}$');
         set(l,'Interpreter','LaTeX');
         set(l,'Box','on')
         set(l,'location','northeast')
         set(l,'FontSize',9)
         set(l,'Orientation','vertical')

% ylim([-90*3,10])
yticks(-180*6:180:180*4)
xlim([1e-3,10])

ylabel('Phase [$^\circ$]')
xlabel('Angular frequency [rad/s]')

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'Nordic5_bode_hydro_wind'),'epsc');
end

