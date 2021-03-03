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

clear all

%% Initialize model
load_data_Nordic5_5000MW_110; 
%load_data_Nordic5_5000MW_110;


model = 'Nordic5_5000MW';

ideal_FCR = false; % hydro-FCR (default)
% ideal_FCR = true; % ideal-FCR through controllable loads
fcr_controller_design

run_load_flow = false; % Rund load flow when switching from 240 to 110 GWs scenario
if run_load_flow
    LF = power_loadflow('-v2',model,'solve');
end

local_freq=0; % 0 -> Use machine speed measurement for frequency dependent loads (needed for linearization)

%% Linear Analysis
% PSSData.Kpss = zeros(5,1); % NO PSS IN LINEARIZATION

linearize_model = true;
if linearize_model
    % Specify the model name
%     model = 'Nordic5_unloaded';
    io = getlinio(model);
    op = operpoint(model);
    sys = linearize(model,io,op);
end

sys = sminreal(sys);
%%  Modal Analyis
[A,B,C,D] = ssdata(sys);

disp('------------- CERTIFY THAT INPUTS/OUTPUTS ARE CORRECT -------------')

idx = 4:8;
B = B(:,4:end); % PSS inputs (d_Vstab)
inputs_pss = sys.InputName(idx)

idx = 6:10;
Cw = C(idx,:); % Machine speed (w)
output_omega = sys.OutputName(idx)

idx = 1:5;
Cpa = C(idx,:); % Accelerating power (Pm-Pe)
output_Pa = sys.OutputName(idx)

% State Index of generators speed
idx_gen = [];
for i = 1:size(Cw,1)
idx_gen = [idx_gen, find(Cw(i,:)==1)]; % index of generators
end

%% Modal Analyis
[E,V,Cmodes,data] = modal_vectors(A,[0,5*2*pi]); %% E = V'*A*inv(V');
mode_idx = [1,3,5,7]; % <<<< Make sure that these are the 4 interarea modes <<<<<<<<<<<<<<<<<<<<<<

if data.damping(1)<0
    disp('OBS! Local bus freqeuncy does not work for the linearization')
end

w_all= [];
e_all=[];
for mode_number = 1:4
    i = mode_idx(mode_number); 

    disp(['Mode', num2str(mode_number),...
          ': s=', num2str(E(i)),...
          ' ; freq=', num2str(data.freq(i)), ' Hz'...
          ' ; damping=' num2str(data.damping(i)*100), ' %']); 


    % w = V(idx_gen,i); % Pick out entries corresponding to rotor states
    w = Cw*V(:,i); % Use output matrix C instead.

    % Rotate vector to align with real axis (usefull for H2 design)
    vpos  = sum(w(real(w)>0));
    vneg  = sum(w(real(w)<=0));
    v2 = [vneg;vpos];
    [~,idx] = max(abs(v2));
    ang = angle(v2(idx));
    %
    w_ = w*exp(-1j*ang);
    if real(w_(1))<0
        ang = ang+pi; % Convention: Generator 1 is in positive real direction 
    end

    V = V*exp(-1j*ang);

    w = Cw*V(:,i);
    w = w/max(abs(w));
    
    w_all = [w_all, w];
    e_all = [e_all, E(i)];
end

%% Plot mode shape
i=1;
w = w_all(:,i);
e = e_all(i);

plotFigure = true;
saveFigure = false;
if plotFigure    
    i=1; w = w_all(:,i); e = e_all(i);
    hfig = plot_compass(w,e);   
    if saveFigure
        path = 'C:\Users\joakbj\Dropbox\KTH\Forskning\gitDocumentation\Nordic5\Figures\Compass_240_with_pss\';
        saveas(hfig,fullfile(path, 'mode_1'),'epsc');
    end
    i=2; w = w_all(:,i); e = e_all(i);
    hfig = plot_compass(w,e); 
    if saveFigure
        saveas(hfig,fullfile(path, 'mode_2'),'epsc');
    end
    i=3; w = w_all(:,i); e = e_all(i);
    hfig = plot_compass(w,e);    
    if saveFigure
        saveas(hfig,fullfile(path, 'mode_3'),'epsc');
    end
    i=4; w = w_all(:,i); e = e_all(i);
    hfig = plot_compass(w,e);  
    if saveFigure
        saveas(hfig,fullfile(path, 'mode_4'),'epsc');
    end
end 

%% Plot residues
res_all = [];
for mode_number = 1:4;
    i = mode_idx(mode_number); 
    disp(['Mode', num2str(i),...
          ': s=', num2str(E(i)),...
          ' ; freq=', num2str(data.freq(i)), ' Hz'...
          ' ; damping=' num2str(data.damping(i)*100), ' %']); 

    % Residue
    V_left = V';
    V_right = inv(V_left);
    v_left = V_left(i,:);
    v_right = V_right(:,i);

    v_obsv = Cpa*v_right;
    w_ctrb = v_left*B;

    residu =  v_obsv.*w_ctrb;
    res = diag(residu); % Local loop
    res_all = [res_all, res];
end

plotFigure = true;
saveFigure = false;
if plotFigure
    i=1; res = res_all(:,i); e = e_all(i);
    hfig = plot_compass(res,e); 
    if saveFigure
%       path = 'C:\Users\joakbj\Dropbox\KTH\Forskning\gitDocumentation\Nordic5';
        saveas(hfig,fullfile(path, 'res_1'),'epsc');
    end
    i=2; res = res_all(:,i); e = e_all(i);
    hfig = plot_compass(res,e);     
    if saveFigure
        saveas(hfig,fullfile(path, 'res_2'),'epsc');
    end
    i=3; res = res_all(:,i); e = e_all(i);
    hfig = plot_compass(res,e);
    if saveFigure
        saveas(hfig,fullfile(path, 'res_3'),'epsc');
    end
    i=4; res = res_all(:,i); e = e_all(i);
    hfig = plot_compass(res,e);      
    if saveFigure
        saveas(hfig,fullfile(path, 'res_4'),'epsc');
    end
end 

%% PSS tuning
machines_left = 1:5;
Kpss_all = [];
disp('------------- PSS TUNING -------------')
for i = 1:4
    ii = i*2; % Oscilatory modes come in pairs
    
    res = res_all(:,i); e = e_all(i);
    [~,m_idx] = max(abs(res));
    if Wkin == 110 & i ==1
        disp('Manual selection of first machine')
        m_idx = 2;
    end

    assert(ismember(m_idx,machines_left));
    
    disp(' ')
    disp(['Tuning PSS ( machine ', num2str(m_idx), ') for mode ', num2str(i),...
      ' : ', num2str(data.freq(ii)), ' Hz (' num2str(data.damping(ii)*100), '%)'])
       
    Fpss = tunePSS(e,res(m_idx))
    Kpss_all{m_idx} = Fpss;
    machines_left(machines_left == m_idx) = []; % One PSS per machine    
end

m_idx = machines_left; % One machine left since 4 modes and 5 machines
m_res_all = res_all(m_idx,:); % Finds most controllable mode from the remaining machine
[~,i] = max(abs(m_res_all));
e = e_all(i);

disp(' ')
disp(['Tuning remaining PSS ( machine ', num2str(m_idx), ') for mode ', num2str(i),...
      ' : ', num2str(data.freq(ii)), ' Hz (' num2str(data.damping(ii)*100), '%)'])
  
Fpss = tunePSS(e,res(m_idx))
Kpss_all{m_idx} = Fpss;  



