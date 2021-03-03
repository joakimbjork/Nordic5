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

function [FCR,FCR_ideal,c_factors,H] = ensemble_fcr_controller(K_ideal,c,GenData,varargin)

% c{1} = c_slow, c{2} = c_fast, cell array of participation factors
% varargin = fcr_mode = 0,1,2 (default 2)
% fcr_mode = 0 -> Ideal fcr at hydro buses
% fcr_mode = 1 -> c_slow and c_fast
% fcr_mode = 2 -> c_clow only

x = length(c);
assert(x<=2, 'FCR only updated for two categories, c_fast and c_slow')

if nargin == 3
    fcr_mode = 2;
    disp('No fcr_mode given, setting fcr_mode = 2')
end

if nargin >= 4
    fcr_mode = varargin{1};
    if nargin == 6
        data_NREL = varargin{2};
        wind_turbines = varargin{3};
    end
end

if fcr_mode == 1 & nargin ~= 6
    disp('No wind data given, setting fcr_mode = 2')
    fcr_mode = 2;
end

c_slow = c{1};
c_slow = c_slow/sum(c_slow);
if x == 2
    c_fast = c{2};
    c_fast = c_fast/sum(c_fast);
else
    c_fast = [];
end

n_hydro = length(c_slow);
n_wind = length(c_fast); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WIND
s = tf('s');

%% FCR for hydro and wind
Hhydro = []; c_hydro = [];
Hwind = []; c_wind = []; 
c_slow_sum = tf(0);
c_fast_sum = tf(0);
% if ismember(fcr_mode, [1,2])
    for i = 1:(n_hydro+n_wind)
        
        % Hydro turbine
        if i <= n_hydro 
            % Assumes that hydro turbines are listed first
            Ty = GenData.Ta(i);
            Tw = GenData.TW(i);
            g0 = GenData.PGEN(i)/GenData.MBASE(i);
            Tz_hydro = Tw*g0; % 1/z_hydro

            Gservo = 1/(s*Ty+1); 
            Ghydro = (1 - s*Tz_hydro)/(1 + s*Tz_hydro/2);
            Hhydro{i} = Ghydro*Gservo*GenData.MBASE(i);

            % Participation factor need to include the RHP zero
            if fcr_mode == 2
                w_x = 1; % Add a pole at w_x times the freqeuncy as the RHP zero. (OBS design parameter choice)
%                 if i == 1
%                     disp('Increasing hydro bandwidth to improve FCR response')
%                 end
            else
                w_x = 1;
            end
            
            c_hydro{i} = c_slow(i)*(1 - s*Tz_hydro)/(1 + s*Tz_hydro/w_x) ; 

            if i == n_hydro              
                for ii = 1:n_hydro
                    c_slow_sum = c_slow_sum + c_hydro{ii};
                end
                if fcr_mode == 1
                    c_fast_sum = 1- c_slow_sum;
                else
                    c_fast_sum = 0;
                end
            end
        
        % Wind turbine
        elseif i > n_hydro 
            ii = i - n_hydro;

            if iscell(data_NREL)
                data = data_NREL{ii};
            else
                data = data_NREL;
            end
            v_wind = data.ctrl.v_wind;

            z_wind = data.ctrl.z_wind;
%             a =v_wind * data.ctrl.C*(data.ctrl.k_stab) / data.ctrl.x_min;
%             p_wind = a-z_wind;
            p_wind = data.ctrl.p_wind;

            Gwind_ = (s-z_wind)/(s+p_wind); 
            Hwind{ii} = Gwind_*wind_turbines.POPT(ii);
%                 Hwind{ii} = Gwind_;

            c_wind{ii} = c_fast_sum*c_fast(ii)*Gwind_ ;
        end       
    end
% end
 
c_fast_sum = 0;
for i = 1:n_wind
    c_fast_sum = c_fast_sum + c_wind{i};
end

% Normalize so that c_hydro+c_wind = 1
c_sum = (c_slow_sum+c_fast_sum);

if fcr_mode ~= 0
    normalized = false;
    c_result = tf(0);
    if isstable(c_sum) & isstable(1/c_sum)
        disp('Normalizing participation factors')
        normalized = true;   
        for i = 1:n_hydro
            c_hydro{i} = (c_hydro{i}/c_sum);
            c_result = c_result + c_hydro{i};
        end
        for i = 1:n_wind
            c_wind{i} = (c_wind{i}/c_sum);
            c_result = c_result + c_wind{i};
        end
        c_result = minreal(c_result);
    else
        c_result = c_sum;
    end
end

FCR = []; FCR_ideal = []; H = [];
if fcr_mode == 0;
    disp('Ideal FCR at hydro power buses')
    for i = 1:(n_hydro+n_wind)
        if i <= n_hydro
            FCR_ideal{i} = K_ideal*c_slow(i)/GenData.MBASE(i);
        end
        FCR{i} = tf(0);
    end
else  
    for i = 1:(n_hydro+n_wind)
        if i <= n_hydro % Hydro turbine
            ii = i;             
            FCR_ideal{i} = tf(0); 
            
            FCR{i} = minreal(K_ideal*c_hydro{ii}/Hhydro{ii}); 
            H{i} = Hhydro{ii};
            if normalized
                if i == 1
                    disp('Hydro and Wind FCR')
                    disp('---------- Make sure that controllers are reduced to correct order ----------')
                end
            % Model reduction, remove unnecessary controller states
%                 FCR{i} = modelred_hsv(FCR{i},3);
            else
                if i == 1
                    disp('Hydro FCR') 
                end
            end
        elseif i > n_hydro % Wind turbine
            if fcr_mode == 1
                ii = i - n_hydro;
                FCR{i} = minreal(K_ideal*c_wind{ii}/Hwind{ii});  
                H{i} = Hwind{ii};
            else
                FCR{i} = tf(0);
            end
            if normalized
            % Model reduction, remove unnecessary controller states
%                 FCR{i} = modelred_hsv(FCR{i},4);
            end
        end
    end 
end
   
% Ouput participation factors
c_factors = [];
c_factors.c_hydro = c_hydro;
c_factors.c_wind = c_wind;

end
