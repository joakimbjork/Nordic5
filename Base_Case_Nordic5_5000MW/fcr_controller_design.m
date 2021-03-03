%% FCR controller design
s = tf('s');

fdcmax = 0.4; % Max allowed steady state freq deviation
d_dim = 1400;
R = (d_dim/fdcmax - load_damping)*fn; % MW/(Hz/50) to provide with FCR


% K_ideal = R/(s*7.2+1); % Ideal FCR
K_ideal = R*(s*6.5+1)/((s*2+1)*(s*17+1)); % Ideal FCR

c_slow = [4230;2500;800]; % Distribution of FCR-N in Saarinen's paper
c_tot = sum(c_slow);
c_slow = c_slow/c_tot;

c_fast = [1]; % Distribution of fast FCR
c_tot = sum(c_fast);
c_fast = c_fast/c_tot;

n_hydro = length(c_slow);
n_wind = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WIND

Hhydro = []; c_hydro = [];
Hwind = []; c_wind = [];
for i = 1:(n_hydro+n_wind)
    if i <= n_hydro % Hydro turbine
        % Assumes that hydro turbines are listed first
        Ty = GenData.Ta(i);
        Tw = GenData.TW(i);
        g0 = GenData.PGEN(i)/GenData.MBASE(i);
        Tz_hydro = Tw*g0; % 1/z_hydro

        Gservo = 1/(s*Ty+1); 
        Ghydro = (1 - s*Tz_hydro)/(1 + s*Tz_hydro/2);
        Hhydro{i} = Ghydro*Gservo*Sgen(i);
        
        % Participation factor need to include the RHP zero
        w_x = 1; % Add a pole at w_x times the freqeuncy as the RHP zero. 
        c_hydro{i} = c_slow(i)*(1 - s*Tz_hydro)/(1 + s*Tz_hydro/w_x); 
        
        if i == n_hydro
            c_slow_sum = tf(0);
            for ii = 1:n_hydro
                c_slow_sum = c_slow_sum + c_hydro{ii};
            end
            c_fast_sum = 1- c_slow_sum;
        end
       
    elseif i > n_hydro % Wind turbine
        ii = i - n_hydro;
        z_wind = data_NREL.ctrl.z_wind;
        a =v_wind * data_NREL.ctrl.C*(data_NREL.ctrl.k_stab) / data_NREL.ctrl.x_min;
        p_wind = a-z_wind;
        Gwind = (s-z_wind)/(s+p_wind);
        Hwind{ii} = Gwind*wind_turbines.POPT(ii);
        
        % Participation factor need to include the RHP zero
        c_wind{ii} = c_fast_sum*c_fast(ii)*Gwind;
    end       
end

c_fast_sum = tf(0);
for i = 1:n_wind
    c_fast_sum = c_fast_sum + c_wind{i};
end

% Normalize so that c_hydro+c_wind = 1
c_sum = minreal(c_slow_sum+c_fast_sum);

if isstable(c_sum) & isstable(1/c_sum)
    disp('Normalizing participation factors')
    c_result = tf(0);
    for i = 1:n_hydro
        c_hydro{i} = c_hydro{i}/c_sum;
        c_result = c_result + c_hydro{i};
    end
    for i = 1:n_wind
        c_wind{i} = c_wind{i}/c_sum;
        c_result = c_result + c_wind{i};
    end
    c_result = minreal(c_result);
else
    c_result = c_sum;
end

FCR = [];
for i = 1:(n_hydro+n_wind)
    if i <= n_hydro % Hydro turbine
        ii = i;
        FCR{i} = minreal(K_ideal*c_hydro{ii}/Hhydro{ii});    
    elseif i > n_hydro % Wind turbine
        ii = i - n_hydro;
        FCR{i} = minreal(K_ideal*c_wind{ii}/Hwind{ii});  
    end
end       

FCR_ideal = [];
if ideal_FCR
    disp('Ideal FCR at hydro power buses')
    for i = 1:(n_hydro+n_wind)
        if i <= n_hydro
            FCR_ideal{i} = K_ideal*c_slow(i)/Sgen(i);
        end
        FCR{i} = tf(0);
    end
else
    disp('Hydro FCR')
    for i = 1:n_hydro
        FCR_ideal{i} = tf(0);
    end
end
