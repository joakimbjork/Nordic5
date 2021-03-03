%% FCR controller design
s = tf('s');


c_tot = sum(c_fast);
c_fast = c_fast/c_tot;

n_hydro = length(c_slow);
n_wind = length(c_fast); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WIND

Hhydro = []; c_hydro = [];
Hwind = []; c_wind = [];
if only_Hydro
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
            w_x = 10; % Add a pole at w_x times the freqeuncy as the RHP zero. 
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
            if iscell(data_NREL)
                data_ = data_NREL{ii};
                v_ = v_wind(ii);
            else
                data_ = data_NREL;
                v_ = v_wind;
            end
            z_wind = data_.ctrl.z_wind;
            a =v_ * data_.ctrl.C*(data_.ctrl.k_stab) / data_.ctrl.x_min;
            p_wind = a-z_wind;
            Gwind_ = (s-z_wind)/(s+p_wind); 
            Hwind{ii} = Gwind_*wind_turbines.POPT(ii);
            c_wind{ii} = 0;
        end       
    end
else
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
            if iscell(data_NREL)
                data_ = data_NREL{ii};
                v_ = v_wind(ii);
            else
                data_ = data_NREL;
                v_ = v_wind;
            end
            z_wind = data_.ctrl.z_wind;
            a =v_ * data_.ctrl.C*(data_.ctrl.k_stab) / data_.ctrl.x_min;
            p_wind = a-z_wind;
            Gwind_ = (s-z_wind)/(s+p_wind); 
            Hwind{ii} = Gwind_*wind_turbines.POPT(ii);
            c_wind{ii} = c_fast_sum*c_fast(ii)*Gwind_;
        end       
    end
end
    
c_fast_sum = tf(0);
for i = 1:n_wind
    c_fast_sum = c_fast_sum + c_wind{i};
end

% Normalize so that c_hydro+c_wind = 1
c_sum = minreal(c_slow_sum+c_fast_sum);

normalized = false;
if isstable(c_sum) & isstable(1/c_sum)
    disp('Normalizing participation factors')
    normalized = true;
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
    if normalized
        disp('Hydro and Wind FCR')
        disp('---------- Make sure that controllers are reduced to correct order ----------')
    else
        disp('Hydro FCR')
    end    
    for i = 1:n_hydro
        FCR_ideal{i} = tf(0);        
        if normalized
            % Model reduction, remove unnecessary controller states
            FCR{i} = modelred_hsv(FCR{i},3);
        end
    end
    for i = (n_hydro+1):(n_hydro+n_wind)
        if normalized
        % Model reduction, remove unnecessary controller states
            FCR{i} = modelred_hsv(FCR{i},4);
        end
       
    end
end
