clc;
clearvars;
close all;


%--------------------------Main settings-----------------------------------
% Low and High pass settinsg
fc_low      = 240;          % Cutoff frequency low in Hz
fc_high     = 4000 * 0.5;   % Cutoff frequency high in Hz (ignore 0.5)
c1_low      = 0.1e-6;       % Choose standard values for C1 and C2 for LPF
c2_low      = 0.68e-6;      
c_high      = 0.01e-6;      % Choose values for all 3 caps on high pass filter

% Band pass settings
fc_low_bp   = 350 * 0.5;    % Choose lower corner frequency for BPF
fc_high_bp  = 2450;         % Choose higher corner frequency for BPF
c1_low_bp   = 0.01e-6;      % Choose standard values for C1 and C2 for LPF
c2_low_bp   = 0.068e-6;     
c_high_bp   = 0.1e-6;       % Choose values for all 3 caps on high pass filter

% General settings
second_order            = 0;
fourth_order            = 1;
resistor_series_name    = 'E24';
capacitor_series_name   = 'E6';     % Pretty redundant as you choose all cap values
K                       = -1;       % Gain


% Fourth Order MFB filter
if fourth_order
    FSF1    = 1.4192;
    Q1      = 0.5219;  % Lowest Q goes first
    FSF2    = 1.5912;
    Q2      = 0.8055;
    order   = 4;

    disp('------------------------------------------------------------')
    disp('-----------------FOURTH ORDER FILTERS-----------------------')
    disp('------------------------------------------------------------')
    disp('--------------------Stage 1 values--------------------------')
    [H_lpf1, H_hpf1] = filter_component_solver(FSF1, Q1, c1_low, c2_low, c_high, K, fc_low, fc_high, order, resistor_series_name, capacitor_series_name);
    disp('--------------------Stage 2 values--------------------------')
    [H_lpf2, H_hpf2] = filter_component_solver(FSF2, Q2, c1_low, c2_low, c_high, K, fc_low, fc_high, order, resistor_series_name, capacitor_series_name);

    disp('------------------------------------------------------------')
    disp('-------------------BANDPASS Values--------------------------')
    disp('------------------------------------------------------------')
    fprintf('K, Q1 and Q2 changed to %f, %f and %f respecitively',  K, Q1, Q2);
    disp('------------------------------------------------------------')
    disp('--------------------Stage 1 values--------------------------')
    [H_lpf1_bp, H_hpf1_bp] = filter_component_solver(FSF1, Q1, c1_low_bp, c2_low_bp, c_high_bp, K, fc_high_bp, fc_low_bp, order, resistor_series_name, capacitor_series_name);
    disp('--------------------Stage 2 values--------------------------')
    [H_lpf2_bp, H_hpf2_bp] = filter_component_solver(FSF2, Q2, c1_low_bp, c2_low_bp, c_high_bp, K, fc_high_bp, fc_low_bp, order, resistor_series_name, capacitor_series_name);
    disp('------------------------------------------------------------')

    H_lpf_4n        = H_lpf1 * H_lpf2;
    H_hpf_4n        = H_hpf1 * H_hpf2;

    H_lpf_4n_bp     = H_lpf1_bp * H_lpf2_bp;
    H_hpf_4n_bp     = H_hpf1_bp * H_hpf2_bp;

    [meas_fc_low]   = plot_bode(H_lpf_4n, "4th Order Low Pass Filter");
    [meas_fc_high]  = plot_bode(H_hpf_4n, "4th Order High Pass Filter");
    [meas_fc_bp]    = plot_bode(H_hpf_4n_bp * H_lpf_4n_bp, "4th Order Band Pass Filter" );
    plot_bode(H_lpf_4n_bp, "4th Order Low Pass Filter for Bandpass");
    plot_bode(H_hpf_4n_bp, "4th Order High Pass Filter for Bandpass");
    plot_bode(H_lpf_4n + H_hpf_4n_bp * H_lpf_4n_bp + H_hpf_4n, "All filters combined" );

    disp('-------------------Important values-------------------------')
    fprintf('Low Pass filter will have a corner frequency at %f\n',  meas_fc_low);
    fprintf('High Pass filter will have a corner frequency at %f\n', meas_fc_high);
    fprintf('Band Pass filter will have a corner frequency at %f\n', meas_fc_bp(1));
    fprintf('Band Pass filter will have a corner frequency at %f\n', meas_fc_bp(2));
    fprintf('Band Pass filter will have a center frequency at %f\n', sqrt(meas_fc_bp(1)*meas_fc_bp(2)));
    disp('------------------------------------------------------------')
end

function [H1, H2] = filter_component_solver(FSF, Q, c1_low, c2_low, c_high, K, fc_low, fc_high, order, resistor_series_name, capacitor_series_name)
    disp('Low pass filter solutions')
    disp('------------------------------------------------------------')
    
    % If you have a cascading network distribute the gain properly
    if order == 4
        K = -1 * sqrt(-1*K);
    end

    % Set up all variables for simultanious equations
    syms m r1 r2 r3
    
    %Solve for m ratio
    n       = c2_low/c1_low;
    eq_m    = Q == sqrt(m * n) / (1 + 2 * m);
    sol_m   = solve(eq_m, m);
    m_value = double(sol_m); % Choose the positive solution
    m       = m_value(1);

    % Solve equations for best solution
    eq1     = K             == -r2 / r1;
    eq2     = r3            == m*r2;
    eq3     = fc_low * FSF  == 1 / (2 * pi * r2 * c1_low * sqrt(m*n));
    sol     = solve([eq1, eq2, eq3], [r1, r2, r3]);
    solr1   = double(sol.r1);
    solr2   = double(sol.r2);
    solr3   = double(sol.r3);

   % Calculate resistor values
    resistor_list_2 = [solr1, solr2, solr3];
    
    % Display resistor values in formatted output
    fprintf('Order %d Resistors (Best solution):\n', order);
    fprintf('R1 = %6.0f Ω\n', resistor_list_2(1));
    fprintf('R2 = %6.0f Ω\n', resistor_list_2(2));
    if length(resistor_list_2) > 2
        fprintf('R3 = %6.0f Ω\n', resistor_list_2(3));
    end
    
    %Make resistor values standard and check that they are real
    H1 =  pick_standard_components(solr1, solr2, solr3, c1_low, c2_low, order, resistor_series_name, 'low_pass');
    
    disp('------------------------------------------------------------')
    disp('High pass filter solutions')
    disp('------------------------------------------------------------')
    
    omega_c = (2*pi*fc_high*FSF); %TODO Theres a 0.5 that needs to be with the fc term here, idk where it is
    a0      = omega_c^2;
    a1      = omega_c / Q;

    eq1     = 1 / (r1*r2*c_high*c_high)                       == a0;
    eq2     = (c_high + c_high + c_high) / (r1*c_high*c_high) == a1;
    sol     = solve([eq1, eq2], [r1, r2]);
    solr1   = double(sol.r1);
    solr2   = double(sol.r2);

    % Display capacitor values
    fprintf('Order %d Capacitors:\n', order);
    fprintf('C1 = %6.2f μF\n', c_high*1e6);
    fprintf('C2 = %6.2f μF\n', c_high*1e6);
    fprintf('C3 = %6.2f μF\n', c_high*1e6);
    
    % Display resistor values
    fprintf('Order %d Resistors:\n', order);
    fprintf('R1 = %6.0f Ω\n', solr1);
    fprintf('R2 = %6.0f Ω\n', solr2)

     %Make capacitor values standard and check that they are real
     H2 = pick_standard_components(c_high, c_high, c_high, solr1, solr2, order, capacitor_series_name, "high_pass");
end

function H = pick_standard_components(a, b, c, d, e, order, resistor_series_name, filter_type)
    % Round component values according to the series, check if it is either
    % a low pass or high pass filter, create their respective transfer
    % functions and return them, all while adding descriptive text

    if isreal(a) && isreal(b) && isreal(c)  % Solutions real?
        % Round values to standard components
        a_std = round63(a, resistor_series_name);
        b_std = round63(b, resistor_series_name);
        c_std = round63(c, resistor_series_name);
        d_std = round63(d, resistor_series_name);
        e_std = round63(e, resistor_series_name);

        if filter_type == "low_pass"
            displ_1 = "Capacitors";
            displ_2 = "Resistors";
            standard_resistor_list  = [a_std, b_std, c_std];
            standard_capacitor_list = [d_std, e_std];
            H = tf_low_pass(a_std, b_std, c_std, d_std, e_std);
        elseif filter_type == "high_pass"
            displ_1 = "Capacitors";
            displ_2 = "Resistors";
            standard_capacitor_list = [a_std, b_std, c_std];
            standard_resistor_list  = [d_std, e_std];
            H = tf_high_pass(a_std, b_std, c_std, d_std, e_std);
        end

        fprintf('Order %d %s (Standard Components):\n', order, displ_1);
        fprintf('C1 = %6.3f μF\n', standard_capacitor_list(1) * 1e6);
        fprintf('C2 = %6.3f μF\n', standard_capacitor_list(2) * 1e6);
        if length(standard_capacitor_list) > 2
            fprintf('C3 = %6.3f μF\n', standard_capacitor_list(3) * 1e6);
        end

        fprintf('Order %d %s (Standard Components):\n', order, displ_2);
        fprintf('R1 = %6.0f Ω\n', standard_resistor_list(1));
        fprintf('R2 = %6.0f Ω\n', standard_resistor_list(2));
        if length(standard_resistor_list) > 2
            fprintf('R3 = %6.0f Ω\n', standard_resistor_list(3));
        end
    else
        standard_resistor_list = 'No standard components available';
        fprintf('%s\n\n', standard_resistor_list);
        error('Cannot plot graphs as no real components exist');
    end
end


function H = tf_low_pass(r1, r2, r3, c1, c2)
    s = tf('s');
    H = (-r2/r1) / ( s^2*(r2*r3*c1*c2) + s*(r3*c1 + r2*c1 + ((r2*r3*c1)/r1  ) ) + 1 );
end

function H  = tf_high_pass(c1, c2, c3, r1, r2)
    s = tf('s');
    H = ((-1*(c1/c2))*s^2) / (s^2 + s*( (c1 + c2 + c3 ) / (r1*c2*c3) ) + 1 / (r1*r2*c2*c3));
end

function [f_3dB] = plot_bode(H, title_text)
    
    % Frequency vector in rad/s
    omega = logspace(log10(10), log10(1e6), 1000); % From 10 rad/s to 10000000 rad/s
    
    % Compute the frequency response
    [mag, phase] = bode(H, omega);

    % Find max in frequency response
    

    
    % Convert omega from rad/s to Hz
    freq_Hz = omega / (2 * pi);
    mag_dB = 20*log10(squeeze(mag));
    figure()
    sgtitle(title_text)

    % Find the maximum magnitude
    max_mag = max(mag_dB);

    % Plot magnitude response
    subplot(2, 1, 1);
    semilogx(freq_Hz, mag_dB);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    ylim([-6, 2]); % Set the lower limit to -60 dB
    grid on;
    
    % Plot phase response
    subplot(2, 1, 2);
    semilogx(freq_Hz, squeeze(phase));
    xlabel('Frequency (Hz)');
    ylabel('Phase (degrees)');
    grid on;

    % Find -3dB points
    if strcmp(title_text, "4th Order High Pass Filter") || strcmp(title_text, "2nd Order High Pass Filter")
        idx_3dB = find(mag_dB >= (max_mag - 3), 1, 'first'); 
        f_3dB = freq_Hz(idx_3dB);
    end
    if strcmp(title_text, "4th Order Low Pass Filter") || strcmp(title_text, "2nd Order Low Pass Filter")
        idx_3dB = find(mag_dB <= (max_mag - 3), 1, 'first'); 
        f_3dB = freq_Hz(idx_3dB);
    end
    if strcmp(title_text, "4th Order Band Pass Filter") || strcmp(title_text, "2nd Order Band Pass Filter")
        idx_3dB_first = find(mag_dB >= (max_mag - 3), 1, 'first'); 
        idx_3dB_last = find(mag_dB(idx_3dB_first:end) <= (max_mag - 3), 1, 'first');
        f_3dB_first = freq_Hz(idx_3dB_first);
        f_3dB_last  = freq_Hz(idx_3dB_last + idx_3dB_first);
        f_3dB = [f_3dB_first, f_3dB_last];
    end
end