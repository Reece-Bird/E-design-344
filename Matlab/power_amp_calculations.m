clc;
clearvars;
close all;

% General Settings
resistor_series_name    = 'E24';
capacitor_series_name   = 'E6';    

% Constants and given infomation

% Assumed/measured Values

% Chosen Values
Is              =  6e-3;
Id              = Is*(22/30);
Bs              = 70;
Bd              = 40;
Imax            = 6.5;
V_eb_be_best    = 0.62;         % Chose 0.62V for design
V_eb_be         = 0.65;         
V_Re            = 1;
Vcc             = 30;
Vee             = -30;

Ta              = 25;
theta_jc_tip41  = 1.92;
theta_ja_tip41  = 62.5;
theta_jc_tip35  = 1;
theta_ja_tip35  = 35.7;
theta_ca_metal  = 1.5;

Pd_q2           = 29.5025 ; %Darlington source 2
Pd_q1           = 0.434893; %Darlington source 1
Pd_q3           = 28.4827 ; %Darlington drain 2
Pd_q4           = 0.599681; %Darlington drain 1
Pd_q7           = 0.209769; %Current source
Pd_q8           = 0.150993; %Current drain
Pd_q9           = 0.228872; %Buffer source 2
Pd_q10          = 0.055025; %Buffer source 1
Pd_q11          = 0.137703; %Buffer drain 2
Pd_q12          = 0.050542; %Buffer drain 1

% Geneal Calculations
P   = 180;
RL  = 8;
V   = sqrt(2*P*RL);
I   = V / RL;

[R1, R2, Re] = current_supply_solver(Is, Bs, V_eb_be_best, resistor_series_name);
[R36, R37, R38] = current_supply_solver(Id, Bd, V_eb_be_best, resistor_series_name);

[Rsense_1, Imax] = output_resistor_solver(Imax, V_eb_be, resistor_series_name);
Rsense_2 = Rsense_1;


T_q2    = find_junction_tempreture(Ta, Pd_q2, theta_jc_tip35 + theta_ca_metal);
T_q1    = find_junction_tempreture(Ta, Pd_q1, theta_ja_tip41);
T_q3    = find_junction_tempreture(Ta, Pd_q3, theta_jc_tip35 + theta_ca_metal);
T_q4    = find_junction_tempreture(Ta, Pd_q4, theta_ja_tip41);
T_q7    = find_junction_tempreture(Ta, Pd_q7, theta_ja_tip41);
T_q8    = find_junction_tempreture(Ta, Pd_q8, theta_ja_tip41);
T_q9    = find_junction_tempreture(Ta, Pd_q9, theta_ja_tip41);
T_q10   = find_junction_tempreture(Ta, Pd_q10, theta_ja_tip41);
T_q11   = find_junction_tempreture(Ta, Pd_q11, theta_ja_tip41);
T_q12   = find_junction_tempreture(Ta, Pd_q12, theta_ja_tip41);



%-----------------------Print Results-----------------------
fprintf(['------------------------------------------------\n' ...
         'Rounded Values\n' ...
         '------------------------------------------------\n' ...
         'R1 = %11.0f Ω\n'...   
         'R2 = %11.0f Ω\n'...  
         'Re = %11.0f Ω\n\n'...
         'R1 = %11.0f Ω\n',...
         'R2 = %11.0f Ω\n'...
         'Re = %11.0f Ω\n\n'...
         'Rsense 1 = %6.2f Ω\n'...  
         'Rsense 2 = %6.2f Ω\n'... 
         'Imax = %10.4f A\n'],... 
          R1, R2, Re, R36, R37, R38, Rsense_1, Rsense_2, Imax)

fprintf('\n\nJunction Temperatures:\n');
fprintf('T_q2:  %.2f°C (Heatsink on) \n', T_q2);
fprintf('T_q1:  %.2f°C \n', T_q1);
fprintf('T_q3:  %.2f°C (Heatsink on) \n', T_q3);
fprintf('T_q4:  %.2f°C \n', T_q4);
fprintf('T_q7:  %.2f°C \n', T_q7);
fprintf('T_q8:  %.2f°C \n', T_q8);
fprintf('T_q9:  %.2f°C \n', T_q9);
fprintf('T_q10: %.2f°C \n', T_q10);
fprintf('T_q11: %.2f°C \n', T_q11);
fprintf('T_q12: %.2f°C \n', T_q12)

% Current Source 
function [R1, R2, Re] = current_supply_solver(Is, Bs, V_eb_be, resistor_series_name)
    Vcc     = 30;
    VRe     = 1;
    Re      = VRe / Is;
    Re      = round63(Re, resistor_series_name);
    VRe     = Is * Re;
    Rth     = 0.1*Bs*Re;
    Vth     = Vcc - VRe - V_eb_be - (Is/Bs)*Rth;
    
    
    syms R1 R2
    eq1     = Vth == (R2 / (R1 + R2)) * Vcc;
    eq2     = Rth == (R1 * R2) / (R1 + R2);
    sol     = solve([eq1, eq2], R1, R2);
    R1      = double(sol.R1);
    R2      = double(sol.R2);

    R1      = round63(R1, resistor_series_name);
    R2      = round63(R2, resistor_series_name);
end

function [R, Imax_new] = output_resistor_solver(Imax, Vbe, resistor_series_name)
    R           = Vbe / Imax;
    R           = round63(R, resistor_series_name);
    Imax_new    = Vbe / R;
end

function [temp] = find_junction_tempreture(Ta, Pd, theta)
    temp = Ta + theta * Pd;
end



