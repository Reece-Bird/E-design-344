clc;
clearvars;
close all;

% General Settings
resistor_series_name    = 'E24';
capacitor_series_name   = 'E6';    

% Constants and given infomation
Vcc     = 30;
Vee     = -30;
Vt      = 26e-3;

% Assumed/measured Values
B       = 100;
Vbe     = 0.7;
Vces    = 0.2;
Rs      = 14;
RL      = 985e3;

% Chosen Values
fc_dominant = 20;
fc_lesser   = fc_dominant / 10;
A           = -28;
Rc          = 9100;
Re          = 600;

    
%-------------------------DC Analyse-------------------------
Rth = 0.1*B*Re;

% Solve for Re1 and Icq by chosing Re and Rc
syms Re1 Icq
gm      = Icq / Vt;
Rpi     = B / gm;       % One unknown, Icq

Rdc     = Rc + Re;
Rac     = para(Rc, RL) + ((B+1)/B)*Re1; % One unknown, Re1
Rib     = (1+B)*Re1 + Rpi; % Two unknowns, Icq and Re1
Rin     = para(Rib, Rth); % Two unknowns, Icq and Re1

% Equations and solving
eq1     = Icq   == ((Vcc - Vee - Vces) )/ (Rdc + Rac); % Two unknowns, Icq and Re1
eq2     = A     == (-1)*B*para(Rc,RL)*(Rin / ((Rin + Rs)*Rib )); % Two unknowns, Icq and Re1
sol     = solve([eq1, eq2], Re1, Icq); % Solve for both unknowns

% Back substitute and save results
Re1     = double(sol.Re1);
Icq     = double(sol.Icq);
gm      = Icq / Vt;
Rpi     = B / gm;
Re2     = Re - Re1;
Rib     = (1+B)*Re1 + Rpi;
Rin     = para(Rib, Rth);
Rout    = Rc;

% Solve eq1 and eq2
R1 = (Rth * (Vcc - Vee)) / (Icq*(Re + Rth/B) + Vbe);
R2 = (Rth *R1) / (R1 - Rth);

% Calculate resistances seen by caps
R_Cin   = Rin + Rs;
R_Cc    = Rc + RL;
R_Ce    = para(Re1, Re2 + Rin);

% Calculate Capacitances
Cin     = solve_capacitance(R_Cin, fc_dominant);   % NB DO NOT USE EMITTER RESISTOR
Cc      = solve_capacitance(R_Cc, fc_lesser);      % HORRIBLE RESULT, SEE PRAC 2 MEMO     
Ce      = solve_capacitance(R_Ce, fc_lesser);      % FROM ELECTRONICS 315

% Make sure gain is correct
A = (-1)*B*para(Rc,RL)*(Rin / (Rib * (Rin + Rs)));

%--------------------Rounded Values--------------------------
Rc_rounded  = round63(Rc, resistor_series_name);
Re1_rounded = round63(Re1, resistor_series_name);
Re2_rounded = round63(Re2, resistor_series_name);
R1_rounded  = round63(R1, resistor_series_name);
R2_rounded  = round63(R2, resistor_series_name);

% Calculate resistances seen by caps with rounded resistors
Rib_rounded  = (1+B)*Re1_rounded + Rpi;
Rin_rounded  = para(Rib_rounded, Rth);
Rout_rounded = Rc_rounded;
Rth_rounded  = para(R1_rounded, R2_rounded);

% Resolve capacitances with rounded resistors
R_Cin_rounded = Rin_rounded + Rs;
R_Cc_rounded  = Rc_rounded + RL;
R_Ce_rounded  = para(Re1_rounded, Re2_rounded + Rin_rounded);


% Resolve and round capacitances
Cin_rounded = round63(solve_capacitance(R_Cin_rounded, fc_dominant), capacitor_series_name);
Cc_rounded  = round63(solve_capacitance(R_Cc_rounded, fc_lesser), capacitor_series_name);
Ce_rounded  = round63(solve_capacitance(R_Ce_rounded, fc_lesser), capacitor_series_name);

%Solve other values
Rdc_rounded = Rc_rounded + Re1_rounded + Re2_rounded;
Rac_rounded = para(Rc_rounded, RL) + ((B+1)/B)*Re1_rounded;
Icq_rounded = ((Vcc - Vee - Vces))/ (Rdc_rounded + Rac_rounded);
A_rounded   = (-1)*B*para(Rc_rounded,RL)*(Rin_rounded / (Rib_rounded * (Rin_rounded + Rs)));

%-----------------------Print Results-----------------------
fprintf(['----------Component Values-------------- \n\n' ...
         'Unrounded Values       |  Rounded Values       \n' ...
         '-----------------------|------------------------\n' ...
         'RC  = %11.0f Ω    |  %11.0f Ω               \n' ...
         'Re1 = %11.0f Ω    |  %11.0f Ω               \n' ...
         'Re2 = %11.0f Ω    |  %11.0f Ω               \n' ...
         'R1  = %11.0f Ω    |  %11.0f Ω               \n' ...
         'R2  = %11.0f Ω    |  %11.0f Ω               \n' ...
         'Cin = %11.2f μF   |  %11.2f μF              \n' ...
         'Cc  = %11.2f μF   |  %11.2f μF              \n' ...
         'Ce  = %11.2f μF   |  %11.2f μF              \n\n'], ...
        Rc, Rc_rounded, Re1, Re1_rounded, Re2, Re2_rounded, R1, R1_rounded,...
        R2, R2_rounded, Cin*1e6, Cin_rounded*1e6, Cc*1e6, Cc_rounded*1e6, ...
        Ce*1e6, Ce_rounded*1e6);

fprintf(['--------Other Useful Values------------ \n\n' ...
         'Unrounded Values       |  Rounded Values       \n' ...
         '-----------------------|------------------------\n' ...
         'Icq = %11.2f mA   |  %11.2f mA              \n' ...
         'A   = %11.3f      |  %11.3f                 \n' ...
         'Rin = %11.0f Ω    |  %11.0f Ω               \n' ...
         'Rout= %11.0f Ω    |  %11.0f Ω               \n\n'], ...
        Icq*1000, Icq_rounded*1000, A, A_rounded, ...
        Rin, Rin_rounded, Rout, Rout_rounded);


fprintf('---------------------------------- \n\n')
fprintf('Note: This spice is used in a sub circuit with no RL or RS\n\n') 
fprintf(['----------SPICE CODE-------------- \n\n' ...
         'QBJT 5 3 6 2N3904 \n\n' ...
         'RRc vcc 5 %.0f \n' ...
         'RRe1 6 8 %.0f \n' ...
         'RRe2 8 vee %.0f \n' ...
         'RR1 vcc 3 %.0f \n' ...
         'RR2 3 vee %.0f \n\n' ...
         'CCin in 3 %.2fu \n' ...
         'CCc 5 out %.2fu \n' ...
         'CCe 8 0 %.2fu \n'], ...
        Rc, Re1, Re2, R1, R2,...
        Cin * 1e6, Cc * 1e6, Ce * 1e6);

fprintf(['\n--------SPICE CODE ROUNDED-------- \n\n' ...
         'QBJT 5 3 6 2N3904 \n\n' ...
         'RRc vcc 5 %.0f \n' ...
         'RRe1 6 8 %.0f \n' ...
         'RRe2 8 vee %.0f \n' ...
         'RR1 vcc 3 %.0f \n' ...
         'RR2 3 vee %.0f \n\n' ...
         'CCin in 3 %.2fu \n' ...
         'CCc 5 out %.2fu \n' ...
         'CCe 8 0 %.2fu \n'], ...
        Rc_rounded, Re1_rounded , Re2_rounded, R1_rounded, R2_rounded,...
        Cin_rounded * 1e6, Cc_rounded * 1e6, Ce_rounded * 1e6);

function Req = para(R1, R2)
    Req = (R1 * R2) / (R1 + R2);
end

function C = solve_capacitance(R, fc)
    C = 1 / (2*pi*R*fc);
end
