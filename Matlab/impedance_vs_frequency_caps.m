% Define component values
L = 1e-3;  % 1 mH
C = 1e-9;  % 100 µF

% Define frequency range (1 Hz to 1 MHz)
f = logspace(0, 6, 500);  % logarithmic scale
omega = 2 * pi * f;       % angular frequency

% Calculate impedance for LC circuit
Z_L = 1i * omega * L;         % Impedance of inductor
Z_C = 1 ./ (1i * omega * C);  % Impedance of capacitor
Z_LC = abs(Z_L + Z_C);        % Total impedance of series LC circuit

% Plot the impedance vs frequency
figure;
semilogx(f, Z_LC, 'r', 'LineWidth', 1.5);
grid on;

% Labels and title
xlabel('Frequency (Hz)');
ylabel('Impedance (Ohms)');
title('Impedance vs Frequency for LC Circuit (1 mH and 100 µF)');

